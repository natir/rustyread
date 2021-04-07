//! Simulate reads

/* mod declaration */
mod description;
mod error;

/* standard use */
use std::io::Write;

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;
use rand::Rng;
use rand::RngCore;
use rand::SeedableRng;
use rayon::prelude::*;

/* local use */
use crate::cli;
use crate::model;
use crate::references::*;
use description::{Description, Origin};

/// main simulate function
#[cfg(not(tarpaulin))]
pub fn simulate(params: cli::simulate::Command) -> Result<()> {
    let mut main_rng = if let Some(seed) = params.seed {
        rand::rngs::StdRng::seed_from_u64(seed)
    } else {
        rand::rngs::StdRng::seed_from_u64(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)?
                .as_secs(),
        )
    };

    log::info!("Start read reference");
    let references = References::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(std::fs::File::open(
            params.reference_path,
        )?)))?
        .0,
    )?;
    log::info!("End read reference");

    log::info!("Start init lenght model");
    let length = model::Length::new(params.length.0 as f64, params.length.1 as f64)?;
    log::info!("End init lenght model");

    log::info!("Start init identity model");
    let identity = model::Identity::new(
        params.identity.0 as f64,
        params.identity.1 as f64,
        params.identity.2 as f64,
    )?;
    log::info!("End init length model");

    log::info!("Start read error model");
    let error = model::Error::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(std::fs::File::open(
            params.error_model,
        )?)))?
        .0,
        &mut main_rng,
    )?;
    let k = error.k();
    log::info!("End read error model");

    log::info!("Start read quality score model");
    let qscore = model::Quality::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(std::fs::File::open(
            params.qscore_model,
        )?)))?
        .0,
    )?;
    log::info!("End read quality score model");

    let total_base = params
        .quantity
        .number_of_base(references.0.iter().map(|x| x.seq.len() as u64).sum());
    log::info!("Total number {}", total_base);
    let base_per_thread = total_base / rayon::current_num_threads() as u64;
    log::info!("Number of base per thread {}", base_per_thread);

    log::info!("Start generate sequences");
    let parameter = (0..rayon::current_num_threads())
        .map(|_| (base_per_thread, main_rng.next_u64()))
        .collect::<Vec<(u64, u64)>>();
    let sequences: Vec<Vec<(Description, Seq, Quality)>> = parameter
        .par_iter()
        .map(|(base, seed)| {
            worker(
                *base,
                &references,
                &length,
                &identity,
                &error,
                &qscore,
                rand::rngs::StdRng::seed_from_u64(*seed),
            )
            .unwrap()
        })
        .collect();
    log::info!("End generate sequences");

    log::info!("Start write sequences");
    let mut output = std::io::BufWriter::new(std::fs::File::create(params.output_path)?);

    for s in sequences {
        for (comment, seq, qual) in s {
            writeln!(
                output,
                "@{} {}\n{}\n+\n{}\n",
                uuid::Uuid::new_v4().to_hyphenated(),
                comment,
                std::str::from_utf8(&seq[k..(seq.len() - k)])?, // begin and end of fragment is just random base
                std::str::from_utf8(&qual[k..seq.len() - k])?
            )?;
        }
    }
    log::info!("End write sequences");

    Ok(())
}

type Seq = Vec<u8>;
type Quality = Vec<u8>;

/// Function realy generate read
#[cfg(not(tarpaulin))]
fn worker(
    target: u64,
    references: &References,
    length_model: &model::Length,
    identity_model: &model::Identity,
    error_model: &model::Error,
    qscore_model: &model::Quality,
    mut rng: rand::rngs::StdRng,
) -> Result<Vec<(Description, Seq, Quality)>> {
    let mut data = Vec::new();

    let k = error_model.k();
    let mut generate = 0;
    let ref_dist = WeightedIndex::new(&references.1).unwrap();

    while generate < target {
        let local_ref = &references.0[ref_dist.sample(&mut rng)];
        let strand = ['+', '-'][rng.gen_range(0..1) as usize];
        let start_pos = rng.gen_range(0..local_ref.seq.len()) as usize;
        let length = length_model.get_length(&mut rng) as usize;
        let identity = identity_model.get_identity(&mut rng);

        let end_pos = if start_pos + length > local_ref.seq.len() {
            local_ref.seq.len() - 1
        } else {
            start_pos + length
        };

        let mut raw_fragment = Vec::with_capacity(end_pos - start_pos + error_model.k() * 2);
        raw_fragment.extend(crate::random_seq(k, &mut rng));
        raw_fragment.extend(&local_ref.seq[start_pos..end_pos]);
        raw_fragment.extend(crate::random_seq(k, &mut rng));

        let err_fragment: Vec<u8> =
            error::add_error(identity, error_model, &raw_fragment, &mut rng);

        let (real_id, quality) =
            generate_quality(&raw_fragment, &err_fragment, qscore_model, &mut rng)?;

        log::debug!("target {} real {}", identity, real_id);

        let ori = Origin::new(
            local_ref.id.clone(),
            strand,
            start_pos,
            end_pos,
            false,
            false,
        );
        let des = Description::new(ori, None, length, real_id * 100.0);

        data.push((des, err_fragment, quality));

        generate += (end_pos - start_pos) as u64;
    }

    Ok(data)
}
/// Generate quality string
fn generate_quality(
    raw: &[u8],
    err: &[u8],
    model: &model::Quality,
    rng: &mut rand::rngs::StdRng,
) -> Result<(f64, Quality)> {
    let mut qual = Vec::with_capacity(err.len());
    let margin = (model.max_k() - 1) / 2;

    let (edit, cigar) = crate::alignment::align(err, raw);

    for i in 0..cigar.len() {
        if cigar[i] == b'D' {
            continue;
        }

        let (start, end) = if i < margin {
            (0, i + i + 1)
        } else if i >= cigar.len() - margin {
            (i - (cigar.len() - i - 1), cigar.len())
        } else {
            (i - margin, i + margin + 1)
        };

        qual.push(model.get_qscore(&cigar[start..end], rng)?);
    }

    Ok((1.0 - (edit as f64 / err.len() as f64), qual))
}

#[cfg(test)]
mod t {
    use super::*;

    fn init() {
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    static MODEL: &[u8] = b"overall;3;1:0.000065,2:0.002683,3:0.011969,4:0.022982,5:0.03079,6:0.036889,7:0.042461,8:0.047495,9:0.051971,10:0.055876,11:0.058956,12:0.0606,13:0.060621,14:0.059446,15:0.056487,16:0.052028,17:0.047379,18:0.043445,19:0.040219,20:0.037767,21:0.035767,22:0.033816,23:0.031369,24:0.027969,25:0.023457,26:0.018417,27:0.008983,28:0.000093,
I;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
X;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
=;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
";

    #[test]
    fn quality() {
        init();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let qual = model::Quality::from_stream(MODEL).unwrap();

        let raw = b"CGACCCTTGCGCAGTCACTGAAAGAGATCAAAGGCGTGCAGTGGATCCCGGACTGGGGCCAGGCGCGTATCGAGTCGATGGTTGCTAACCGTCCTGACTGGTGTATCGCAAAAT";
        let err = b"CGACCCTTGACGCAGTCACTGGAAGAAGTCAAAGGCTGCAGTGGATCCCGGATGGGGCCAGGGTATGAGTCGATGGTTGCTAACCGCCTGACTGGTGTATCGCAAAAT";

        let (identity, qual) = generate_quality(raw, err, &qual, &mut rng).unwrap();

        assert_eq!(0.8981481481481481, identity);

        assert_eq!(err.len(), qual.len());
        assert_eq!(b",,-*%+/2'#5,*',($,%.2,,&*+0(*.+353.(-*&+9-+%'*+72,01-(**((-+7'+,'*(,.&+-3,++5%+/)**/.-&*/5&2,00)'03)')5/)'2'".to_vec(), qual);
    }
}
