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
