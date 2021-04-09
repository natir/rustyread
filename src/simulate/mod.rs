//! Simulate reads

/* mod declaration */
mod description;
mod error;
mod quality;

/* standard use */
use std::io::Write;

/* crate use */
use anyhow::Result;
use rand::Rng;
use rand::RngCore;
use rand::SeedableRng;
use rayon::prelude::*;

/* local use */
use crate::cli;
use crate::model;
use crate::references::*;
use description::{Description, Origin};

#[cfg(not(tarpaulin_include))]
/// main simulate function
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

    log::info!("Start init adapter model");
    let adapter = model::Adapter::new(
        params.start_adapter_seq.as_bytes().to_vec(),
        params.end_adapter_seq.as_bytes().to_vec(),
        params.start_adapter.0 as f64,
        params.start_adapter.1 as f64,
        params.end_adapter.0 as f64,
        params.end_adapter.1 as f64,
    )?;
    log::info!("End init adapter model");

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

    let total_base = params.quantity.number_of_base(
        references
            .sequences
            .iter()
            .map(|x| x.seq.len() as u64)
            .sum(),
    );
    log::info!("Target number of base {}", total_base);

    log::info!("Start generate sequences");
    let sequences: Vec<(Description, Seq, Quality)> =
        LenIdSeed::new(total_base, &length, &identity, &mut main_rng)
            .par_bridge()
            .map(|(len, id, seed)| {
                generate_read(
                    &references,
                    len as usize,
                    id,
                    &adapter,
                    &error,
                    &qscore,
                    rand::rngs::StdRng::seed_from_u64(seed),
                )
                .unwrap()
            })
            .collect();
    log::info!("End generate sequences");

    log::info!("Start write sequences");
    let mut output = std::io::BufWriter::new(std::fs::File::create(params.output_path)?);

    for (comment, seq, qual) in sequences {
        writeln!(
            output,
            "@{} {}\n{}\n+\n{}\n",
            uuid::Uuid::new_v4().to_hyphenated(),
            comment,
            std::str::from_utf8(&seq[k..(seq.len() - k)])?, // begin and end of fragment is just random base
            std::str::from_utf8(&qual[k..seq.len() - k])?
        )?;
    }
    log::info!("End write sequences");

    Ok(())
}

type Seq = Vec<u8>;
type Quality = Vec<u8>;

#[cfg(not(tarpaulin_include))]
/// Function realy generate read
fn generate_read(
    references: &References,
    length: usize,
    identity: f64,
    adapter_model: &model::Adapter,
    error_model: &model::Error,
    qscore_model: &model::Quality,
    mut rng: rand::rngs::StdRng,
) -> Result<(Description, Seq, Quality)> {
    let k = error_model.k();

    // Generate fragment
    let (ref_fragment, origin) = get_fragment(length, references, &mut rng);
    let start_adapter = adapter_model.get_start(&mut rng);
    let end_adapter = adapter_model.get_end(&mut rng);

    let mut raw_fragment = Vec::with_capacity(
        ref_fragment.len() + error_model.k() * 2 + start_adapter.len() + end_adapter.len(),
    );
    raw_fragment.extend(crate::random_seq(k, &mut rng));
    raw_fragment.extend(&start_adapter);
    raw_fragment.extend(&ref_fragment);
    raw_fragment.extend(&end_adapter);
    raw_fragment.extend(crate::random_seq(k, &mut rng));

    // Add error in fragment and produce quality
    let (err_fragment, diffpos) = error::add_error(identity, error_model, &raw_fragment, &mut rng);

    let (real_id, mut quality) = quality::generate_quality(
        &raw_fragment,
        &err_fragment,
        qscore_model,
        diffpos,
        &mut rng,
    )?;

    if quality.len() != err_fragment.len() {
        log::warn!("read and quality string have different length, if you use seed please send all run information to author.");
        quality.resize(err_fragment.len(), b'!');
    }

    // Generate information on read
    let des = Description::new(origin, None, length, real_id * 100.0);

    Ok((des, err_fragment, quality))
}

fn get_fragment(
    length: usize,
    references: &References,
    rng: &mut rand::rngs::StdRng,
) -> (Vec<u8>, Origin) {
    get_ref_fragment(length, references, rng)
}

fn get_ref_fragment(
    length: usize,
    references: &References,
    rng: &mut rand::rngs::StdRng,
) -> (Vec<u8>, Origin) {
    let (id, local_ref, strand) = &references.get_reference(rng);
    let start_pos = rng.gen_range(0..local_ref.len()) as usize;

    let end_pos = if start_pos + length > local_ref.len() {
        local_ref.len() - 1
    } else {
        start_pos + length
    };

    (
        local_ref[start_pos..end_pos].to_vec(),
        Origin::new(id.clone(), *strand, start_pos, end_pos, false, false),
    )
}

/// An iterator produce a length and u64 seed, until sum of length not reach target
struct LenIdSeed<'a> {
    target: u64,
    length_model: &'a model::Length,
    identity_model: &'a model::Identity,
    rng: &'a mut rand::rngs::StdRng,
}

impl<'a> LenIdSeed<'a> {
    /// Create a new LenIdSeed
    pub fn new(
        target: u64,
        length_model: &'a model::Length,
        identity_model: &'a model::Identity,
        rng: &'a mut rand::rngs::StdRng,
    ) -> Self {
        Self {
            target,
            length_model,
            identity_model,
            rng,
        }
    }
}

impl<'a> Iterator for LenIdSeed<'a> {
    type Item = (u64, f64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.target == 0 {
            return None;
        }

        let length = self.length_model.get_length(self.rng);

        if length > self.target {
            self.target = 0;
        } else {
            self.target -= length;
        }

        Some((
            length,
            self.identity_model.get_identity(self.rng),
            self.rng.next_u64(),
        ))
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn length_iterator() {
        let length = model::Length::new(100.0, 5.0).unwrap();
        let identity = model::Identity::new(90.0, 100.0, 5.0).unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let it = LenIdSeed::new(1_000, &length, &identity, &mut rng);

        let lengths: Vec<(u64, f64, u64)> = it.collect();

        assert_eq!(
            vec![
                (100, 0.8766417197970516, 633513173585076202),
                (99, 0.812764901371403, 59990589126097438),
                (109, 0.915656282506864, 9648369018563374588),
                (96, 0.8970375784344082, 11924907057293251871),
                (104, 0.8990886407136385, 6475006809388010320),
                (99, 0.9411813514218517, 12075250104723286995),
                (99, 0.8741520990047399, 10644393278569127094),
                (99, 0.9104462472912016, 11204491199501125651),
                (99, 0.9485885683427931, 9064963017249372811),
                (104, 0.8928093807556134, 14957698519116010673)
            ],
            lengths
        );
    }
}
