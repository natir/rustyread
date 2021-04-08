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
        LenSeed::new(total_base, &mut main_rng, &length)
            .par_bridge()
            .map(|(length, seed)| {
                generate_read(
                    &references,
                    length as usize,
                    &identity,
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
    identity_model: &model::Identity,
    error_model: &model::Error,
    qscore_model: &model::Quality,
    mut rng: rand::rngs::StdRng,
) -> Result<(Description, Seq, Quality)> {
    let k = error_model.k();

    let (id, local_ref, strand) = &references.get_reference(&mut rng);
    let start_pos = rng.gen_range(0..local_ref.len()) as usize;
    let identity = identity_model.get_identity(&mut rng);

    let end_pos = if start_pos + length > local_ref.len() {
        local_ref.len() - 1
    } else {
        start_pos + length
    };

    let mut raw_fragment = Vec::with_capacity(end_pos - start_pos + error_model.k() * 2);
    raw_fragment.extend(crate::random_seq(k, &mut rng));
    raw_fragment.extend(&local_ref[start_pos..end_pos]);
    raw_fragment.extend(crate::random_seq(k, &mut rng));

    let (err_fragment, diffpos) = error::add_error(identity, error_model, &raw_fragment, &mut rng);

    let (real_id, mut quality) = quality::generate_quality(
        &raw_fragment,
        &err_fragment,
        qscore_model,
        diffpos,
        &mut rng,
    )?;

    let ori = Origin::new(id.clone(), *strand, start_pos, end_pos, false, false);
    let des = Description::new(ori, None, length, real_id * 100.0);

    // TODO remove that
    quality.resize(err_fragment.len(), b'!');

    Ok((des, err_fragment, quality))
}

/// An iterator produce a length and u64 seed, until sum of length not reach target
struct LenSeed<'a> {
    target: u64,
    rng: &'a mut rand::rngs::StdRng,
    length_model: &'a model::Length,
}

impl<'a> LenSeed<'a> {
    /// Create a new LenSeed
    pub fn new(
        target: u64,
        rng: &'a mut rand::rngs::StdRng,
        length_model: &'a model::Length,
    ) -> Self {
        Self {
            target,
            rng,
            length_model,
        }
    }
}

impl<'a> Iterator for LenSeed<'a> {
    type Item = (u64, u64);

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

        Some((length, self.rng.next_u64()))
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn length_iterator() {
        let model = model::Length::new(100.0, 5.0).unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let it = LenSeed::new(1_000, &mut rng, &model);

        let lengths: Vec<(u64, u64)> = it.collect();

        assert_eq!(
            vec![
                (100, 11740708795755607249),
                (99, 7654602743214997928),
                (101, 2421668070752551774),
                (88, 9336807752121363465),
                (98, 9648369018563374588),
                (96, 9578448464351515635),
                (96, 15602080788219557311),
                (100, 1079037894117179173),
                (97, 14653814560646992085),
                (97, 12075250104723286995),
                (99, 17007309290515904425)
            ],
            lengths
        );
    }
}
