//! Simulate reads

/* mod declaration */
mod description;

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
                std::str::from_utf8(&seq)?,
                std::str::from_utf8(&qual)?
            )?;
        }
    }
    log::info!("End write sequences");

    Ok(())
}

type Seq = Box<[u8]>;
type Quality = Box<[u8]>;

/// Function realy generate read
fn worker(
    target: u64,
    references: &References,
    length_model: &model::Length,
    identity_model: &model::Identity,
    error_model: &model::Error,
    qscore_model: &model::Quality,
    mut rng: rand::rngs::StdRng,
) -> Vec<(Description, Seq, Quality)> {
    let mut data = Vec::new();

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

        let raw_fragment = local_ref.seq[start_pos..end_pos]
            .to_vec()
            .into_boxed_slice();

        let ori = Origin::new(
            local_ref.id.clone(),
            strand,
            start_pos,
            end_pos,
            false,
            false,
        );
        let des = Description::new(ori, None, length, identity);

        data.push((des, raw_fragment, vec![65].into_boxed_slice()));

        generate += (end_pos - start_pos) as u64;
    }

    data
}
