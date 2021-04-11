//! Simulate reads

/* mod declaration */
mod description;
mod error;
mod quality;

/* standard use */
use std::io::Write;

/* crate use */
use anyhow::{Context, Result};
use rand::Rng;
use rand::RngCore;
use rand::SeedableRng;
use rayon::prelude::*;

/* local use */
use crate::cli;
use crate::model;
use crate::references::*;
use description::{Description, Origin};

/* constant definition */
const CHIMERA_START_ADAPTER_CHANCE: f64 = 0.25;
const CHIMERA_END_ADAPTER_CHANCE: f64 = 0.25;

#[cfg(not(tarpaulin_include))]
/// main simulate function
pub fn simulate(params: cli::simulate::Command) -> Result<()> {
    let mut main_rng = if let Some(seed) = params.seed {
        rand::rngs::StdRng::seed_from_u64(seed)
    } else {
        rand::rngs::StdRng::seed_from_u64(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .with_context(|| "Get seed for rng")?
                .as_secs(),
        )
    };

    log::info!("Start read reference");
    let references = References::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(
            std::fs::File::open(params.reference_path).with_context(|| "Read reference file")?,
        )))
        .with_context(|| "Read reference file niffler")?
        .0,
    )?;
    log::info!("End read reference");

    log::info!("Start init lenght model");
    let length = model::Length::new(params.length.0 as f64, params.length.1 as f64)
        .with_context(|| "Init length model")?;
    log::info!("End init lenght model");

    log::info!("Start init identity model");
    let identity = model::Identity::new(
        params.identity.0 as f64,
        params.identity.1 as f64,
        params.identity.2 as f64,
    )
    .with_context(|| "Init identity model")?;
    log::info!("End init length model");

    log::info!("Start init adapter model");
    let adapter = model::Adapter::new(
        params.start_adapter_seq.as_bytes().to_vec(),
        params.end_adapter_seq.as_bytes().to_vec(),
        params.start_adapter.0 as f64,
        params.start_adapter.1 as f64,
        params.end_adapter.0 as f64,
        params.end_adapter.1 as f64,
    )
    .with_context(|| "Init adapter model")?;
    log::info!("End init adapter model");

    log::info!("Start init glitches model");
    let glitches = model::Glitch::new(
        params.glitches.0 as f64,
        params.glitches.1 as f64,
        params.glitches.2 as f64,
    )
    .with_context(|| "Init glitches model")?;
    log::info!("End init glitches model");

    log::info!("Start read error model");
    let error_path = crate::cli::simulate::found_model(params.error_model, "error".to_string())
        .with_context(|| "Get path of error model")?;
    let error = model::Error::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(
            std::fs::File::open(error_path).with_context(|| "Open error model")?,
        )))
        .with_context(|| "Open error model")?
        .0,
        &mut main_rng,
    )
    .with_context(|| "Init error model")?;
    let k = error.k();
    log::info!("End read error model");

    log::info!("Start read quality score model");
    let qscore_path = crate::cli::simulate::found_model(params.qscore_model, "qscore".to_string())
        .with_context(|| "Get path of qscore model")?;
    let qscore = model::Quality::from_stream(
        niffler::get_reader(Box::new(std::io::BufReader::new(
            std::fs::File::open(qscore_path).with_context(|| "Open qscore model")?,
        )))
        .with_context(|| "Open qscore model")?
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

    let junk_rate = params.junk / 100.0;
    let random_rate = params.random / 100.0;
    let chimera_rate = params.chimera / 100.0;
    log::info!("Start generate sequences");
    let sequences: Vec<(Description, Seq, Quality)> =
        LenIdSeed::new(total_base, chimera_rate, &length, &identity, &mut main_rng)
            .par_bridge()
            .map(|(len, len2, id, seed)| {
                generate_read(
                    &references,
                    (len, len2),
                    id,
                    junk_rate,
                    random_rate,
                    &adapter,
                    &error,
                    &glitches,
                    &qscore,
                    rand::rngs::StdRng::seed_from_u64(seed),
                )
                .unwrap()
            })
            .collect();
    log::info!("End generate sequences");

    log::info!("Start write sequences");
    let mut output: std::io::BufWriter<Box<dyn std::io::Write>> =
        if let Some(output_path) = params.output_path {
            std::io::BufWriter::new(Box::new(
                std::fs::File::create(output_path).with_context(|| "Open output file")?,
            ))
        } else {
            std::io::BufWriter::new(Box::new(std::io::stdout()))
        };

    for (comment, seq, qual) in sequences {
        writeln!(
            output,
            "@{} {}\n{}\n+\n{}",
            uuid::Uuid::new_v3(
                &uuid::Uuid::NAMESPACE_X500,
                &main_rng.gen::<u128>().to_be_bytes()
            )
            .to_hyphenated(),
            comment,
            std::str::from_utf8(&seq[k..(seq.len() - k)])
                .with_context(|| "Write read in output file")?, // begin and end of fragment is just random base
            std::str::from_utf8(&qual[k..seq.len() - k])
                .with_context(|| "Write read in output file")?
        )
        .with_context(|| "Write read in output file")?;
    }
    log::info!("End write sequences");

    Ok(())
}

type Seq = Vec<u8>;
type Quality = Vec<u8>;

#[allow(clippy::too_many_arguments)]
#[cfg(not(tarpaulin_include))]
/// Function realy generate read
fn generate_read(
    references: &References,
    length: (usize, Option<usize>),
    identity: f64,
    junk_rate: f64,
    random_rate: f64,
    adapter_model: &model::Adapter,
    error_model: &model::Error,
    glitch_model: &model::Glitch,
    qscore_model: &model::Quality,
    mut rng: rand::rngs::StdRng,
) -> Result<(Description, Seq, Quality)> {
    let k = error_model.k();

    // Estimate size of final fragment all edit is consider as insertion -> it's overestimation
    let mut estimate_length = 2 * k + adapter_model.max_len() + length.0 + length.1.unwrap_or(0);
    estimate_length += error::number_of_edit(identity, estimate_length).round() as usize;

    // Generate fragment
    let mut raw_fragment = Vec::with_capacity(estimate_length);
    raw_fragment.extend(crate::random_seq(k, &mut rng));

    let start_adapter = adapter_model.get_start(&mut rng);
    raw_fragment.extend(&start_adapter);

    let (ref_fragment, origin) =
        get_fragment(length.0, junk_rate, random_rate, references, &mut rng);
    raw_fragment.extend(&ref_fragment);

    // Add chimeric part
    let chimera_origin = if let Some(length2) = length.1 {
        if rng.gen_bool(CHIMERA_END_ADAPTER_CHANCE) {
            raw_fragment.extend(adapter_model.get_end(&mut rng));
        }
        if rng.gen_bool(CHIMERA_START_ADAPTER_CHANCE) {
            raw_fragment.extend(adapter_model.get_start(&mut rng));
        }

        let (other_fragment, other_origin) =
            get_fragment(length2, junk_rate, random_rate, references, &mut rng);
        raw_fragment.extend(other_fragment);

        Some(other_origin)
    } else {
        None
    };

    let end_adapter = adapter_model.get_end(&mut rng);
    raw_fragment.extend(&end_adapter);

    raw_fragment.extend(crate::random_seq(k, &mut rng));

    // Add error in fragment and produce quality
    let (err_fragment, diffpos) =
        error::add_error(identity, &raw_fragment, error_model, glitch_model, &mut rng);

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
    let des = Description::new(
        origin,
        chimera_origin,
        raw_fragment.len() - (k * 2),
        real_id * 100.0,
    );

    Ok((des, err_fragment, quality))
}

fn get_fragment(
    length: usize,
    junk_rate: f64,
    random_rate: f64,
    references: &References,
    rng: &mut rand::rngs::StdRng,
) -> (Vec<u8>, Origin) {
    if rng.gen_bool(junk_rate) {
        get_junk(length, rng)
    } else if rng.gen_bool(random_rate) {
        get_random(length, rng)
    } else {
        get_ref_fragment(length, references, rng)
    }
}

fn get_junk(length: usize, rng: &mut rand::rngs::StdRng) -> (Vec<u8>, Origin) {
    let small_seq = crate::random_seq(rng.gen_range(1..=5), rng);

    (
        small_seq.repeat(length / small_seq.len()),
        Origin::junk(small_seq.len() * (length / small_seq.len())),
    )
}

fn get_random(length: usize, rng: &mut rand::rngs::StdRng) -> (Vec<u8>, Origin) {
    (crate::random_seq(length, rng), Origin::random(length))
}

fn get_ref_fragment(
    length: usize,
    references: &References,
    rng: &mut rand::rngs::StdRng,
) -> (Vec<u8>, Origin) {
    let (id, local_ref, strand, circular) = &references.get_reference(rng);

    let start_pos = rng.gen_range(0..local_ref.len()) as usize;
    let mut end_pos = start_pos + length;

    let mut seq = Vec::new();

    if length > local_ref.len() && !circular {
        seq.reserve(local_ref.len());
        seq.extend(&local_ref[..]);

        return (
            seq,
            Origin::new(id.clone(), *strand, 0, local_ref.len() - 1, false, false),
        );
    }

    if *circular {
        seq.reserve(length);
        seq.extend(&local_ref[start_pos..]);

        if length - start_pos > local_ref.len() {
            end_pos = start_pos - 1;
        } else {
            end_pos = length - start_pos;
        }
        seq.extend(&local_ref[..=end_pos]);
    } else {
        if end_pos >= local_ref.len() {
            end_pos = local_ref.len() - 1;
        }
        seq.reserve(end_pos - start_pos);
        seq.extend(&local_ref[start_pos..=end_pos]);
    }

    (
        seq,
        Origin::new(id.clone(), *strand, start_pos, end_pos, false, false),
    )
}

/// An iterator produce a length and u64 seed, until sum of length not reach target
struct LenIdSeed<'a> {
    target: u64,
    chimera_rate: f64,
    length_model: &'a model::Length,
    identity_model: &'a model::Identity,
    rng: &'a mut rand::rngs::StdRng,
}

impl<'a> LenIdSeed<'a> {
    /// Create a new LenIdSeed
    pub fn new(
        target: u64,
        chimera_rate: f64,
        length_model: &'a model::Length,
        identity_model: &'a model::Identity,
        rng: &'a mut rand::rngs::StdRng,
    ) -> Self {
        Self {
            target,
            chimera_rate,
            length_model,
            identity_model,
            rng,
        }
    }
}

impl<'a> Iterator for LenIdSeed<'a> {
    type Item = (usize, Option<usize>, f64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.target == 0 {
            return None;
        }

        let length = self.length_model.get_length(self.rng) as usize;

        let length2 = if self.rng.gen_bool(self.chimera_rate) {
            Some(self.length_model.get_length(self.rng) as usize)
        } else {
            None
        };

        let tt_length = (length + length2.unwrap_or(0)) as u64;

        if tt_length > self.target {
            self.target = 0;
        } else {
            self.target -= tt_length;
        }

        Some((
            length,
            length2,
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

        let it = LenIdSeed::new(1_000, 0.1, &length, &identity, &mut rng);

        let lengths: Vec<(usize, Option<usize>, f64, u64)> = it.collect();

        assert_eq!(
            vec![
                (100, None, 0.9135596717952422, 7654602743214997928),
                (101, None, 0.8990623479644879, 2598777197943013338),
                (100, Some(100), 0.8746476160545346, 9181438499313657906),
                (100, None, 0.9077129516303598, 3754557543903678404),
                (99, None, 0.9387255684776452, 7196201290701505999),
                (93, None, 0.8821726775189294, 520811042500337745),
                (93, None, 0.8928093807556134, 14957698519116010673),
                (102, None, 0.9152291848490597, 4532762459406505246),
                (102, None, 0.946352414331155, 9855488741258874048),
                (95, None, 0.8890601689755928, 1094680473675979714),
            ],
            lengths
        );
    }

    #[test]
    fn junk_seq() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let (seq, ori) = get_junk(100, &mut rng);
        assert_eq!(b"AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT".to_vec(), seq);
        assert_eq!(ori, Origin::junk(100));
    }

    #[test]
    fn random_seq() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let (seq, ori) = get_random(100, &mut rng);

        assert_eq!(b"TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCG".to_vec(), seq);
        assert_eq!(ori, Origin::random(100));
    }

    static FASTA: &'static [u8] = b">random_seq_0
TCCTAACGTG
>random_seq_1 depth=1.5
TCACGATTAC
>random_seq_2 circular=true
CCTATCCGAT
>random_seq_3
TGCAAGATCA
>random_seq_4
TAGCCGTGGT
>random_seq_5
CGCTTTGTGA
>random_seq_6 circular=true depth=1
CACATGGGCG
>random_seq_7 depth=2.0 circular=false
ATCTAATGCG
>random_seq_8 depth=0.5 circular=true
CGGAACTCAG
>random_seq_9
TCCCGCTGTC
";

    #[test]
    fn read_reference() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();

        let frags: Vec<(Vec<u8>, Origin)> = (0..10)
            .map(|_| get_ref_fragment(100, &refs, &mut rng))
            .collect();

        assert_eq!(
            vec![
                (
                    b"CGCTTTGTGA".to_vec(),
                    Origin {
                        ref_id: "random_seq_5".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"TGATCTTGCA".to_vec(),
                    Origin {
                        ref_id: "random_seq_3".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"TCACGATTAC".to_vec(),
                    Origin {
                        ref_id: "random_seq_1".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"TCACAAAGCG".to_vec(),
                    Origin {
                        ref_id: "random_seq_5".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![65, 84, 67, 84, 65, 65, 84, 71, 67, 71],
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"CGCTTTGTGA".to_vec(),
                    Origin {
                        ref_id: "random_seq_5".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![67, 71, 67, 65, 84, 84, 65, 71, 65, 84],
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"TGCAAGATCA".to_vec(),
                    Origin {
                        ref_id: "random_seq_3".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"GTAATCGTGA".to_vec(),
                    Origin {
                        ref_id: "random_seq_1".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        junk: false,
                        random: false
                    }
                ),
                (
                    b"CGCTGAGTTC".to_vec(),
                    Origin {
                        ref_id: "random_seq_8".to_string(),
                        strand: '-',
                        start: 8,
                        end: 7,
                        junk: false,
                        random: false
                    }
                )
            ],
            frags
        );
    }

    static FASTA_LARGE: &'static [u8] = b">random_seq_0
TCCTAACGTG
TCACGATTAC
CCTATCCGAT
TGCAAGATCA
TAGCCGTGGT
CGCTTTGTGA
CACATGGGCG
ATCTAATGCG
CGGAACTCAG
TCCCGCTGTC
";

    #[test]
    fn read_larger_reference() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA_LARGE)).unwrap();

        let frags: Vec<(Vec<u8>, Origin)> = (0..5)
            .map(|_| get_ref_fragment(18, &refs, &mut rng))
            .collect();

        assert_eq!(
            vec![
                (
                    vec![67],
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '+',
                        start: 99,
                        end: 99,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![
                        65, 71, 67, 67, 71, 84, 71, 71, 84, 67, 71, 67, 84, 84, 84, 71, 84, 71, 65
                    ],
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '+',
                        start: 41,
                        end: 59,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![71, 65, 65, 67, 84, 67, 65, 71, 84, 67, 67, 67, 71, 67, 84, 71, 84, 67],
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '+',
                        start: 82,
                        end: 99,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![
                        71, 84, 67, 65, 67, 65, 65, 65, 71, 67, 71, 65, 67, 67, 65, 67, 71, 71, 67
                    ],
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '-',
                        start: 39,
                        end: 57,
                        junk: false,
                        random: false
                    }
                ),
                (
                    vec![
                        65, 67, 65, 71, 67, 71, 71, 71, 65, 67, 84, 71, 65, 71, 84, 84, 67, 67, 71
                    ],
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '-',
                        start: 1,
                        end: 19,
                        junk: false,
                        random: false
                    }
                )
            ],
            frags
        );
    }
}
