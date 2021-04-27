//! A collections of sequence to store reference sequence

/* standard use */
use std::str::FromStr;

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;

/* local use */
use crate::model;

type Text = Box<[u8]>;

/// Store a reference sequence
#[derive(Debug, PartialEq)]
pub struct Reference {
    pub id: String,
    pub seq: Text,
    pub revcomp: Text,
    pub circular: bool,
}

impl Reference {
    /// Build a new refrence
    pub fn new(id: String, seq: Text, circular: bool) -> Self {
        let revcomp = bio::alphabets::dna::revcomp(&seq[..]).into_boxed_slice();

        Self {
            id,
            seq,
            revcomp,
            circular,
        }
    }
}

/// A collections of sequence
pub struct References {
    pub sequences: Vec<Reference>,
    pub dist: rand::distributions::WeightedIndex<f64>,
}

impl References {
    /// Read a collection of sequence in fasta format from an input stream.
    ///
    /// If sequence contains 'circular=true' in his description sequence is consider as circular
    /// If sequence contains 'depth=(\[\\d.\]+)' number of reads from this sequence is multiply by float value in capturing group
    pub fn from_stream<R>(input: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        let (seqs, prob) = References::read_reference(input)?;

        Ok(Self {
            sequences: seqs,
            dist: rand::distributions::WeightedIndex::new(prob)?,
        })
    }

    /// Same as from_stream but small sequence have increase weighted to fix bias.
    pub fn from_stream_adjusted_weigth<R>(
        input: R,
        small_plasmid_bias: bool,
        length_model: &model::Length,
        rng: &mut rand::rngs::StdRng,
    ) -> Result<Self>
    where
        R: std::io::Read,
    {
        let (seqs, mut prob) = References::read_reference(input)?;

        prob = References::adjust_depth(&seqs, prob, small_plasmid_bias, length_model, rng)?;

        Ok(Self {
            sequences: seqs,
            dist: rand::distributions::WeightedIndex::new(prob)?,
        })
    }

    /// Randomly get a reference index and strand according to depth
    pub fn choose_reference<R>(&self, rng: &mut R) -> (usize, char)
    where
        R: rand::Rng,
    {
        match ['+', '-'][rng.gen_range(0..=1) as usize] {
            '+' => (self.dist.sample(rng), '+'),
            '-' => (self.dist.sample(rng), '-'),
            _ => unreachable!(),
        }
    }

    /// Adjust depth of reference to fix bias in small sequence representation
    fn adjust_depth(
        sequences: &[Reference],
        mut weight: Vec<f64>,
        small_plasmid_bias: bool,
        model: &model::Length,
        rng: &mut rand::rngs::StdRng,
    ) -> Result<Vec<f64>> {
        let lengths: Vec<u64> = (0..100_000).map(|_| model.get_length(rng)).collect();
        let total: u64 = lengths.iter().sum();

        for (i, reference) in sequences.iter().enumerate() {
            if !small_plasmid_bias && reference.circular {
                let passing: u64 = lengths
                    .iter()
                    .filter(|x| x <= &&(reference.seq.len() as u64))
                    .sum();
                if passing == 0 {
                    anyhow::bail!(crate::error::Cli::SmallPlasmidBias);
                }

                weight[i] *= total as f64 / passing as f64;
            }

            if !reference.circular {
                let passing: u64 = lengths
                    .iter()
                    .map(|x| u64::min(*x, reference.seq.len() as u64))
                    .sum();

                weight[i] *= total as f64 / passing as f64;
            }
        }

        Ok(weight)
    }

    /// Read reference from stream
    fn read_reference<R>(input: R) -> Result<(Vec<Reference>, Vec<f64>)>
    where
        R: std::io::Read,
    {
        let mut me_seq = Vec::new();
        let mut me_pro = Vec::new();
        let mut records = bio::io::fasta::Reader::new(input).records();

        let weight_re = regex::Regex::new(r"depth=([\d.]+)").unwrap(); // we ignore result this regex is static

        while let Some(Ok(record)) = records.next() {
            let weight = if let Some(d) = record.desc() {
                if let Some(c) = weight_re.captures(d) {
                    f64::from_str(&c[1])?
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let circular = if let Some(d) = record.desc() {
                d.contains("circular=true")
            } else {
                false
            };

            me_seq.push(Reference::new(
                record.id().into(),
                record.seq().into(),
                circular,
            ));
            me_pro.push(weight * record.seq().len() as f64);
        }

        Ok((me_seq, me_pro))
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

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
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();

        assert_eq!(
            refs.sequences,
            vec![
                Reference {
                    id: "random_seq_0".to_string(),
                    seq: Box::new([84, 67, 67, 84, 65, 65, 67, 71, 84, 71]),
                    revcomp: Box::new([67, 65, 67, 71, 84, 84, 65, 71, 71, 65]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_1".to_string(),
                    seq: Box::new([84, 67, 65, 67, 71, 65, 84, 84, 65, 67]),
                    revcomp: Box::new([71, 84, 65, 65, 84, 67, 71, 84, 71, 65]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_2".to_string(),
                    seq: Box::new([67, 67, 84, 65, 84, 67, 67, 71, 65, 84]),
                    revcomp: Box::new([65, 84, 67, 71, 71, 65, 84, 65, 71, 71]),
                    circular: true,
                },
                Reference {
                    id: "random_seq_3".to_string(),
                    seq: Box::new([84, 71, 67, 65, 65, 71, 65, 84, 67, 65]),
                    revcomp: Box::new([84, 71, 65, 84, 67, 84, 84, 71, 67, 65]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_4".to_string(),
                    seq: Box::new([84, 65, 71, 67, 67, 71, 84, 71, 71, 84]),
                    revcomp: Box::new([65, 67, 67, 65, 67, 71, 71, 67, 84, 65]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_5".to_string(),
                    seq: Box::new([67, 71, 67, 84, 84, 84, 71, 84, 71, 65]),
                    revcomp: Box::new([84, 67, 65, 67, 65, 65, 65, 71, 67, 71]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_6".to_string(),
                    seq: Box::new([67, 65, 67, 65, 84, 71, 71, 71, 67, 71]),
                    revcomp: Box::new([67, 71, 67, 67, 67, 65, 84, 71, 84, 71]),
                    circular: true,
                },
                Reference {
                    id: "random_seq_7".to_string(),
                    seq: Box::new([65, 84, 67, 84, 65, 65, 84, 71, 67, 71]),
                    revcomp: Box::new([67, 71, 67, 65, 84, 84, 65, 71, 65, 84]),
                    circular: false,
                },
                Reference {
                    id: "random_seq_8".to_string(),
                    seq: Box::new([67, 71, 71, 65, 65, 67, 84, 67, 65, 71]),
                    revcomp: Box::new([67, 84, 71, 65, 71, 84, 84, 67, 67, 71]),
                    circular: true,
                },
                Reference {
                    id: "random_seq_9".to_string(),
                    seq: Box::new([84, 67, 67, 67, 71, 67, 84, 71, 84, 67]),
                    revcomp: Box::new([71, 65, 67, 65, 71, 67, 71, 71, 71, 65]),
                    circular: false,
                }
            ]
        );
    }

    #[test]
    fn get_reference() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();

        let seqs: Vec<(usize, char)> = (0..10).map(|_| refs.choose_reference(&mut rng)).collect();

        assert_eq!(
            vec![
                (2, '+'),
                (6, '-'),
                (4, '+'),
                (7, '-'),
                (4, '-'),
                (9, '+'),
                (8, '-'),
                (7, '-'),
                (1, '-'),
                (1, '-')
            ],
            seqs
        );
    }

    #[test]
    fn read_reference_adjust() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let length = model::Length::new(10.0, 8.0).unwrap();
        let refs = References::from_stream_adjusted_weigth(
            std::io::Cursor::new(FASTA),
            false,
            &length,
            &mut rng,
        )
        .unwrap();

        let seqs: Vec<(usize, char)> = (0..10).map(|_| refs.choose_reference(&mut rng)).collect();

        assert_eq!(
            vec![
                (6, '+'),
                (9, '+'),
                (2, '+'),
                (3, '-'),
                (1, '+'),
                (4, '-'),
                (6, '+'),
                (6, '+'),
                (1, '+'),
                (7, '+')
            ],
            seqs
        );
    }
}
