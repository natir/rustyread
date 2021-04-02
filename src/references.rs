//! A collections of sequence to store reference sequence

/* standard use */
use std::str::FromStr;

/* crate use */
use anyhow::Result;

/* local use */

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
pub type References = (Vec<Reference>, Vec<f64>);

pub trait AbsReferences {
    fn from_stream<R>(input: R) -> Result<Self>
    where
        R: std::io::BufRead,
        Self: Sized;
}

impl<'a> AbsReferences for References {
    /// Read a collection of sequence in fasta format from an input stream.
    ///
    /// If sequence contains 'circular=true' in his description sequence is consider as circular
    /// If sequence contains 'depth=(\[\\d.\]+)' number of reads from this sequence is multiply by float value in capturing group
    fn from_stream<R>(input: R) -> Result<Self>
    where
        R: std::io::BufRead,
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
            me_pro.push(weight);
        }

        Ok((me_seq, me_pro))
    }
}

#[cfg(test)]
mod t {
    use super::*;

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
>random_seq_7 depth=50.0 circular=false
ATCTAATGCG
>random_seq_8 depth=0.5 circular=true
CGGAACTCAG
>random_seq_9
TCCCGCTGTC
";

    #[test]
    fn read_reference() {
        let (me, prob) = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();

        assert_eq!(
            me,
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
        assert_eq!(
            prob,
            vec![1.0, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 50.0, 0.5, 1.0]
        );
    }
}
