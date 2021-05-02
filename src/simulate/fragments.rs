//! Generate fragments

/* standard use */

/* crate use */

/* local use */
use crate::model;
use crate::references::*;
use crate::simulate::description::{Description, Origin, ReadType};

/// An iterator produce fragment, a description and a seed
pub struct Fragments<'a, R>
where
    R: rand::Rng,
{
    target: u64,
    junk_rate: f64,
    random_rate: f64,
    chimera_rate: f64,
    references: &'a References,
    length_model: &'a model::Length,
    identity_model: &'a model::Identity,
    rng: &'a mut R,
}

impl<'a, R> Fragments<'a, R>
where
    R: rand::Rng,
{
    /// Create a new Fragments
    pub fn new(
        target: u64,
        rates: (f64, f64, f64),
        references: &'a References,
        length_model: &'a model::Length,
        identity_model: &'a model::Identity,
        rng: &'a mut R,
    ) -> Self
    where
        R: rand::Rng,
    {
        Self {
            target,
            junk_rate: rates.0 / 100.0,
            random_rate: rates.1 / 100.0,
            chimera_rate: rates.2 / 100.0,
            references,
            length_model,
            identity_model,
            rng,
        }
    }

    /// Get the read type
    pub fn get_read_type(&mut self) -> ReadType {
        if self.rng.gen_bool(self.junk_rate) {
            ReadType::Junk
        } else if self.rng.gen_bool(self.random_rate) {
            ReadType::Random
        } else {
            ReadType::Real
        }
    }

    /// Return true fragment must be a chimera
    pub fn is_chimera(&mut self) -> bool {
        self.rng.gen_bool(self.chimera_rate)
    }

    /// Produce a fragment
    pub fn generate_fragment(&mut self) -> (usize, usize, Origin) {
        let read_type = self.get_read_type();
        let length = self.length_model.get_length(self.rng) as usize;

        match read_type {
            ReadType::Real => {
                let (mut ref_index, mut strand) = self.references.choose_reference(self.rng);
                let mut reference = &self.references.sequences[ref_index];

                while read_type == ReadType::Real
                    && !fragment_is_possible(length, reference.seq.len(), reference.circular)
                {
                    let (r, s) = self.references.choose_reference(self.rng);
                    ref_index = r;
                    strand = s;
                    reference = &self.references.sequences[ref_index];
                }

                let try_begin = self.rng.gen_range(0..reference.seq.len()) as usize;
                let (begin, end, real_length) = if try_begin + length < reference.seq.len() {
                    (try_begin, try_begin + length, length)
                } else if reference.circular {
                    (
                        try_begin,
                        length - (reference.seq.len() - try_begin),
                        length,
                    )
                } else {
                    (0, reference.seq.len() - 1, reference.seq.len())
                };

                (
                    ref_index,
                    real_length,
                    Origin::reference(reference.id.clone(), strand, begin, end),
                )
            }
            ReadType::Junk => (0, length, Origin::junk(length)),
            ReadType::Random => (0, length, Origin::random(length)),
        }
    }
}

fn fragment_is_possible(frag_len: usize, ref_len: usize, circular: bool) -> bool {
    if frag_len >= ref_len {
        !circular
    } else {
        true
    }
}

impl<'a, R> Iterator for Fragments<'a, R>
where
    R: rand::Rng,
{
    // ((ref_index, chimeric_ref_index), Description, Seed)
    type Item = (usize, usize, Description, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.target == 0 {
            return None;
        }

        let (first_index, first_length, first_origin) = self.generate_fragment();
        let (second_index, second_length, second_origin) = if self.is_chimera() {
            let tmp = self.generate_fragment();
            (tmp.0, tmp.1, Some(tmp.2))
        } else {
            (0, 0, None)
        };

        let tt_length = first_length + second_length;

        if tt_length as u64 > self.target {
            self.target = 0;
        } else {
            self.target -= tt_length as u64;
        }

        Some((
            first_index,
            second_index,
            Description::new(
                first_origin,
                second_origin,
                tt_length,
                self.identity_model.get_identity(self.rng),
            ),
            self.rng.next_u64(),
        ))
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
>random_seq_10
TCCTAACGTGTCACGATTACCCTATCCGATTGCAAGATCATAGCCGTGGTCGCTTTGTGACACATGGGCGATCTAATGCGCGGAACTCAGTCCCGCTGTC
";

    #[test]
    fn read_type() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();
        let length = model::Length::new(8.0, 2.0).unwrap();
        let identity = model::Identity::new(85.0, 95.0, 5.0).unwrap();
        let mut fragments = Fragments::new(
            10_000,
            (20.0, 20.0, 20.0),
            &refs,
            &length,
            &identity,
            &mut rng,
        );

        let read_type: Vec<ReadType> = (0..10).map(|_| fragments.get_read_type()).collect();
        assert_eq!(
            vec![
                ReadType::Real,
                ReadType::Real,
                ReadType::Junk,
                ReadType::Real,
                ReadType::Random,
                ReadType::Junk,
                ReadType::Real,
                ReadType::Random,
                ReadType::Real,
                ReadType::Junk
            ],
            read_type
        );
    }

    #[test]
    fn chimera() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();
        let length = model::Length::new(8.0, 2.0).unwrap();
        let identity = model::Identity::new(85.0, 95.0, 5.0).unwrap();
        let mut fragments = Fragments::new(
            10_000,
            (20.0, 20.0, 20.0),
            &refs,
            &length,
            &identity,
            &mut rng,
        );

        let chimera: Vec<bool> = (0..10).map(|_| fragments.is_chimera()).collect();
        assert_eq!(
            vec![false, false, false, false, true, false, false, false, true, true],
            chimera
        );
    }

    #[test]
    fn fragments() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();
        let length = model::Length::new(8.0, 2.0).unwrap();
        let identity = model::Identity::new(85.0, 95.0, 5.0).unwrap();
        let mut fragments = Fragments::new(
            10_000,
            (20.0, 20.0, 20.0),
            &refs,
            &length,
            &identity,
            &mut rng,
        );

        let frags: Vec<(usize, usize, Origin)> =
            (0..20).map(|_| fragments.generate_fragment()).collect();
        assert_eq!(
            vec![
                (
                    7,
                    10,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    4,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 4,
                        read_type: ReadType::Random
                    }
                ),
                (
                    0,
                    10,
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    9,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Junk
                    }
                ),
                (
                    10,
                    4,
                    Origin {
                        ref_id: "random_seq_10".to_string(),
                        strand: '+',
                        start: 16,
                        end: 20,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    17,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 17,
                        read_type: ReadType::Random
                    }
                ),
                (
                    1,
                    7,
                    Origin {
                        ref_id: "random_seq_1".to_string(),
                        strand: '-',
                        start: 0,
                        end: 7,
                        read_type: ReadType::Real
                    }
                ),
                (
                    7,
                    10,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    7,
                    10,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    8,
                    8,
                    Origin {
                        ref_id: "random_seq_8".to_string(),
                        strand: '+',
                        start: 5,
                        end: 3,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    6,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 6,
                        read_type: ReadType::Junk
                    }
                ),
                (
                    10,
                    6,
                    Origin {
                        ref_id: "random_seq_10".to_string(),
                        strand: '-',
                        start: 48,
                        end: 54,
                        read_type: ReadType::Real
                    }
                ),
                (
                    9,
                    10,
                    Origin {
                        ref_id: "random_seq_9".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    10,
                    9,
                    Origin {
                        ref_id: "random_seq_10".to_string(),
                        strand: '-',
                        start: 25,
                        end: 34,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    7,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 7,
                        read_type: ReadType::Junk
                    }
                ),
                (
                    0,
                    8,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 8,
                        read_type: ReadType::Random
                    }
                ),
                (
                    5,
                    10,
                    Origin {
                        ref_id: "random_seq_5".to_string(),
                        strand: '-',
                        start: 0,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    7,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 7,
                        read_type: ReadType::Random
                    }
                ),
                (
                    3,
                    6,
                    Origin {
                        ref_id: "random_seq_3".to_string(),
                        strand: '-',
                        start: 0,
                        end: 6,
                        read_type: ReadType::Real
                    }
                ),
                (
                    0,
                    8,
                    Origin {
                        ref_id: "".to_string(),
                        strand: '*',
                        start: 0,
                        end: 8,
                        read_type: ReadType::Random
                    }
                )
            ],
            frags
        );
    }

    #[test]
    fn iter() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let refs = References::from_stream(std::io::Cursor::new(FASTA)).unwrap();
        let length = model::Length::new(8.0, 2.0).unwrap();
        let identity = model::Identity::new(85.0, 95.0, 5.0).unwrap();
        let fragments = Fragments::new(
            10_000,
            (20.0, 20.0, 20.0),
            &refs,
            &length,
            &identity,
            &mut rng,
        );

        let frags: Vec<(usize, usize, Description, u64)> = fragments.take(10).collect();

        println!("{:?}", frags);
        assert_eq!(
            vec![
                (
                    7,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_7".to_string(),
                            strand: '+',
                            start: 0,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 10,
                        identity: 0.9023903395427547
                    },
                    17195042692806716983
                ),
                (
                    0,
                    6,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_0".to_string(),
                            strand: '+',
                            start: 0,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: Some(Origin {
                            ref_id: "random_seq_6".to_string(),
                            strand: '+',
                            start: 4,
                            end: 2,
                            read_type: ReadType::Real
                        }),
                        length: 18,
                        identity: 0.785919024034962
                    },
                    7410303534117827570
                ),
                (
                    10,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_10".to_string(),
                            strand: '-',
                            start: 35,
                            end: 44,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 9,
                        identity: 0.8336097597069272
                    },
                    657338316926129147
                ),
                (
                    7,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_7".to_string(),
                            strand: '+',
                            start: 0,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 10,
                        identity: 0.7943651602000301
                    },
                    10605392195150115091
                ),
                (
                    10,
                    10,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_10".to_string(),
                            strand: '+',
                            start: 0,
                            end: 99,
                            read_type: ReadType::Real
                        },
                        chimera: Some(Origin {
                            ref_id: "random_seq_10".to_string(),
                            strand: '+',
                            start: 68,
                            end: 74,
                            read_type: ReadType::Real
                        }),
                        length: 106,
                        identity: 0.9166196996085733
                    },
                    11312190434313393638
                ),
                (
                    10,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_10".to_string(),
                            strand: '+',
                            start: 25,
                            end: 33,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 8,
                        identity: 0.8409338668084709
                    },
                    5274222100112014305
                ),
                (
                    7,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_7".to_string(),
                            strand: '+',
                            start: 0,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 10,
                        identity: 0.9103369460151146
                    },
                    10567391463651436578
                ),
                (
                    10,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_10".to_string(),
                            strand: '+',
                            start: 81,
                            end: 88,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 7,
                        identity: 0.8210852839903914
                    },
                    12595372283568864177
                ),
                (
                    0,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "".to_string(),
                            strand: '*',
                            start: 0,
                            end: 8,
                            read_type: ReadType::Random
                        },
                        chimera: None,
                        length: 8,
                        identity: 0.8383956529757561
                    },
                    14078074552533106200
                ),
                (
                    0,
                    6,
                    Description {
                        origin: Origin {
                            ref_id: "".to_string(),
                            strand: '*',
                            start: 0,
                            end: 4,
                            read_type: ReadType::Random
                        },
                        chimera: Some(Origin {
                            ref_id: "random_seq_6".to_string(),
                            strand: '-',
                            start: 2,
                            end: 0,
                            read_type: ReadType::Real
                        }),
                        length: 12,
                        identity: 0.8815059110082734
                    },
                    14485571221210959617
                )
            ],
            frags
        );
    }
}
