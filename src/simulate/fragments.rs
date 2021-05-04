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

                while !fragment_is_possible(length, reference.seq.len(), reference.circular) {
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
                    (
                        try_begin,
                        reference.seq.len() - 1,
                        reference.seq.len() - try_begin,
                    )
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
    use std::io::Seek;
    use std::io::Write;

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
                    3,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 7,
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
                    5,
                    Origin {
                        ref_id: "random_seq_0".to_string(),
                        strand: '+',
                        start: 5,
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
                    2,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 8,
                        end: 9,
                        read_type: ReadType::Real
                    }
                ),
                (
                    7,
                    8,
                    Origin {
                        ref_id: "random_seq_7".to_string(),
                        strand: '+',
                        start: 2,
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
                    4,
                    Origin {
                        ref_id: "random_seq_9".to_string(),
                        strand: '-',
                        start: 6,
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
                    9,
                    Origin {
                        ref_id: "random_seq_5".to_string(),
                        strand: '-',
                        start: 1,
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

        assert_eq!(
            vec![
                (
                    7,
                    0,
                    Description {
                        origin: Origin {
                            ref_id: "random_seq_7".to_string(),
                            strand: '+',
                            start: 7,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 3,
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
                            start: 5,
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
                        length: 13,
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
                            start: 8,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 2,
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
                            start: 96,
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
                        length: 10,
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
                            start: 7,
                            end: 9,
                            read_type: ReadType::Real
                        },
                        chimera: None,
                        length: 3,
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

    #[test]
    fn shape() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let length = model::Length::new(18200.0, 15500.0).unwrap();

        let data: Vec<f64> = (0..100_000)
            .map(|_| length.get_length(&mut rng) as f64)
            .collect();

        let sum = data.iter().sum::<f64>();
        let avg = sum / data.len() as f64;
        let std =
            (data.iter().map(|x| (x - avg).powf(2.0)).sum::<f64>() / data.len() as f64).sqrt();

        assert_eq!(sum, 1823347447.0);
        assert_eq!(avg, 18233.47447);
        assert_eq!(std, 15450.618044403509);

        let identity = model::Identity::new(87.0, 100.0, 9.0).unwrap();

        let data: Vec<f64> = (0..100_000)
            .map(|_| identity.get_identity(&mut rng))
            .collect();

        let sum = data.iter().sum::<f64>();
        let avg = sum / data.len() as f64;
        let std =
            (data.iter().map(|x| (x - avg).powf(2.0)).sum::<f64>() / data.len() as f64).sqrt();

        assert_eq!(sum, 86957.66076320085);
        assert_eq!(avg, 0.8695766076320085);
        assert_eq!(std, 0.09016183161726986);

        let mut file: std::io::Cursor<Vec<u8>> = std::io::Cursor::new(Vec::new());

        writeln!(
            file,
            ">5_000_000\n{}",
            String::from_utf8(crate::random_seq(5_000_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">2_500\n{}",
            String::from_utf8(crate::random_seq(2_500, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">500\n{}",
            String::from_utf8(crate::random_seq(500, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">140_000\n{}",
            String::from_utf8(crate::random_seq(140_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">80_000\n{}",
            String::from_utf8(crate::random_seq(80_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">75_000\n{}",
            String::from_utf8(crate::random_seq(75_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">65_000\n{}",
            String::from_utf8(crate::random_seq(65_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">8_000\n{}",
            String::from_utf8(crate::random_seq(8_000, &mut rng)).unwrap()
        )
        .unwrap();
        writeln!(
            file,
            ">7_000\n{}",
            String::from_utf8(crate::random_seq(7_000, &mut rng)).unwrap()
        )
        .unwrap();

        file.seek(std::io::SeekFrom::Start(0)).unwrap();

        let refs = References::from_stream_adjusted_weight(file, false, &length, &mut rng).unwrap();

        let fragments = Fragments::new(
            refs.sequences
                .iter()
                .map(|x| x.seq.len() as u64)
                .sum::<u64>()
                * 50,
            (1.0, 1.0, 1.0),
            &refs,
            &length,
            &identity,
            &mut rng,
        );

        let mut lengths = vec![Vec::new(); refs.sequences.len()];
        let mut identitys = Vec::new();

        let mut type_count = [0; 4];
        for (ori1, ori2, des, _) in fragments {
            match des.origin.read_type {
                ReadType::Real => type_count[0] += 1,
                ReadType::Junk => type_count[1] += 1,
                ReadType::Random => type_count[2] += 1,
            }

            lengths[ori1].push((des.origin.end - des.origin.start) as f64);

            if let Some(chimera) = des.chimera {
                type_count[3] += 1;
                lengths[ori2].push((chimera.end - chimera.start) as f64);
            }

            identitys.push(des.identity);
        }

        assert_eq!(type_count, [14710, 155, 158, 130]);

        let read_per_ref: Vec<usize> = lengths.iter().map(|x| x.len()).collect();

        //5_000_000, 2_500, 500, 140_000, 80_000, 75_000, 65_000, 8_000, 7_000
        assert_eq!(
            read_per_ref,
            vec![13949, 46, 51, 367, 229, 211, 182, 63, 55]
        );

        assert_eq!(
            lengths[0].iter().cloned().sum::<f64>() / lengths[0].len() as f64,
            18149.65775324396
        );
        assert_eq!(
            lengths[1].iter().cloned().sum::<f64>() / lengths[1].len() as f64,
            1135.4782608695652
        );
        assert_eq!(
            lengths[2].iter().cloned().sum::<f64>() / lengths[2].len() as f64,
            256.6862745098039
        );
        assert_eq!(
            lengths[3].iter().cloned().sum::<f64>() / lengths[3].len() as f64,
            16951.59673024523
        );
        assert_eq!(
            lengths[4].iter().cloned().sum::<f64>() / lengths[4].len() as f64,
            15845.423580786026
        );
        assert_eq!(
            lengths[5].iter().cloned().sum::<f64>() / lengths[5].len() as f64,
            13841.938388625593
        );
        assert_eq!(
            lengths[6].iter().cloned().sum::<f64>() / lengths[6].len() as f64,
            13571.34065934066
        );
        assert_eq!(
            lengths[7].iter().cloned().sum::<f64>() / lengths[7].len() as f64,
            3849.095238095238
        );
        assert_eq!(
            lengths[8].iter().cloned().sum::<f64>() / lengths[8].len() as f64,
            3376.0363636363636
        );

        let sum_len: f64 = lengths.iter().map(|x| x.iter().sum::<f64>()).sum::<f64>();
        let avg_len: f64 = sum_len / lengths.len() as f64;
        let std_len: f64 = (lengths
            .iter()
            .map(|x| x.iter().map(|x| (x - avg_len).powf(2.0)).sum::<f64>())
            .sum::<f64>()
            / lengths.len() as f64)
            .sqrt();

        assert_eq!(sum_len, 268903545.0);
        assert_eq!(avg_len, 29878171.666666668);
        assert_eq!(std_len, 1225248306.06916);

        let sum_id: f64 = identitys.iter().sum::<f64>();
        let avg_id: f64 = sum_id / identitys.len() as f64;
        let std_id: f64 = (identitys
            .iter()
            .map(|x| (x - avg_id).powf(2.0))
            .sum::<f64>()
            / identitys.len() as f64)
            .sqrt();

        assert_eq!(sum_id, 13077.741035790727);
        assert_eq!(avg_id, 0.870514613312303);
        assert_eq!(std_id, 0.09010848672092044);
    }
}
