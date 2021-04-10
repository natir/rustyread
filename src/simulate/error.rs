//! Add error on reads

/* standard use */

/* crate use */
use rand::Rng;

/* local use */
use crate::model;

type Seq = Vec<u8>;

/// From identity an seq length compute number of edit
pub fn number_of_edit(target: f64, length: usize) -> f64 {
    (1.0 - target) * length as f64
}

/// Apply error on read
pub fn add_error(
    target: f64,
    seq: &[u8],
    error_model: &model::Error,
    glitch_model: &model::Glitch,
    rng: &mut rand::rngs::StdRng,
) -> (Seq, DiffPos) {
    let k = error_model.k();
    let mut changes = generate_change(seq, number_of_edit(target, seq.len()), error_model, rng);
    changes.extend(generate_glitches(seq.len(), glitch_model, rng));

    changes.sort_by_cached_key(|x| x.begin);

    let asm_changes = assemble_change(changes, k);

    (
        apply_changes(seq, &asm_changes),
        DiffPos::from_change(&asm_changes),
    )
}

/// Generate glitches
pub fn generate_glitches(
    length: usize,
    model: &model::Glitch,
    rng: &mut rand::rngs::StdRng,
) -> Vec<Change> {
    let mut glitches = Vec::new();

    let mut position = 0;

    while let Some((begin, end, seq)) = model.get_glitch(rng) {
        position += begin;
        if position > length || position + end > length {
            break;
        }
        glitches.push(Change::new(position, end - begin, &seq));
    }

    glitches
}

/// Generate change have to apply on read
pub fn generate_change(
    seq: &[u8],
    mut target: f64,
    model: &model::Error,
    rng: &mut rand::rngs::StdRng,
) -> Vec<Change> {
    let mut changes = Vec::new();
    let k = model.k();

    // generate change
    while target > 0.0 {
        let pos = rng.gen_range(0..(seq.len() - k));

        let (kmer, edit) = model.add_errors_to_kmer(&seq[pos..pos + k], rng);
        if edit == 0 {
            continue;
        }
        changes.push(Change::new(pos, k, &kmer));
        target -= edit as f64;
    }

    changes
}

#[derive(Debug, Clone, PartialEq)]
pub struct Change {
    pub begin: usize,
    pub end_raw: usize,
    pub end_err: usize,
    pub err_seq: Vec<u8>,
}

impl Change {
    pub fn new(begin: usize, k: usize, err_seq: &[u8]) -> Self {
        Self {
            begin,
            end_raw: begin + k,
            end_err: begin + err_seq.len(),
            err_seq: err_seq.to_owned(),
        }
    }
}

/// Assemble overlapping change
pub fn assemble_change(old: Vec<Change>, k: usize) -> Vec<Change> {
    let mut changes = Vec::new();

    log::trace!("RAW CHANGE {:?}", old);
    let mut prev: Option<Change> = None;
    for change in old.iter() {
        log::trace!("TOP change {:?} prev {:?}", change, prev);
        if let Some(ref mut p) = prev {
            if change.begin < p.end_err {
                p.end_raw = change.begin + k;

                let ovl = p.end_err - change.begin;
                if change.err_seq.len() > ovl {
                    p.err_seq.extend(&change.err_seq[ovl..]);
                }
                p.end_err = p.begin + p.err_seq.len();
                log::trace!("OVL ERR {:?}", prev);
            } else if change.begin < p.end_raw {
                p.end_raw = change.begin + k;

                p.err_seq.extend(&change.err_seq);
                p.end_err = p.begin + p.err_seq.len();
                log::trace!("OVL RAW {:?}", prev)
            } else {
                changes.push(p.clone());
                prev = Some(change.clone());
                log::trace!("NEW PREV {:?}", prev);
            }
        } else {
            prev = Some(change.clone());
            log::trace!("new prev {:?}", prev);
        }
    }

    if let Some(p) = prev {
        changes.push(p);
    }

    changes
}

#[derive(Debug)]
pub struct DiffPos {
    pub raw: Vec<usize>,
    pub err: Vec<usize>,
}

impl DiffPos {
    pub fn from_change(changes: &[Change]) -> Self {
        let mut raw = Vec::new();
        let mut err = Vec::new();
        let mut err_len = 0;

        log::trace!("CHANGES {:?}", changes);
        for change in changes {
            let match_dist = if !raw.is_empty() {
                log::trace!("begin {} raw[-1] {}", change.begin, raw[raw.len() - 1]);
                change.begin - raw[raw.len() - 1]
            } else {
                change.begin
            };

            raw.push(change.begin);
            raw.push(change.end_raw);

            err_len += match_dist;
            err.push(err_len);
            err_len += change.err_seq.len();
            err.push(err_len);
        }

        Self { raw, err }
    }
}

/// Apply changes on sequencs
pub fn apply_changes(seq: &[u8], changes: &[Change]) -> Seq {
    let mut err = Vec::with_capacity(seq.len());
    let mut pos_in_raw = 0;

    for change in changes {
        err.extend(&seq[pos_in_raw..change.begin]);
        err.extend(&change.err_seq);
        pos_in_raw = change.end_raw;
    }
    err.extend(&seq[pos_in_raw..]);
    err
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    fn init() {
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    #[test]
    fn nb_of_edit() {
        assert_eq!(5.000000000000004, number_of_edit(0.95, 100));
        assert_eq!(131.99999999999997, number_of_edit(0.89, 1200));
    }

    #[test]
    fn get_glitches() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = model::Glitch::new(1_000.0, 25.0, 25.0).unwrap();

        let glitches = generate_glitches(10_000, &model, &mut rng);

        assert_eq!(
            vec![
                Change {
                    begin: 381,
                    end_raw: 391,
                    end_err: 404,
                    err_seq: vec![
                        65, 84, 84, 65, 84, 65, 71, 84, 65, 67, 71, 71, 84, 65, 84, 65, 71, 84, 71,
                        71, 84, 84, 65
                    ]
                },
                Change {
                    begin: 2553,
                    end_raw: 2555,
                    end_err: 2573,
                    err_seq: vec![
                        71, 67, 67, 84, 65, 65, 71, 84, 71, 71, 67, 71, 67, 67, 67, 71, 84, 84, 71,
                        84
                    ]
                },
                Change {
                    begin: 2607,
                    end_raw: 2670,
                    end_err: 2634,
                    err_seq: vec![
                        65, 71, 71, 65, 65, 84, 67, 67, 65, 67, 84, 84, 65, 84, 65, 84, 65, 65, 67,
                        65, 67, 65, 71, 71, 84, 65, 84
                    ]
                },
                Change {
                    begin: 2675,
                    end_raw: 2733,
                    end_err: 2676,
                    err_seq: vec![67]
                },
                Change {
                    begin: 2842,
                    end_raw: 2897,
                    end_err: 2877,
                    err_seq: vec![
                        71, 71, 67, 65, 84, 71, 67, 67, 84, 65, 84, 65, 84, 84, 67, 84, 65, 84, 71,
                        65, 67, 65, 71, 67, 65, 71, 71, 65, 84, 84, 65, 84, 71, 71, 65
                    ]
                },
                Change {
                    begin: 3936,
                    end_raw: 3938,
                    end_err: 3940,
                    err_seq: vec![71, 67, 84, 67]
                },
                Change {
                    begin: 4324,
                    end_raw: 4379,
                    end_err: 4332,
                    err_seq: vec![65, 67, 71, 84, 84, 84, 71, 71]
                },
                Change {
                    begin: 5139,
                    end_raw: 5196,
                    end_err: 5160,
                    err_seq: vec![
                        67, 67, 67, 71, 84, 65, 71, 67, 65, 67, 71, 65, 67, 67, 71, 71, 67, 84, 65,
                        84, 71
                    ]
                },
                Change {
                    begin: 5164,
                    end_raw: 5168,
                    end_err: 5189,
                    err_seq: vec![
                        84, 84, 84, 84, 67, 84, 84, 71, 71, 65, 67, 65, 84, 65, 71, 84, 84, 84, 67,
                        71, 84, 67, 67, 65, 67
                    ]
                },
                Change {
                    begin: 5313,
                    end_raw: 5358,
                    end_err: 5370,
                    err_seq: vec![
                        84, 65, 67, 65, 65, 71, 71, 65, 67, 71, 67, 84, 84, 71, 71, 71, 65, 65, 84,
                        65, 71, 71, 71, 67, 65, 71, 67, 71, 71, 65, 71, 84, 84, 65, 84, 67, 71, 84,
                        71, 84, 65, 67, 67, 84, 67, 67, 84, 65, 71, 67, 84, 84, 84, 84, 65, 71, 84
                    ]
                },
                Change {
                    begin: 5912,
                    end_raw: 5955,
                    end_err: 5933,
                    err_seq: vec![
                        65, 67, 65, 71, 84, 71, 84, 65, 65, 67, 65, 84, 84, 71, 71, 71, 65, 67, 71,
                        67, 84
                    ]
                },
                Change {
                    begin: 6209,
                    end_raw: 6240,
                    end_err: 6215,
                    err_seq: vec![67, 71, 67, 67, 71, 71]
                },
                Change {
                    begin: 6774,
                    end_raw: 6777,
                    end_err: 6806,
                    err_seq: vec![
                        84, 84, 67, 67, 84, 84, 71, 65, 67, 84, 65, 84, 65, 67, 67, 71, 65, 84, 67,
                        71, 84, 71, 71, 65, 71, 84, 84, 67, 65, 84, 71, 67
                    ]
                },
                Change {
                    begin: 7774,
                    end_raw: 7775,
                    end_err: 7789,
                    err_seq: vec![67, 67, 84, 67, 65, 71, 67, 71, 84, 84, 67, 84, 67, 71, 71]
                }
            ],
            glitches
        );
    }

    #[test]
    fn get_changes() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let seq = crate::random_seq(100, &mut rng);
        let model = crate::model::Error::from_stream(
            std::fs::File::open("tests/data/error_model_seed42_100.csv").unwrap(),
            &mut rng,
        )
        .unwrap();

        let changes = generate_change(&seq, 10.0, &model, &mut rng);

        assert_eq!(
            vec![
                Change {
                    begin: 31,
                    end_raw: 38,
                    end_err: 39,
                    err_seq: vec![71, 71, 84, 65, 71, 67, 67, 84]
                },
                Change {
                    begin: 41,
                    end_raw: 48,
                    end_err: 47,
                    err_seq: vec![84, 71, 67, 71, 67, 67]
                },
                Change {
                    begin: 50,
                    end_raw: 57,
                    end_err: 56,
                    err_seq: vec![84, 84, 71, 84, 71, 65]
                },
                Change {
                    begin: 50,
                    end_raw: 57,
                    end_err: 58,
                    err_seq: vec![84, 84, 71, 84, 65, 71, 71, 65]
                },
                Change {
                    begin: 90,
                    end_raw: 97,
                    end_err: 97,
                    err_seq: vec![65, 67, 65, 71, 67, 65, 84]
                },
                Change {
                    begin: 71,
                    end_raw: 78,
                    end_err: 77,
                    err_seq: vec![84, 65, 65, 67, 67, 65]
                },
                Change {
                    begin: 21,
                    end_raw: 28,
                    end_err: 29,
                    err_seq: vec![84, 65, 71, 84, 84, 84, 65, 67]
                },
                Change {
                    begin: 83,
                    end_raw: 90,
                    end_err: 89,
                    err_seq: vec![65, 65, 84, 67, 71, 71]
                },
                Change {
                    begin: 66,
                    end_raw: 73,
                    end_err: 73,
                    err_seq: vec![84, 84, 71, 84, 65, 84, 65]
                }
            ],
            changes
        );
    }

    #[test]
    fn asm_change() {
        init();

        let mut changes = vec![
            Change::new(1, 7, &[65, 65, 84, 65, 65, 67, 67, 65]),
            Change::new(2, 7, &[65, 84, 65, 71, 67, 65, 67]),
            Change::new(3, 7, &[84, 65, 71, 67, 65, 65, 84]),
        ];

        let mut asm = assemble_change(changes, 7);

        assert_eq!(
            vec![Change {
                begin: 1,
                end_raw: 10,
                end_err: 10,
                err_seq: vec![65, 65, 84, 65, 65, 67, 67, 65, 84]
            }],
            asm
        );

        changes = vec![
            Change::new(5, 7, &[65, 71, 65, 65, 71, 71]),
            Change::new(5, 7, &[65, 71, 65, 71]),
            Change::new(6, 7, &[71, 65, 71]),
            Change::new(7, 7, &[65, 65, 71, 71, 71, 67]),
            Change::new(7, 7, &[65, 71, 71, 71, 71, 67]),
        ];

        asm = assemble_change(changes, 7);

        assert_eq!(
            vec![Change {
                begin: 5,
                end_raw: 14,
                end_err: 13,
                err_seq: vec![65, 71, 65, 65, 71, 71, 71, 67],
            }],
            asm
        );
    }

    #[test]
    fn asm_change_hard() {
        init();

        let changes = vec![
            Change::new(72, 7, &[67, 71, 67, 84, 71, 71]),
            Change::new(76, 7, &[84, 65, 71, 67]),
            Change::new(83, 7, &[65, 65, 67, 84, 84, 67, 71, 67]),
            Change::new(88, 7, &[84, 67, 67, 84, 84, 67, 67]),
            Change::new(88, 7, &[84, 67, 67, 84, 84, 67, 67]),
            Change::new(90, 7, &[67, 84, 84, 67, 71, 71]),
            Change::new(92, 7, &[84, 67, 65, 71, 71, 67]),
        ];

        let t_asm = vec![
            Change {
                begin: 72,
                end_raw: 83,
                end_err: 80,
                err_seq: vec![67, 71, 67, 84, 71, 71, 71, 67],
            },
            Change {
                begin: 83,
                end_raw: 99,
                end_err: 98,
                err_seq: vec![65, 65, 67, 84, 84, 67, 71, 67, 84, 84, 67, 67, 71, 71, 67],
            },
        ];

        assert_eq!(t_asm, assemble_change(changes, 7));
    }

    #[test]
    fn match_pos() {
        let changes = vec![
            Change::new(21, 7, &[84, 65, 71, 84, 84, 84, 65, 67]),
            Change::new(31, 7, &[71, 71, 84, 65, 71, 67, 67, 84]),
            Change::new(41, 7, &[84, 71, 67, 71, 67, 67]),
            Change::new(50, 7, &[84, 84, 71, 84, 65, 71, 71, 65]),
            Change::new(50, 7, &[84, 84, 71, 84, 71, 65]),
            Change::new(66, 7, &[84, 84, 71, 84, 65, 84, 65]),
            Change::new(71, 7, &[84, 65, 65, 67, 67, 65]),
            Change::new(83, 7, &[65, 65, 84, 67, 71, 71]),
            Change::new(90, 7, &[65, 67, 65, 71, 67, 65, 84]),
        ];

        let match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        let raws: Vec<(usize, usize)> = match_pos
            .raw
            .chunks_exact(2)
            .map(|x| (x[0], x[1]))
            .collect();
        let errs: Vec<(usize, usize)> = match_pos
            .err
            .chunks_exact(2)
            .map(|x| (x[0], x[1]))
            .collect();

        assert_eq!(
            vec![
                (21, 28),
                (31, 38),
                (41, 48),
                (50, 57),
                (66, 78),
                (83, 90),
                (90, 97)
            ],
            raws
        );
        assert_eq!(
            vec![
                (21, 29),
                (32, 40),
                (43, 49),
                (51, 59),
                (68, 79),
                (84, 90),
                (90, 97)
            ],
            errs
        );
    }

    #[test]
    fn match_pos_simple() {
        let mut changes = vec![
            Change::new(3, 7, &vec![84; 8]),
            Change::new(24, 7, &vec![71; 6]),
            Change::new(34, 7, &vec![71; 7]),
        ];

        let mut match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![3, 10, 24, 31, 34, 41], match_pos.raw);
        assert_eq!(vec![3, 11, 25, 31, 34, 41], match_pos.err);

        changes = vec![
            Change::new(56, 7, &vec![71; 6]),
            Change::new(60, 7, &vec![71; 7]),
        ];

        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![56, 67], match_pos.raw);
        assert_eq!(vec![56, 67], match_pos.err);

        changes = vec![
            Change::new(80, 7, &vec![71; 8]),
            Change::new(85, 7, &vec![71; 7]),
        ];

        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![80, 92], match_pos.raw);
        assert_eq!(vec![80, 92], match_pos.err);
    }

    #[test]
    fn many_change_at_same_pos() {
        init();

        let mut changes = vec![
            Change::new(1, 7, &[65, 65, 84, 65, 65, 67, 67, 65]),
            Change::new(2, 7, &[65, 84, 65, 71, 67, 65, 67]),
            Change::new(3, 7, &[84, 65, 71, 67, 65, 65, 84]),
        ];

        let mut match_pos = DiffPos::from_change(&assemble_change(changes, 7));
        assert_eq!(vec![1, 10], match_pos.raw);
        assert_eq!(vec![1, 10], match_pos.err);

        changes = vec![
            Change::new(5, 7, &[65, 71, 65, 65, 71, 71]),
            Change::new(5, 7, &[65, 71, 65, 71]),
            Change::new(6, 7, &[71, 65, 71]),
            Change::new(7, 7, &[65, 65, 71, 71, 71, 67]),
            Change::new(7, 7, &[65, 71, 71, 71, 71, 67]),
        ];
        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![5, 14], match_pos.raw);
        assert_eq!(vec![5, 13], match_pos.err);
    }

    #[test]
    fn apply_changes_() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let raw = crate::random_seq(100, &mut rng);

        let changes = vec![
            Change::new(14, 7, &[71, 71, 84, 65, 84, 71, 65, 71]),
            Change::new(21, 7, &[84, 65, 71, 84, 84, 84, 65, 67]),
            Change::new(30, 7, &[84, 71, 84, 65, 67, 67, 67]),
            Change::new(31, 7, &[71, 71, 84, 65, 71, 67, 67, 84]),
            Change::new(41, 7, &[84, 71, 67, 71, 67, 67]),
            Change::new(44, 7, &[67, 71, 67, 67, 71, 84]),
            Change::new(49, 7, &[71, 84, 84, 71, 84, 71, 65, 71]),
            Change::new(50, 7, &[84, 84, 71, 84, 65, 71, 71, 65]),
            Change::new(50, 7, &[84, 84, 71, 84, 71, 65]),
            Change::new(52, 7, &[71, 84, 65, 71, 65, 71]),
            Change::new(54, 7, &[65, 71, 65, 71, 65, 65]),
            Change::new(66, 7, &[84, 84, 71, 84, 65, 84, 65]),
            Change::new(71, 7, &[84, 65, 65, 67, 67, 65]),
            Change::new(79, 7, &[71, 84, 65, 84, 65, 71, 84]),
            Change::new(82, 7, &[84, 65, 84, 67, 67, 71]),
            Change::new(83, 7, &[65, 65, 84, 67, 71, 71]),
            Change::new(87, 7, &[67, 71, 71, 65, 67, 65, 71]),
            Change::new(90, 7, &[65, 67, 65, 71, 67, 65, 84]),
            Change::new(91, 7, &[67, 71, 71, 67, 65, 67, 84, 71]),
        ];

        let err = apply_changes(&raw, &assemble_change(changes, 7));

        assert_eq!(b"TTAGATTATAGTACGGTATGAGAGTTTACTATGTACCCCTAAGTGCGCCCGTTTGTGAGAAATCCACTTGTATAACCAGGTATAGTCGGGACAGCATTGCG".to_vec(), err);
    }
}
