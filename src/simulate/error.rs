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
    rng: &mut rand::rngs::StdRng,
) -> (Seq, DiffPos) {
    let k = error_model.k();
    let changes = generate_change(seq, number_of_edit(target, seq.len()), error_model, rng);

    let asm_changes = assemble_change(changes, k);

    (
        apply_changes(seq, &asm_changes),
        DiffPos::from_change(&asm_changes),
    )
}

/// Generate change have to apply on read
pub fn generate_change(
    seq: &[u8],
    mut target: f64,
    model: &model::Error,
    rng: &mut rand::rngs::StdRng,
) -> Vec<(usize, Vec<u8>)> {
    let mut changes: Vec<(usize, Vec<u8>)> = Vec::new();
    let k = model.k();

    // generate change
    while target > 0.0 {
        let pos = rng.gen_range(0..(seq.len() - k));

        let (kmer, edit) = model.add_errors_to_kmer(&seq[pos..pos + k], rng);
        if edit == 0 {
            continue;
        }
        changes.push((pos, kmer));
        target -= edit as f64;
    }

    changes.sort();

    changes
}

#[derive(Debug, Clone, PartialEq)]
pub struct AsmChange {
    pub begin: usize,
    pub end_raw: usize,
    pub end_err: usize,
    pub err_seq: Vec<u8>,
}

impl AsmChange {
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
pub fn assemble_change(old: Vec<(usize, Vec<u8>)>, k: usize) -> Vec<AsmChange> {
    let mut changes = Vec::new();

    log::trace!("RAW CHANGE {:?}", old);
    let mut prev: Option<AsmChange> = None;
    for change in old.iter().map(|x| AsmChange::new(x.0, k, &x.1)) {
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

                p.err_seq.extend(change.err_seq);
                p.end_err = p.begin + p.err_seq.len();
                log::trace!("OVL RAW {:?}", prev)
            } else {
                changes.push(p.clone());
                prev = Some(change);
                log::trace!("NEW PREV {:?}", prev);
            }
        } else {
            prev = Some(change);
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
    pub fn from_change(changes: &[AsmChange]) -> Self {
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
pub fn apply_changes(seq: &[u8], changes: &[AsmChange]) -> Seq {
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
    fn get_changes() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let seq = crate::random_seq(100, &mut rng);
        let model = crate::model::Error::from_stream(
            std::fs::File::open("tests/data/error_model_seed42_100.csv").unwrap(),
            &mut rng,
        )
        .unwrap();

        let changes = generate_change(&seq, 20.0, &model, &mut rng);

        assert_eq!(
            vec![
                (14, [71, 71, 84, 65, 84, 71, 65, 71].to_vec()),
                (21, [84, 65, 71, 84, 84, 84, 65, 67].to_vec()),
                (30, [84, 71, 84, 65, 67, 67, 67].to_vec()),
                (31, [71, 71, 84, 65, 71, 67, 67, 84].to_vec()),
                (41, [84, 71, 67, 71, 67, 67].to_vec()),
                (44, [67, 71, 67, 67, 71, 84].to_vec()),
                (49, [71, 84, 84, 71, 84, 71, 65, 71].to_vec()),
                (50, [84, 84, 71, 84, 65, 71, 71, 65].to_vec()),
                (50, [84, 84, 71, 84, 71, 65].to_vec()),
                (52, [71, 84, 65, 71, 65, 71].to_vec()),
                (54, [65, 71, 65, 71, 65, 65].to_vec()),
                (66, [84, 84, 71, 84, 65, 84, 65].to_vec()),
                (71, [84, 65, 65, 67, 67, 65].to_vec()),
                (79, [71, 84, 65, 84, 65, 71, 84].to_vec()),
                (82, [84, 65, 84, 67, 67, 71].to_vec()),
                (83, [65, 65, 84, 67, 71, 71].to_vec()),
                (87, [67, 71, 71, 65, 67, 65, 71].to_vec()),
                (90, [65, 67, 65, 71, 67, 65, 84].to_vec()),
                (91, [67, 71, 71, 67, 65, 67, 84, 71].to_vec())
            ],
            changes
        );
    }

    #[test]
    fn asm_change() {
        init();

        let mut changes = vec![
            (1, vec![65, 65, 84, 65, 65, 67, 67, 65]),
            (2, vec![65, 84, 65, 71, 67, 65, 67]),
            (3, vec![84, 65, 71, 67, 65, 65, 84]),
        ];

        let mut asm = assemble_change(changes, 7);

        assert_eq!(
            vec![AsmChange {
                begin: 1,
                end_raw: 10,
                end_err: 10,
                err_seq: vec![65, 65, 84, 65, 65, 67, 67, 65, 84]
            }],
            asm
        );

        changes = vec![
            (5, vec![65, 71, 65, 65, 71, 71]),
            (5, vec![65, 71, 65, 71]),
            (6, vec![71, 65, 71]),
            (7, vec![65, 65, 71, 71, 71, 67]),
            (7, vec![65, 71, 71, 71, 71, 67]),
        ];

        asm = assemble_change(changes, 7);

        assert_eq!(
            vec![AsmChange {
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
            (72, vec![67, 71, 67, 84, 71, 71]),
            (76, vec![84, 65, 71, 67]),
            (83, vec![65, 65, 67, 84, 84, 67, 71, 67]),
            (88, vec![84, 67, 67, 84, 84, 67, 67]),
            (88, vec![84, 67, 67, 84, 84, 67, 67]),
            (90, vec![67, 84, 84, 67, 71, 71]),
            (92, vec![84, 67, 65, 71, 71, 67]),
        ];

        let t_asm = vec![
            AsmChange {
                begin: 72,
                end_raw: 83,
                end_err: 80,
                err_seq: vec![67, 71, 67, 84, 71, 71, 71, 67],
            },
            AsmChange {
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
            (21, [84, 65, 71, 84, 84, 84, 65, 67].to_vec()),
            (31, [71, 71, 84, 65, 71, 67, 67, 84].to_vec()),
            (41, [84, 71, 67, 71, 67, 67].to_vec()),
            (50, [84, 84, 71, 84, 65, 71, 71, 65].to_vec()),
            (50, [84, 84, 71, 84, 71, 65].to_vec()),
            (66, [84, 84, 71, 84, 65, 84, 65].to_vec()),
            (71, [84, 65, 65, 67, 67, 65].to_vec()),
            (83, [65, 65, 84, 67, 71, 71].to_vec()),
            (90, [65, 67, 65, 71, 67, 65, 84].to_vec()),
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
        let mut changes = vec![(3, vec![84; 8]), (24, vec![71; 6]), (34, vec![71; 7])];

        let mut match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![3, 10, 24, 31, 34, 41], match_pos.raw);
        assert_eq!(vec![3, 11, 25, 31, 34, 41], match_pos.err);

        changes = vec![(56, vec![71; 6]), (60, vec![71; 7])];

        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![56, 67], match_pos.raw);
        assert_eq!(vec![56, 67], match_pos.err);

        changes = vec![(80, vec![71; 8]), (85, vec![71; 7])];

        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![80, 92], match_pos.raw);
        assert_eq!(vec![80, 92], match_pos.err);
    }

    #[test]
    fn many_change_at_same_pos() {
        init();

        let mut changes = vec![
            (1, vec![65, 65, 84, 65, 65, 67, 67, 65]),
            (2, vec![65, 84, 65, 71, 67, 65, 67]),
            (3, vec![84, 65, 71, 67, 65, 65, 84]),
        ];

        let mut match_pos = DiffPos::from_change(&assemble_change(changes, 7));
        assert_eq!(vec![1, 10], match_pos.raw);
        assert_eq!(vec![1, 10], match_pos.err);

        changes = vec![
            (5, vec![65, 71, 65, 65, 71, 71]),
            (5, vec![65, 71, 65, 71]),
            (6, vec![71, 65, 71]),
            (7, vec![65, 65, 71, 71, 71, 67]),
            (7, vec![65, 71, 71, 71, 71, 67]),
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
            (14, [71, 71, 84, 65, 84, 71, 65, 71].to_vec()),
            (21, [84, 65, 71, 84, 84, 84, 65, 67].to_vec()),
            (30, [84, 71, 84, 65, 67, 67, 67].to_vec()),
            (31, [71, 71, 84, 65, 71, 67, 67, 84].to_vec()),
            (41, [84, 71, 67, 71, 67, 67].to_vec()),
            (44, [67, 71, 67, 67, 71, 84].to_vec()),
            (49, [71, 84, 84, 71, 84, 71, 65, 71].to_vec()),
            (50, [84, 84, 71, 84, 65, 71, 71, 65].to_vec()),
            (50, [84, 84, 71, 84, 71, 65].to_vec()),
            (52, [71, 84, 65, 71, 65, 71].to_vec()),
            (54, [65, 71, 65, 71, 65, 65].to_vec()),
            (66, [84, 84, 71, 84, 65, 84, 65].to_vec()),
            (71, [84, 65, 65, 67, 67, 65].to_vec()),
            (79, [71, 84, 65, 84, 65, 71, 84].to_vec()),
            (82, [84, 65, 84, 67, 67, 71].to_vec()),
            (83, [65, 65, 84, 67, 71, 71].to_vec()),
            (87, [67, 71, 71, 65, 67, 65, 71].to_vec()),
            (90, [65, 67, 65, 71, 67, 65, 84].to_vec()),
            (91, [67, 71, 71, 67, 65, 67, 84, 71].to_vec()),
        ];

        let err = apply_changes(&raw, &assemble_change(changes, 7));

        assert_eq!(b"TTAGATTATAGTACGGTATGAGAGTTTACTATGTACCCCTAAGTGCGCCCGTTTGTGAGAAATCCACTTGTATAACCAGGTATAGTCGGGACAGCATTGCG".to_vec(), err);
    }
}
