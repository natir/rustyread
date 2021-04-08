//! Add error on reads

/* standard use */

/* crate use */
use rand::Rng;

/* local use */
use crate::model;

type Seq = Vec<u8>;

/// Apply error on read
pub fn add_error(
    target: f64,
    error_model: &model::Error,
    seq: &[u8],
    rng: &mut rand::rngs::StdRng,
) -> (Seq, DiffPos) {
    let k = error_model.k();
    let changes = generate_change(seq, (1.0 - target) * seq.len() as f64, error_model, rng);

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

/// Assemble overlapping change
pub fn assemble_change(old: Vec<(usize, Vec<u8>)>, k: usize) -> Vec<(usize, usize, Vec<u8>)> {
    let mut changes = Vec::new();

    log::debug!("RAW CHANGE {:?}", old);
    let mut prev: Option<(usize, usize, Vec<u8>)> = None;
    for change in old {
        log::debug!("TOP CHANGE {:?} {:?}", change, prev);
        if let Some(ref mut p) = prev {
            if change.0 > (p.0 + p.2.len()) {
                changes.push(p.clone());
                p.0 = change.0;
                p.1 = change.0 + 7;
                p.2 = change.1.clone();
                log::debug!("end of ovl {:?}", changes);
            } else {
                log::debug!("change {:?}, prev {:?}", change, p);
                if change.0 + change.1.len() <= p.1 {
                    continue;
                } else {
                    let ovl = p.0 + p.2.len() - change.0;
                    log::debug!("p {:?} change {:?} ovl {}", p, change, ovl);
                    if ovl < change.1.len() {
                        p.2.extend(&change.1[ovl..]);
                    }

                    if ovl < 7 {
                        p.1 += 7 - ovl;
                    } else {
                        p.1 += ovl - 7;
                    }
                }
            }
        } else {
            prev = Some((change.0, change.0 + k, change.1.clone()));
            log::debug!("new prev {:?}", prev);
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
    pub fn from_change(changes: &[(usize, usize, Vec<u8>)]) -> Self {
        let mut raw = Vec::new();
        let mut err = Vec::new();
        let mut err_len = 0;

        log::debug!("CHANGES {:?}", changes);
        for (begin, end, change) in changes {
            let match_dist = if !raw.is_empty() {
                log::debug!("begin {} raw[-1] {}", begin, raw[raw.len() - 1]);
                begin - raw[raw.len() - 1]
            } else {
                *begin
            };

            raw.push(*begin);
            raw.push(*end);

            err_len += match_dist;
            err.push(err_len);
            err_len += change.len();
            err.push(err_len);
        }

        Self { raw, err }
    }
}

/// Apply changes on sequencs
pub fn apply_changes(seq: &[u8], changes: &[(usize, usize, Vec<u8>)]) -> Seq {
    let mut err = Vec::with_capacity(seq.len());
    let mut pos_in_raw = 0;

    for change in changes {
        err.extend(&seq[pos_in_raw..change.0]);
        err.extend(&change.2);
        pos_in_raw = change.1;
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

        assert_eq!(vec![(1, 9, vec![65, 65, 84, 65, 65, 67, 67, 65, 84])], asm);

        changes = vec![
            (5, vec![65, 71, 65, 65, 71, 71]),
            (5, vec![65, 71, 65, 71]),
            (6, vec![71, 65, 71]),
            (7, vec![65, 65, 71, 71, 71, 67]),
            (7, vec![65, 71, 71, 71, 71, 67]),
        ];

        asm = assemble_change(changes, 7);

        assert_eq!(vec![(5, 15, vec![65, 71, 65, 65, 71, 71, 71, 67])], asm);
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
            (72, 84, vec![67, 71, 67, 84, 71, 71, 71, 67]),
            (
                83,
                99,
                vec![65, 65, 67, 84, 84, 67, 71, 67, 84, 84, 67, 67, 71, 71, 67],
            ),
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

        assert_eq!(vec![56, 68], match_pos.raw);
        assert_eq!(vec![56, 67], match_pos.err);

        changes = vec![(80, vec![71; 8]), (85, vec![71; 7])];

        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![80, 91], match_pos.raw);
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
        assert_eq!(vec![1, 9], match_pos.raw);
        assert_eq!(vec![1, 10], match_pos.err);

        changes = vec![
            (5, vec![65, 71, 65, 65, 71, 71]),
            (5, vec![65, 71, 65, 71]),
            (6, vec![71, 65, 71]),
            (7, vec![65, 65, 71, 71, 71, 67]),
            (7, vec![65, 71, 71, 71, 71, 67]),
        ];
        match_pos = DiffPos::from_change(&assemble_change(changes, 7));

        assert_eq!(vec![5, 15], match_pos.raw);
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
        assert_eq!(b"TTAGATTATAGTACGGTATGAGAGTTTACCTATGTACCCCTAAGTGCGCCCGTTTGTGAGGAACCACTTGTATAACCAGGTATAGTCGGGACAGCATTGG".to_vec(), err);
    }
}
