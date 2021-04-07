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
) -> Seq {
    let changes = generate_change(seq, (1.0 - target) * seq.len() as f64, error_model, rng);

    apply_changes(seq, changes, error_model.k())
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

/// Apply changes on sequencs
pub fn apply_changes(seq: &[u8], changes: Vec<(usize, Vec<u8>)>, k: usize) -> Seq {
    let mut err = Vec::with_capacity(seq.len());
    let mut pos_in_raw = 0;

    for change in changes {
        if change.0 < pos_in_raw {
            let remove = pos_in_raw - change.0;
            err.drain((err.len() - remove)..);
            pos_in_raw -= remove;
        }

        err.extend(&seq[pos_in_raw..change.0]);
        err.extend(change.1);

        pos_in_raw = change.0 + k;
    }

    err.extend(&seq[pos_in_raw..]);

    err
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

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

        let err = apply_changes(&raw, changes, 7);

        assert_eq!(b"TTAGATTATAGTACGGTATGAGTAGTTTACTATGGTAGCCTAAGTGCGCCGTTTGAGAGAATCCACTTGTATAACCAGGTAAATCGGACGGCACTGCG".to_vec(), err);
    }
}
