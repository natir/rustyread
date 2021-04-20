//! Add error on reads

/* standard use */

/* crate use */
use rand::Rng;

/* local use */
use crate::model;

pub type Seq = Vec<u8>;
pub type Cigar = Vec<u8>;

#[derive(Debug, Clone, PartialEq)]
/// A struct to represent a change in a sequence
pub struct Change {
    begin: usize,
    end: usize,
    seq: Seq,
    cigar: Cigar,
    edit_distance: u64,
}

impl Change {
    /// Build a new Change from a begin, end, a change sequence and original sequence
    pub fn from_seq(begin: usize, end: usize, seq: Seq, original: &[u8]) -> Self {
        let mut new = Self {
            begin,
            end,
            seq,
            cigar: b"".to_vec(),
            edit_distance: 0,
        };

        new.update_align(original);

        new
    }

    /// Try to merge two Change.
    /// This function assume other.begin is greater than self.begin
    ///
    /// Return an None if self and other Change can't be merge
    /// Return Some(i64) change in edit distance
    pub fn merge(&mut self, other: &Change, original: &[u8]) -> Option<i64> {
        if self.contain(other) {
            return None;
        }

        if other.begin() < self.end_err() {
            let ovl = self.end_err() - other.begin();
            if other.seq().len() > ovl {
                self.seq.extend(&other.seq()[ovl..]);
            }

            self.end = self.end_raw() + other.end_raw() - other.begin() - ovl;
        } else if other.begin() < self.end_raw() {
            self.seq.extend(&other.seq()[..]);

            self.end = other.begin() + other.seq().len();
        } else {
            return None;
        }

        let prev_ed = self.edit_distance;
        self.update_align(original);

        Some(self.edit_distance as i64 - prev_ed as i64)
    }

    /// Get begin of Change
    pub fn begin(&self) -> usize {
        self.begin
    }

    /// Get end of change in raw sequence
    pub fn end_raw(&self) -> usize {
        self.end
    }

    /// Get end of change in erroneous sequence
    pub fn end_err(&self) -> usize {
        self.begin + self.seq.len()
    }

    /// Get erroneous sequence
    pub fn seq(&self) -> &Seq {
        &self.seq
    }

    /// Get cigar string
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Get edit distance
    pub fn edit(&self) -> u64 {
        self.edit_distance
    }

    /// Generate key associate to this change
    pub fn key(&self) -> ChangeKey {
        ChangeKey {
            begin: self.begin(),
            end_raw: self.end_raw(),
            end_err: self.end_err(),
        }
    }

    #[allow(clippy::suspicious_operation_groupings)]
    /// Return true if other is contain in self
    /// This function assume other.begin is greater than self.begin
    pub fn contain(&self, other: &Change) -> bool {
        other.end_raw() < self.end_err()
            || other.end_raw() < self.end_raw()
            || other.end_err() < self.end_err()
            || other.end_err() < self.end_err()
    }

    /// Update alignment info of error
    fn update_align(&mut self, original: &[u8]) {
        let (edit_distance, cigar) = if self.end_raw() > original.len() {
            crate::alignment::align(&self.seq, &original[self.begin()..])
        } else {
            crate::alignment::align(&self.seq, &original[self.begin()..self.end_raw()])
        };

        self.cigar = cigar.into_vec();
        self.edit_distance = edit_distance;
    }
}

#[derive(Debug, Clone)]
/// A struct to store Key parameter of change
pub struct ChangeKey {
    begin: usize,
    end_raw: usize,
    end_err: usize,
}

impl Eq for ChangeKey {}

impl std::cmp::PartialOrd for ChangeKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.eq(other) {
            Some(std::cmp::Ordering::Equal)
        } else {
            self.begin.partial_cmp(&other.begin)
        }
    }
}

impl std::cmp::Ord for ChangeKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.eq(other) {
            std::cmp::Ordering::Equal
        } else {
            self.begin.cmp(&other.begin)
        }
    }
}

/// WARNING this implementation don't do what you think
/// If two ChangeKey share an overlap in raw or in erroneous PartialEq return True
impl std::cmp::PartialEq for ChangeKey {
    fn eq(&self, other: &Self) -> bool {
        let (first, second) = if self.begin < other.begin {
            (&self, &other)
        } else {
            (&other, &self)
        };

        second.begin < first.end_err || second.begin < first.end_raw
    }
}

pub type Changes = std::collections::BTreeMap<ChangeKey, Change>;

/// From identity an seq length compute number of edit
pub fn number_of_edit(target: f64, length: usize) -> f64 {
    (1.0 - target) * length as f64
}

/// Apply error on read
pub fn add_error(
    identity: f64,
    seq: &[u8],
    error_model: &model::Error,
    glitch_model: &model::Glitch,
    rng: &mut rand::rngs::StdRng,
) -> (Seq, Cigar, f64) {
    let k = error_model.k();
    let target = number_of_edit(identity, seq.len());

    let mut raw_changes = get_glitches(seq, glitch_model, rng);
    get_error(k, target, seq, &mut raw_changes, error_model, rng);
    let (changes, real_edit) = asm_change(raw_changes, seq);

    let mut pos_in_raw = 0;
    let mut err: Vec<u8> =
        Vec::with_capacity(seq.len() + number_of_edit(target, seq.len()) as usize);
    let mut cig: Vec<u8> =
        Vec::with_capacity(seq.len() + number_of_edit(target, seq.len()) as usize);
    for change in changes.values() {
        if change.begin() < pos_in_raw {
            continue;
        }

        err.extend(&seq[pos_in_raw..change.begin()]);
        cig.extend(std::iter::repeat(b'=').take(change.begin() - pos_in_raw));

        err.extend(change.seq());
        cig.extend(change.cigar());

        pos_in_raw = change.end_raw();
    }

    err.extend(&seq[pos_in_raw..]);
    cig.extend(std::iter::repeat(b'=').take(seq.len() - pos_in_raw));

    (err, cig, (1.0 - (real_edit / seq.len() as f64)))
}

/// Merge changes they share overlap
pub fn asm_change(mut inputs: Changes, seq: &[u8]) -> (Changes, f64) {
    let mut real_edit = 0.0;
    let mut changes = Changes::new();
    let mut prev_wrapp: Option<Change> = None;

    for change in inputs.values_mut() {
        if let Some(ref mut prev) = prev_wrapp {
            if prev.contain(change) {
                continue;
            }

            let (first, second) = if prev.begin() < change.begin() {
                (prev, change)
            } else {
                (change, prev)
            };

            if first.merge(second, &seq).is_none() {
                real_edit += first.edit() as f64;
                changes.insert(first.key(), first.clone());
                prev_wrapp = Some(second.clone());
            }
        } else {
            prev_wrapp = Some(change.clone());
        }
    }
    if let Some(prev) = prev_wrapp {
        real_edit += prev.edit() as f64;
        changes.insert(prev.key(), prev);
    }

    (changes, real_edit)
}

/// Create Change correspond to glitches
pub fn get_glitches(raw: &[u8], model: &model::Glitch, rng: &mut rand::rngs::StdRng) -> Changes {
    let mut changes = Changes::new();
    let mut position = 0;

    while let Some((begin, end, seq)) = model.get_glitch(rng) {
        position += begin;
        if position > raw.len() || position + end > raw.len() {
            break;
        }

        let change = Change::from_seq(position + begin, position + end, seq, raw);
        changes.insert(change.key(), change);
    }

    changes
}

/// Add Change correspond to error in changes
pub fn get_error(
    k: usize,
    target: f64,
    seq: &[u8],
    changes: &mut Changes,
    error_model: &model::Error,
    rng: &mut rand::rngs::StdRng,
) {
    let mut sum_of_edit = 0.0;

    while sum_of_edit < target {
        let pos = rng.gen_range(0..(seq.len() - k));

        let (kmer, edit) = error_model.add_errors_to_kmer(&seq[pos..pos + k], rng);
        if edit == 0 {
            continue;
        }

        let mut new_change = Change::from_seq(pos, pos + k, kmer, seq);
        let key = new_change.key();

        if let Some(mut old_change) = changes.remove(&key) {
            let (first, second) = if old_change.begin() < new_change.begin() {
                (&mut old_change, &mut new_change)
            } else {
                (&mut new_change, &mut old_change)
            };

            if let Some(edit_change) = first.merge(second, seq) {
                sum_of_edit += edit_change as f64;
            }
            changes.insert(first.key(), first.clone());
        } else {
            sum_of_edit += new_change.edit() as f64;
            changes.insert(new_change.key(), new_change);
        }
    }
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
    fn build_change() {
        init();

        let raw = b"GGCGGTGTCCTTGCTAACAT";
        let mut change = Change::from_seq(3, 10, b"GGTAACC".to_vec(), raw);
        assert_eq!(change.begin(), 3);
        assert_eq!(change.end_raw(), 10);
        assert_eq!(change.end_err(), 10);
        assert_eq!(change.seq(), &b"GGTAACC".to_vec());
        assert_eq!(change.cigar(), &b"===XX==".to_vec());
        assert_eq!(change.edit(), 2);

        change = Change::from_seq(6, 13, b"GTTTG".to_vec(), raw);
        assert_eq!(change.begin(), 6);
        assert_eq!(change.end_raw(), 13);
        assert_eq!(change.end_err(), 11);
        assert_eq!(change.seq(), &b"GTTTG".to_vec());
        assert_eq!(change.cigar(), &b"==DD===".to_vec());
        assert_eq!(change.edit(), 2);

        change = Change::from_seq(6, 13, b"GTCAACTTG".to_vec(), raw);
        assert_eq!(change.begin(), 6);
        assert_eq!(change.end_raw(), 13);
        assert_eq!(change.end_err(), 15);
        assert_eq!(change.seq(), &b"GTCAACTTG".to_vec());
        assert_eq!(change.cigar(), &b"===II====".to_vec());
        assert_eq!(change.edit(), 2);
    }

    #[test]
    fn merge_change() {
        init();

        let raw = b"TCAGGAAGATCCACA";
        // TCAGGAAGATCCACA
        //    GGCCGAT
        //       AGA-CCA
        //    GGCCGAT-CA
        let mut first = Change::from_seq(3, 10, b"GGCCGAT".to_vec(), raw);
        let second = Change::from_seq(6, 13, b"AGACCA".to_vec(), raw);

        first.merge(&second, raw);

        assert_eq!(first.begin(), 3);
        assert_eq!(first.end_raw(), 13);
        assert_eq!(first.end_err(), 12);
        assert_eq!(first.seq(), &b"GGCCGATCA".to_vec());
        assert_eq!(first.cigar(), &b"==XX===D==".to_vec());
        assert_eq!(first.edit(), 3);

        let rax = b"GTGACCGTTTCTC";
        //   GA---TTGTAACTC
        // GTGACCGTT-T--CTC
        //   GA---TT
        //       GTAACTC
        let mut first = Change::from_seq(2, 9, b"GATT".to_vec(), rax);
        let second = Change::from_seq(6, 13, b"GTAACTC".to_vec(), rax);

        first.merge(&second, rax);

        assert_eq!(first.begin(), 2);
        assert_eq!(first.end_raw(), 13);
        assert_eq!(first.end_err(), 13);
        assert_eq!(first.seq(), &b"GATTGTAACTC".to_vec());
        assert_eq!(first.cigar(), &b"==DDD==I=II===".to_vec());
        assert_eq!(first.edit(), 6);

        let ray = b"TTTGCACTTGTTGCCACAGG";
        // TTTGCACTTGTTGCCACAGG
        //   TGCCCTT
        //            TGCCTCA
        let mut first = Change::from_seq(2, 9, b"TGCCCTT".to_vec(), ray);
        let second = Change::from_seq(11, 18, b"TGCCTCA".to_vec(), ray);

        first.merge(&second, ray);

        assert_eq!(first.begin(), 2);
        assert_eq!(first.end_raw(), 9);
        assert_eq!(first.end_err(), 9);
        assert_eq!(first.seq(), &b"TGCCCTT".to_vec());
        assert_eq!(first.cigar(), &b"===X===".to_vec());
        assert_eq!(first.edit(), 1);

        assert_eq!(second.begin(), 11);
        assert_eq!(second.end_raw(), 18);
        assert_eq!(second.end_err(), 18);
        assert_eq!(second.seq(), &b"TGCCTCA".to_vec());
        assert_eq!(second.cigar(), &b"====X==".to_vec());
        assert_eq!(second.edit(), 1);
    }

    #[test]
    fn merge_containment() {
        let raw = b"GCTATGCGTCAG";

        let mut change = Change {
            begin: 2,
            end: 9,
            seq: b"TATGCGTC".to_vec(),
            cigar: b"===I====".to_vec(),
            edit_distance: 1,
        };
        let glitche = Change {
            begin: 2,
            end: 2,
            seq: b"T".to_vec(),
            cigar: b"I".to_vec(),
            edit_distance: 1,
        };
        let del = Change {
            begin: 4,
            end: 11,
            seq: b"TGCG".to_vec(),
            cigar: b"====DDD".to_vec(),
            edit_distance: 3,
        };

        change.merge(&glitche, raw);

        assert_eq!(change.begin(), 2);
        assert_eq!(change.end_raw(), 9);
        assert_eq!(change.end_err(), 10);
        assert_eq!(change.seq(), &b"TATGCGTC".to_vec());
        assert_eq!(change.cigar(), &b"===I====".to_vec());
        assert_eq!(change.edit(), 1);

        assert!(change.contain(&del));
        change.merge(&del, raw);

        assert_eq!(change.begin(), 2);
        assert_eq!(change.end_raw(), 9);
        assert_eq!(change.end_err(), 10);
        assert_eq!(change.seq(), &b"TATGCGTC".to_vec());
        assert_eq!(change.cigar(), &b"===I====".to_vec());
        assert_eq!(change.edit(), 1);
    }

    #[test]
    fn key() {
        let ray = b"TTTGCACTTGTTGCCACAGG";
        let first = Change::from_seq(2, 9, b"TGCCCTT".to_vec(), ray);
        let second = Change::from_seq(11, 18, b"TGCCTA".to_vec(), ray);

        assert_eq!(
            first.key(),
            ChangeKey {
                begin: 2,
                end_raw: 9,
                end_err: 9
            }
        );
        assert_eq!(
            second.key(),
            ChangeKey {
                begin: 11,
                end_raw: 18,
                end_err: 17
            }
        );
    }

    #[test]
    fn key_cmp() {
        let raw = b"TCAGGAAGATCCACA";
        let first = Change::from_seq(3, 10, b"GGCCGAT".to_vec(), raw);
        let second = Change::from_seq(6, 13, b"AGACCA".to_vec(), raw);

        assert!(first.key() == second.key());

        let rax = b"GTGACCGTTTCTC";
        let first = Change::from_seq(2, 9, b"GATT".to_vec(), rax);
        let second = Change::from_seq(6, 13, b"GTAACTC".to_vec(), rax);

        assert!(first.key() == second.key());

        let ray = b"TTTGCACTTGTTGCCACAGG";
        let first = Change::from_seq(2, 9, b"TGCCCTT".to_vec(), ray);
        let second = Change::from_seq(11, 18, b"TGCCTCA".to_vec(), ray);

        assert!(first.key() < second.key());
        assert!(first.key() != second.key());
    }

    #[test]
    fn identity2edit() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let edits: Vec<f64> = (0..10).map(|_| number_of_edit(rng.gen(), 100)).collect();

        assert_eq!(
            edits,
            [
                47.34425909972262,
                45.72747900968561,
                36.35349008561051,
                59.409824176922335,
                96.56571820450438,
                58.50431538146399,
                26.257557227560657,
                15.074839555058372,
                86.87211108325653,
                99.67479036470401
            ]
        );
    }

    #[test]
    fn glitches() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = model::Glitch::new(1_000.0, 25.0, 25.0).unwrap();

        let glitches = get_glitches(&crate::random_seq(10_000, &mut rng), &model, &mut rng);

        assert_eq!(
            vec![
                ChangeKey {
                    begin: 1560,
                    end_raw: 1561,
                    end_err: 1613
                },
                ChangeKey {
                    begin: 2320,
                    end_raw: 2339,
                    end_err: 2389
                },
                ChangeKey {
                    begin: 4040,
                    end_raw: 4043,
                    end_err: 4047
                },
                ChangeKey {
                    begin: 4550,
                    end_raw: 4558,
                    end_err: 4550
                },
                ChangeKey {
                    begin: 5155,
                    end_raw: 5179,
                    end_err: 5190
                },
                ChangeKey {
                    begin: 5269,
                    end_raw: 5274,
                    end_err: 5297
                },
                ChangeKey {
                    begin: 6193,
                    end_raw: 6244,
                    end_err: 6215
                }
            ],
            glitches.keys().cloned().collect::<Vec<ChangeKey>>()
        );

        assert_eq!(
            vec![
                Change {
                    begin: 1560,
                    end: 1561,
                    seq: b"GTTACATGCGATACACCATGTGATGCGGTTGCCTCCAATCCCGGTAATATAAC".to_vec(),
                    cigar: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII=I".to_vec(),
                    edit_distance: 52
                },
                Change {
                    begin: 2320,
                    end: 2339,
                    seq: b"TGCAACTATTACACGTTTGATTTGACCACCCCACTCCGCCGCTCCCAATTGTAACTTCGATCCTCGTAC"
                        .to_vec(),
                    cigar: b"I=IIIIIIIIII=IIIII=IIII=IIIIII==I==IIIIIII=III==IIII===I==I=IIIIII=I="
                        .to_vec(),
                    edit_distance: 50
                },
                Change {
                    begin: 4040,
                    end: 4043,
                    seq: b"ATCATCG".to_vec(),
                    cigar: b"III==I=".to_vec(),
                    edit_distance: 4
                },
                Change {
                    begin: 4550,
                    end: 4558,
                    seq: b"".to_vec(),
                    cigar: b"DDDDDDDD".to_vec(),
                    edit_distance: 8
                },
                Change {
                    begin: 5155,
                    end: 5179,
                    seq: b"GCCCGGGCACCGCCACTACCTATCGTCGCTCGGGT".to_vec(),
                    cigar: b"IIIIIIIX=IX=II=DDX==D==D==I=I=III=X=D==I".to_vec(),
                    edit_distance: 25
                },
                Change {
                    begin: 5269,
                    end: 5274,
                    seq: b"ATCATCGACGATGAGAGATAAGGGCACG".to_vec(),
                    cigar: b"IIII=IIIII==IIIIII=IIIIIII=I".to_vec(),
                    edit_distance: 23
                },
                Change {
                    begin: 6193,
                    end: 6244,
                    seq: b"ACGACGCAGTCAGAGTCCTCGC".to_vec(),
                    cigar: b"I=DD==DDD==D=DX=DDDDD=DDD=DDDDDDX=D==D==D=DD=D=DDD=X".to_vec(),
                    edit_distance: 34
                }
            ],
            glitches.values().cloned().collect::<Vec<Change>>()
        );
    }

    #[test]
    fn errors() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = model::Error::random(7);

        let mut errors = Changes::new();
        get_error(
            7,
            20.0,
            &crate::random_seq(100, &mut rng),
            &mut errors,
            &model,
            &mut rng,
        );

        assert_eq!(
            vec![
                ChangeKey {
                    begin: 4,
                    end_raw: 17,
                    end_err: 20
                },
                ChangeKey {
                    begin: 21,
                    end_raw: 39,
                    end_err: 38
                },
                ChangeKey {
                    begin: 41,
                    end_raw: 48,
                    end_err: 48
                },
                ChangeKey {
                    begin: 52,
                    end_raw: 60,
                    end_err: 59
                },
                ChangeKey {
                    begin: 65,
                    end_raw: 73,
                    end_err: 73
                },
                ChangeKey {
                    begin: 75,
                    end_raw: 86,
                    end_err: 87
                },
                ChangeKey {
                    begin: 88,
                    end_raw: 98,
                    end_err: 99
                }
            ],
            errors.keys().cloned().collect::<Vec<ChangeKey>>()
        );

        assert_eq!(
            vec![
                Change {
                    begin: 4,
                    end: 17,
                    seq: b"ATTATGAGTACGTAAT".to_vec(),
                    cigar: b"=====I======IIX=".to_vec(),
                    edit_distance: 4
                },
                Change {
                    begin: 21,
                    end: 39,
                    seq: b"TGTTACACTATGTCCTA".to_vec(),
                    cigar: b"=D===II=======DD====".to_vec(),
                    edit_distance: 5
                },
                Change {
                    begin: 41,
                    end: 48,
                    seq: b"TGGTGCC".to_vec(),
                    cigar: b"===X===".to_vec(),
                    edit_distance: 1
                },
                Change {
                    begin: 52,
                    end: 60,
                    seq: b"GTGGAGG".to_vec(),
                    cigar: b"==X====D".to_vec(),
                    edit_distance: 2
                },
                Change {
                    begin: 65,
                    end: 73,
                    seq: b"CTTATATT".to_vec(),
                    cigar: b"=======X".to_vec(),
                    edit_distance: 1
                },
                Change {
                    begin: 75,
                    end: 86,
                    seq: b"ACGGGTAATAAT".to_vec(),
                    cigar: b"==X===I=====".to_vec(),
                    edit_distance: 2
                },
                Change {
                    begin: 88,
                    end: 98,
                    seq: b"GGAGGGCCATG".to_vec(),
                    cigar: b"===X==I====".to_vec(),
                    edit_distance: 2
                }
            ],
            errors.values().cloned().collect::<Vec<Change>>()
        );
    }

    #[test]
    fn assemble_change() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(6);
        let raw = crate::random_seq(100, &mut rng);
        let e_model = model::Error::random(7);

        let mut changes = Changes::new();
        get_error(7, 25.0, &raw, &mut changes, &e_model, &mut rng);

        let (asm, edit) = asm_change(changes, &raw);

        assert_eq!(
            vec![
                ChangeKey {
                    begin: 1,
                    end_raw: 12,
                    end_err: 11
                },
                ChangeKey {
                    begin: 20,
                    end_raw: 34,
                    end_err: 34
                },
                ChangeKey {
                    begin: 34,
                    end_raw: 46,
                    end_err: 45
                },
                ChangeKey {
                    begin: 70,
                    end_raw: 86,
                    end_err: 86
                }
            ],
            asm.keys().cloned().collect::<Vec<ChangeKey>>()
        );

        assert_eq!(
            vec![
                Change {
                    begin: 1,
                    end: 12,
                    seq: vec![67, 71, 67, 84, 84, 71, 71, 67, 67, 67],
                    cigar: vec![88, 61, 61, 61, 61, 61, 61, 68, 61, 61, 88],
                    edit_distance: 3
                },
                Change {
                    begin: 20,
                    end: 34,
                    seq: vec![65, 67, 67, 67, 67, 84, 65, 67, 65, 65, 65, 84, 67, 84],
                    cigar: vec![61, 61, 61, 61, 61, 61, 88, 61, 73, 61, 61, 61, 61, 61, 68],
                    edit_distance: 3
                },
                Change {
                    begin: 34,
                    end: 46,
                    seq: vec![65, 71, 71, 67, 84, 71, 71, 84, 65, 67, 84],
                    cigar: vec![68, 61, 61, 61, 61, 61, 73, 61, 61, 61, 61, 61, 68],
                    edit_distance: 3
                },
                Change {
                    begin: 70,
                    end: 86,
                    seq: vec![84, 84, 67, 67, 67, 67, 71, 65, 65, 67, 84, 65, 71, 84, 67, 84],
                    cigar: vec![61, 61, 61, 61, 88, 61, 61, 61, 61, 73, 61, 61, 61, 61, 61, 68, 61],
                    edit_distance: 3
                }
            ],
            asm.values().cloned().collect::<Vec<Change>>()
        );

        assert_eq!(12.0, edit);
    }

    #[test]
    fn all() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let e_model = model::Error::random(7);
        let g_model = model::Glitch::new(40.0, 5.0, 5.0).unwrap();

        let raw = crate::random_seq(100, &mut rng);

        let (err, cigar, edit) = add_error(0.95, &raw, &e_model, &g_model, &mut rng);

        assert_eq!(b"TTAGATTATAGTACGGTATAGTGGTTCTCTAGTAGCCTCCCGTTGTAGAGGAATCCACTTAATATAACACAGGTATAATCAGGACGGCATTCG".to_vec(), err);
        assert_eq!(b"==========================D==X=I======DDDD=DDDD======================I===================X=========X==".to_vec(), cigar);
        assert_eq!(0.86, edit);

        let raw = crate::random_seq(100, &mut rng);

        let (err, cigar, edit) = add_error(0.75, &raw, &e_model, &g_model, &mut rng);

        assert_eq!(b"TCTGGAAAGATGGTGCTCTAGTACGTTTGGGGCACCCCGTAGGCACACCGGCATGTGGTTTTCCTCGAACAATTTCGTCCCGATTATATACAGGGGACGC".to_vec(), err);
        assert_eq!(b"=X===I===============D=============X======I====D======D====I=====I==X=X===D=D========D===I========IX======".to_vec(), cigar);
        assert_eq!(0.8, edit);
    }
}
