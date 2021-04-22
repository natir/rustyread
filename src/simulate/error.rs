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
    /// Build a new Change warning CIGAR String and edit_distance isn't compute use update_align before
    pub fn new(begin: usize, end: usize, seq: Seq) -> Self {
        Self {
            begin,
            end,
            seq,
            cigar: b"".to_vec(),
            edit_distance: 0,
        }
    }

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
    pub fn edit(&self) -> f64 {
        self.edit_distance as f64
    }

    #[allow(clippy::suspicious_operation_groupings)]
    /// Return true if other is contain in self
    /// This function assume other.begin is greater than self.begin
    pub fn contain(&self, other: &Change) -> bool {
        other.end_raw() <= self.end_err()
            || other.end_raw() <= self.end_raw()
            || other.end_err() <= self.end_err()
            || other.end_err() <= self.end_err()
    }

    /// Return True if self and Change overlap
    /// This function assume other.begin is greater than self.begin
    pub fn overlap(&self, other: &Change) -> bool {
        if self.contain(other) {
            false
        } else {
            other.begin() < self.end_err() || other.begin() < self.end_raw()
        }
    }

    /// Try to merge two Change.
    /// This function assume other.begin is greater than self.begin
    ///
    /// Return an None if self and other Change can't be merge
    /// Return Some(i64) change in edit distance
    pub fn merge(&mut self, other: &Change, original: &[u8]) -> Option<i64> {
        if !self.overlap(other) {
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
        }

        let prev_ed = self.edit();
        self.update_align(original);

        Some(self.edit() as i64 - prev_ed as i64)
    }

    /// Update alignment info of error
    pub fn update_align(&mut self, original: &[u8]) {
        let (edit_distance, cigar) = if self.end_raw() > original.len() {
            crate::alignment::align(&self.seq, &original[self.begin()..])
        } else {
            crate::alignment::align(&self.seq, &original[self.begin()..self.end_raw()])
        };

        self.cigar = cigar.into_vec();
        self.edit_distance = edit_distance;
    }
}

pub type Changes = Vec<Change>;

trait AbsChanges {
    fn add_change(&mut self, change: Change, raw: &[u8]);

    fn insert_change(&mut self, change: Change, idx: usize);
}

impl AbsChanges for Changes {
    fn add_change(&mut self, mut change: Change, raw: &[u8]) {
        //println!("{:?}", change);

        if self.is_empty() {
            self.push(change);
            return;
        }

        match self.binary_search_by_key(&change.begin(), |x| x.begin()) {
            Ok(idx) => {
                if self[idx].overlap(&change) {
                    self[idx].merge(&change, raw);
                    if idx + 1 < self.len() && self[idx].overlap(&self[idx + 1]) {
                        let next = self.remove(idx + 1);
                        self[idx].merge(&next, raw);
                    }
                }
            }
            Err(idx) => {
                //println!("ERR idx {} len {}", idx, self.len());
                if idx == 0 {
                    if idx + 1 < self.len() {
                        if self[idx + 1].overlap(&change) {
                            self[idx + 1].merge(&change, raw);
                        } else if change.overlap(&self[idx]) {
                            change.merge(&self[idx], raw);
                            self[idx] = change;
                        } else {
                            change.update_align(&raw);
                            self.insert_change(change, idx);
                        }
                    } else {
                        change.update_align(&raw);
                        self.insert_change(change, idx);
                    }
                } else if idx > 0
                    && self[idx - 1].begin() <= change.begin()
                    && self[idx - 1].contain(&change)
                {
                    //println!("contain");
                } else if idx > 0 && self[idx - 1].overlap(&change) {
                    self[idx - 1].merge(&change, raw);
                    if idx < self.len() && self[idx - 1].overlap(&self[idx]) {
                        let next = self.remove(idx);
                        self[idx - 1].merge(&next, raw);
                    }
                } else if idx < self.len()
                    && change.begin() <= self[idx].begin()
                    && change.contain(&self[idx])
                {
                    //println!("contain");
                } else if idx < self.len() && change.overlap(&self[idx]) {
                    change.merge(&self[idx], raw);
                    self[idx] = change;
                    if self[idx].overlap(&self[idx]) {
                        let next = self.remove(idx);
                        self[idx].merge(&next, raw);
                    }
                } else {
                    change.update_align(&raw);
                    self.insert_change(change, idx);
                }
            }
        }
    }

    fn insert_change(&mut self, change: Change, idx: usize) {
        self.insert(idx, change);
    }
}

/// From identity an seq length compute number of edit
pub fn number_of_edit(target: f64, length: usize) -> f64 {
    (1.0 - target) * length as f64
}

/// Apply error on read
pub fn sequence(
    identity: f64,
    seq: &[u8],
    error_model: &model::Error,
    glitch_model: &model::Glitch,
    rng: &mut rand::rngs::StdRng,
) -> (Seq, Cigar, f64) {
    let k = error_model.k();
    let target = number_of_edit(identity, seq.len());

    let mut changes = Changes::with_capacity(target.round() as usize);

    add_glitches(seq, &mut changes, glitch_model, rng);
    add_error(k, target, seq, &mut changes, error_model, rng);

    let mut real_edit = 0.0;
    let mut pos_in_raw = 0;
    let mut err: Vec<u8> =
        Vec::with_capacity(seq.len() + number_of_edit(target, seq.len()) as usize);
    let mut cig: Vec<u8> =
        Vec::with_capacity(seq.len() + number_of_edit(target, seq.len()) as usize);
    for change in changes {
        if change.begin() < pos_in_raw {
            continue;
        }

        real_edit += change.edit() as f64;

        err.extend(&seq[pos_in_raw..change.begin()]);
        cig.extend(std::iter::repeat(b'=').take(change.begin() - pos_in_raw));

        err.extend(change.seq());
        cig.extend(change.cigar());

        pos_in_raw = change.end_raw();
    }

    if pos_in_raw < seq.len() {
        err.extend(&seq[pos_in_raw..]);
        cig.extend(std::iter::repeat(b'=').take(seq.len() - pos_in_raw));
    }

    log::debug!(
        "identtity {} real {}",
        identity,
        1.0 - (real_edit / seq.len() as f64)
    );
    (err, cig, (1.0 - (real_edit / seq.len() as f64)))
}

/// Create Change correspond to glitches
pub fn add_glitches(
    raw: &[u8],
    changes: &mut Changes,
    model: &model::Glitch,
    rng: &mut rand::rngs::StdRng,
) {
    let mut position = 0;

    while let Some((begin, end, seq)) = model.get_glitch(rng) {
        position += begin;
        if position > raw.len() || position + end > raw.len() {
            break;
        }

        changes.push(Change::from_seq(
            position,
            position + (end - begin),
            seq,
            raw,
        ));
    }
}

/// Add Change correspond to error in changes
pub fn add_error(
    k: usize,
    target: f64,
    raw: &[u8],
    changes: &mut Changes,
    error_model: &model::Error,
    rng: &mut rand::rngs::StdRng,
) {
    let mut sum_of_edit = 0.0;

    while sum_of_edit < target {
        let pos = rng.gen_range(0..(raw.len() - k));

        let (kmer, edit) = error_model.add_errors_to_kmer(&raw[pos..pos + k], rng);
        if edit == 0 {
            continue;
        }

        changes.add_change(Change::new(pos, pos + k, kmer), raw);
        sum_of_edit = changes.iter().map(|x| x.edit()).sum();
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
        assert_eq!(change.edit(), 2.0);

        change = Change::from_seq(6, 13, b"GTTTG".to_vec(), raw);
        assert_eq!(change.begin(), 6);
        assert_eq!(change.end_raw(), 13);
        assert_eq!(change.end_err(), 11);
        assert_eq!(change.seq(), &b"GTTTG".to_vec());
        assert_eq!(change.cigar(), &b"==DD===".to_vec());
        assert_eq!(change.edit(), 2.0);

        change = Change::from_seq(6, 13, b"GTCAACTTG".to_vec(), raw);
        assert_eq!(change.begin(), 6);
        assert_eq!(change.end_raw(), 13);
        assert_eq!(change.end_err(), 15);
        assert_eq!(change.seq(), &b"GTCAACTTG".to_vec());
        assert_eq!(change.cigar(), &b"===II====".to_vec());
        assert_eq!(change.edit(), 2.0);
    }

    #[test]
    fn contain() {
        init();

        let raw = b"TCAGGAAGATCCACA";
        let first = Change::from_seq(3, 15, b"GGCCGAT".to_vec(), raw);
        let second = Change::from_seq(6, 13, b"AGACCA".to_vec(), raw);

        assert_eq!(first.contain(&second), true);

        let ray = b"TTTGCACTTGTTGCCACAGG";
        let first = Change::from_seq(2, 9, b"TGCCCTT".to_vec(), ray);
        let second = Change::from_seq(11, 18, b"TGCCTCA".to_vec(), ray);

        assert_eq!(first.contain(&second), false);
    }

    #[test]
    fn overlap() {
        init();

        let raw = b"TCAGGAAGATCCACA";
        let first = Change::from_seq(3, 10, b"GGCCGAT".to_vec(), raw);
        let second = Change::from_seq(6, 13, b"AGACCA".to_vec(), raw);

        assert_eq!(first.overlap(&second), true);

        let ray = b"TTTGCACTTGTTGCCACAGG";
        let first = Change::from_seq(2, 9, b"TGCCCTT".to_vec(), ray);
        let second = Change::from_seq(11, 18, b"TGCCTCA".to_vec(), ray);

        assert_eq!(first.overlap(&second), false);
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
        assert_eq!(first.edit(), 3.0);

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
        assert_eq!(first.edit(), 6.0);

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
        assert_eq!(first.edit(), 1.0);

        assert_eq!(second.begin(), 11);
        assert_eq!(second.end_raw(), 18);
        assert_eq!(second.end_err(), 18);
        assert_eq!(second.seq(), &b"TGCCTCA".to_vec());
        assert_eq!(second.cigar(), &b"====X==".to_vec());
        assert_eq!(second.edit(), 1.0);
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
        assert_eq!(change.edit(), 1.0);

        assert!(change.contain(&del));
        change.merge(&del, raw);

        assert_eq!(change.begin(), 2);
        assert_eq!(change.end_raw(), 9);
        assert_eq!(change.end_err(), 10);
        assert_eq!(change.seq(), &b"TATGCGTC".to_vec());
        assert_eq!(change.cigar(), &b"===I====".to_vec());
        assert_eq!(change.edit(), 1.0);
    }

    #[test]
    fn changes() {
        let raw = vec![b'X'; 150];

        let mut input = vec![
            Change {
                begin: 37,
                end: 50,
                seq: vec![b'A'; 13],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 90,
                end: 100,
                seq: vec![b'B'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 10,
                end: 30,
                seq: vec![b'C'; 20],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 60,
                end: 70,
                seq: vec![b'D'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 25,
                end: 35,
                seq: vec![b'E'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 30,
                end: 40,
                seq: vec![b'F'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 90,
                end: 95,
                seq: vec![b'G'; 5],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 105,
                end: 120,
                seq: vec![b'H'; 15],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 90,
                end: 110,
                seq: vec![b'I'; 20],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 140,
                end: 150,
                seq: vec![b'J'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
            Change {
                begin: 135,
                end: 145,
                seq: vec![b'K'; 10],
                cigar: vec![b'I'; 0],
                edit_distance: 12,
            },
        ];
        let mut changes = Changes::new();

        for change in input.drain(..) {
            changes.add_change(change, &raw);
        }

        assert_eq!(
            vec![
                Change {
                    begin: 10,
                    end: 50,
                    seq: b"CCCCCCCCCCCCCCCCCCCCEEEEEFFFFFAAAAAAAAAA".to_vec(),
                    cigar: b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX".to_vec(),
                    edit_distance: 40
                },
                Change {
                    begin: 60,
                    end: 70,
                    seq: b"DDDDDDDDDD".to_vec(),
                    cigar: b"XXXXXXXXXX".to_vec(),
                    edit_distance: 10,
                },
                Change {
                    begin: 90,
                    end: 120,
                    seq: b"BBBBBBBBBBIIIIIIIIIIHHHHHHHHHH".to_vec(),
                    cigar: b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX".to_vec(),
                    edit_distance: 30,
                },
                Change {
                    begin: 135,
                    end: 150,
                    seq: b"KKKKKKKKKKJJJJJ".to_vec(),
                    cigar: b"XXXXXXXXXXXXXXX".to_vec(),
                    edit_distance: 15
                },
            ],
            changes
        );
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

        let mut glitches = Changes::new();
        add_glitches(
            &crate::random_seq(10_000, &mut rng),
            &mut glitches,
            &model,
            &mut rng,
        );

        assert_eq!(
            vec![
                Change {
                    begin: 1160,
                    end: 1179,
                    seq: b"TGCAACTATTACACGTTTGATTTGACCACCCCACTCCGCCGCTCCCAATTGTAACTTCGATCCTCGTAC"
                        .to_vec(),
                    cigar: b"II=II=III=IIII=II==IIII=III=III==IIIIIII=IIII=III=I=I=IIII=IIIIII==I="
                        .to_vec(),
                    edit_distance: 50
                },
                Change {
                    begin: 1360,
                    end: 1361,
                    seq: b"GTTACATGCGATACACCATGTGATGCGGTTGCCTCCAATCCCGGTAATATAAC".to_vec(),
                    cigar: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII=IIIIIIIII".to_vec(),
                    edit_distance: 52
                },
                Change {
                    begin: 2700,
                    end: 2703,
                    seq: b"ATCATCG".to_vec(),
                    cigar: b"X=I=III".to_vec(),
                    edit_distance: 5
                },
                Change {
                    begin: 3625,
                    end: 3633,
                    seq: b"".to_vec(),
                    cigar: b"DDDDDDDD".to_vec(),
                    edit_distance: 8
                },
                Change {
                    begin: 4447,
                    end: 4452,
                    seq: b"ATCATCGACGATGAGAGATAAGGGCACG".to_vec(),
                    cigar: b"IIIIIIIIIIIII=I=I==IIIIIII=I".to_vec(),
                    edit_distance: 23
                },
                Change {
                    begin: 4801,
                    end: 4825,
                    seq: b"GCCCGGGCACCGCCACTACCTATCGTCGCTCGGGT".to_vec(),
                    cigar: b"IIIXXX=X=XXX===II=D=X=I=X=II=I==II=I".to_vec(),
                    edit_distance: 22
                },
                Change {
                    begin: 5497,
                    end: 5548,
                    seq: b"ACGACGCAGTCAGAGTCCTCGC".to_vec(),
                    cigar: b"DDDDDDDDDDD==DDDD=DDD=DDDD===D=D==D==I=X=DDX==D=DD=X".to_vec(),
                    edit_distance: 34
                }
            ],
            glitches
        );
    }

    #[test]
    fn errors() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = model::Error::random(7);

        let mut errors = Changes::new();
        add_error(
            7,
            20.0,
            &crate::random_seq(100, &mut rng),
            &mut errors,
            &model,
            &mut rng,
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
                    end: 59,
                    seq: b"GTGGAGG".to_vec(),
                    cigar: b"==X====".to_vec(),
                    edit_distance: 1
                },
                Change {
                    begin: 60,
                    end: 73,
                    seq: b"AATCCACTATATTC".to_vec(),
                    cigar: b"I======D====I=X".to_vec(),
                    edit_distance: 4
                },
                Change {
                    begin: 75,
                    end: 97,
                    seq: b"ACGGGTAATAATGGGGACGGCATG".to_vec(),
                    cigar: b"==X===I=====XX=========I".to_vec(),
                    edit_distance: 5
                },
            ],
            errors
        );
    }

    #[test]
    fn all() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let e_model = model::Error::random(7);
        let g_model = model::Glitch::new(40.0, 5.0, 5.0).unwrap();

        let raw = crate::random_seq(100, &mut rng);

        let (err, cigar, edit) = sequence(0.95, &raw, &e_model, &g_model, &mut rng);

        assert_eq!(b"TTAGATTATAGTACGGTACTGGTTCTTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCAGGACGGCATGCG".to_vec(), err);
        assert_eq!(b"==================DDX=====D=I===========================================================X============".to_vec(), cigar);
        assert_eq!(0.94, edit);

        let raw = crate::random_seq(100, &mut rng);

        let (err, cigar, edit) = sequence(0.75, &raw, &e_model, &g_model, &mut rng);

        assert_eq!(b"GCAGGGATTAGAGAAGTGGTCTCCTTAGATACGTTTGGGCATACCTTAAAACGACCCCGGCTAGTGGTTTTTTTGTGACATACTTCGCGTCTCACGATTAACA".to_vec(), err);
        assert_eq!(b"===I======D=I====D====D==I=I===========D=====I=D==X==XX====II=======D==I=====X===I======X==IX====I======D==D===".to_vec(), cigar);
        assert_eq!(0.75, edit);
    }
}
