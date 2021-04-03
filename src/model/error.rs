//! Model to add error sequence

/* standard use */
use std::str::FromStr;

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;

/* local use */
use crate::error::Model;

type Kmer = Vec<u8>;
type KmerEdit = (Kmer, u64);

/// Struct to load and apply error model
pub struct Error {
    kmer2alts_edit_prob: rustc_hash::FxHashMap<Kmer, (Vec<KmerEdit>, Vec<f64>)>,
}

impl Error {
    /// Load model from an stdin
    pub fn from_stream<R, A>(input: R, rng: &mut A) -> Result<Self>
    where
        R: std::io::Read,
        A: rand::Rng,
    {
        let mut data = rustc_hash::FxHashMap::default();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b';')
            .has_headers(false)
            .flexible(true)
            .from_reader(input);
        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            let mut alts: Vec<(Kmer, u64)> = Vec::new();
            let mut prob = Vec::new();

            for row in record.iter() {
                if row.is_empty() {
                    continue;
                }

                let mut info = row.split(',');

                let kmer: Vec<u8> = info.next().ok_or(Model::ErrorParsing)?.bytes().collect();
                let p = f64::from_str(info.next().ok_or(Model::ErrorParsing)?)?;

                let dist = if !alts.is_empty() {
                    crate::alignment::edit_distance(&alts[0].0, &kmer)
                } else {
                    0
                };

                alts.push((kmer, dist));
                prob.push(p);
            }

            let sum_of_prob: f64 = prob.iter().sum();
            if sum_of_prob < 1.0 {
                alts.push((random_error(&alts[0].0, rng), 1));
                prob.push(1.0 - sum_of_prob);
            }
            data.insert(alts[0].0.clone(), (alts, prob));
        }

        Ok(Self {
            kmer2alts_edit_prob: data,
        })
    }

    /// Add error to a kmer
    pub fn add_errors_to_kmer<R>(&self, kmer: &[u8], rng: &mut R) -> (Kmer, u64)
    where
        R: rand::Rng,
    {
        if let Some(values) = self.kmer2alts_edit_prob.get(kmer) {
            //values -> (Vec<(kmer, edit_distance)>, Vec<weight>)
            let dist = WeightedIndex::new(&values.1).unwrap();
            values.0[dist.sample(rng)].clone()
        } else {
            (kmer.to_vec(), 0)
        }
    }
}

/// Add a single random error in a kmer
pub fn random_error<R>(kmer: &[u8], rng: &mut R) -> Kmer
where
    R: rand::Rng,
{
    let error_pos = rng.gen_range(0..kmer.len());

    let mut new_kmer: Vec<u8> = Vec::with_capacity(kmer.len() + 1);
    new_kmer.extend(&kmer[0..error_pos]);

    // 1 -> substitution, 2 -> insertion, 3 -> deletion
    match rng.gen_range(1..=3) {
        1 => new_kmer.push(crate::random_base_diff(kmer[error_pos], rng)),
        2 => {new_kmer.push(crate::random_base(rng)); new_kmer.push(kmer[error_pos])},
        3 => (),
        _ => unreachable!("contact author with you command line dataset and this message model::Error::random_error"),
    }
    new_kmer.extend(kmer[error_pos + 1..].iter());

    new_kmer
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn random_error_() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        assert_eq!(b"AAAGAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAAAAAG", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAAAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAAATAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"ATAAAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAAAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAAAAAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
        assert_eq!(b"AAGAAAAA", &random_error(b"AAAAAAA", &mut rng)[..]);
    }

    static MODEL: &[u8] = b"ACAGTTG,0.25;ACGGTTG,0.25;ACAGG,0.25;";

    #[test]
    fn read_error_model() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Error::from_stream(MODEL, &mut rng).unwrap();

        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACGGTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACAGTTG".to_vec(), 0)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACAGG".to_vec(), 2)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACGGTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACATTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACAGG".to_vec(), 2)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACATTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACAGG".to_vec(), 2)
        );
    }
}
