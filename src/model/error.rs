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
type KmerEditWeight = (Vec<KmerEdit>, Vec<f64>);

/// Struct to load and apply error model
pub struct Error {
    length: usize,
    kmer2alts_edit_prob: Option<rustc_hash::FxHashMap<Kmer, KmerEditWeight>>,
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
        let mut kmer_length = 0;

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

            let key = alts.remove(0);
            prob.remove(0);

            kmer_length = key.0.len();
            data.insert(key.0.clone(), (alts, prob));
        }

        Ok(Self {
            length: kmer_length,
            kmer2alts_edit_prob: Some(data),
        })
    }

    /// Setup a random error model
    pub fn random(k: usize) -> Self {
        Self {
            length: k,
            kmer2alts_edit_prob: None,
        }
    }

    /// Add error to a kmer
    pub fn add_errors_to_kmer<R>(&self, kmer: &[u8], rng: &mut R) -> (Kmer, u64)
    where
        R: rand::Rng,
    {
        if let Some(data) = &self.kmer2alts_edit_prob {
            if let Some(values) = data.get(kmer) {
                let dist = WeightedIndex::new(&values.1).unwrap();
                values.0[dist.sample(rng)].clone()
            } else {
                (random_error(kmer, rng), 1)
            }
        } else {
            (random_error(kmer, rng), 1)
        }
    }

    /// Kmer length of model
    pub fn k(&self) -> usize {
        self.length
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
            (b"ACAGG".to_vec(), 2)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACGGTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACAGG".to_vec(), 2)
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
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACATTTG".to_vec(), 1)
        );
        assert_eq!(
            model.add_errors_to_kmer(b"ACAGTTG", &mut rng),
            (b"ACATTTG".to_vec(), 1)
        );
    }

    #[test]
    fn random() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Error::random(7);

        let kmers: Vec<(Vec<u8>, u64)> = (0..10)
            .map(|_| model.add_errors_to_kmer(b"ACAGTTG", &mut rng))
            .collect();

        assert_eq!(
            vec![
                (b"ACATTTG".to_vec(), 1),
                (b"ACGTTG".to_vec(), 1),
                (b"CAGTTG".to_vec(), 1),
                (b"ACAGTTTG".to_vec(), 1),
                (b"ATAGTTG".to_vec(), 1),
                (b"CAGTTG".to_vec(), 1),
                (b"ACAAGTTG".to_vec(), 1),
                (b"ACGAGTTG".to_vec(), 1),
                (b"AGCAGTTG".to_vec(), 1),
                (b"AAGTTG".to_vec(), 1)
            ],
            kmers
        );
    }
}
