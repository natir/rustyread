//! Model to generate qscore

/* standard use */
use std::str::FromStr;

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;

/* local use */
use crate::error::Model;

type Cigar = Vec<u8>;
type Scores = Vec<u8>;
type Weights = Vec<f64>;

/// Struct to load and apply quality model
pub struct Quality {
    max_k: usize,
    cigar2score_weight: rustc_hash::FxHashMap<Cigar, (Scores, Weights)>,
}

impl Quality {
    /// Load model from an stdin
    pub fn from_stream<RNG>(input: RNG) -> Result<Self>
    where
        RNG: std::io::Read,
    {
        let mut data = rustc_hash::FxHashMap::default();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b';')
            .has_headers(false)
            .flexible(true)
            .from_reader(input);
        let mut records = reader.records().skip(1);

        let mut kmer_length = 0;

        while let Some(Ok(record)) = records.next() {
            let mut rows = record.iter();
            let cigar: Vec<u8> = rows.next().ok_or(Model::QualityParsing)?.bytes().collect();
            let _cigar_count = rows.next();
            let scores_weight = rows.next().ok_or(Model::QualityParsing)?.split(',');

            let mut scores = Vec::new();
            let mut weights = Vec::new();

            for row in scores_weight {
                if row.is_empty() {
                    continue;
                }

                let mut info = row.split(':');
                scores.push(u8::from_str(info.next().ok_or(Model::ErrorParsing)?)?);
                weights.push(f64::from_str(info.next().ok_or(Model::ErrorParsing)?)?);
            }

            if cigar.len() > kmer_length {
                kmer_length = cigar.len();
            }
            data.insert(cigar, (scores, weights));
        }

        if !data.contains_key(&vec![b'='])
            || !data.contains_key(&vec![b'X'])
            || !data.contains_key(&vec![b'I'])
        {
            Err(anyhow::Error::new(Model::QualityNotMinimalCigarString))
        } else {
            Ok(Self {
                max_k: kmer_length,
                cigar2score_weight: data,
            })
        }
    }

    /// Build a random quality score model
    pub fn random() -> Self {
        let mut data = rustc_hash::FxHashMap::default();

        data.insert(b"=".to_vec(), ((1..=20).collect(), vec![1.0; 20]));
        data.insert(b"X".to_vec(), ((1..=20).collect(), vec![1.0; 20]));
        data.insert(b"I".to_vec(), ((1..=20).collect(), vec![1.0; 20]));

        Self {
            max_k: 1,
            cigar2score_weight: data,
        }
    }

    /// Build an ideal quality score model
    pub fn ideal() -> Self {
        let mut data = rustc_hash::FxHashMap::default();

        data.insert(b"X".to_vec(), ((1..=3).collect(), vec![1.0; 3]));
        data.insert(b"I".to_vec(), ((1..=3).collect(), vec![1.0; 3]));

        data.insert(b"=".to_vec(), ((4..=7).collect(), vec![1.0; 4]));
        data.insert(b"===".to_vec(), ((8..=20).collect(), vec![1.0; 13]));
        data.insert(b"=====".to_vec(), ((21..=30).collect(), vec![1.0; 10]));
        data.insert(b"=======".to_vec(), ((31..=40).collect(), vec![1.0; 10]));
        data.insert(b"=========".to_vec(), ((41..=50).collect(), vec![1.0; 10]));

        Self {
            max_k: 9,
            cigar2score_weight: data,
        }
    }

    /// Generate error associate to a cigar string with odd length
    pub fn get_qscore<RNG>(&self, cigar: &[u8], rng: &mut RNG) -> Result<u8>
    where
        RNG: rand::Rng,
    {
        let mut c = cigar;

        while !c.is_empty() {
            if c.len() % 2 == 0 {
                anyhow::bail!(Model::QualityCigarLenNotOdd);
            }

            if let Some((scores, weights)) = self.cigar2score_weight.get(c) {
                let dist = WeightedIndex::new(weights).unwrap();
                return Ok(scores[dist.sample(rng)] + 33);
            } else {
                c = &c[1..c.len() - 1];
            }
        }

        Err(anyhow::Error::new(Model::QualityNotMinimalCigarString))
    }

    /// Kmer length of model
    pub fn max_k(&self) -> usize {
        self.max_k
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    static MODEL: &[u8] = b"overall;3;1:0.000065,2:0.002683,3:0.011969,4:0.022982,5:0.03079,6:0.036889,7:0.042461,8:0.047495,9:0.051971,10:0.055876,11:0.058956,12:0.0606,13:0.060621,14:0.059446,15:0.056487,16:0.052028,17:0.047379,18:0.043445,19:0.040219,20:0.037767,21:0.035767,22:0.033816,23:0.031369,24:0.027969,25:0.023457,26:0.018417,27:0.008983,28:0.000093,
I;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
X;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
=;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
";

    #[test]
    fn read_quality_model() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Quality::from_stream(MODEL).unwrap();

        let m: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"=", &mut rng).unwrap())
            .collect();
        assert_eq!(
            m,
            vec![
                44, 44, 45, 42, 37, 43, 47, 50, 39, 35, 53, 44, 42, 39, 44, 40, 36, 44, 37, 46, 50,
                44, 44, 38, 42, 43, 48, 40, 42, 46, 43, 51, 53, 51, 46, 40, 45, 42, 38, 43, 57, 45,
                43, 37, 39, 42, 43, 55, 50, 44, 48, 49, 45, 40, 42, 42, 40, 40, 45, 43, 55, 39, 43,
                44, 39, 42, 40, 44, 46, 38, 43, 45, 51, 44, 43, 43, 53, 37, 43, 47, 41, 42, 42, 47,
                46, 45, 38, 42, 47, 53, 38, 50, 44, 48, 48, 41, 39, 48, 51, 41
            ]
        );

        let i: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"I", &mut rng).unwrap())
            .collect();
        assert_eq!(
            i,
            vec![
                39, 41, 53, 47, 41, 39, 50, 39, 42, 42, 45, 51, 39, 44, 43, 48, 44, 51, 39, 50, 40,
                55, 57, 39, 56, 51, 40, 43, 38, 46, 48, 39, 42, 43, 49, 40, 38, 38, 44, 39, 41, 46,
                41, 39, 39, 47, 43, 41, 44, 45, 42, 51, 43, 37, 45, 36, 46, 39, 44, 36, 40, 46, 52,
                43, 43, 38, 46, 41, 38, 54, 41, 40, 35, 44, 40, 50, 50, 43, 42, 55, 45, 44, 38, 43,
                47, 56, 44, 47, 40, 38, 35, 47, 45, 41, 47, 41, 47, 56, 48, 54
            ]
        );

        let x: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"X", &mut rng).unwrap())
            .collect();
        assert_eq!(
            x,
            vec![
                49, 47, 38, 39, 46, 51, 40, 42, 43, 49, 48, 41, 42, 38, 44, 48, 46, 41, 51, 51, 47,
                41, 38, 39, 46, 46, 47, 49, 51, 41, 38, 43, 45, 50, 41, 42, 41, 44, 42, 37, 44, 37,
                44, 44, 42, 51, 45, 46, 43, 46, 51, 42, 44, 47, 55, 40, 40, 39, 40, 45, 47, 48, 39,
                39, 43, 46, 41, 47, 39, 47, 43, 42, 50, 54, 41, 57, 40, 38, 48, 51, 37, 39, 43, 45,
                45, 38, 44, 48, 47, 46, 50, 42, 47, 37, 40, 42, 38, 50, 51, 37
            ]
        );

        assert!(model.get_qscore(b"==", &mut rng).is_err());
        assert!(model.get_qscore(b"===", &mut rng).is_ok());
        assert!(model.get_qscore(b"=I=", &mut rng).is_ok());
        assert!(model.get_qscore(b"=X=", &mut rng).is_ok());

        assert!(model.get_qscore(b"", &mut rng).is_err());
        assert!(model.get_qscore(b"bepo", &mut rng).is_err());
    }

    #[test]
    fn random() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Quality::random();

        let m: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"=", &mut rng).unwrap())
            .collect();
        assert_eq!(
            m,
            vec![
                44, 44, 46, 42, 34, 42, 48, 50, 36, 34, 52, 44, 41, 36, 44, 38, 34, 44, 35, 46, 50,
                43, 44, 35, 41, 42, 49, 38, 42, 47, 42, 51, 52, 51, 46, 38, 45, 41, 35, 42, 53, 46,
                43, 34, 37, 41, 43, 53, 51, 44, 49, 50, 45, 38, 41, 41, 38, 38, 46, 43, 53, 37, 43,
                44, 36, 40, 38, 45, 47, 35, 43, 46, 51, 45, 43, 42, 52, 34, 43, 49, 40, 40, 41, 48,
                47, 45, 35, 41, 48, 52, 35, 50, 45, 49, 49, 40, 37, 49, 51, 40
            ]
        );

        let i: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"I", &mut rng).unwrap())
            .collect();
        assert_eq!(
            i,
            vec![
                36, 39, 52, 48, 39, 36, 50, 37, 42, 41, 45, 51, 36, 45, 43, 49, 44, 51, 36, 50, 37,
                53, 53, 36, 53, 51, 38, 43, 35, 47, 49, 36, 41, 43, 50, 38, 35, 35, 44, 36, 39, 47,
                40, 36, 37, 48, 43, 39, 45, 45, 42, 51, 43, 34, 46, 34, 47, 36, 44, 34, 38, 47, 52,
                43, 42, 35, 47, 39, 35, 52, 40, 37, 34, 43, 38, 51, 51, 42, 41, 53, 46, 44, 35, 43,
                48, 53, 44, 48, 37, 35, 34, 48, 46, 39, 48, 39, 48, 53, 49, 53
            ]
        );

        let x: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"X", &mut rng).unwrap())
            .collect();
        assert_eq!(
            x,
            vec![
                50, 48, 35, 36, 47, 51, 38, 40, 43, 50, 49, 39, 42, 35, 45, 49, 47, 39, 51, 51, 48,
                40, 36, 36, 47, 47, 48, 50, 51, 39, 35, 42, 45, 51, 39, 41, 40, 45, 41, 34, 44, 34,
                44, 45, 41, 51, 46, 47, 43, 46, 51, 41, 44, 48, 53, 37, 38, 36, 38, 46, 48, 49, 36,
                36, 43, 47, 40, 48, 36, 48, 43, 41, 51, 52, 39, 53, 38, 35, 49, 51, 34, 36, 43, 45,
                46, 35, 44, 49, 48, 47, 51, 41, 48, 34, 38, 41, 35, 51, 51, 34
            ]
        );

        assert!(model.get_qscore(b"==", &mut rng).is_err());
        assert!(model.get_qscore(b"===", &mut rng).is_ok());
        assert!(model.get_qscore(b"=I=", &mut rng).is_ok());
        assert!(model.get_qscore(b"=X=", &mut rng).is_ok());

        assert!(model.get_qscore(b"", &mut rng).is_err());
        assert!(model.get_qscore(b"bepo", &mut rng).is_err());
    }

    #[test]
    fn ideal() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Quality::ideal();

        let m: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"=", &mut rng).unwrap())
            .collect();
        assert_eq!(
            m,
            vec![
                39, 39, 39, 38, 37, 38, 39, 40, 37, 37, 40, 39, 38, 37, 39, 37, 37, 39, 37, 39, 40,
                38, 39, 37, 38, 38, 40, 37, 38, 39, 38, 40, 40, 40, 39, 37, 39, 38, 37, 38, 40, 39,
                38, 37, 37, 38, 38, 40, 40, 39, 40, 40, 39, 37, 38, 38, 37, 37, 39, 38, 40, 37, 38,
                39, 37, 38, 37, 39, 39, 37, 38, 39, 40, 39, 38, 38, 40, 37, 38, 40, 38, 38, 38, 39,
                39, 39, 37, 38, 39, 40, 37, 40, 39, 40, 40, 38, 37, 40, 40, 38
            ]
        );

        let i: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"I", &mut rng).unwrap())
            .collect();
        assert_eq!(
            i,
            vec![
                34, 34, 36, 36, 34, 34, 36, 34, 35, 35, 35, 36, 34, 35, 35, 36, 35, 36, 34, 36, 34,
                36, 36, 34, 36, 36, 34, 35, 34, 35, 36, 34, 35, 35, 36, 34, 34, 34, 35, 34, 34, 36,
                34, 34, 34, 36, 35, 34, 35, 35, 35, 36, 35, 34, 35, 34, 35, 34, 35, 34, 34, 36, 36,
                35, 35, 34, 36, 34, 34, 36, 34, 34, 34, 35, 34, 36, 36, 35, 35, 36, 35, 35, 34, 35,
                36, 36, 35, 36, 34, 34, 34, 36, 35, 34, 36, 34, 36, 36, 36, 36
            ]
        );

        let x: Vec<u8> = (0..100)
            .map(|_| model.get_qscore(b"X", &mut rng).unwrap())
            .collect();
        assert_eq!(
            x,
            vec![
                36, 36, 34, 34, 36, 36, 34, 35, 35, 36, 36, 34, 35, 34, 35, 36, 35, 34, 36, 36, 36,
                34, 34, 34, 36, 36, 36, 36, 36, 34, 34, 35, 35, 36, 34, 35, 34, 35, 35, 34, 35, 34,
                35, 35, 35, 36, 35, 36, 35, 35, 36, 35, 35, 36, 36, 34, 34, 34, 34, 35, 36, 36, 34,
                34, 35, 36, 34, 36, 34, 36, 35, 35, 36, 36, 34, 36, 34, 34, 36, 36, 34, 34, 35, 35,
                35, 34, 35, 36, 36, 36, 36, 35, 36, 34, 34, 35, 34, 36, 36, 34
            ]
        );

        assert!(model.get_qscore(b"==", &mut rng).is_err());
        assert!(model.get_qscore(b"===", &mut rng).is_ok());
        assert!(model.get_qscore(b"=I=", &mut rng).is_ok());
        assert!(model.get_qscore(b"=X=", &mut rng).is_ok());

        assert!(model.get_qscore(b"", &mut rng).is_err());
        assert!(model.get_qscore(b"bepo", &mut rng).is_err());
    }
}
