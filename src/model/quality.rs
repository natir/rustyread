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
    cigar2score_weight: rustc_hash::FxHashMap<Cigar, (Scores, Weights)>,
}

impl Quality {
    /// Load model from an stdin
    pub fn from_stream<R>(input: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        let mut data = rustc_hash::FxHashMap::default();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b';')
            .has_headers(false)
            .flexible(true)
            .from_reader(input);
        let mut records = reader.records().skip(1);

        while let Some(Ok(record)) = records.next() {
            let mut rows = record.iter();
            let cigar = rows.next().ok_or(Model::QualityParsing)?.bytes().collect();
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

            data.insert(cigar, (scores, weights));
        }

        if !data.contains_key(&vec![b'='])
            || !data.contains_key(&vec![b'X'])
            || !data.contains_key(&vec![b'I'])
        {
            Err(anyhow::Error::new(Model::QualityNotMinimalCigarString))
        } else {
            Ok(Self {
                cigar2score_weight: data,
            })
        }
    }

    /// Generate error associate to a cigar string with odd length
    pub fn get_qscore<R>(&self, cigar: &[u8], rng: &mut R) -> Result<u8>
    where
        R: rand::Rng,
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
}
