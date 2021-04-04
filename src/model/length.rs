//! Model to get length of reads

/* standard use */

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;

/* local use */
use crate::error::Model;

/// Struct to generate length of fragment
pub struct Length {
    mean: f64,
    stdev: f64,
    dist: rand_distr::Gamma<f64>,
}

impl Length {
    /// Create model from parameter
    pub fn new(mean: f64, stdev: f64) -> Result<Length> {
        if mean <= 0.0 || stdev <= 0.0 {
            anyhow::bail!(Model::LengthParamMustBeUpperThan0);
        }

        let k = mean.powf(2.0) / stdev.powf(2.0);
        let t = stdev.powf(2.0) / mean;

        Ok(Self {
            mean,
            stdev,
            dist: rand_distr::Gamma::new(k, t)?,
        })
    }

    /// Get length from model
    pub fn get_length<R>(&self, rng: &mut R) -> u64
    where
        R: rand::Rng,
    {
        if self.stdev == 0.0 {
            self.mean.round() as u64
        } else {
            self.dist.sample(rng).round() as u64
        }
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn create_dist() {
        assert!(Length::new(0.0, 0.0).is_err());
        assert!(Length::new(0.0, 1.0).is_err());
        assert!(Length::new(-10.0, -0.8).is_err());

        assert!(Length::new(1.0, 1.0).is_ok());
    }

    #[test]
    fn get_value() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let dist = Length::new(50.0, 20.0).unwrap();

        let data: Vec<u64> = (0..100).map(|_| dist.get_length(&mut rng)).collect();

        assert_eq!(
            data,
            vec![
                49, 53, 36, 51, 26, 93, 40, 48, 31, 35, 64, 48, 37, 59, 45, 45, 71, 53, 51, 24, 73,
                45, 24, 47, 86, 61, 55, 44, 31, 56, 72, 44, 31, 28, 59, 46, 87, 46, 80, 47, 38, 36,
                72, 25, 56, 40, 54, 68, 29, 69, 40, 84, 44, 100, 43, 59, 31, 44, 49, 22, 30, 83,
                106, 42, 29, 51, 40, 96, 37, 48, 34, 37, 27, 44, 52, 43, 47, 55, 53, 49, 32, 104,
                44, 62, 39, 32, 77, 42, 54, 31, 66, 49, 28, 16, 37, 40, 67, 56, 55, 23
            ]
        )
    }
}
