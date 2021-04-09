//! Model to add glitches sequence

/* standard use */

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;

/* local use */

/// Struct to generate glitches in fragment
pub struct Glitches {
    distance: Option<rand_distr::Geometric>,
    length_size: rand_distr::Geometric,
    length_skip: rand_distr::Geometric,
}

impl Glitches {
    /// Create model from parameter
    pub fn new(rate: f64, size: f64, skip: f64) -> Result<Glitches> {
        let distance = if rate != 0.0 {
            if rate > 1.0 {
                Some(rand_distr::Geometric::new(1.0 / rate)?)
            } else {
                Some(rand_distr::Geometric::new(1.0)?)
            }
        } else {
            None
        };

        let length_size = if size > 1.0 {
            rand_distr::Geometric::new(1.0 / size)?
        } else {
            rand_distr::Geometric::new(1.0)?
        };

        let length_skip = if skip > 1.0 {
            rand_distr::Geometric::new(1.0 / skip)?
        } else {
            rand_distr::Geometric::new(1.0)?
        };

        Ok(Self {
            distance,
            length_size,
            length_skip,
        })
    }

    /// Get glitch
    pub fn get_glitch<R>(&self, rng: &mut R) -> Option<(usize, usize, Vec<u8>)>
    where
        R: rand::Rng,
    {
        if let Some(dist) = self.distance {
            let begin = dist.sample(rng) as usize;
            let end = begin + self.length_skip.sample(rng) as usize;
            let seq = crate::random_seq(self.length_size.sample(rng) as usize, rng);

            Some((begin, end, seq))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn create_dist() {
        assert!(Glitches::new(10_000.0, 25.0, 25.0).is_ok());
        assert!(Glitches::new(10_000.0, 25.0, 0.0).is_ok());
        assert!(Glitches::new(10_000.0, 0.0, 25.0).is_ok());
        assert!(Glitches::new(0.0, 25.0, 25.0).is_ok());
        assert!(Glitches::new(0.0, 0.0, 0.0).is_ok());
    }

    #[test]
    fn get_value() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Glitches::new(10_000.0, 25.0, 25.0).unwrap();

        let data: Vec<(usize, usize, Vec<u8>)> = (0..20)
            .map(|_| model.get_glitch(&mut rng).unwrap())
            .collect();

        assert_eq!(
            vec![
                (
                    4477,
                    4487,
                    vec![
                        65, 84, 84, 65, 84, 65, 71, 84, 65, 67, 71, 71, 84, 65, 84, 65, 71, 84, 71,
                        71, 84, 84, 65
                    ]
                ),
                (
                    23676,
                    23678,
                    vec![
                        71, 67, 67, 84, 65, 65, 71, 84, 71, 71, 67, 71, 67, 67, 67, 71, 84, 84, 71,
                        84
                    ]
                ),
                (
                    4735,
                    4762,
                    vec![
                        71, 65, 65, 84, 67, 67, 65, 67, 84, 84, 65, 84, 65, 84, 65, 65, 67, 65, 67,
                        65, 71, 71, 84, 65, 84, 65, 65, 84, 67
                    ]
                ),
                (
                    727,
                    737,
                    vec![
                        67, 71, 71, 67, 65, 84, 71, 67, 71, 67, 65, 71, 71, 67, 65, 84, 71, 67, 67,
                        84, 65, 84, 65, 84, 84, 67, 84, 65, 84, 71, 65
                    ]
                ),
                (
                    16090,
                    16100,
                    vec![
                        65, 71, 71, 65, 84, 84, 65, 84, 71, 71, 65, 65, 71, 65, 84, 71, 71, 84, 71,
                        67, 84, 67, 84, 65, 71, 65, 84, 65, 67, 71, 84
                    ]
                ),
                (
                    5935,
                    5992,
                    vec![
                        67, 67, 67, 71, 84, 65, 71, 67, 65, 67, 71, 65, 67, 67, 71, 71, 67, 84, 65,
                        84, 71
                    ]
                ),
                (
                    25,
                    29,
                    vec![
                        84, 84, 84, 84, 67, 84, 84, 71, 71, 65, 67, 65, 84, 65, 71, 84, 84, 84, 67,
                        71, 84, 67, 67, 65, 67
                    ]
                ),
                (
                    3221,
                    3266,
                    vec![
                        84, 65, 67, 65, 65, 71, 71, 65, 67, 71, 67, 84, 84, 71, 71, 71, 65, 65, 84,
                        65, 71, 71, 71, 67, 65, 71, 67, 71, 71, 65, 71, 84, 84, 65, 84, 67, 71, 84,
                        71, 84, 65, 67, 67, 84, 67, 67, 84, 65, 71, 67, 84, 84, 84, 84, 65, 71, 84
                    ]
                ),
                (
                    6743,
                    6786,
                    vec![
                        65, 67, 65, 71, 84, 71, 84, 65, 65, 67, 65, 84, 84, 71, 71, 71, 65, 67, 71,
                        67, 84
                    ]
                ),
                (7465, 7496, vec![67, 71, 67, 67, 71, 71]),
                (15218, 15232, vec![]),
                (726, 773, vec![67, 84, 65, 84, 65, 67, 67]),
                (5263, 5291, vec![67, 71, 84, 71]),
                (
                    2313,
                    2340,
                    vec![67, 71, 67, 71, 67, 71, 71, 65, 84, 67, 67, 67, 84, 67, 65]
                ),
                (
                    1517,
                    1521,
                    vec![84, 67, 71, 71, 71, 65, 65, 71, 67, 71, 67, 71, 65, 65, 67]
                ),
                (
                    8099,
                    8108,
                    vec![
                        67, 84, 84, 65, 84, 65, 67, 84, 65, 65, 84, 84, 67, 67, 65, 67, 71, 67, 65,
                        65, 84, 71, 84, 65, 67, 84, 67, 71, 67, 84, 84, 65, 67, 71, 65, 84, 84, 71,
                        67, 65, 65, 84, 84, 84, 71, 67, 65, 65, 65, 84
                    ]
                ),
                (
                    7455,
                    7506,
                    vec![65, 84, 67, 71, 84, 67, 67, 67, 84, 84, 67, 65, 84, 65]
                ),
                (
                    1452,
                    1470,
                    vec![71, 71, 65, 71, 67, 84, 84, 71, 65, 84, 67, 67, 84, 71, 65, 65, 84, 71]
                ),
                (
                    10384,
                    10391,
                    vec![
                        84, 67, 67, 71, 67, 71, 71, 67, 65, 84, 71, 71, 67, 84, 65, 65, 71, 84, 65,
                        67, 67, 65, 67, 67, 71, 84, 71, 71, 65
                    ]
                ),
                (
                    4776,
                    4853,
                    vec![
                        65, 71, 84, 67, 71, 84, 67, 84, 67, 84, 84, 67, 65, 84, 65, 67, 84, 71, 84
                    ]
                )
            ],
            data
        );

        assert!(Glitches::new(0.0, 25.0, 25.0)
            .unwrap()
            .get_glitch(&mut rng)
            .is_none());
    }
}
