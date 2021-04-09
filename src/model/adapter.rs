//! Model to get sequence adapter

/* standard use */

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;

/* local use */

/// Struct to get sequence adapter
pub struct Adapter {
    start: Vec<u8>,
    end: Vec<u8>,
    start_rate: f64,
    end_rate: f64,
    start_dist: Option<rand_distr::Beta<f64>>,
    end_dist: Option<rand_distr::Beta<f64>>,
}

impl Adapter {
    /// Create model from parameter
    pub fn new(
        start: Vec<u8>,
        end: Vec<u8>,
        mut start_rate: f64,
        mut start_amount: f64,
        mut end_rate: f64,
        mut end_amount: f64,
    ) -> Result<Self> {
        let start_dist = if start.is_empty() || start_rate == 0.0 || start_amount == 0.0 {
            None
        } else {
            start_amount /= 100.0;
            Some(rand_distr::Beta::new(
                2.0 * start_amount,
                2.0 - 2.0 * (start_amount / 100.0),
            )?)
        };

        let end_dist = if end.is_empty() || end_rate == 0.0 || end_amount == 0.0 {
            None
        } else {
            end_amount /= 100.0;
            Some(rand_distr::Beta::new(
                2.0 * end_amount,
                2.0 - 2.0 * end_amount,
            )?)
        };

        start_rate /= 100.0;
        end_rate /= 100.0;

        Ok(Self {
            start,
            end,
            start_rate,
            end_rate,
            start_dist,
            end_dist,
        })
    }

    pub fn get_start<R>(&self, rng: &mut R) -> Vec<u8>
    where
        R: rand::Rng,
    {
        if rng.gen_bool(self.start_rate) {
            if let Some(dist) = self.start_dist {
                self.start[0..(self.start.len() as f64 * dist.sample(rng)).round() as usize]
                    .to_vec()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    }

    pub fn get_end<R>(&self, rng: &mut R) -> Vec<u8>
    where
        R: rand::Rng,
    {
        if rng.gen_bool(self.end_rate) {
            if let Some(dist) = self.end_dist {
                self.end[0..(self.end.len() as f64 * dist.sample(rng)).round() as usize].to_vec()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    }
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn create_model() {
        assert!(Adapter::new(b"AA".to_vec(), b"AA".to_vec(), 0.0, 0.0, 0.0, 0.0).is_ok());
        assert!(Adapter::new(b"".to_vec(), b"".to_vec(), 90.0, 60.0, 90.0, 60.0).is_ok());
    }

    #[test]
    fn generate() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let model = Adapter::new(
            b"CGATCACAAA".to_vec(),
            b"CAACATATGT".to_vec(),
            90.0,
            60.0,
            50.0,
            20.0,
        )
        .unwrap();

        let starts: Vec<Vec<u8>> = (0..20).map(|_| model.get_start(&mut rng)).collect();
        let ends: Vec<Vec<u8>> = (0..20).map(|_| model.get_end(&mut rng)).collect();

        assert_eq!(
            vec![
                b"CGAT".to_vec(),
                b"".to_vec(),
                b"CGATCAC".to_vec(),
                b"CGATCACAA".to_vec(),
                b"C".to_vec(),
                b"".to_vec(),
                b"CGATC".to_vec(),
                b"CGAT".to_vec(),
                b"CGA".to_vec(),
                b"CGA".to_vec(),
                b"CGATCACA".to_vec(),
                b"CGAT".to_vec(),
                b"CGA".to_vec(),
                b"CGA".to_vec(),
                b"CGA".to_vec(),
                b"".to_vec(),
                b"CGAT".to_vec(),
                b"CGAT".to_vec(),
                b"CGA".to_vec(),
                b"CGATC".to_vec(),
            ],
            starts
        );

        assert_eq!(vec![
            b"".to_vec(),
            b"CAA".to_vec(),
            b"CAACAT".to_vec(),
            b"".to_vec(),
            b"".to_vec(),
            b"CA".to_vec(),
            b"".to_vec(),
            b"".to_vec(),
            b"CAA".to_vec(),
            b"CA".to_vec(),
            b"CAACAT".to_vec(),
            b"".to_vec(),
            b"".to_vec(),
            b"".to_vec(),
            b"CAAC".to_vec(),
            b"".to_vec(),
            b"CAACAT".to_vec(),
            b"".to_vec(),
            b"".to_vec(),
            b"CAACATAT".to_vec(),
        ], ends);
    }
}
