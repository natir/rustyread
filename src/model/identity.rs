//! Model to get read identity

/* standard use */

/* crate use */
use anyhow::Result;
use rand::distributions::Distribution;

/* local use */
use crate::error::Model;

/// Struct to generate length of fragment
#[derive(Debug)]
pub struct Identity {
    mean: f64,
    max: f64,
    dist: Option<rand_distr::Beta<f64>>,
}

impl Identity {
    /// Create model from parameter
    pub fn new(mut mean: f64, mut max: f64, mut stdev: f64) -> Result<Identity> {
        if mean <= 0.0 || stdev <= 0.0 {
            anyhow::bail!(Model::IdentityParamMustBeUpperThan0);
        }

        mean /= 100.0;
        max /= 100.0;
        stdev /= 100.0;

        let dist = if (mean - max).abs() > f64::EPSILON {
            let beta_a = (((1.0 - (mean / max)) / ((stdev / max).powf(2.0))) - (max / mean))
                * ((mean / max).powf(2.0));
            let beta_b = beta_a * ((max / mean) - 1.0);

            Some(rand_distr::Beta::new(beta_a, beta_b)?)
        } else {
            None
        };

        Ok(Self { mean, max, dist })
    }

    /// Get identity from model
    pub fn get_identity<RNG>(&self, rng: &mut RNG) -> f64
    where
        RNG: rand::Rng,
    {
        if let Some(dist) = self.dist {
            self.max * dist.sample(rng)
        } else {
            self.mean
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
    fn create_dist() {
        init();

        assert!(Identity::new(0.0, 0.0, 0.0).is_err());
        assert!(Identity::new(0.0, 1.0, 0.0).is_err());
        assert!(Identity::new(-10.0, 0.0, -0.8).is_err());
        assert!(Identity::new(-10.0, 0.0, 0.8).is_err());
        assert!(Identity::new(10.0, 0.0, -0.8).is_err());

        assert!(Identity::new(1.0, 1.0, 1.0).is_ok());
    }

    #[test]
    fn get_value() {
        init();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let mut dist = Identity::new(85.0, 85.0, 5.0).unwrap();

        assert!((0.85 - dist.get_identity(&mut rng)).abs() < f64::EPSILON);

        dist = Identity::new(85.0, 95.0, 5.0).unwrap();
        let data: Vec<f64> = (0..100).map(|_| dist.get_identity(&mut rng)).collect();

        assert!(
            ((data.iter().sum::<f64>() / data.len() as f64) - 0.8489702298458458).abs()
                < f64::EPSILON
        );
        assert_eq!(
            data,
            vec![
                0.8458900869695859,
                0.8266538854981239,
                0.9236269584128707,
                0.8034865877150157,
                0.9023903395427547,
                0.8656678044869676,
                0.8464437950555987,
                0.9327644014072227,
                0.9189364199899893,
                0.7647425652996294,
                0.8490884148493464,
                0.8710502522772242,
                0.785919024034962,
                0.8641472012042868,
                0.8610523699746363,
                0.8241679295609741,
                0.8375776316664659,
                0.9148647985266006,
                0.8569910245307752,
                0.8986572917056497,
                0.8512969128721771,
                0.7620161191597191,
                0.7943651602000301,
                0.8379412623444956,
                0.8666881255552256,
                0.8876123818467793,
                0.8270479168255866,
                0.855173685107204,
                0.905951658892641,
                0.8854418993366224,
                0.8168036832036264,
                0.8532959404680499,
                0.7517399836760904,
                0.8521340032515394,
                0.702418906086419,
                0.8518751488823618,
                0.8751966637130258,
                0.8667176493376229,
                0.8165377637852074,
                0.9120448459919036,
                0.8097332846232775,
                0.9103369460151146,
                0.8382941343361601,
                0.7958910184445049,
                0.8966955191574987,
                0.7486294996900733,
                0.9003303734288046,
                0.8823350902699747,
                0.7658694316773411,
                0.8631809573545917,
                0.8355254661126231,
                0.9006486696475493,
                0.8539656186330764,
                0.843608901838129,
                0.9048342559019861,
                0.895771639500239,
                0.8915642215659686,
                0.909214662051522,
                0.7924956254488512,
                0.8652285631020808,
                0.7714131628716486,
                0.9167083604775889,
                0.8486951888053259,
                0.8826009315975513,
                0.8758660085107146,
                0.8962395016444161,
                0.8558438393764741,
                0.8416411315009604,
                0.8627243745111065,
                0.8513177703149141,
                0.8305559162542443,
                0.8208532280198885,
                0.847292343720847,
                0.8893556893040299,
                0.8582042658094021,
                0.8174828761929837,
                0.8751514428106903,
                0.886755890812422,
                0.7602763904321761,
                0.8655242896932386,
                0.8275229788135905,
                0.9127426034333019,
                0.8456011346179603,
                0.892100004954406,
                0.8324030684955885,
                0.8037140594351367,
                0.7750889163453305,
                0.9170595098908257,
                0.8144668273044616,
                0.8903911403485442,
                0.8532951362787395,
                0.7958103985550533,
                0.8638181439336684,
                0.8387222639155002,
                0.8225538335662824,
                0.800535759072639,
                0.9075943692038912,
                0.8179780300228882,
                0.809036484166593,
                0.7515133475251731,
            ]
        )
    }
}
