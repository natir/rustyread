/* crate use */

/* module declaration */
pub mod alignment;
pub mod cli;
pub mod error;
pub mod model;
pub mod references;
pub mod simulate;

/* constant definition */
const CHIMERA_START_ADAPTER_CHANCE: f64 = 0.25;
const CHIMERA_END_ADAPTER_CHANCE: f64 = 0.25;

const NUCS: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Get a random base
pub fn random_base<RNG>(rng: &mut RNG) -> u8
where
    RNG: rand::Rng,
{
    NUCS[rng.gen_range(0..=3)]
}

/// Get a random base diffrent than nuc
pub fn random_base_diff<RNG>(nuc: u8, rng: &mut RNG) -> u8
where
    RNG: rand::Rng,
{
    loop {
        let idx = rng.gen_range(0..=3);

        if NUCS[idx] != nuc {
            return NUCS[idx];
        }
    }
}

/// Get random sequences
pub fn random_seq<RNG>(length: usize, rng: &mut RNG) -> Vec<u8>
where
    RNG: rand::Rng,
{
    (0..length).map(|_| random_base(rng)).collect()
}

#[cfg(test)]
mod t {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn random_base_() {
        let truth: Vec<u8> = [
            b'T', b'T', b'A', b'G', b'A', b'T', b'T', b'A', b'T', b'A', b'G', b'T', b'A', b'C',
            b'G', b'G', b'T', b'A', b'T', b'A',
        ]
        .to_vec();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let data: Vec<u8> = (0..20).map(|_| random_base(&mut rng)).collect();

        assert_eq!(truth, data);
    }

    #[test]
    fn random_base_diff_() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let mut truth: Vec<u8> = [
            b'T', b'T', b'G', b'T', b'T', b'T', b'G', b'T', b'C', b'G', b'G', b'T', b'T', b'G',
            b'T', b'G', b'G', b'T', b'T', b'C',
        ]
        .to_vec();
        let mut data: Vec<u8> = (0..20).map(|_| random_base_diff(b'A', &mut rng)).collect();

        assert_eq!(truth, data);

        truth = [
            b'T', b'A', b'T', b'G', b'T', b'A', b'G', b'T', b'A', b'A', b'G', b'T', b'G', b'G',
            b'G', b'G', b'T', b'T', b'G', b'T',
        ]
        .to_vec();
        data = (0..20).map(|_| random_base_diff(b'C', &mut rng)).collect();
        assert_eq!(truth, data);

        truth = [
            b'A', b'G', b'A', b'G', b'G', b'A', b'A', b'C', b'C', b'A', b'C', b'A', b'A', b'A',
            b'A', b'C', b'A', b'C', b'A', b'G',
        ]
        .to_vec();
        data = (0..20).map(|_| random_base_diff(b'T', &mut rng)).collect();
        assert_eq!(truth, data);

        truth = [
            b'T', b'A', b'T', b'A', b'A', b'T', b'C', b'C', b'A', b'C', b'C', b'A', b'T', b'C',
            b'C', b'A', b'C', b'A', b'T', b'C',
        ]
        .to_vec();
        data = (0..20).map(|_| random_base_diff(b'G', &mut rng)).collect();
        assert_eq!(truth, data);
    }

    #[test]
    fn random_seq_() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        assert_eq!(b"TTAGATTATAGTACGGTATA".to_vec(), random_seq(20, &mut rng));
    }
}
