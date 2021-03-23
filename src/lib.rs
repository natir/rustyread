/* crate use */

/* module declaration */
pub mod alignment;
pub mod cli;
pub mod error;
pub mod model;
pub mod references;

const NUCS: [u8; 4] = [b'A', b'C', b'T', b'G'];

pub fn random_base<R>(rng: &mut R) -> u8
where
    R: rand::Rng,
{
    NUCS[rng.gen_range(0..=3)]
}

pub fn random_base_diff<R>(nuc: u8, rng: &mut R) -> u8
where
    R: rand::Rng,
{
    loop {
        let idx = rng.gen_range(0..=3);

        if NUCS[idx] != nuc {
            return NUCS[idx];
        }
    }
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
}
