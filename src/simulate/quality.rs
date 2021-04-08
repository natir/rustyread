//! Generate quality for simulate

/* standard use */

/* crate use */
use anyhow::Result;

/* local use */
use crate::model;
use crate::simulate::error;

type Quality = Vec<u8>;

/// Generate quality string
pub fn generate_quality(
    raw: &[u8],
    err: &[u8],
    model: &model::Quality,
    diffs: error::DiffPos,
    rng: &mut rand::rngs::StdRng,
) -> Result<(f64, Quality)> {
    let mut qual = Vec::with_capacity(err.len());
    let margin = (model.max_k() - 1) / 2;

    let (edit, cigar) = rebuild_cigar(raw, err, diffs);

    for i in 0..cigar.len() {
        if cigar[i] == b'D' {
            continue;
        }

        let (start, end) = if i < margin {
            (0, i + i + 1)
        } else if i >= cigar.len() - margin {
            (i - (cigar.len() - i - 1), cigar.len())
        } else {
            (i - margin, i + margin + 1)
        };

        qual.push(model.get_qscore(&cigar[start..end], rng)?);
    }

    Ok((1.0 - (edit as f64 / err.len() as f64), qual))
}

/// Build cigar string by align only changed position
pub fn rebuild_cigar(raw: &[u8], err: &[u8], diffs: error::DiffPos) -> (usize, Vec<u8>) {
    let mut edit: usize = 0;
    let mut cigar = Vec::with_capacity(err.len());
    let mut prev_e = 0;

    //crate::alignment::align(err, raw);
    for (r, e) in diffs.raw.chunks_exact(2).zip(diffs.err.chunks_exact(2)) {
        if e[0] > prev_e {
            cigar.extend((0..(e[0] - prev_e)).map(|_| b'='));
        }

        let (ed, c) = if e[1] >= err.len() {
            log::warn!("err.len {} e[1] {}", err.len(), e[1]);
            prev_e = err.len();
            crate::alignment::align(&err[e[0]..], &raw[r[0]..r[1]])
        } else {
            prev_e = e[1];
            crate::alignment::align(&err[e[0]..e[1]], &raw[r[0]..r[1]])
        };
        cigar.extend(&c[..]);

        edit += ed;
    }

    cigar.extend((0..(err.len() - prev_e)).map(|_| b'='));

    (edit, cigar)
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

    static MODEL: &[u8] = b"overall;3;1:0.000065,2:0.002683,3:0.011969,4:0.022982,5:0.03079,6:0.036889,7:0.042461,8:0.047495,9:0.051971,10:0.055876,11:0.058956,12:0.0606,13:0.060621,14:0.059446,15:0.056487,16:0.052028,17:0.047379,18:0.043445,19:0.040219,20:0.037767,21:0.035767,22:0.033816,23:0.031369,24:0.027969,25:0.023457,26:0.018417,27:0.008983,28:0.000093,
I;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
X;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
=;1;1:0.000076,2:0.00327,3:0.014147,4:0.034226,5:0.053392,6:0.066246,7:0.078339,8:0.078491,9:0.085032,10:0.082446,11:0.078111,12:0.066854,13:0.063356,14:0.053772,15:0.040006,16:0.036431,17:0.033541,18:0.027989,19:0.020992,20:0.019851,21:0.019319,22:0.015364,23:0.012017,24:0.008442,25:0.005476,26:0.002282,27:0.000532,
";

    #[test]
    fn quality() {
        init();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let qual = model::Quality::from_stream(MODEL).unwrap();

        let raw = b"TCAGCCACACATCCAGCCCCGTCTCCATACGCTTAATGGTGTAGCTAATGGCGGAAGTGGTTAAACCCAACTCTTCTGCGGCTTTACTGAAGCTGCCAAAACGCGCAGTCCATG";
        let err = b"TCAGCCACACTATCCAGCCCGTCTCCATACGCTTAATGGTGTGCTAATGGCGGAAGTGGTTAAACCCAGCTCTTCTGCGGCTTTGCTGAAACTGCCAAAAACGCAGTCCATG";

        let diffpos = error::DiffPos {
            raw: vec![7, 18, 41, 48, 67, 74, 80, 93, 93, 100, 100, 107],
            err: vec![7, 18, 41, 47, 66, 73, 79, 92, 92, 100, 100, 105],
        };

        let (identity, qual) = generate_quality(raw, err, &qual, diffpos, &mut rng).unwrap();

        assert_eq!(0.9196428571428571, identity);

        assert_eq!(err.len(), qual.len());
        assert_eq!(b",,-*%+/2'#5,*',($,%.2,,&*+0(*.+353.(-*&+9-+%'*+72,01-(**((-+7'+,'*(,.&+-3,++5%+/)**/.-&*/5&2,00)'03)')5/)'2'**-3".to_vec(), qual);
    }

    #[test]
    fn reconstruct_cigar() {
        let raw = b"TTTGTTCTGCCATCGGCCCTTACTGCGTGCCGGTGGTTAACCTCGAGGCGAACGTCGATCAACTGAACGTCAACATGGTCACCTGCGGCGGCCAGGCCACCATTCCACCATATT";
        let err = b"TTTGTTCTGGCCATCGGCCCTTACTGCGTGCCGGTGGTTAACCTCGAGGCGAACGTCGATCAACTGAACGTCACATGGTCACCTCGCGGCGGCCAGGCCACCATTCCACATATT";

        let diffpos = error::DiffPos {
            raw: vec![6, 13, 66, 73, 83, 90, 104, 111],
            err: vec![6, 14, 67, 73, 83, 91, 105, 111],
        };

        let (edit, cigar) = rebuild_cigar(raw, err, diffpos);
        let (t_e, t_c) = crate::alignment::align(err, raw);

        assert_eq!(edit, t_e);
        assert_eq!(cigar, t_c.to_vec());
    }
}
