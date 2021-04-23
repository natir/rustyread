//! Generate quality for simulate

/* standard use */

/* crate use */
use anyhow::Result;

/* local use */
use crate::model;

type Quality = Vec<u8>;

/// Generate quality string
pub fn generate_quality(
    cigar: &[u8],
    model: &model::Quality,
    rng: &mut rand::rngs::StdRng,
) -> Result<Quality> {
    let mut qual = Vec::with_capacity(cigar.len());
    let margin = (model.max_k() - 1) / 2;

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

    Ok(qual)
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

        let (_, cigar) = crate::alignment::align(err, raw);

        let qual = generate_quality(&cigar, &qual, &mut rng).unwrap();

        assert_eq!(err.len(), qual.len());
        assert_eq!(b",,-*%+/2'#5,*',($,%.2,,&*+0(*.+353.(-*&+9-+%'*+72,01-(**((-+7'+,'*(,.&+-3,++5%+/)**/.-&*/5&2,00)'03)')5/)'2'**-3".to_vec(), qual);
    }
}
