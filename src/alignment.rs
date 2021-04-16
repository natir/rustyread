//! Function to align sequence

/* standard use */

/* crate use */

/* local use */

pub fn edit_distance(seq1: &[u8], seq2: &[u8]) -> u64 {
    bio::alignment::distance::levenshtein(seq1, seq2) as u64
}

pub fn align(query: &[u8], target: &[u8]) -> (u64, Box<[u8]>) {
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(
        target.len(),
        query.len(),
        0,
        0,
        |a: u8, b: u8| if a == b { 1i32 } else { 0i32 },
    );

    let alignment = aligner.global(query, target);

    let mut cigar = Vec::with_capacity(query.len().max(target.len()));

    for op in alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Xclip(len) => cigar.extend(vec![b'I'; len]),
            bio::alignment::AlignmentOperation::Yclip(len) => cigar.extend(vec![b'I'; len]),
            bio::alignment::AlignmentOperation::Subst => cigar.push(b'X'),
            bio::alignment::AlignmentOperation::Ins => cigar.push(b'I'),
            bio::alignment::AlignmentOperation::Del => cigar.push(b'D'),
            bio::alignment::AlignmentOperation::Match => cigar.push(b'='),
        }
    }

    (
        cigar.iter().filter(|x| **x != b'=').count() as u64,
        cigar.into_boxed_slice(),
    )
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn compute_edit() {
        assert_eq!(edit_distance(b"GCCTACGCAA", b"GCCTACCCAA"), 1);
        assert_eq!(edit_distance(b"GTACTGTCGG", b"GTACTCGG"), 2);
        assert_eq!(edit_distance(b"TGTGCAAGCG", b"TGTGAAAACG"), 2);
        assert_eq!(edit_distance(b"CTTGTACTAT", b"CTTGAACTAT"), 1);
        assert_eq!(edit_distance(b"TATCCCCTAA", b"TATGGGGTAA"), 4);
    }

    #[test]
    fn align_edit() {
        assert_eq!(
            align(b"GCCTACGCAA", b"GCCTACCCAA"),
            (1, b"======X===".to_vec().into_boxed_slice())
        );

        assert_eq!(
            align(b"GTACTGTCGG", b"GTACTCGG"),
            (2, b"====II====".to_vec().into_boxed_slice())
        );

        assert_eq!(
            align(b"TGTGCAAGCG", b"TGTGAAAACG"),
            (2, b"====X==X==".to_vec().into_boxed_slice())
        );

        assert_eq!(
            align(b"CTTGTACTAT", b"CTTGACTAT"),
            (1, b"====I=====".to_vec().into_boxed_slice())
        );

        assert_eq!(
            align(b"TATTAA", b"TATCCCCTAA"),
            (4, b"===DDDD===".to_vec().into_boxed_slice())
        );

        assert_eq!(
            align(b"CAAG", b"CAGCAAGGCC"),
            (6, b"DDD===D=DD".to_vec().into_boxed_slice())
        );
    }
}
