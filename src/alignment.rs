//! Function to align sequence

/* standard use */

/* crate use */

/* local use */

pub fn edit_distance(seq1: &[u8], seq2: &[u8]) -> u64 {
    bio::alignment::distance::levenshtein(seq1, seq2) as u64
}

pub fn align(query: &[u8], target: &[u8]) -> (usize, Box<[u8]>) {
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
            bio::alignment::AlignmentOperation::Xclip(len) => cigar.extend(vec![b'X'; len]),
            bio::alignment::AlignmentOperation::Yclip(len) => cigar.extend(vec![b'X'; len]),
            bio::alignment::AlignmentOperation::Subst => cigar.push(b'X'),
            bio::alignment::AlignmentOperation::Ins => cigar.push(b'I'),
            bio::alignment::AlignmentOperation::Del => cigar.push(b'D'),
            bio::alignment::AlignmentOperation::Match => cigar.push(b'='),
        }
    }

    (
        cigar.iter().filter(|x| **x != b'=').count(),
        cigar.into_boxed_slice(),
    )
}