//! Function to align sequence

/* standard use */

/* crate use */

/* local use */

pub fn edit_distance(seq1: &[u8], seq2: &[u8]) -> u64 {
    bio::alignment::distance::levenshtein(seq1, seq2) as u64
}

pub fn identity(perfect: &[u8], error: &[u8]) -> f64 {
    edit_distance(perfect, error) as f64 / perfect.len() as f64
}
