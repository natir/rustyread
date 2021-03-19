//! A collections of sequence to store reference sequence

/// A collections of sequence
pub type References = rustc_hash::FxHashMap<Box<str>, Box<[u8]>>;

pub trait AbsReferences {
    fn from_stream<R>(input: R) -> Self
    where
        R: std::io::BufRead;
}

impl<'a> AbsReferences for References {
    fn from_stream<R>(input: R) -> Self
    where
        R: std::io::BufRead,
    {
        let mut me = References::default();
        let mut records = bio::io::fasta::Reader::new(input).records();

        while let Some(Ok(record)) = records.next() {
            me.insert(record.id().into(), record.seq().into());
        }

        me
    }
}

#[cfg(test)]
mod t {
    use super::*;

    static FASTA: &'static [u8] = b">random_seq_0
TCCTAACGTG
>random_seq_1
TCACGATTAC
>random_seq_2
CCTATCCGAT
>random_seq_3
TGCAAGATCA
>random_seq_4
TAGCCGTGGT
>random_seq_5
CGCTTTGTGA
>random_seq_6
CACATGGGCG
>random_seq_7
ATCTAATGCG
>random_seq_8
CGGAACTCAG
>random_seq_9
TCCCGCTGTC
";

    #[test]
    fn read_reference() {
        let me = References::from_stream(std::io::Cursor::new(FASTA));

        assert_eq!(
            me.get("random_seq_0"),
            Some(&b"TCCTAACGTG".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_1"),
            Some(&b"TCACGATTAC".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_2"),
            Some(&b"CCTATCCGAT".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_3"),
            Some(&b"TGCAAGATCA".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_4"),
            Some(&b"TAGCCGTGGT".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_5"),
            Some(&b"CGCTTTGTGA".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_6"),
            Some(&b"CACATGGGCG".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_7"),
            Some(&b"ATCTAATGCG".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_8"),
            Some(&b"CGGAACTCAG".to_vec().into_boxed_slice())
        );
        assert_eq!(
            me.get("random_seq_9"),
            Some(&b"TCCCGCTGTC".to_vec().into_boxed_slice())
        );

        assert_eq!(me.get("bépo"), None);
    }
}
