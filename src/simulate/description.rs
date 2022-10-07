//! Manage read description

/* standard use */

/* crate use */

/* local use */

/// An enum to represent type of read
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadType {
    Real,
    Junk,
    Random,
}

/// Store information about origin of read
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Origin {
    pub ref_id: String,
    pub strand: char,
    pub start: usize,
    pub end: usize,
    pub read_type: ReadType,
}

impl Origin {
    pub fn reference(ref_id: String, strand: char, start: usize, end: usize) -> Self {
        Origin {
            ref_id,
            strand,
            start,
            end,
            read_type: ReadType::Real,
        }
    }

    pub fn junk(length: usize) -> Self {
        Origin {
            ref_id: "".to_string(),
            strand: '*',
            start: 0,
            end: length,
            read_type: ReadType::Junk,
        }
    }

    pub fn random(length: usize) -> Self {
        Origin {
            ref_id: "".to_string(),
            strand: '*',
            start: 0,
            end: length,
            read_type: ReadType::Random,
        }
    }
}

impl std::fmt::Display for Origin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.read_type {
            ReadType::Junk => write!(f, "junk_seq"),
            ReadType::Random => write!(f, "random_seq"),
            ReadType::Real => write!(
                f,
                "{},{}strand,{}-{}",
                self.ref_id, self.strand, self.start, self.end
            ),
        }
    }
}

/// Store information about read
#[derive(Debug, Clone, PartialEq)]
pub struct Description {
    pub origin: Origin,
    pub chimera: Option<Origin>,
    pub length: usize,
    pub identity: f64,
}

impl Description {
    pub fn new(origin: Origin, chimera: Option<Origin>, length: usize, identity: f64) -> Self {
        Description {
            origin,
            chimera,
            length,
            identity,
        }
    }
}

impl std::fmt::Display for Description {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if let Some(chimera) = &self.chimera {
            write!(f, "{} chimera {} ", self.origin, chimera)?;
        } else {
            write!(f, "{} ", self.origin)?;
        }

        write!(
            f,
            "length={} error-free_length={} read_identity={}%",
            self.origin.end - self.origin.start,
            self.length,
            self.identity
        )
    }
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn origin() {
        let mut test = Origin::reference("bépo".to_string(), '+', 100, 400);

        assert_eq!("bépo,+strand,100-400", format!("{}", test));

        test.read_type = ReadType::Junk;

        assert_eq!("junk_seq", format!("{}", test));

        test.read_type = ReadType::Random;

        assert_eq!("random_seq", format!("{}", test));

        assert_eq!("junk_seq", format!("{}", Origin::junk(100)));

        assert_eq!("random_seq", format!("{}", Origin::random(100)));
    }

    #[test]
    fn description() {
        let ori = Origin::reference("bépo".to_string(), '+', 100, 400);

        let mut des = Description::new(ori.clone(), None, 301, 99.99);

        assert_eq!(
            "bépo,+strand,100-400 length=300 error-free_length=301 read_identity=99.99%",
            format!("{}", des)
        );

        des.chimera = Some(ori);

        assert_eq!("bépo,+strand,100-400 chimera bépo,+strand,100-400 length=300 error-free_length=301 read_identity=99.99%", format!("{}", des));

        des.origin.read_type = ReadType::Junk;

        assert_eq!("junk_seq chimera bépo,+strand,100-400 length=300 error-free_length=301 read_identity=99.99%", format!("{}", des));

        des.chimera = None;

        assert_eq!(
            "junk_seq length=300 error-free_length=301 read_identity=99.99%",
            format!("{}", des)
        );

        des.origin.read_type = ReadType::Random;

        assert_eq!(
            "random_seq length=300 error-free_length=301 read_identity=99.99%",
            format!("{}", des)
        );
    }
}
