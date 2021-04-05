//! Manage read description

/* standard use */

/* crate use */

/* local use */

/// Store information about origin of read
#[derive(Debug, Clone)]
pub struct Origin {
    pub ref_id: String,
    pub strand: char,
    pub start: usize,
    pub end: usize,
    pub junk: bool,
    pub random: bool,
}

impl Origin {
    pub fn new(
        ref_id: String,
        strand: char,
        start: usize,
        end: usize,
        junk: bool,
        random: bool,
    ) -> Self {
        Origin {
            ref_id,
            strand,
            start,
            end,
            junk,
            random,
        }
    }
}

impl std::fmt::Display for Origin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.junk {
            write!(f, "junk_seq")
        } else if self.random {
            write!(f, "random_seq")
        } else {
            write!(
                f,
                "{},{}strand,{}-{}",
                self.ref_id, self.strand, self.start, self.end
            )
        }
    }
}

/// Store information about read
#[derive(Debug, Clone)]
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
            self.length,
            self.origin.end - self.origin.start,
            self.identity
        )
    }
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn origin() {
        let mut test = Origin::new("bépo".to_string(), '+', 100, 400, false, false);

        assert_eq!("bépo,+strand,100-400", format!("{}", test));

        test.junk = true;

        assert_eq!("junk_seq", format!("{}", test));

        test.junk = false;
        test.random = true;

        assert_eq!("random_seq", format!("{}", test));
    }

    #[test]
    fn description() {
        let ori = Origin::new("bépo".to_string(), '+', 100, 400, false, false);

        let mut des = Description::new(ori.clone(), None, 301, 99.99);

        assert_eq!(
            "bépo,+strand,100-400 length=301 error-free_length=300 read_identity=99.99%",
            format!("{}", des)
        );

        des.chimera = Some(ori);

        assert_eq!("bépo,+strand,100-400 chimera bépo,+strand,100-400 length=301 error-free_length=300 read_identity=99.99%", format!("{}", des));

        des.origin.junk = true;

        assert_eq!("junk_seq chimera bépo,+strand,100-400 length=301 error-free_length=300 read_identity=99.99%", format!("{}", des));

        des.chimera = None;

        assert_eq!(
            "junk_seq length=301 error-free_length=300 read_identity=99.99%",
            format!("{}", des)
        );

        des.origin.junk = false;
        des.origin.random = true;

        assert_eq!(
            "random_seq length=301 error-free_length=300 read_identity=99.99%",
            format!("{}", des)
        );
    }
}
