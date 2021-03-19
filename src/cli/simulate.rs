//! All stuff relate to simulate subcommand

/// Store quantity as coverage of number of base
#[derive(Debug, PartialEq)]
pub struct Quantity {
    coverage: u64,
    base: Option<u64>,
}

impl std::str::FromStr for Quantity {
    type Err = crate::error::Cli;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.as_bytes() {
            &[rest @ .., b'x'] => match std::str::from_utf8(rest) {
                Ok(number) => match u64::from_str(number) {
                    Ok(n) => Ok(Quantity {
                        coverage: n,
                        base: None,
                    }),
                    Err(_) => Err(crate::error::Cli::CantParssQuantity),
                },
                Err(_) => Err(crate::error::Cli::CantParssQuantity),
            },
            &[rest @ .., b'G'] => match std::str::from_utf8(rest) {
                Ok(number) => match u64::from_str(number) {
                    Ok(n) => Ok(Quantity {
                        coverage: 0,
                        base: Some(n * 1_000_000_000),
                    }),
                    Err(_) => Err(crate::error::Cli::CantParssQuantity),
                },
                Err(_) => Err(crate::error::Cli::CantParssQuantity),
            },
            &[rest @ .., b'M'] => match std::str::from_utf8(rest) {
                Ok(number) => match u64::from_str(number) {
                    Ok(n) => Ok(Quantity {
                        coverage: 0,
                        base: Some(n * 1_000_000),
                    }),
                    Err(_) => Err(crate::error::Cli::CantParssQuantity),
                },
                Err(_) => Err(crate::error::Cli::CantParssQuantity),
            },
            &[rest @ .., b'K'] => match std::str::from_utf8(rest) {
                Ok(number) => match u64::from_str(number) {
                    Ok(n) => Ok(Quantity {
                        coverage: 0,
                        base: Some(n * 1_000),
                    }),
                    Err(_) => Err(crate::error::Cli::CantParssQuantity),
                },
                Err(_) => Err(crate::error::Cli::CantParssQuantity),
            },
            rest => match std::str::from_utf8(rest) {
                Ok(number) => match u64::from_str(number) {
                    Ok(n) => Ok(Quantity {
                        coverage: 0,
                        base: Some(n),
                    }),
                    Err(_) => Err(crate::error::Cli::CantParssQuantity),
                },
                Err(_) => Err(crate::error::Cli::CantParssQuantity),
            },
        }
    }
}

impl Quantity {
    /// Convert Quantity in a number of base, if base is set return just number of base else return genome_length times coverage
    pub fn number_of_base(&self, genome_length: u64) -> u64 {
        match self.base {
            Some(n) => n,
            None => genome_length * self.coverage,
        }
    }
}

/// Struct use to parse simulate subcommand argument
#[derive(clap::Clap, Debug)]
#[clap(about = "Generate fake long read")]
pub struct Command {
    /// Path to reference sequence in fasta format
    #[clap(
        long = "reference",
        about = "Reference fasta (can be gzipped, bzip2ped, xzped)"
    )]
    pub reference_path: String,

    /// Quantity of base badread have to generate
    #[clap(
        long = "quantity",
        about = "Either an absolute value (e.g. 250M) or a relative depth (e.g. 25x)"
    )]
    pub quantity: Quantity,
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::str::FromStr;

    #[test]
    fn parse_quantity() {
        assert_eq!(
            Quantity::from_str("50").unwrap(),
            Quantity {
                coverage: 0,
                base: Some(50)
            }
        );
        assert_eq!(
            Quantity::from_str("50K").unwrap(),
            Quantity {
                coverage: 0,
                base: Some(50_000)
            }
        );
        assert_eq!(
            Quantity::from_str("50M").unwrap(),
            Quantity {
                coverage: 0,
                base: Some(50_000_000)
            }
        );
        assert_eq!(
            Quantity::from_str("50G").unwrap(),
            Quantity {
                coverage: 0,
                base: Some(50_000_000_000)
            }
        );
        assert_eq!(
            Quantity::from_str("50x").unwrap(),
            Quantity {
                coverage: 50,
                base: None
            }
        );

        assert!(Quantity::from_str("").is_err());
        assert!(Quantity::from_str("50z").is_err());
        assert!(Quantity::from_str("bépo50").is_err());
        assert!(Quantity::from_str("bépo50K").is_err());
        assert!(Quantity::from_str("bépo50M").is_err());
        assert!(Quantity::from_str("bépo50G").is_err());
        assert!(Quantity::from_str("bépo50x").is_err());
    }

    #[test]
    fn quantity2number_of_base() {
        assert_eq!(Quantity::from_str("50").unwrap().number_of_base(50), 50);
        assert_eq!(
            Quantity::from_str("50K").unwrap().number_of_base(50),
            50_000
        );
        assert_eq!(
            Quantity::from_str("50M").unwrap().number_of_base(50),
            50_000_000
        );
        assert_eq!(
            Quantity::from_str("50G").unwrap().number_of_base(50),
            50_000_000_000
        );
        assert_eq!(
            Quantity::from_str("50x").unwrap().number_of_base(50),
	    2500
        );
    }
}
