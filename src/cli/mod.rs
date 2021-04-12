//! All stuff relate to command line

/* module declaration */
pub mod simulate;
// pub(crate) mod error_model;
// pub(crate) mod qscore_model;
// pub(crate) mod plot;

#[derive(clap::Clap, Debug)]
#[clap(
    name = "rustyread",
    version = "0.2 Electabuzz",
    author = "Pierre Marijon <pierre.marijon@hhu.de>",
    about = "A long read simulator based on badread idea"
)]
pub struct Command {
    /// Subcommand call
    #[clap(subcommand)]
    pub subcmd: SubCommand,

    /// Number of thread rustyread can use
    #[clap(
        short = 't',
        long = "threads",
        about = "Number of thread use by rustyread, 0 use all avaible core, default value 0"
    )]
    pub threads: Option<usize>,

    /// Verbosity level
    #[clap(
        short = 'v',
        long = "verbosity",
        parse(from_occurrences),
        about = "verbosity level also control by environment variable RUSTYREAD_LOG if flag is set RUSTYREAD_LOG value is ignored"
    )]
    pub verbosity: i8,
}

#[derive(clap::Clap, Debug)]
pub enum SubCommand {
    Simulate(simulate::Command),
    // ErrorModel(error_model::Command),
    // QScoreModel(qscore_model::Command),
    // Plot(plot::Command),
}

/// Convert verbosity level (number of v) is log::Level
pub fn i82level(level: i8) -> Option<log::Level> {
    match level {
        std::i8::MIN..=0 => None,
        1 => Some(log::Level::Error),
        2 => Some(log::Level::Warn),
        3 => Some(log::Level::Info),
        4 => Some(log::Level::Debug),
        5..=std::i8::MAX => Some(log::Level::Trace),
    }
}

/// set number of global rayon thread pool
pub fn set_nb_threads(nb_threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(nb_threads)
        .build_global()
        .unwrap();
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn loglevel() {
        assert_eq!(i82level(i8::MIN), None);
        assert_eq!(i82level(-3), None);
        assert_eq!(i82level(1), Some(log::Level::Error));
        assert_eq!(i82level(2), Some(log::Level::Warn));
        assert_eq!(i82level(3), Some(log::Level::Info));
        assert_eq!(i82level(4), Some(log::Level::Debug));
        assert_eq!(i82level(5), Some(log::Level::Trace));
        assert_eq!(i82level(i8::MAX), Some(log::Level::Trace));
    }

    #[test]
    fn change_number_of_thread() {
        set_nb_threads(16);
        assert_eq!(rayon::current_num_threads(), 16);
    }
}
