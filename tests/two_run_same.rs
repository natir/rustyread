/* std use */
use std::io::Read;
use std::process::{Command, Stdio};

/* crate use */
use bio::io::fastq::Reader;

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    fn diff_unorder(first: &str, second: &str) {
        let first_file = Reader::new(std::io::BufReader::new(
            std::fs::File::open(first).expect(&format!("Impossible to open {}", first)),
        ));

        let mut first_id = std::collections::HashSet::new();
        let mut first_desc = std::collections::HashSet::new();
        let mut first_seq = std::collections::HashSet::new();
        let mut first_qual = std::collections::HashSet::new();

        for res in first_file.records() {
            let line = res.unwrap();
            first_id.insert(format!("{}", line.id()));
            first_desc.insert(format!("{}", line.desc().unwrap()));
            first_seq.insert(format!("{}", std::str::from_utf8(line.seq()).unwrap()));
            first_qual.insert(format!("{}", std::str::from_utf8(line.qual()).unwrap()));
        }

        let second_file = Reader::new(std::io::BufReader::new(
            std::fs::File::open(second).expect(&format!("Impossible to open {}", second)),
        ));

        let mut second_id = std::collections::HashSet::new();
        let mut second_desc = std::collections::HashSet::new();
        let mut second_seq = std::collections::HashSet::new();
        let mut second_qual = std::collections::HashSet::new();

        for res in second_file.records() {
            let line = res.unwrap();
            second_id.insert(format!("{}", line.id()));
            second_desc.insert(format!("{}", line.desc().unwrap()));
            second_seq.insert(format!("{}", std::str::from_utf8(line.seq()).unwrap()));
            second_qual.insert(format!("{}", std::str::from_utf8(line.qual()).unwrap()));
        }

        assert_eq!(first_id, second_id);
        assert_eq!(first_desc, second_desc);
        assert_eq!(first_seq, second_seq);
        assert_eq!(first_qual, second_qual);
    }

    fn run_process(path: &str, params: &[&str]) {
        let mut child = Command::new(path)
            .args(params)
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Couldn't create subprocess");

        if !child.wait().expect("Error durring run").success() {
            let mut stdout = String::new();
            let mut stderr = String::new();

            child.stdout.unwrap().read_to_string(&mut stdout).unwrap();
            child.stderr.unwrap().read_to_string(&mut stderr).unwrap();

            println!("stdout: {}", stdout);
            println!("stderr: {}", stderr);
            panic!();
        }
    }

    #[test]
    fn two_run_same_result() {
        init();

        run_process(
            "./target/debug/rustyread",
            &[
                "-vvvvv",
                "--threads",
                "2",
                "simulate",
                "--reference",
                "./tests/data/ref_100000.fasta",
                "--quantity",
                "5x",
                "--seed",
                "42",
                "--error_model",
                "random",
                "--qscore_model",
                "random",
                "--output",
                "./tests/run1.fastq",
            ],
        );

        run_process(
            "./target/debug/rustyread",
            &[
                "-vvvvv",
                "--threads",
                "2",
                "simulate",
                "--reference",
                "./tests/data/ref_100000.fasta",
                "--quantity",
                "5x",
                "--seed",
                "42",
                "--error_model",
                "random",
                "--qscore_model",
                "random",
                "--output",
                "./tests/run2.fastq",
            ],
        );

        diff_unorder("./tests/run1.fastq", "./tests/run2.fastq")
    }
}
