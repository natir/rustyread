mod common;

#[test]
fn two_run_same_result() {
    common::init();

    common::run_process(
        "./target/debug/rustyread",
        &[
            "-vvvvv",
            "--threads",
            "2",
            "simulate",
            "--reference",
            "./tests/data/ref_100000.fasta",
            "--quantity",
            "50000",
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

    common::run_process(
        "./target/debug/rustyread",
        &[
            "-vvvvv",
            "--threads",
            "2",
            "simulate",
            "--reference",
            "./tests/data/ref_100000.fasta",
            "--quantity",
            "50000",
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

    common::diff_unorder_fasta("./tests/run1.fastq", "./tests/run2.fastq")
}
