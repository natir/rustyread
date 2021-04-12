![Test](https://github.com/natir/rustyread/workflows/Test/badge.svg)
![Lints](https://github.com/natir/rustyread/workflows/Lints/badge.svg)
![MRV](https://github.com/natir/rustyread/workflows/MRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/rustyread/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/rustyread)
[![Documentation](https://github.com/natir/rustyread/workflows/Documentation/badge.svg)](https://natir.github.io/rustyread/rustyread)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/rustyread/blob/master/LICENSE)

<p align="center"><img src="images/logo.svg" alt="Rustyread" width="75%"></p>

Rustyread is a drop in replacement of `badread simulate`. Rustyread is very heavily inspired by [badread](https://github.com/rrwick/Badread), it reuses the same error and quality model file. But Rustyreads is multi-threaded and benefits from other optimizations.

- [Usage](#usage)
- [Installation](#installation)
- [Minimum supported Rust version](#minimum-supported-rust-version)
- [Difference with badread](#difference-with-badread)


**WARNING**:
- Rustyread has not yet been evaluated or even compared to *any* other long read generators
- Rustyread is tested only on Linux
- Rustyread is still in developpement many thing can change or be break


## Usage

If previously you called badread like this:

```
badread simulate --reference {reference path} --quantity {quantity} > {reads}.fastq
```

you can now replace badread by rustyread:

```
rustyread simulate --reference {reference path} --quantity {quantity} > {reads}.fastq
```

But by default rustyread use all avaible core you can control it with option `threads`:

```
rustyread --theads {number of thread} simulate --reference {reference path} --quantity {quantity} > {reads}.fastq
```

If you have `badread` installed in your python `sys.path` rustyread can found error and quality model automatically, but you can still use `--error_model` and `--qscore_model` option.

### Full usage

```
rustyread 0.2 Electabuzz
Pierre Marijon <pierre.marijon@hhu.de>
A long read simulator based on badread idea

USAGE:
    rustyread [FLAGS] [OPTIONS] <SUBCOMMAND>

FLAGS:
    -h, --help         Prints help information
    -v, --verbosity    verbosity level also control by environment variable RUSTYREAD_LOG if flag is
                       set RUSTYREAD_LOG value is ignored
    -V, --version      Prints version information

OPTIONS:
    -t, --threads <threads>    Number of thread use by rustyread, 0 use all avaible core, default
                               value 0

SUBCOMMANDS:
    help        Prints this message or the help of the given subcommand(s)
    simulate    Generate fake long read
```

```
rustyread-simulate
Generate fake long read

USAGE:
    rustyread simulate [FLAGS] [OPTIONS] --reference <reference-path> --quantity <quantity>

FLAGS:
    -h, --help                  Prints help information
        --small_plasmid_bias    If set, then small circular plasmids are lost when the fragment
                                length is too high (default: small plasmids are included regardless
                                of fragment length)
    -V, --version               Prints version information

OPTIONS:
        --chimera <chimera>
            Percentage at which separate fragments join together [default: 1]

        --end_adapter <end-adapter>
            Adapter parameters for read ends (rate and amount) [default: 50,20]

        --end_adapter_seq <end-adapter-seq>
            Adapter parameters for read ends [default: GCAATACGTAACTGAACGAAGT]

        --error_model <error-model>
            Path to an error model file [default: nanopore2020]

        --glitches <glitches>
            Read glitch parameters (rate, size and skip) [default: 10000,25,25]

        --identity <identity>
            Sequencing identity distribution (mean, max and stdev) [default: 85,95,5]

        --junk_reads <junk>
            This percentage of reads wil be low complexity junk [default: 1]

        --length <length>
            Fragment length distribution (mean and stdev) [default: 15000,13000]

        --output <output-path>                     Where read is write
        --qscore_model <qscore-model>
            Path to an quality score model file [default: nanopore2020]

        --quantity <quantity>
            Either an absolute value (e.g. 250M) or a relative depth (e.g. 25x)

        --random_reads <random>
            This percentage of reads wil be random sequence [default: 1]

        --reference <reference-path>               Reference fasta (can be gzipped, bzip2ped, xzped)
        --seed <seed>
            Random number generator seed for deterministic output (default: different output each
            time)

        --start_adapter <start-adapter>
            Adapter parameters for read starts (rate and amount) [default: 90,60]

        --start_adapter_seq <start-adapter-seq>
            Adapter parameters for read starts [default: AATGTACTTCGTTCAGTTACGTATTGCT]
```

## Installation

### Bioconda

If you haven't bioconda setup follow [this instruction](https://bioconda.github.io/user/install.html)

```
conda|mamba install rustyread
```

### With rust environment

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

#### With cargo

```
cargo install --git https://github.com/natir/rustyread.git
```

### From source

```
git clone https://github.com/natir/rustyread.git
cd rustyread
git checkout 0.2
cargo install --path .
```


## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.

## Difference with badread

- option `small_plasmid_bias` is silently ignored but small plasmid is 'sequence'
