![Test](https://github.com/natir/rustyread/workflows/Test/badge.svg)
![Lints](https://github.com/natir/rustyread/workflows/Lints/badge.svg)
![MRV](https://github.com/natir/rustyread/workflows/MRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/rustyread/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/rustyread)
[![Documentation](https://github.com/natir/rustyread/workflows/Documentation/badge.svg)](https://natir.github.io/rustyread/rustyread)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/rustyread/blob/master/LICENSE)

# Rustyread, a rewrite of Badread in rust 🧬 💻

- [Instalation](#instalation)
- [Usage](#usage)
- [Minimum supported Rust version](#minimum-supported-rust-version)

## Instalation

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

### With cargo

Recommended solution.

```
cargo install --git https://github.com/natir/rustyread.git
```

### With source

```
git clone https://github.com/natir/rustyread.git
cd rustyread
cargo install --path .
```

## Usage

WIP

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.

## Difference with badread

### Error model

- If sum of alternative kmer probability is lower than 1.0, python generate a random error on fly, we build a static random error durring error model generation
