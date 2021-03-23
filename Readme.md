[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/badread-rs/blob/master/LICENSE)
![CI](https://github.com/natir/badread-rs/workflows/CI/badge.svg)
[![Documentation](https://github.com/natir/badread-rs/workflows/Documentation/badge.svg)](https://natir.github.io/badread-rs/badread_rs)
[![CodeCov](https://codecov.io/gh/natir/badread-rs/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/badread-rs)

# Badread-rs, a rewrite of Badread in rust 🧬 💻

- [Instalation](#instalation)
- [Usage](#usage)
- [Minimum supported Rust version](#minimum-supported-rust-version)

## Instalation

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

### With cargo

Recommended solution.

```
cargo install --git https://github.com/natir/badread-rs.git
```

### With source

```
git clone https://github.com/natir/badread-rs.git
cd badread-rs
cargo install --path .
```

## Usage

WIP

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.

## Difference with python badread

### Error model

- If some of alternative kmer is lower than 1.0, python generate a random error on fly, we build a static random error durring error model generation
