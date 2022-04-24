# FeOs - gc-PC-SAFT

[![crate](https://img.shields.io/crates/v/feos-gc-pcsaft.svg)](https://crates.io/crates/feos-gc-pcsaft)
[![documentation](https://docs.rs/feos-gc-pcsaft/badge.svg)](https://docs.rs/feos-gc-pcsaft)
[![documentation](https://img.shields.io/badge/docs-github--pages-blue)](https://feos-org.github.io/feos/)

Implementation of the (heterosegmented) group contribution PC-SAFT equation of state[^gross2003][^sauer2014] and corresponding Helmholtz energy functional[^mairhofer2018][^rehner2021] within the FeOs project. This project contains a Rust implementation as well as bindings to Python.

## Usage in Python

If you want to use `feos-gc-pcsaft` in Python, take a look at the [`feos`-repository](https://github.com/feos-org/feos). `feos` contains multiple equation of state implementations in a single, easy-to-use Python package.

## FeOs

> FeOs is a framework for equations of state and classical density function theory

You can learn more about the principles behind `FeOs` [here](https://feos-org.github.io/feos/).

## Installation

Add this to your `Cargo.toml`

```toml
[dependencies]
feos-gc-pcsaft = "0.1"
```

## Test building python wheel

From within a Python virtual environment with `maturin` installed, type

```
maturin build --release --out dist --no-sdist -m build_wheel/Cargo.toml
```

[^gross2003]: [J. Gross, O. Spuhl, F. Tumakaka and G. Sadowski (2003). *Industrial & Engineering Chemistry Research*, 42(6), 1266-1274.](https://doi.org/10.1021/ie020509y)
[^sauer2014]: [E. Sauer, M. Stavrou and J. Gross (2014). *Industrial & Engineering Chemistry Research*, 53(38), 14854-14864.](https://doi.org/10.1021/ie502203w)
[^mairhofer2018]: [J. Mairhofer, B. Xiao and J. Gross (2018). *Fluid Phase Equilibria*, 472, 117-127.](https://doi.org/10.1016/j.fluid.2018.05.016)
[^rehner2021]: [P. Rehner, B. Bursik and J. Gross (2021). *Industrial & Engineering Chemistry Research*, 60(19), 7111-7123.](https://doi.org/10.1021/acs.iecr.1c00169)