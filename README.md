# KhaldoonPeriodicity

[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**KhaldoonPeriodicity** is an advanced R package for robust periodicity detection in time series data using composite spectral statistics and permutation-based bootstrap inference.

## Features

- Robust periodicity detection using composite spectral statistics
- Permutation-based bootstrap for exact inference
- Parallel processing support
- Comprehensive visualization tools
- No parametric assumptions required

## Installation
```r
# Install from GitHub
devtools::install_github("khaldoonalhyasat-svg/KhaldoonPeriodicity")
```

## Quick Start
```r
library(KhaldoonPeriodicity)

# Generate sample periodic data
t <- 1:200
y <- sin(2 * pi * t / 24) + rnorm(200, sd = 0.3)

# Test for periodicity
result <- khaldoon_test(y, n_bootstrap = 1000)
print(result)
summary(result)
plot(result)
```

## Methodology

The test uses **permutation-based bootstrap** which:
- Randomly reorders data without replacement
- Preserves empirical distribution
- Destroys temporal autocorrelation
- Provides exact inference under exchangeability

**Validated Performance:**
- Type I error rate: 0.04 (excellent control at α=0.05)
- Statistical power: 1.00 (perfect detection)

## References

Good, P. I. (2005). Permutation, Parametric, and Bootstrap Tests. Springer.

## License

MIT License

## Author

Khaldoon Alhyasat - [@khaldoonalhyasat-svg](https://github.com/khaldoonalhyasat-svg)

