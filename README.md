
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Lilace

<!-- badges: start -->
<!-- badges: end -->
<p align="left">
<img src="man/figures/lilace-logo.png" width="150">
</p>

**Lilace** is an R package for scoring FACS-based DMS experiments with
uncertainty quantification. It takes in a negative control group
(usually synonymous variants) and scores each variant relative to the
negative control group.

## Installation

Lilace relies on [cmdstanr](https://mc-stan.org/cmdstanr/), which should
be properly installed first.

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 2)
```

The compiler requirements can be seen at
[stan-dev](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms#supported-cpp-versions-and-compilers).
If you run into issues with installation, please ensure your gcc version
is \> 5.

Then, install Lilace from GitHub.

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("pimentellab/lilace")
library(lilace)
```

## How to use Lilace

An introductory vignette can be found
[here](https://github.com/jermoef/lilace/tree/main/vignettes). If you
run into any problems, please submit an issue on github or email
<jfreudenberg@ucla.edu>.
