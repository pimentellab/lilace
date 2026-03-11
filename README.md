
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

### R package installation

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
is \> 5. Additionally, if you run into issues with the TBB library, try
downgrading CmdStan to version 2.33.1. If it is not possible to install
CmdStan, our variational implementation of Lilace is developed in Python
and is available in a [separate
repo](https://github.com/jermoef/lilace-VI).

Then, install Lilace from GitHub.

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("pimentellab/lilace")
library(lilace)
```

### Docker installation

If you prefer to use docker or run into issues with the regular
installation, a docker image is available

    docker pull jfreudenberg/lilace

To connect the container to an interactive command line environment run

    docker container run -it lilace bash

To instead launch Rstudio in the container, you can specify a port and
run

    docker run -p 8888:8787 -e PASSWORD=<password> lilace

Then go to <http://localhost:8888/> and use username “rstudio” and your
input password to login. From there, you can call `library(lilace)` and
check that you can run the intro vignette code.

If you would like to build the container yourself or make adjustments,
the Dockerfile is available and can be built with

    docker buildx build --platform linux/amd64,linux/arm64 -t lilace .

## How to use Lilace

An introductory vignette can be found
[here](https://pimentellab.com/lilace/articles/intro.html). Variational
Lilace can be found [here](https://github.com/jermoef/lilace-VI). If you
run into any problems, please submit an issue on github or email
<jfreudenberg@ucla.edu>.

## Citing Lilace

If you use Lilace for your research, please cite our paper:

> Freudenberg, J., Rao, J., Howard, M.K. *et al.* Accurate variant
> effect estimation in FACS-based deep mutational scanning data with
> Lilace. *Genome Biol* **27**, 48 (2026). DOI:
> [10.1186/s13059-026-03934-1](https://doi.org/10.1186/s13059-026-03934-1)
