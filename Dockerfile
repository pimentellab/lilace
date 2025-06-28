# get R
FROM rocker/rstudio:4.4

# Install system dependencies
RUN apt-get update && apt-get install -y git cmake g++ && apt-get clean

RUN mkdir -p /opt/cmdstan

# install R packages
ARG WHEN
RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/${WHEN}'))"
# RUN R -e "install.packages('tidyverse', dependencies = TRUE);     if (!library(tidyverse, logical.return=T)) quit(status=10)" \
RUN R -e "install.packages('cowplot', dependencies = TRUE);   if (!library(cowplot, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('ggpubr', dependencies = TRUE); if (!library(ggpubr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('colorspace', dependencies = TRUE);      if (!library(colorspace, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('dplyr', dependencies = TRUE);      if (!library(dplyr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('ggplot2', dependencies = TRUE);   if (!library(ggplot2, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('remotes', dependencies = TRUE); if (!library(remotes, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('readr', dependencies = TRUE);      if (!library(readr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('stringr', dependencies = TRUE);   if (!library(stringr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('tidyr', dependencies = TRUE); if (!library(tidyr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('utils', dependencies = TRUE);      if (!library(utils, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('yaml', dependencies = TRUE); if (!library(yaml, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('knitr', dependencies = TRUE);      if (!library(knitr, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('rmarkdown', dependencies = TRUE); if (!library(rmarkdown, logical.return=T)) quit(status=10)" \
&& R -e "install.packages('posterior', dependencies = TRUE);      if (!library(posterior, logical.return=T)) quit(status=10)"

# install cmdstanr and CmdStan
# RUN R -e "install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))"
RUN R -e "remotes::install_github('stan-dev/cmdstanr@v0.9.0')"
RUN R -e "library(cmdstanr); install_cmdstan(dir='/opt/cmdstan', cores=4, version='2.32.2')"
ENV CMDSTAN="/opt/cmdstan/cmdstan-2.32.2"
RUN R -e "library(cmdstanr); print(cmdstan_version())"

# install Lilace
RUN R -e "remotes::install_github('pimentellab/lilace')"