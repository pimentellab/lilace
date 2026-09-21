# Introduction to Lilace

## Overview

Lilace is a tool for scoring FACS-based DMS experiments with uncertainty
quantification. It takes in a negative control group (usually synonymous
variants) and scores each variant relative to the negative control
group.

The standard workflow is:

1.  Import data into Lilace format
2.  Normalize to cell sorting percentages (if available)
3.  Run Lilace. On a full length protein, Lilace usually takes a few
    hours to run (see paper supplement for more details)
4.  Analyze Lilace output

This notebook shows an example run of Lilace from installation to
analysis on a toy dataset. We provide functions to load counts into
Lilace’s input format (step 1), normalize counts to cell sorting
proportions (step 2), run Lilace (step 3), and plot the results (step
4). By default, Lilace incorporates variance in the negative control
scores into effect size uncertainty (standard errors) and makes use of
similar effects at a position to improve estimation. Both of these
behaviors can optionally be turned off when running Lilace.

## Installation

Lilace is built using `stan` and interacts with `CmdStan` using the
package `cmdstanr`. This requires C++ compilation. The compiler
requirements can be seen at
[stan-dev](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms#supported-cpp-versions-and-compilers).
If you run into issues with installation, please ensure your gcc version
is \> 5.

First, we install `cmdstanr`. We recommend running this in a fresh R
session or restarting your current session.

``` r

install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

Then, we use `cmdstanr` to install `CmdStan`. This requires a working
C++ toolchain and compiler.

``` r

library(cmdstanr)
install_cmdstan(cores = 2) # number of cores to use for installation
```

We can check the `CmdStan` version to verify correct installation.

``` r

cmdstan_version()
```

Now, we can install Lilace from github.

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("pimentellab/lilace")
library(lilace)
```

## Load input data

Next, we load in the toy dataset, which includes the first 20 positions
of a GPR68 screen from
(<https://www.cell.com/cell/fulltext/S0092-8674(24)01373-4>). This data
is provided with a standard installation of Lilace.

``` r

load(system.file("extdata", "gpr68_ph55_20pos.RData", package="lilace"))
head(dataset)
```

    ## # A tibble: 6 × 14
    ##   hgvs_exp  hgvs  position wildtype mutation type    c_0   c_1   c_2   c_3 exp  
    ##   <chr>     <chr>    <dbl> <chr>    <chr>    <chr> <dbl> <dbl> <dbl> <dbl> <chr>
    ## 1 p.(G2A)_… p.(G…        2 G        A        miss…   118   146   215   260 ph55 
    ## 2 p.(G2A)_… p.(G…        2 G        A        miss…   146   108   127   284 ph55 
    ## 3 p.(G2A)_… p.(G…        2 G        A        miss…   144    55   106   134 ph55 
    ## 4 p.(G2C)_… p.(G…        2 G        C        miss…   148   101   384   507 ph55 
    ## 5 p.(G2C)_… p.(G…        2 G        C        miss…   203   201   307   616 ph55 
    ## 6 p.(G2C)_… p.(G…        2 G        C        miss…   236   175   262   547 ph55 
    ## # ℹ 3 more variables: rep <chr>, n_counts <dbl>, total_counts <dbl>

We then put the data into a Lilace object using the
[`lilace_from_counts()`](http://pimentellab.com/lilace/reference/lilace_from_counts.md)
function, which requires:

1.  variant identifier (e.g. hgvs), ensure this is unique for each
    variant
2.  mutation type (e.g. synonymous, missense)–this will be used to
    identify the negative controls
3.  residue position
4.  replicate information
5.  the FACS bin counts
6.  any metadata you wish to carry over into the final output.

``` r

# load from components
lilace_obj <- lilace_from_counts(variant_id=dataset$hgvs, mutation_type=dataset$type, 
                      position=dataset$position, replicate=dataset$rep, 
                      counts=dataset %>% select(starts_with("c_")), 
                      metadata=dataset %>% select(wildtype, mutation, exp))
head(lilace_obj$data)
```

    ## # A tibble: 6 × 13
    ##   variant wildtype mutation exp   type    position rep     c_0   c_1   c_2   c_3
    ##   <chr>   <chr>    <chr>    <chr> <chr>      <dbl> <chr> <dbl> <dbl> <dbl> <dbl>
    ## 1 p.(G2A) G        A        ph55  missen…        2 R1      118   146   215   260
    ## 2 p.(G2A) G        A        ph55  missen…        2 R2      146   108   127   284
    ## 3 p.(G2A) G        A        ph55  missen…        2 R3      144    55   106   134
    ## 4 p.(G2C) G        C        ph55  missen…        2 R1      148   101   384   507
    ## 5 p.(G2C) G        C        ph55  missen…        2 R2      203   201   307   616
    ## 6 p.(G2C) G        C        ph55  missen…        2 R3      236   175   262   547
    ## # ℹ 2 more variables: n_counts <dbl>, total_counts <dbl>

We also provide the functions
[`lilace_from_enrich()`](http://pimentellab.com/lilace/reference/lilace_from_enrich.md),
which can load counts from enrich2’s output format, and
[`lilace_from_files()`](http://pimentellab.com/lilace/reference/lilace_from_files.md),
which can load counts from separate bin and replicate count files, but
they may not generalize to all dataset formats. See the end of this
vignette for an example call of each.

## Normalize to cell sorting proportions

FACS gating is often set to contain equal proportions of the overall
cell fluorescence distribution, but differences in bin PCR amplification
can lead to the over and under representation of certain bins in the
read counts. If this is the case, the data can be normalized to the
sorting proportions using the
[`lilace_sorting_normalize()`](http://pimentellab.com/lilace/reference/lilace_sorting_normalize.md)
function, which takes in a Lilace object and the proportions to
normalize to.

![Example sort proportions.](gpr68_ph55_sort.png)

Example sort proportions.

In this case the sort proportions don’t add up to 1, but that’s okay,
they are scaled to add up to 1 internally.

``` r

# normalize to sorting proportions
GPR_FACS_trace_percent_55 <- c(0.2184, 0.1972, 0.1869, 0.1236)

lilace_obj <- lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, rep_specific=F)
```

    ## Warning in lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, :
    ## Sorting proportions sum up to 0.7261 instead of 1.

``` r

head(lilace_obj$normalized_data)
```

    ## # A tibble: 6 × 13
    ##   variant wildtype mutation exp   type    position rep     c_0   c_1   c_2   c_3
    ##   <chr>   <chr>    <chr>    <chr> <chr>      <dbl> <chr> <dbl> <dbl> <dbl> <dbl>
    ## 1 p.(G2A) G        A        ph55  missen…        2 R1      177   164   199   160
    ## 2 p.(G2A) G        A        ph55  missen…        2 R2      192   131   148   150
    ## 3 p.(G2A) G        A        ph55  missen…        2 R3      199    69    92    86
    ## 4 p.(G2C) G        C        ph55  missen…        2 R1      222   113   355   312
    ## 5 p.(G2C) G        C        ph55  missen…        2 R2      267   243   356   324
    ## 6 p.(G2C) G        C        ph55  missen…        2 R3      327   218   227   351
    ## # ℹ 2 more variables: n_counts <dbl>, total_counts <dbl>

If we wish to normalize to replicate-specific sorting proportions, we
can format it like the following. The replicate labels in the list (R1,
R2, R3) should correspond to the labels used in the Lilace input.

``` r

# normalize to sorting proportions
GPR_FACS_trace_percent_55 <- list(R1=c(0.2184, 0.1972, 0.1869, 0.1236),
                                  R2=c(0.2184, 0.1972, 0.1869, 0.1236),
                                  R3=c(0.2184, 0.1972, 0.1869, 0.1236))

lilace_obj <- lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, rep_specific=T)
```

    ## Warning in lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, :
    ## Sorting proportions sum up to 0.7261 instead of 1.
    ## Warning in lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, :
    ## Sorting proportions sum up to 0.7261 instead of 1.
    ## Warning in lilace_sorting_normalize(lilace_obj, GPR_FACS_trace_percent_55, :
    ## Sorting proportions sum up to 0.7261 instead of 1.

``` r

head(lilace_obj$normalized_data)
```

    ## # A tibble: 6 × 13
    ##   variant wildtype mutation exp   type    position rep     c_0   c_1   c_2   c_3
    ##   <chr>   <chr>    <chr>    <chr> <chr>      <dbl> <chr> <dbl> <dbl> <dbl> <dbl>
    ## 1 p.(G2A) G        A        ph55  missen…        2 R1      177   164   199   160
    ## 2 p.(G2A) G        A        ph55  missen…        2 R2      192   131   148   150
    ## 3 p.(G2A) G        A        ph55  missen…        2 R3      199    69    92    86
    ## 4 p.(G2C) G        C        ph55  missen…        2 R1      222   113   355   312
    ## 5 p.(G2C) G        C        ph55  missen…        2 R2      267   243   356   324
    ## 6 p.(G2C) G        C        ph55  missen…        2 R3      327   218   227   351
    ## # ℹ 2 more variables: n_counts <dbl>, total_counts <dbl>

## Run Lilace

Once, we are ready to run Lilace, we first specify an output directory
where the output and model fitting logs will be saved.

``` r

output_dir <- "output"
```

Then, we use the function
[`lilace_fit_model()`](http://pimentellab.com/lilace/reference/lilace_fit_model.md)
to run the main Lilace algorithm. This should only take a few minutes on
this truncated dataset, but usually takes a few hours on a full length
protein.

By default, the negative controls will be identified by variants labeled
as `synonymous` in the `lilace$data$type` column. This can be changed by
setting the `control_label` argument to a different value in the column.
To disable the negative control-based bias correction, set
`control_correction` to `FALSE`. To disable position-level grouping, set
`use_positions` to `FALSE`. We provide example outputs for both of these
options towards the end of the vignette.

``` r

lilace_obj <- lilace_fit_model(lilace_obj, output_dir, control_label="synonymous", 
                               control_correction=T, use_positions=T, pseudocount=T,
                               n_parallel_chains=4, seed=1)
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 3 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 2 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 1 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 4 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 2 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 1 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 4 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 4 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 91.6 seconds.
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 94.7 seconds.
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 126.2 seconds.
    ## Chain 4 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 152.8 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 116.3 seconds.
    ## Total execution time: 152.9 seconds.

Once it’s done running, Lilace saves the scores for each variant to
`output/lilace_output/variant_scores.tsv`. The scores dataframe is also
saved in `lilace_obj$scores`. The scores are also combined with the
input dataframe at `lilace_obj$fitted_data`.

``` r

scores <- read_tsv("output/lilace_output/variant_scores.tsv", show_col_types = FALSE)
head(scores)
```

    ## # A tibble: 6 × 12
    ##   variant        type  wildtype mutation exp   position  effect effect_se   lfsr
    ##   <chr>          <chr> <chr>    <chr>    <chr>    <dbl>   <dbl>     <dbl>  <dbl>
    ## 1 p.(G2del)      dele… G        del0     ph55         2 -0.0719     0.285 0.374 
    ## 2 p.(G2_N3del)   dele… G        del1     ph55         2  0.349      0.285 0.0978
    ## 3 p.(G2_I4del)   dele… G        del2     ph55         2  0.228      0.285 0.205 
    ## 4 p.(G2_N3insG)  inse… G        ins0     ph55         2  0.281      0.281 0.151 
    ## 5 p.(G2_N3insGS) inse… G        ins1     ph55         2  0.104      0.282 0.367 
    ## 6 p.(G2A)        miss… G        A        ph55         2  0.104      0.282 0.362 
    ## # ℹ 3 more variables: pos_mean <dbl>, pos_sd <dbl>, discovery05 <dbl>

The `effect` column represents the estimated effect size / score for a
variant, while `effect_se` gives the standard error. `lfsr` is the local
false sign rate, which is a Bayesian measure of significance analogous
to a p-value. For example, a default way to call discoveries would be to
set a threshold of `lfsr < 0.05`. We encode this in the `discovery05`
column where a `-1` indicates a significant leftward shift in
fluorescence, `+1` a significant rightward shift, and `0` is not
significant. `pos_mean` and `pos_sd` are the estimated position-level
mean and variance parameters–these parameters are only included in the
final score table if Lilace was run with `use_positions=T`.

## Plot results

To analyze the output of the model, we can first plot the score
distributions of the different mutation types using
[`lilace_score_density()`](http://pimentellab.com/lilace/reference/lilace_score_density.md)
and ensure they align with expectations.

``` r

lilace_score_density(scores, output_dir, score.col="effect", 
                     name="score_histogram", hist=T, scale.free=T)
```

![](intro_files/figure-html/unnamed-chunk-13-1.png)

We observe a substantial LOF mode for indels, a smaller one for missense
mutations, and a distribution around zero for synonymous mutations.

Next, we plot the estimated scores on a heatmap.

``` r

lilace_score_heatmap(scores, output_dir, score.col="effect", name="score_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-14-1.png)

As well as the discoveries, at a threshold of `lfsr < 0.05`. Again, -1
indicates a significant leftward shift in fluorescence, +1 a rightward
shift, and 0 is not significant.

``` r

lilace_score_heatmap(scores, output_dir, score.col="discovery05", name="discovery05_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-15-1.png)

From here on, further analysis can be done, such as mapping to protein
structure and qualitatively analyzing the results.

## Example runs without negative control correction or position hierarchy

#### Without the negative control-based bias correction

Without the negative control-based bias correction, the effect size
standard errors are smaller and more effects are called significant (at
the cost of more false positives). You may want to run this in cases
where variance in negative controls do not provide a good estimate of
experimental uncertainty (e.g. only a single WT variant), but you still
want Lilace scores for each variant. In this case, Lilace uncertainty
estimates (effect standard errors) do not incorporate any notion of
negative control variance, which may result in more false positives.
(note: a control label column is still required as a reference set to
score against)

``` r

# specify outpur dir
output_dir <- "output_no_nc_correction"
# run the model on data
fitted_data <- lilace_fit_model(lilace_obj, output_dir, control_label="synonymous", 
                                control_correction=F, seed=1)
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 3 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 2 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 1 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 4 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 1 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 2 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 4 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 4 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 88.4 seconds.
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 93.4 seconds.
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 123.2 seconds.
    ## Chain 4 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 150.7 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 114.0 seconds.
    ## Total execution time: 150.9 seconds.

``` r

# read scores file
scores <- read_tsv("output_no_nc_correction/lilace_output/variant_scores.tsv", show_col_types = FALSE)
print(head(scores))
```

    ## # A tibble: 6 × 12
    ##   variant       type  wildtype mutation exp   position  effect effect_se    lfsr
    ##   <chr>         <chr> <chr>    <chr>    <chr>    <dbl>   <dbl>     <dbl>   <dbl>
    ## 1 p.(G2del)     dele… G        del0     ph55         2 -0.0642     0.125 0.313  
    ## 2 p.(G2_N3del)  dele… G        del1     ph55         2  0.357      0.125 0.00225
    ## 3 p.(G2_I4del)  dele… G        del2     ph55         2  0.228      0.110 0.02   
    ## 4 p.(G2_N3insG) inse… G        ins0     ph55         2  0.279      0.114 0.005  
    ## 5 p.(G2_N3insG… inse… G        ins1     ph55         2  0.105      0.105 0.158  
    ## 6 p.(G2A)       miss… G        A        ph55         2  0.111      0.121 0.188  
    ## # ℹ 3 more variables: pos_mean <dbl>, pos_sd <dbl>, discovery05 <dbl>

``` r

lilace_score_heatmap(scores, output_dir, score.col="effect", name="score_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-16-1.png)

``` r

lilace_score_heatmap(scores, output_dir, score.col="discovery05", name="discovery05_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-16-2.png)

#### Without the position hierarchy

Without using position effects to inform variant effect estimates, each
effect size will be completely independent of the others at its
position. You may want to run this if the data has more than a single
residue mutation at a time or complete position independence is desired
for downstream analysis (e.g. clinical classification that uses multiple
effects at a position as a line of evidence).

``` r

# specify outpur dir
output_dir <- "output_nopos"
# run the model on data
fitted_data <- lilace_fit_model(lilace_obj, output_dir, control_label="synonymous", 
                                use_positions=F, seed=1)
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
    ## Chain 4 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 1 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 3 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 2 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 3 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 1 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 2 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 1 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 4 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 1 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 3 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 2 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 92.1 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 93.6 seconds.
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 93.9 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 94.1 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 93.4 seconds.
    ## Total execution time: 94.3 seconds.

``` r

# read scores file
scores <- read_tsv("output_nopos/lilace_output/variant_scores.tsv", show_col_types = FALSE)
print(head(scores))
```

    ## # A tibble: 6 × 10
    ##   variant        type  wildtype mutation exp   position  effect effect_se   lfsr
    ##   <chr>          <chr> <chr>    <chr>    <chr>    <dbl>   <dbl>     <dbl>  <dbl>
    ## 1 p.(G2del)      dele… G        del0     ph55         2 -0.201      0.295 0.232 
    ## 2 p.(G2_N3del)   dele… G        del1     ph55         2  0.449      0.298 0.0575
    ## 3 p.(G2_I4del)   dele… G        del2     ph55         2  0.247      0.292 0.191 
    ## 4 p.(G2_N3insG)  inse… G        ins0     ph55         2  0.323      0.295 0.128 
    ## 5 p.(G2_N3insGS) inse… G        ins1     ph55         2  0.0826     0.288 0.398 
    ## 6 p.(G2A)        miss… G        A        ph55         2  0.0663     0.297 0.420 
    ## # ℹ 1 more variable: discovery05 <dbl>

``` r

lilace_score_heatmap(scores, output_dir, score.col="effect", name="score_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-17-1.png)

``` r

lilace_score_heatmap(scores, output_dir, score.col="discovery05", name="discovery05_heatmap", 
                     x.text=4, seq.text=1.5, y.text=3)
```

![](intro_files/figure-html/unnamed-chunk-17-2.png)

## View posterior samples

For convenience and advanced analysis, the posterior samples are saved
in `lilace_output/posterior_samples.rds`. They can be analyzed further
with the `posterior` package. Here is an example of how to access them:

``` r

library(posterior)
draws <- readRDS("output/lilace_output/posterior_samples.rds")
posterior_samples <- posterior::as_draws_rvars(draws)
head(posterior::draws_of(posterior_samples$a))
```

    ##        [,1]       [,2]      [,3]
    ## 1 0.0186883 0.00572697 0.0245361
    ## 2 0.0157604 0.00556233 0.0201306
    ## 3 0.0214507 0.01091660 0.0230276
    ## 4 0.0171415 0.00890432 0.0241538
    ## 5 0.0203703 0.00793405 0.0195903
    ## 6 0.0131362 0.00361930 0.0232617

## Format input data from Enrich2 output format

This function expects data in the following format created through
processing counts with Enrich2. ![Enrich2 counts example
format](enrich_format.png)

``` r

enrich_file <- system.file("extdata", "kir21_enrich_format.tsv", package="lilace")
lilace_obj_from_enrich <- lilace_from_enrich(enrich_file, pheno="abundance") # filtering to just the abundance phenotype
```

    ## Formatting from enrich2 output. Synonymous mutations will be used for negative control-based bias correction. If this is not desirable, please input using the format_data() function.

``` r

head(lilace_obj_from_enrich$data)
```

    ## # A tibble: 6 × 13
    ##   variant  position wildtype mutation type     c_0   c_1   c_2   c_3 exp   rep  
    ##   <chr>       <dbl> <chr>    <chr>    <chr>  <dbl> <dbl> <dbl> <dbl> <chr> <chr>
    ## 1 p.(A22A)       22 A        A        synon…   116   181   176   278 abun… R1   
    ## 2 p.(A22A)       22 A        A        synon…   221   230   420   397 abun… R2   
    ## 3 p.(A22A)       22 A        A        synon…    94   133   119   257 abun… R4   
    ## 4 p.(A22A)       22 A        A        synon…   207   210   385   412 abun… R5   
    ## 5 p.(A22C)       22 A        C        misse…    60   187    96   115 abun… R1   
    ## 6 p.(A22C)       22 A        C        misse…   102   148   104    60 abun… R2   
    ## # ℹ 2 more variables: n_counts <dbl>, total_counts <dbl>

## Format input from separate gating files

Here, we provide an example of formatting from multiple input files
where the counts for each bin and replicate are in different files (such
as from the Dumpling processing pipeline). The counts will be joined by
shared column names, so each file must have at least a variant id
column, a position column, and a mutation type column. For example the
counts from bin 1 in replicate 1 could look like: ![Separate file counts
example format](sep_file_counts.png)

``` r

# input files in list format 
# (can add more replicates or bins--make sure they are in order from left to right)
file_list <- list(R1=list(
  system.file("extdata", "multi_input/R1_A.tsv", package="lilace"), # rep 1 bin 1
  system.file("extdata", "multi_input/R1_B.tsv", package="lilace"), # rep 1 bin 2
  system.file("extdata", "multi_input/R1_C.tsv", package="lilace"), # rep 1 bin 3
  system.file("extdata", "multi_input/R1_D.tsv", package="lilace")  # rep 1 bin 4 
  ),
                  R2=list(
  system.file("extdata", "multi_input/R2_A.tsv", package="lilace"), # rep 2 bin 1 
  system.file("extdata", "multi_input/R2_B.tsv", package="lilace"), # rep 2 bin 2
  system.file("extdata", "multi_input/R2_C.tsv", package="lilace"), # rep 2 bin 3
  system.file("extdata", "multi_input/R2_D.tsv", package="lilace")) # rep 2 bin 4
  ) 
lilace_obj_from_multi_file <- lilace_from_files(file_list, 
                                                variant_id_col="hgvs",
                                                position_col="position",
                                                mutation_type_col="type",
                                                count_col="count", # column containing count is called count
                                                delim="\t") 

head(lilace_obj_from_multi_file$data)
```

    ## # A tibble: 6 × 12
    ##   variant wildtype mutation type     position rep     c_0   c_1   c_2   c_3
    ##   <chr>   <chr>    <chr>    <chr>       <int> <chr> <int> <int> <int> <int>
    ## 1 p.(G2A) G        A        missense        2 R1      118   146   215   260
    ## 2 p.(G2A) G        A        missense        2 R2      146   108   127   284
    ## 3 p.(G2C) G        C        missense        2 R1      148   101   384   507
    ## 4 p.(G2C) G        C        missense        2 R2      203   201   307   616
    ## 5 p.(G2D) G        D        missense        2 R1      201   185   212   328
    ## 6 p.(G2D) G        D        missense        2 R2      176   142   119   343
    ## # ℹ 2 more variables: n_counts <dbl>, total_counts <dbl>
