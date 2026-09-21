# Run Lilace

Runs Lilace on given data. Lilace will run on \$normalized_data if it
exists, otherwise it will use \$data.

## Usage

``` r
lilace_fit_model(
  lilace_obj,
  output_dir,
  control_label = "synonymous",
  control_correction = TRUE,
  use_positions = TRUE,
  pseudocount = TRUE,
  seed = NULL,
  min_total_counts = 15,
  n_parallel_chains = 4,
  nc_seed = NULL
)
```

## Arguments

- lilace_obj:

  initialized lilace object

- output_dir:

  output directory to write scores and sampling logs to

- control_label:

  label from \$data\$type column to use as negative controls to score
  against

- control_correction:

  boolean for whether to use negative control scores as bias correction

- use_positions:

  boolean for whether to use position hierarchy to improve estimation

- pseudocount:

  boolean for whether to add pseudocount (+1 to all counts) for fitting
  model

- seed:

  random seed for sampling process to get exactly reproducible results.
  A NULL value indicates no fixed seed.

- min_total_counts:

  minimum total counts for a variant–anything less will be filtered out

- n_parallel_chains:

  number of chains to run in parallel

- nc_seed:

  negative control seed (default: same as seed param)

## Value

lilace object with \$scores and \$fitted_data entries

## Examples

``` r
if (FALSE) { # \dontrun{
lilace_obj <- lilace_fit_model(lilace_obj, output_dir, control_label="synonymous",
                               control_correction=T, use_positions=T, pseudocount=T)
} # }
```
