# summarize posterior samples (intended for internal use)

summarize posterior samples (intended for internal use)

## Usage

``` r
.summarize_posteriors(
  fit,
  data,
  input,
  control_correction,
  control_label,
  use_positions,
  nc_seed = NULL
)
```

## Arguments

- fit:

  stan output fit object

- data:

  data

- input:

  stan input list

- control_correction:

  logical for whether to do negative control correction

- control_label:

  label for negative control variant types

- use_positions:

  logical for whether positions used in running Lilace

- nc_seed:

  negative control correction seed
