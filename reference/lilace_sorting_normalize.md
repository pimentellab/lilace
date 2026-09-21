# Normalize data to cell sorting proportions

Scales the read counts in a bin to match the sort proportions given.

## Usage

``` r
lilace_sorting_normalize(lilace_obj, sort_prop, rep_specific)
```

## Arguments

- lilace_obj:

  a lilace obj, such as created using lilace_from_counts()

- sort_prop:

  a list of cell sorting bin proportions to normalize If replicate
  specific, input as list(\<rep_label\>=c(0.25, 0.25, 0.25, 0.25)) where
  \<rep_label\> corresponds to the lilace object replicate labels If not
  replicate specific, input as vector.

- rep_specific:

  boolean indicating whether sort_prop is replicate specific or not

## Value

lilace object with \$normalized_data element

## Examples

``` r
# Example 1 (not replicate specific)
if (FALSE) { # \dontrun{
sort_prop <- c(0.2184, 0.1972, 0.1869, 0.1236)
lilace_obj <- lilace_sorting_normalize(lilace_obj, sort_prop, rep_specific=F)
} # }

# Example 2 (replicate specific)
if (FALSE) { # \dontrun{
sort_prop=list(R1=c(0.2184, 0.1972, 0.1869, 0.1236),
               R2=c(0.25, 0.25, 0.25, 0.25),
               R3=c(0.3, 0.2, 0.3, 0.2))
lilace_obj <- lilace_sorting_normalize(lilace_obj, sort_prop, rep_specific=T)
} # }
```
