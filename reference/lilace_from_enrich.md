# Create Lilace object from enrich processed counts output

Create Lilace object from enrich processed counts output

## Usage

``` r
lilace_from_enrich(file, pheno = "abundance")
```

## Arguments

- file:

  enrich counts file (see end of intro vignette for example input
  format)

- pheno:

  name of phenotype ("condition" row in counts table)

## Value

a lilace object, with data stored in \$data
