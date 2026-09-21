# Create Lilace object from counts

Formats required information for Lilace to run into a Lilace object. See
intro vignette for example usage.

## Usage

``` r
lilace_from_counts(
  variant_id,
  mutation_type,
  position,
  replicate,
  counts,
  metadata
)
```

## Arguments

- variant_id:

  a length V vector of unique variant identifiers such as hgvs
  nomenclature

- mutation_type:

  a vector of mutation types (e.g. synonymous, missense, etc). The label
  used to identify negative control label should be included in this.

- position:

  a vector of residue positions

- replicate:

  a vector of replicate ids

- counts:

  a (V x K) matrix of bin counts, where K is the number of bins

- metadata:

  a matrix or dataframe of any additional variant information to be
  incorporated into the final lilace output (e.g. wildtype residue)

## Value

a lilace object, with data stored in \$data
