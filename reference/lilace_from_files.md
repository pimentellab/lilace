# Create Lilace object from separate count files for each bin and replicate

Create Lilace object from separate count files for each bin and
replicate

## Usage

``` r
lilace_from_files(
  file_list,
  variant_id_col = "hgvs",
  position_col = "position",
  mutation_type_col = "type",
  count_col = "count",
  delim = "\t"
)
```

## Arguments

- file_list:

  a list containing a separate list of bin count files for each
  replicate (see end of intro vignette for example)

- variant_id_col:

  column name for variant id

- position_col:

  column name for position information

- mutation_type_col:

  column name for mutation type information

- count_col:

  column name for counts

- delim:

  the field separator character (see ?readr::read.delim)

## Value

a lilace object, with data stored in \$data
