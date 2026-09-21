# Generate score histograms/density plots by mutation type

generates a density plot for visualizing the distribution of scores
across different mutation types (same as Rosace scoreDensity function)

## Usage

``` r
lilace_score_density(
  data,
  savedir,
  type.col = "type",
  score.col = "mean",
  hist = FALSE,
  nbins = 30,
  c.fill = c("#FF7575", "lightgreen", "#7298BF"),
  alpha = 0.5,
  x.text = 10,
  y.text = 10,
  scale.free = FALSE,
  space.free = FALSE,
  ht = 10,
  wd = 8,
  name = "DensityPlot",
  savepdf = TRUE,
  savesvg = FALSE,
  show = TRUE
)
```

## Arguments

- data:

  Scores data frame. Expected to have columns containing information
  about position, control amino acid, mutated amino acid, mutation type,
  and score.

- savedir:

  Character string specifying the directory to save plots.

- type.col:

  Column name in \`data\` for mutation types. Default is \`type\`.

- score.col:

  Column name in \`data\` for mutation scores. Default is \`mean\`.

- hist:

  Logical indicating whether to plot count histogram or density. Default
  is \`FALSE\`.

- nbins:

  Numeric value specifying the number of bins for the histogram. Default
  is \`30\`.

- c.fill:

  Vector indicating the fill color for mutation types.

- alpha:

  Numeric value between 0-1 indicating the fill transparency. Default is
  \`0.5\`

- x.text:

  Numeric value for x-axis text size. Default is \`10\`.

- y.text:

  Numeric value for x-axis text size. Default is \`10\`.

- scale.free:

  Logical indicating whether to make score range proportional. Default
  is \`FALSE\`.

- space.free:

  Logical indicating whether to make panel heights variable. Default is
  \`FALSE\`.

- ht:

  Numeric value for the height of the saved plot. Default is \`10\`.

- wd:

  Numeric value for the width of the saved plot. Default is \`8\`.

- name:

  Character string specifying the base name of the saved file.

- savepdf:

  Logical indicating whether to also save a PDF version of the plot.
  Default is \`TRUE\`.

- savesvg:

  logical indicating whether to save as svg

- show:

  Logical indicating whether or not to display the plot in the viewer.
  Default is \`TRUE\`.

## Value

NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
lilace_score_density(scores, output_dir, score.col="effect",
                     name="score_histogram", hist=T, scale.free=T)
} # }
```
