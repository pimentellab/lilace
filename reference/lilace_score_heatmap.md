# Generate score heatmap

generates an amino acid by position heatmap for input dataset and given
column (derived from Rosace scoreHeatmap function)

## Usage

``` r
lilace_score_heatmap(
  data,
  savedir,
  ctrl.name = "synonymous",
  pos.col = "position",
  wt.col = "wildtype",
  mut.col = "mutation",
  type.col = "type",
  score.col = "mean",
  aa.order = NA,
  npos = 100,
  ncol = 1,
  pos.step = 5,
  x.text = 6,
  y.text = 3,
  seq.text = 1.1,
  c.pallete = "RdBu",
  c.scale = list(),
  ht = 11,
  wd = 8.5,
  name = "Heatmap",
  savepdf = TRUE,
  savesvg = FALSE,
  show = FALSE,
  factor_score = FALSE,
  discovery_score = FALSE,
  compare_score = FALSE,
  category_score = FALSE,
  cat_name = "discovery"
)
```

## Arguments

- data:

  Scores data frame. Expected to have columns containing information
  about position, control amino acid, mutated amino acid, mutation type,
  and score.

- savedir:

  Character string specifying the directory to save plots.

- ctrl.name:

  String specifying the control mutation name. Default is
  \`synonymous\`.

- pos.col:

  Column name in \`data\` for mutation positions. Default is
  \`position\`.

- wt.col:

  Column name in \`data\` for wildtype amino acids. Default is
  \`wildtype\`.

- mut.col:

  Column name in \`data\` for mutated amino acids. Default is
  \`mutation\`.

- type.col:

  Column name in \`data\` for mutation types. Default is \`type\`.

- score.col:

  Column name in \`data\` for mutation scores. Default is \`mean\`.

- aa.order:

  Character vector defining the order of amino acid mutations in the
  y-axis. Default to using mutations from \`data\` in alphabetical
  order.

- npos:

  Integer specifying the number of positions per subplot. Default is
  \`100\`.

- ncol:

  Integer specifying the number of columns of subplots. Default is
  \`1\`.

- pos.step:

  Integer specifying the steps between x-axis labels. Default is \`5\`.

- x.text:

  Numeric value for x-axis text size. Default is \`6\`.

- y.text:

  Numeric value for y-axis text size. Default is \`3\`.

- seq.text:

  Numeric value for wildtype sequence text size. Default is \`1.1\`

- c.pallete:

  Character string or vector defining the color palette. Default is
  \`'RdBu'\`.

- c.scale:

  List of parameters definning the score color scale.

- ht:

  Numeric value for the height of the saved plot. Default is \`11\`.

- wd:

  Numeric value for the width of the saved plot. Default is \`8.5\`.

- name:

  Character string specifying the base name of the saved file.

- savepdf:

  Logical indicating whether to also save a PDF version of the plot.
  Default is \`TRUE\`.

- savesvg:

  Logical indicating whether to also save a SVG version of the plot.
  Default is \`TRUE\`.

- show:

  Logical indicating whether or not to display the plot in the viewer.
  Default is \`FALSE\`.

- factor_score:

  boolean for plotting factor variable with up to 6 levels

- discovery_score:

  boolean for plotting vector with "LOF", "Not significant", and "GOF"
  labels

- compare_score:

  boolean for plotting labels with format "\<+/-/=\>\<GOF/LOF\>" or
  "agree" to compare discoveries

- category_score:

  boolean for plotting score column as factor

- cat_name:

  name of category plotted in factor_score, discovery_score, or
  compare_score

## Examples

``` r
if (FALSE) { # \dontrun{
lilace_score_heatmap(scores, output_dir, score.col="effect", name="score_heatmap",
                     x.text=4, seq.text=1.5, y.text=3)
} # }
```
