# Functions for summarizing and plotting Lilace output
#' @import dplyr
#' @import ggplot2
#' @import colorspace
#' @import ggpubr
#' @importFrom cowplot theme_cowplot
#' @import utils
#' @importFrom stats na.omit
#'
NULL



#' Generate score heatmap
#'
#' @description
#' generates an amino acid by position heatmap for input dataset and given column (derived from Rosace scoreHeatmap function)
#'
#' @param data Scores data frame. Expected to have columns containing information about
#' position, control amino acid, mutated amino acid, mutation type, and score.
#' @param savedir Character string specifying the directory to save plots.
#' @param ctrl.name String specifying the control mutation name. Default is `synonymous`.
#' @param pos.col Column name in `data` for mutation positions. Default is `position`.
#' @param wt.col Column name in `data` for wildtype amino acids. Default is `wildtype`.
#' @param mut.col Column name in `data` for mutated amino acids. Default is `mutation`.
#' @param type.col Column name in `data` for mutation types. Default is `type`.
#' @param score.col Column name in `data` for mutation scores. Default is `mean`.
#' @param aa.order Character vector defining the order of amino acid mutations in the y-axis.
#' Default to using mutations from `data` in alphabetical order.
#' @param npos Integer specifying the number of positions per subplot. Default is `100`.
#' @param ncol Integer specifying the number of columns of subplots. Default is `1`.
#' @param pos.step Integer specifying the steps between x-axis labels. Default is `5`.
#' @param x.text Numeric value for x-axis text size. Default is `6`.
#' @param y.text Numeric value for y-axis text size. Default is `3`.
#' @param seq.text Numeric value for wildtype sequence text size. Default is `1.1`
#' @param c.pallete Character string or vector defining the color palette. Default is `'RdBu'`.
#' @param c.scale List of parameters definning the score color scale.
#' @param ht Numeric value for the height of the saved plot. Default is `11`.
#' @param wd Numeric value for the width of the saved plot. Default is `8.5`.
#' @param name Character string specifying the base name of the saved file.
#' @param savepdf Logical indicating whether to also save a PDF version of the plot.
#' Default is `TRUE`.
#' @param savesvg Logical indicating whether to also save a SVG version of the plot.
#' Default is `TRUE`.
#' @param show Logical indicating whether or not to display the plot in the viewer.
#' Default is `FALSE`.
#' @param factor_score boolean for plotting factor variable with up to 6 levels
#' @param discovery_score boolean for plotting vector with "LOF", "Not significant", and "GOF" labels
#' @param compare_score boolean for plotting labels with format "<+/-/=><GOF/LOF>" or "agree" to compare discoveries
#' @param category_score boolean for plotting score column as factor
#' @param cat_name name of category plotted in factor_score, discovery_score, or compare_score
#' @return NULL
#' @examples
#' \dontrun{
#' lilace_score_heatmap(scores, output_dir, score.col="effect", name="score_heatmap",
#'                      x.text=4, seq.text=1.5, y.text=3)
#' }
#'
#' @export
lilace_score_heatmap <- function(data,
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
                         c.pallete = 'RdBu',
                         c.scale = list(),
                         ht = 11,
                         wd = 8.5,
                         name = "Heatmap",
                         savepdf = TRUE,
                         savesvg = FALSE,
                         show = FALSE,
                         factor_score =FALSE,
                         discovery_score=FALSE,
                         compare_score=FALSE,
                         category_score=FALSE,
                         cat_name="discovery"
){

  # determine certain plot properties
  if (is.na(aa.order)) {
    aa.order = unique(na.omit(data[[mut.col]]))
  }
  pos.order <- levels(factor(data[[pos.col]]))
  # pos.order <- min(data[[pos.col]]):max(data[[pos.col]])
  npanel <- length(pos.order) %/% npos +1
  nrow <- npanel / ncol
  c.default <- list(palette = c.pallete, mid = 0, rev=TRUE, na.value = '#E5E5E5')
  c.args <- modifyList(c.default, c.scale)

  # parse the positions
  starts <- seq(1, length(pos.order), by = npos)
  if (length(starts) > 1) {
    ends <- c((starts - 1)[2:length(starts)], length(pos.order))
  } else {
    ends <- length(pos.order)
  }

  plot_list <- lapply(1:length(starts), function(i) {

    start <- starts[i]
    end <- ends[i]
    legend = ifelse(i == length(starts), 'right', 'none')

    # subset data for positions in range
    sub_data <- data[data[[pos.col]] %in% pos.order[start:end], ]
    sub_pos.order <- levels(factor(sub_data[[pos.col]]))

    # create subplot heatmaps
    if (factor_score) {
      colorscale <- c("1" = "#5778a4", "2" = "#e49444", "3" = "#6a9f58", "4" = "#967662", "5" = "#a87c9f", "6" = "#d62727")
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]),
                                       y = factor(.data[[mut.col]], levels = aa.order), fill = factor(.data[[score.col]]))) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (discovery_score) {
      colorscale <- c("LOF" = "#e49444", "Not significant" = "lightgray", "GOF" = "#5778a4", "4" = "#967662", "5" = "#a87c9f", "6" = "#d62727")
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]),
                                       y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (compare_score) {
      # colorscale <- c("both" = "darkgray", "Lilace" = "navy", "original approach" = "maroon", "neither" = "lightgray", "5" = "#a87c9f", "6" = "#d62727")
      colorscale <- c("agree" = "darkgray", "+GOF" = "blue", "-GOF" = "navy", "-LOF" = "maroon", "+LOF" = "red", "=LOF" = "pink", "=GOF" = "lightblue")

      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]),
                                       y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (category_score) {
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]),
                                       y = factor(.data[[mut.col]], levels = aa.order), fill = as.factor(.data[[score.col]]))) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_brewer(palette = "Set3")
    }
    else {
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]),
                                       y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        do.call(scale_fill_continuous_divergingx, c.args)
    }
    p <- p +
      scale_color_manual(values = c(NA,'darkgreen'), name = ctrl.name, guide = "none") +
      scale_x_discrete(breaks = sub_pos.order[seq(1, length(sub_pos.order), by = pos.step)]) +
      coord_fixed(ratio = 1, clip = "off") +
      theme_cowplot() +
      theme(
        plot.margin = unit(c(1, 0.25, 0.25, 0.25), "lines"),
        axis.title.y = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = x.text),
        axis.text.y = element_text(size = y.text),
        axis.ticks = element_blank(),
        line = element_blank()
      ) +
      geom_text(data = sub_data,
                aes(x = factor(.data[[pos.col]], levels=sort(pos.order)), label = factor(.data[[wt.col]]), y = Inf),
                vjust = -1, check_overlap = TRUE, size = seq.text) +
      labs(y = "Mutation", x = "Position")
    return(p)
  })

  p_all <- do.call(ggarrange, c(plot_list, list(nrow = nrow, ncol = ncol, common.legend=T, legend="right")))

  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }

  ggsave(file.path(savedir, paste0(name,".png")), plot = p_all, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p_all, height = ht, width = wd)
  }
  if (savesvg) {
    ggsave(file.path(savedir, paste0(name,".svg")), plot = p_all, height = ht, width = wd)
  }

  if (show) {
    cat(paste("Showing the first", as.character(npos) , "positions.",
              "Full figure can be found in the saved directory."))
    print(plot_list[[1]])
  }
  return(p_all)
}

#' Generate score histograms/density plots by mutation type
#'
#' @description
#' generates a density plot for visualizing the distribution of scores across different mutation types (same as Rosace scoreDensity function)
#' @param data Scores data frame. Expected to have columns containing information about
#' position, control amino acid, mutated amino acid, mutation type, and score.
#' @param savedir Character string specifying the directory to save plots.
#' @param type.col Column name in `data` for mutation types. Default is `type`.
#' @param score.col Column name in `data` for mutation scores. Default is `mean`.
#' @param hist Logical indicating whether to plot count histogram or density.
#' Default is `FALSE`.
#' @param nbins Numeric value specifying the number of bins for the histogram.
#' Default is `30`.
#' @param c.fill Vector indicating the fill color for mutation types.
#' @param alpha Numeric value between 0-1 indicating the fill transparency.
#' Default is `0.5`
#' @param x.text Numeric value for x-axis text size. Default is `10`.
#' @param y.text Numeric value for x-axis text size. Default is `10`.
#' @param scale.free Logical indicating whether to make score range proportional.
#' Default is `FALSE`.
#' @param space.free Logical indicating whether to make panel heights variable.
#' Default is `FALSE`.
#' @param ht Numeric value for the height of the saved plot. Default is `10`.
#' @param wd Numeric value for the width of the saved plot. Default is `8`.
#' @param name Character string specifying the base name of the saved file.
#' @param savepdf Logical indicating whether to also save a PDF version of the plot.
#' Default is `TRUE`.
#' @param savesvg logical indicating whether to save as svg
#' @param show Logical indicating whether or not to display the plot in the viewer.
#' Default is `TRUE`.
#' @return NULL.
#' @examples
#' \dontrun{
#' lilace_score_density(scores, output_dir, score.col="effect",
#'                      name="score_histogram", hist=T, scale.free=T)
#' }
#' @export
lilace_score_density <- function(data,
                         savedir,
                         type.col = "type",
                         score.col = "mean",
                         hist = FALSE,
                         nbins = 30,
                         c.fill = c('#FF7575', 'lightgreen', "#7298BF"),
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
){

  sc <- ifelse(scale.free, "free_y", "fixed")
  sp <- ifelse(space.free, "free_y", "fixed")
  if (length(c.fill) != length(levels(factor(data[[type.col]])))) {
    c.fill = c('#FF7575', "orange", 'lightgreen', "#7298BF", "#8828a8")
    warning("Length of color vector does not match the number of mutation types.")
  }

  p <- ggplot(data, aes(x = .data[[score.col]], fill = .data[[type.col]])) +
    facet_grid(.data[[type.col]] ~ ., scales = sc, space = sp) +
    theme_minimal() +
    scale_fill_manual(values = c.fill) +
    theme(
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(size = x.text),
      axis.text.y = element_text(size = y.text)
    )

  if (hist) {
    p <- p + geom_histogram(aes(alpha = 0.5), bins = nbins, position="dodge") +
      labs(x = "Score", y = "Count") +
      scale_alpha(guide = "none")
  }
  else{
    p <- p + geom_density(alpha = alpha) + labs(x = "Score", y = "Density")
  }

  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }

  ggsave(file.path(savedir, paste0(name,".png")), plot = p, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p, height = ht, width = wd)
  }
  if (savesvg) {
    ggsave(file.path(savedir, paste0(name,".svg")), plot = p, height = ht, width = wd)
  }

  if (show) {
    print(p)
  }
}




