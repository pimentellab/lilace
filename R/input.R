# input processing functions to setup for lilace
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#'
NULL

#' Create Lilace object from counts
#'
#' @description
#' Formats required information for Lilace to run into a Lilace object. See intro vignette for example usage.
#'
#' @param variant_id a length V vector of unique variant identifiers such as hgvs nomenclature
#' @param mutation_type a vector of mutation types (e.g. synonymous, missense, etc).
#'  The label used to identify negative control label should be included in this.
#' @param position a vector of residue positions
#' @param replicate a vector of replicate ids
#' @param counts a (V x K) matrix of bin counts, where K is the number of bins
#' @param metadata a matrix or dataframe of any additional variant information to be incorporated into the final lilace output (e.g. wildtype residue)
#' @returns a lilace object, with data stored in $data
#' @export
lilace_from_counts <- function(variant_id, mutation_type, position, replicate, counts, metadata) {
  colnames(counts) <- paste0("c_", 1:ncol(counts)-1)
  # TODO: add input checks
  # all have same length as a variant_id
  data <- data.frame(variant=variant_id, metadata, type=mutation_type, position=position, rep=replicate, counts)
  # add count info
  data$n_counts <- apply(data %>% dplyr::select(starts_with("c_")), 1, sum)
  data <- data %>% group_by(.data$variant) %>% mutate(total_counts=sum(.data$n_counts))
  data <- data %>% arrange(.data$position, .data$mutation)
  # create lilace object
  lilace_obj <- list(data=data)
  lilace_obj$metadata_cols <- colnames(metadata)
  return(lilace_obj)
}


#' Create Lilace object from enrich processed counts output
#'
#' @param file enrich counts file (see end of intro vignette for example input format)
#' @param pheno name of phenotype ("condition" row in counts table)
#' @returns a lilace object, with data stored in $data
#' @export
lilace_from_enrich <- function(file, pheno="abundance") {
  message("Formatting from enrich2 output. Synonymous mutations will be used for negative control-based bias correction. If this is not desirable, please input using the format_data() function.")
  data_raw <- read_tsv(file, skip = 4, col_names = FALSE, show_col_types = FALSE)
  title <- read_tsv(file, n_max = 3, col_names = FALSE, show_col_types = FALSE)
  colnames(data_raw) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(data_raw)[1] <- "hgvs"
  rm(title)

  data <- data_raw %>%
    pivot_longer(!.data$hgvs) %>%
    separate("name", into = c("exp", "rep", "time"), sep = "[.]") %>%
    separate("rep", into = c("exp2", "rep"), sep = -2) %>%
    dplyr::select(-.data$exp2) %>%
    pivot_wider(names_from = .data$time, values_from = .data$value)

  data <- data %>% filter(.data$hgvs != "_wt")

  data <- data %>%
    mutate(
      hgvs2 = substr(.data$hgvs, 4, nchar(.data$hgvs) - 1),
      position = gsub("\\_.*", "", .data$hgvs2),
      position = as.numeric(gsub("[[:alpha:]]", "", .data$position)),
      wildtype = substr(.data$hgvs2, 1, 1),
      hgvs2 = substr(.data$hgvs2, 2, nchar(.data$hgvs2))
    ) %>%
    separate("hgvs2", sep = "[_]", into = c("a1", "a2"), remove = FALSE) %>%
    mutate(mutation = ifelse(is.na(.data$a2), gsub("[[:digit:]]", "", .data$hgvs2), NA)) %>%
    mutate(a1 =  as.numeric(gsub("[[:alpha:]]", "", .data$a1)),
           a2 = as.numeric(gsub("[[:alpha:]]", "", .data$a2)),
           a = .data$a2 - .data$a1, a = ifelse(is.na(.data$a), 0, .data$a)) %>%
    mutate(mutation = ifelse(str_detect(.data$hgvs, "del"), paste("del", .data$a, sep = ""), .data$mutation)) %>%
    dplyr::select(-.data$a, -.data$a1, -.data$a2) %>%
    mutate(ins = gsub(".*ins", "", .data$hgvs2),
           mutation = ifelse(str_detect(.data$hgvs, "ins"), paste("ins", nchar(.data$ins)-1, sep = ""), .data$mutation)) %>%
    dplyr::select(-.data$ins, -.data$hgvs2)


  func_map <- function(wt, mut) {
    if (wt == mut) {
      return("synonymous")
    } else if (str_detect(mut, "del")) {
      return("deletion")
    } else if (str_detect(mut, "ins")) {
      return("insertion")
    } else {
      return("missense")
    }
  }

  data <- data %>%
    rowwise() %>%
    mutate(type = func_map(.data$wildtype, .data$mutation)) %>%
    ungroup() %>%
    arrange(.data$position) %>%
    dplyr::select(1, all_of(8:11), all_of(4:7), all_of(2:3))
  colnames(data)[colnames(data) == "hgvs"] <- "variant"
  # add count info
  data <- data[data$exp==pheno,]
  data$n_counts <- apply(data %>% dplyr::select(starts_with("c_")), 1, sum)
  data <- data %>% group_by(.data$variant) %>% mutate(total_counts=sum(.data$n_counts))
  # create lilace object
  lilace_obj <- list(data=data)
  return(lilace_obj)
}

#' Create Lilace object from separate count files for each bin and replicate
#'
#' @param file_list a list containing a separate list of bin count files for each replicate  (see end of intro vignette for example)
#' @param variant_id_col column name for variant id
#' @param position_col column name for position information
#' @param mutation_type_col column name for mutation type information
#' @param count_col column name for counts
#' @param delim the field separator character (see ?readr::read.delim)
#' @returns a lilace object, with data stored in $data
#' @export
lilace_from_files <- function(file_list, variant_id_col="hgvs", position_col="position", mutation_type_col="type", count_col="count", delim="\t") {
  data <- c()
  for (rep in 1:length(file_list)) {
    rep_list <- file_list[[rep]]
    rep_df <- c()
    n_bins <- length(rep_list)
    for (bin in 1:n_bins) {
      file <- rep_list[[bin]]
      counts <- read.delim(file)
      counts[[paste0("c_", bin-1)]] <- counts[[count_col]]
      counts[[count_col]] <- NULL
      if (is.null(rep_df)) {
        rep_df <- counts
      } else {
        rep_df <- merge(rep_df, counts)
      }
    }
    rep_df$rep <- paste0("R", rep)
    data <- rbind(data, rep_df)
  }
  lilace_obj <- lilace_from_counts(variant_id=data[[variant_id_col]],
                                   mutation_type=data[[mutation_type_col]],
                                   position=data[[position_col]],
                                   replicate=data$rep,
                                   counts=data %>% dplyr::select(starts_with("c_")),
                                   metadata=data %>% dplyr::select(-c(variant_id_col, mutation_type_col, position_col, rep, starts_with("c_"))))
}

#' Normalize data to cell sorting proportions
#'
#' @description
#' Scales the read counts in a bin to match the sort proportions given.
#'
#' @param lilace_obj a lilace obj, such as created using lilace_from_counts()
#' @param sort_prop a list of cell sorting bin proportions to normalize
#'  If replicate specific, input as list(<rep_label>=c(0.25, 0.25, 0.25, 0.25)) where <rep_label> corresponds to the lilace object replicate labels
#'  If not replicate specific, input as vector.
#' @param rep_specific boolean indicating whether sort_prop is replicate specific or not
#' @returns lilace object with $normalized_data element
#' @examples
#' # Example 1 (not replicate specific)
#' \dontrun{
#' sort_prop <- c(0.2184, 0.1972, 0.1869, 0.1236)
#' lilace_obj <- lilace_sorting_normalize(lilace_obj, sort_prop, rep_specific=F)
#' }
#'
#' # Example 2 (replicate specific)
#' \dontrun{
#' sort_prop=list(R1=c(0.2184, 0.1972, 0.1869, 0.1236),
#'                R2=c(0.25, 0.25, 0.25, 0.25),
#'                R3=c(0.3, 0.2, 0.3, 0.2))
#' lilace_obj <- lilace_sorting_normalize(lilace_obj, sort_prop, rep_specific=T)
#' }
#' @export
lilace_sorting_normalize <- function(lilace_obj, sort_prop, rep_specific) {
  data <- lilace_obj$data %>% ungroup()
  rep_set <- unique(data$rep)
  # check if rep specific that sort prop labels correspond to lilace_obj$data labels
  if (rep_specific) {
    for (rep in rep_set) {
      if (!(rep %in% names(sort_prop))) {
        stop(paste0("Replicate ", rep, " not in given sorting proportions."))
      }
    }
  }
  # check actually proportions
  if (!rep_specific && sum(sort_prop) != 1) {
    warning(paste0("Sorting proportions sum up to ", sum(sort_prop), " instead of 1."))
    sort_prop <- sort_prop / sum(sort_prop) # normalize to actual percents
  }

  # normalization
  K <- ncol(data %>% dplyr::select(starts_with("c_")))
  for (rep in rep_set) {
    if (rep_specific && sum(sort_prop[[rep]]) != 1) {
      warning(paste0("Sorting proportions sum up to ", sum(sort_prop[[rep]]), " instead of 1."))
      sort_prop[[rep]] <- sort_prop[[rep]] / sum(sort_prop[[rep]]) # normalize to actual percents
    }
    rep_total <- sum(data[data$rep==rep,paste0("c_", (1:K)-1)])
    for (i in 1:K) {
      bin <- paste0("c_", 0:(K-1))[i]
      if (rep_specific) {
        rescaled_bin_total <- (sort_prop[[rep]] / sum(sort_prop[[rep]])) * rep_total
      } else {
        rescaled_bin_total <- sort_prop * rep_total
      }
      data[data$rep==rep,][[bin]] <- ceiling(data[data$rep==rep,][[bin]] / sum(data[data$rep==rep,][[bin]]) * rescaled_bin_total[i])
    }
  }
  # recalculate count info
  data$n_counts <- apply(data %>% dplyr::select(starts_with("c_")), 1, sum)
  data <- data %>% group_by(.data$variant) %>% mutate(total_counts=sum(.data$n_counts))
  lilace_obj$normalized_data <- data
  return(lilace_obj)
}






