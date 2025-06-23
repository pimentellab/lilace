# Functions to run Lilace
#' @import dplyr
#' @import tidyr
#' @import cmdstanr

#' @importFrom yaml write_yaml
#' @importFrom readr read_tsv
#' @importFrom readr write_tsv
#' @importFrom stringr str_detect
NULL


#' filter minimum total counts and zero count rows
#' @param data data
#' @param min_total_counts minimum total counts for variant
#' @export
.filter_data_low_counts <- function(data, min_total_counts=15) {
  # remove 0 count rows and filter to min total counts
  n_zero_removed <- nrow(data) - nrow(data[data$n_counts > 0,])
  if (n_zero_removed > 0) {
    message(paste("Removing", nrow(data) - nrow(data[data$n_counts > 0,]), "zero count rows"))
    data <- data[data$n_counts > 0,]
  }
  n_less_15 <- sum(data$total_counts<=min_total_counts)
  if (n_less_15 > 0) {
    message(paste("Removing", n_less_15, "observations due to less than", min_total_counts, "total counts in the variant."))
    data <- data[data$total_counts>=min_total_counts,]
  }
  return(data)
}

#' format model input (intended for internal use)
#' @param data data
#' @param control_label control label
#' @param pseudocount logical for pseudocount
#' @export
.generate_model_input <- function(data, control_label, pseudocount=F) {
  # check controls
  if (!(control_label %in% data$type)) {
    stop(paste0("No negative controls with label ", control_label, " in mutation type argument.
                  You can specify negative control label with control_label=<label>."))
  }
  if (sum(data$type == control_label) < 10) {
    warning("Less than 10 negative controls present. Negative control-based bias correction may not perform optimally.")
  }
  for (col in c("variant", "type", "position", "rep", "n_counts", "total_counts")) {
    if (!(col %in% colnames(data))) {
      stop(paste0("Missing column: ", col))
    }
  }

  # add pseudocount
  if (pseudocount) {
    data <- data %>% mutate(across(starts_with("c_"), ~ .x + 1))
  }

  # create stan input
  V <- length(unique(data$variant))
  N <- nrow(data)
  nMAPv <- as.numeric(factor(data$variant))
  n_counts <- as.numeric(data$n_counts)
  n_syn_counts <- as.numeric(data[data$type=="synonymous",]$n_counts)
  distinct_df <- (distinct(data.frame(variant = data$variant, nMAPv = nMAPv,
                                      total_counts = data$total_counts, type=data$type)) %>%
                    arrange(nMAPv))

  vMAPp <- distinct(data.frame(variant = data$variant, nMAPv = nMAPv, type = data$type, position = data$position)) %>%
    arrange(nMAPv)
  vMAPp$position <- as.numeric(vMAPp$position)

  # group positions if needed
  positions <- sort(unique(vMAPp$position))
  for (i in 1:length(positions)) {
    if (i != length(positions)) {
      if (sum(vMAPp$position == positions[i]) < 10) {
        message("Grouping positions ", positions[i], " with ", positions[i+1], " due to less than 10 variants at position.")
        vMAPp[vMAPp$position == positions[i],]$position <- positions[i+1]
      }
    }
  }
  # ensure end also has 10 variants
  while (sum(vMAPp$position == positions[i]) < 10) {
    message("Grouping positions ", positions[i], " with ", positions[i-1], " due to less than 10 variants at position.")
    vMAPp[vMAPp$position == positions[i],]$position <- positions[i-1]
    i <- i - 1
  }

  # encode negative control positions
  vMAPp[vMAPp$type=="synonymous",]$position <- 0 # encode synonymous as diff position
  S <- sum(vMAPp$type=="synonymous")
  vMAPs <- rep(-1, nrow(vMAPp)) # -1 for non syn, o.w. each syn gets own index
  vMAPs[vMAPp$type=="synonymous"] <- 1:S

  vMAPp <- as.numeric(as.factor(vMAPp$position))
  P <- length(unique(vMAPp))

  R <- length(unique(data$rep))
  nMAPr <- as.numeric(factor(data$rep))
  sMAPr <- nMAPr[data$type=="synonymous"]
  N_syn <- length(sMAPr)

  y <- as.matrix(data %>% dplyr::select(starts_with("c_")))
  y_syn <- y[data$type=="synonymous",]
  K <- ncol(y)

  input <- list(V = V, S = S, N = N, N_syn=N_syn, P=P, R=R, nMAPv = nMAPv, vMAPs = vMAPs, nMAPr = nMAPr, vMAPp = vMAPp, sMAPr=sMAPr,
                n_counts=n_counts, n_syn_counts=n_syn_counts, K = K, y = y, y_syn = y_syn)
  return(input)
}

#' compute negative control bias correction (intended for internal use)
#' @param mu mu posterior samples
#' @param data data
#' @param input stan model input parameters
#' @param control_label control label
#' @export
.correct_negative_control <- function(mu, data, input, control_label) {
  syn_mu <- as.matrix(mu[,unique(input$nMAPv[data$type==control_label])])
  syn_indices <- t(replicate(nrow(mu), sample.int(ncol(syn_mu), ncol(mu), replace=T)))
  mu2 <- matrix(NA, nrow=nrow(mu), ncol=ncol(mu))
  for (i in 1:nrow(mu)) {
    mu2[i,] <- as.numeric(mu[i,]) - as.numeric(syn_mu[i, syn_indices[i,]])
  }
  colnames(mu2) <- colnames(mu)
  return(mu2)
}

.get_mu_summary <- function(mu, V) {
  effect <- apply(mu, 2, mean)
  effect_se <- apply(mu, 2, stats::sd)
  lfsr <- apply(mu, 2, function(samples) {
    p_less <- sum(samples <= 0) / length(samples)
    p_great <- sum(samples >= 0) / length(samples)
    min(c(p_less, p_great))
  })
  mu_stats <- data.frame(t(rbind(effect, effect_se, lfsr)))
  mu_stats$param <- rownames(mu_stats)
  mu_stats <- mu_stats %>% filter(substr(.data$param, 1, 2) == "mu")
  mu_stats$map <- 1:V
  return(mu_stats)
}

.get_pos_summary <- function(fit, P) {
  theta <- fit$draws(variables = "theta", format = "df")
  sigma <- fit$draws(variables = "sigma", format = "df")
  pos_mean <- apply(theta, 2, mean)
  pos_sd <- apply(theta, 2, stats::sd)
  sigma_mean <- apply(sigma, 2, mean)
  sigma_sd <- apply(sigma, 2, stats::sd)
  pos_lfsr <- apply(theta, 2, function(samples) {
    p_less <- sum(samples <= 0) / length(samples)
    p_great <- sum(samples >= 0) / length(samples)
    min(c(p_less, p_great))
  })
  pos_stats <- data.frame(t(rbind(pos_mean, pos_sd)))
  pos_stats$pos_param <- rownames(pos_stats)
  pos_stats <- pos_stats %>% filter(substr(.data$pos_param, 1, 5) == "theta")
  pos_stats$p_map <- 2:P
  return(pos_stats)
}

#' get phi parameter summary (intended for internal use)
#' @param fit stan output fit object
#' @export
.get_phi_summary <- function(fit) {
  a <- fit$draws(variables = "a", format = "matrix")
  b <- fit$draws(variables = "b", format = "matrix")
  a_mean <- apply(a, 2, mean)
  a_sd <- apply(a, 2, stats::sd)
  b_mean <- apply(b, 2, mean)
  b_sd <- apply(b, 2, stats::sd)
  phi_stats <- data.frame(t(rbind(a_mean, a_sd, b_mean, b_sd)))
  phi_stats$param <- rownames(phi_stats)
  phi_stats <- phi_stats %>% filter(substr(.data$param, 1, 1) == "a" | substr(.data$param, 1, 1) == "b")
  return(phi_stats)
}

#' summarize posterior samples (intended for internal use)
#' @param fit stan output fit object
#' @param data data
#' @param input stan input list
#' @param control_correction logical for whether to do negative control correction
#' @param control_label label for negative control variant types
#' @param use_positions logical for whether positions used in running Lilace
#' @export
.summarize_posteriors <- function(fit, data, input, control_correction, control_label, use_positions) {
  # get effect size summary
  mu <- fit$draws(variables = "mu", format = "matrix")
  # apply negative control correction
  if (control_correction) {
    mu <- .correct_negative_control(mu, data, input, control_label)
  }
  mu_stats <- .get_mu_summary(mu, input$V)
  # get pos summary
  if (use_positions) {
    pos_stats <- .get_pos_summary(fit, input$P)
    mu_stats$p_map <- input$vMAPp
    mu_stats <- distinct(mu_stats %>% left_join(pos_stats))
  }
  # join with data
  n_pos <- data.frame(variant = data$variant,
                      rep = data$rep,
                      map = input$nMAPv)
  n_pos <- distinct(n_pos %>% left_join(mu_stats))
  data <- data %>% left_join(n_pos)
  return(data)
}

#' get variant scores dataframe
#' @param lilace_obj fitted Lilace object
#' @param use_positions whether Lilace was run with use_positions=T
#' @export
.get_score_df <- function(lilace_obj, use_positions) {
  if (use_positions) {
    score_df <- lilace_obj$fitted_data %>% select(c(.data$variant, .data$type,
                                                    lilace_obj$metadata_cols, .data$position,
                                                    .data$effect, .data$effect_se, .data$lfsr,
                                                    .data$pos_mean, .data$pos_sd))
  } else {
    score_df <- lilace_obj$fitted_data %>% select(c(.data$variant, .data$type, lilace_obj$metadata_cols,
                                                    .data$position, .data$effect,
                                                    .data$effect_se, .data$lfsr))
  }
  score_df <- distinct(score_df)
  score_df$position <- as.numeric(score_df$position)
  score_df$discovery05 <- sign(score_df$effect) * as.numeric(score_df$lfsr < 0.05)
  score_df <- score_df[order(score_df$position, score_df$type),]
  lilace_obj$scores <- score_df
  return(lilace_obj)
}



#' Run Lilace
#'
#' @description
#' Runs Lilace on given data. Lilace will run on $normalized_data if it exists, otherwise it will use $data.
#'
#' @param lilace_obj initialized lilace object
#' @param output_dir output directory to write scores and sampling logs to
#' @param control_label label from $data$type column to use as negative controls to score against
#' @param control_correction boolean for whether to use negative control scores as bias correction
#' @param use_positions boolean for whether to use position hierarchy to improve estimation
#' @param pseudocount boolean for whether to add pseudocount (+1 to all counts) for fitting model
#' @param seed random seed for sampling process to get exactly reproducible results. A NULL value indicates no fixed seed.
#' @param min_total_counts minimum total counts for a variant--anything less will be filtered out
#' @param n_parallel_chains number of chains to run in parallel
#' @returns lilace object with $scores and $fitted_data entries
#' @examples
#' \dontrun{
#' lilace_obj <- lilace_fit_model(lilace_obj, output_dir, control_label="synonymous",
#'                                control_correction=T, use_positions=T, pseudocount=T)
#' }
#' @export
lilace_fit_model <- function(lilace_obj, output_dir, control_label="synonymous", control_correction=TRUE,
                             use_positions=TRUE, pseudocount=TRUE, seed=NULL, min_total_counts=15, n_parallel_chains=4) {
  if (!is.null(lilace_obj$normalized_data)) {
    data <- lilace_obj$normalized_data %>% ungroup()
    message("Running Lilace on sorting normalized counts")
  } else {
    data <- lilace_obj$data %>% ungroup()
    message("Running Lilace on input counts")
  }
  if (!dir.exists(output_dir)) {
    warning(paste0("Output directory ", output_dir, " does not exist. Creating...\n"))
    dir.create(output_dir, recursive=T)
  }
  output_folder <- gsub("//", "/", file.path(output_dir, "lilace_output"))
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  data <- .filter_data_low_counts(data, min_total_counts)

  input <- .generate_model_input(data, control_label, pseudocount)
  if (use_positions) {
    model_file <- system.file("stan", "lilace.stan", package="lilace")
    init <- list(list(sigma=rep(1, input$P-1)), list(sigma=rep(1, input$P-1)), list(sigma=rep(1, input$P-1)), list(sigma=rep(1, input$P-1)))
  } else {
    model_file <- system.file("stan", "lilace_nopos.stan", package="lilace")
    init <- list(list(sigma=rep(1, input$V)), list(sigma=rep(1, input$V)), list(sigma=rep(1, input$V)), list(sigma=rep(1, input$V)))
  }
  mod <- cmdstan_model(model_file)
  start <- Sys.time()
  fit <- mod$sample(
    data = input,
    chains = 4,
    parallel_chains = n_parallel_chains,
    refresh = 250,
    seed = seed,
    init=init,
    show_messages=TRUE,
    show_exceptions=FALSE
  )
  end <- Sys.time()
  sink(file.path(output_folder, "sampling.log"))
  print(fit$output())
  sink()
  # save posterior samples to output_dir
  samples <- fit$draws(format = "df")

  saveRDS(samples, file = file.path(output_folder, "posterior_samples.RData"))
  sampling_diagnostics <- fit$diagnostic_summary()
  difftime <- end - start
  runtime <- paste(difftime, attr(difftime, "units"))
  sampling_diagnostics$runtime <- runtime
  lilace_obj$sampling_diagnostics <- sampling_diagnostics
  write_yaml(sampling_diagnostics, file.path(output_folder, "/sampling_diagnostics.yaml"))

  # integrate scores with data and save to file
  fit_df <- .summarize_posteriors(fit, data, input, control_correction, control_label, use_positions)
  lilace_obj$fitted_data <- fit_df
  lilace_obj <- .get_score_df(lilace_obj, use_positions)
  write_tsv(lilace_obj$scores, file.path(output_folder, "/variant_scores.tsv"))

  closeAllConnections()
  return(lilace_obj)
}
