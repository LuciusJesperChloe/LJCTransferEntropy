#' Improved Transfer Entropy Analysis with Rolling Window
#'
#' @param x Time series data for variable x.
#' @param y Time series data for variable y.
#' @param lx Embedding dimension for x.
#' @param ly Embedding dimension for y.
#' @param q Quantile for Transfer Entropy.
#' @param entropy Type of entropy to use.
#' @param shuffles Number of shuffles for significance testing.
#' @param type Type of binning.
#' @param quantiles Quantiles for binning.
#' @param bins Number of bins.
#' @param limits Custom limits for binning.
#' @param nboot Number of bootstrap iterations.
#' @param burn Burn-in period for bootstrap.
#' @param quiet Suppress output.
#' @param seed Random seed for reproducibility.
#' @param na.rm Remove missing values.
#' @param eff_entropy Use efficient entropy calculation.
#' @param cr Criterion for stationarity test.
#'
#' @return A data frame containing Transfer Entropy results for each window.
#' @export

LJCTransferEntropy <- function(x,
                               y,
                               lx = 1,
                               ly = 1,
                               q = 0.1,
                               entropy = "Shannon",
                               shuffles = 100,
                               type = "quantiles",
                               quantiles = c(5, 95),
                               bins = NULL,
                               limits = NULL,
                               nboot = 300,
                               burn = 50,
                               quiet = NULL,
                               seed = NULL,
                               na.rm = TRUE,
                               eff_entropy = TRUE,
                               cr = "AIC") {

  # Prepare data
  start_year <- Data$Year[1]
  end_year <- Data$Year[nrow(Data)]
  end_window_size <- (end_year - start_year + 1)
  years <- start_year:end_year

  # Set parameters for incrementing window sizes
  start_window_size <- 10
  rolling_step_size <- 1

  # Initialize a list to store results for each window size
  results_list <- list()

  # Initialize progress bar
  total_iterations <- sum(end_window_size - start_window_size + 1) # Total iterations across all windows
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

  iteration_count <- 0

  # Loop through different window sizes
  for (window_size in start_window_size:end_window_size) {
    temp_results <- data.frame(Start_Year = numeric(),
                               End_Year = numeric(),
                               TE_x_to_y = numeric(),
                               TE_y_to_x = numeric(),
                               Pval_x_to_y = numeric(),
                               Pval_y_to_x = numeric(),
                               Stat_x = numeric(),
                               Stat_y = numeric(),
                               Window_size = numeric())

    for (i in seq(1, length(years) - window_size + 1, by = rolling_step_size)) {
      windowed_x <- x[i:(i + window_size - 1)]
      windowed_y <- y[i:(i + window_size - 1)]

      set.seed(seed)

      # Perform Transfer Entropy calculation
      te_results <- transfer_entropy(windowed_x, windowed_y, lx, ly, q, entropy, shuffles, type, quantiles, bins, limits, nboot, burn, quiet, seed, na.rm)

      TE_x_to_y <- te_results$coef["X->Y", "te"]
      TE_y_to_x <- te_results$coef["Y->X", "te"]

      # Perform ADF Test for Stationarity
      adf_x <- CADFtest(windowed_x, type = "drift", criterion = cr)
      adf_y <- CADFtest(windowed_y, type = "drift", criterion = cr)

      temp_results <- rbind(temp_results, data.frame(
        Window_size = window_size,
        Start_Year = years[i],
        End_Year = years[i + window_size - 1],
        TE_x_to_y,
        TE_y_to_x,
        Pval_x_to_y = te_results$coef["X->Y", "p-value"],
        Pval_y_to_x = te_results$coef["Y->X", "p-value"],
        Stat_x = adf_x$p.value,
        Stat_y = adf_y$p.value
      ))

      # Update progress bar
      iteration_count <- iteration_count + 1
      setTxtProgressBar(pb, iteration_count)
    }

    results_list[[as.character(window_size)]] <- temp_results
  }

  # Close progress bar
  close(pb)

  # Combine all results into one data frame
  final_output <- do.call(rbind, results_list)

  # return
  return(final_output)
}
