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
LJCTransferEntropy = function(x,
                              y,
                              Time,
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
                              cr = "AIC",
                              start_window_size = 10,
                              rolling_step_size = 1){

  # Check if Time is numeric
  if (!is.numeric(Time)) {
    # Replace Time with a numeric sequence
    Time = seq_along(x)
    message("Time series was not numeric; replaced with a sequence starting from 1.")
  }

  if (length(Time) != length(x) | length(Time) != length(y) | length(x) != length(y)) {
    stop("The lengths of x, y and Time must match!")
  }

  # Load necessary packages
  if (!requireNamespace("RTransferEntropy")) {
    install.packages("RTransferEntropy")
  }
  if (!requireNamespace("CADFtest")) {
    install.packages("CADFtest")
  }
  if (!requireNamespace("future")) {
    install.packages("future")
  }
  if (!requireNamespace("ggplot2")) {
    install.packages("ggplot2")
  }

  library(CADFtest)
  library(RTransferEntropy)
  library(future)
  library(ggplot2)
  library(gridExtra)
  library(grid)

  # Enable parallel processing
  plan(multisession)

  print("Analyzing...")

  # Prepare data
  start_time = Time[1]
  end_time = Time[nrow(Data)]
  end_window_size = (end_time - start_time + 1)
  times = start_time:end_time

  # Initialize a list to store results for each window size
  results_list = list()


  # Loop through different window sizes
  for (window_size in start_window_size:end_window_size) {
    temp_results = data.frame(Start_Time = numeric(),
                              End_Time = numeric(),
                              TE_x_to_y = numeric(),
                              TE_y_to_x = numeric(),
                              Pval_x_to_y = numeric(),
                              Pval_y_to_x = numeric(),
                              Stat_x = numeric(),
                              Stat_y = numeric(),
                              Window_size = numeric())

    for (i in seq(1, length(times) - window_size + 1, by = rolling_step_size)) {
      windowed_x = x[i:(i + window_size - 1)]
      windowed_y = y[i:(i + window_size - 1)]

      if (na.rm) {
        windowed_x = windowed_x[!is.na(windowed_x)]
        windowed_y = windowed_y[!is.na(windowed_y)]
      }

      if (length(windowed_x) < window_size || length(windowed_y) < window_size) {
        next
      }

      set.seed(seed)

      # Perform Transfer Entropy calculation
      tryCatch(
        expr = {
          if(eff_entropy){
            te_results = transfer_entropy(windowed_x, windowed_y, lx, ly, q, entropy, shuffles, type, quantiles, bins, limits, nboot, burn, quiet, seed, na.rm)
            TE_x_to_y = calc_ete(windowed_x, windowed_y, lx, ly, q, entropy, shuffles, type, quantiles, bins, limits, burn, seed, na.rm)
            TE_y_to_x = calc_ete(windowed_y, windowed_x, lx, ly, q, entropy, shuffles, type, quantiles, bins, limits, burn, seed, na.rm)

          }else{
            te_results = transfer_entropy(windowed_x, windowed_y, lx, ly, q, entropy, shuffles, type, quantiles, bins, limits, nboot, burn, quiet, seed, na.rm)
            TE_x_to_y = te_results$coef["X->Y", "te"]
            TE_y_to_x = te_results$coef["Y->X", "te"]
          }
        },
        error = function(e) {
          message(paste("Error in Transfer Entropy calculation:", e))
          next
        }
      )

      # Perform ADF Test for Stationarity
      tryCatch(
        expr = {
          adf_x = CADFtest(windowed_x, type = "drift", criterion = cr)
          adf_y = CADFtest(windowed_y, type = "drift", criterion = cr)
        },
        error = function(e) {
          message(paste("Error in ADF test:", e))
          next
        }
      )

      temp_results = rbind(temp_results, data.frame(
        Window_size = window_size,
        Start_Time = times[i],
        End_Time = times[i + window_size - 1],
        TE_x_to_y,
        TE_y_to_x,
        Pval_x_to_y = te_results$coef["X->Y", "p-value"],
        Pval_y_to_x = te_results$coef["Y->X", "p-value"],
        Stat_x = adf_x$p.value,
        Stat_y = adf_y$p.value
      ))
    }

    results_list[[as.character(window_size)]] = temp_results
  }

  # Combine all results into one data frame
  final_output = do.call(rbind, results_list)

  # Save results to a CSV file on the desktop without NA values appearing in Excel
  write.table(final_output, file = "C:/Users/User/Desktop/LJCTransferEntropy.csv", row.names = FALSE, col.names = TRUE, sep = ",", na = "")

  # Calculate average TE values for each time point across all windows
  avg_results = data.frame(
    Time = numeric(),
    Window_size = numeric(),
    Avg_TE_x_to_y = numeric(),
    Avg_TE_y_to_x = numeric(),
    Avg_Pval_x_to_y = numeric(),
    Avg_Pval_y_to_x = numeric(),
    Avg_Stat_x = numeric(),
    Avg_Stat_y = numeric()
  )

  for (window_size in start_window_size:length(Time)) {
    for (time_point in Time) {
      window_indices = which(final_output$Start_Time <= time_point & final_output$End_Time >= time_point & final_output$Window_size == window_size)

      if (length(window_indices) > 0) {
        avg_TE_x_to_y = mean(final_output$TE_x_to_y[window_indices])
        avg_TE_y_to_x = mean(final_output$TE_y_to_x[window_indices])
        avg_Pval_x_to_y = mean(final_output$Pval_x_to_y[window_indices])
        avg_Pval_y_to_x = mean(final_output$Pval_y_to_x[window_indices])
        avg_Stat_x = mean(final_output$Stat_x[window_indices])
        avg_Stat_y = mean(final_output$Stat_y[window_indices])

        avg_results = rbind(avg_results, data.frame(
          Time = time_point,
          Window_size = window_size,
          Avg_TE_x_to_y = avg_TE_x_to_y,
          Avg_TE_y_to_x = avg_TE_y_to_x,
          Avg_Pval_x_to_y = avg_Pval_x_to_y,
          Avg_Pval_y_to_x = avg_Pval_y_to_x,
          Avg_Stat_x = avg_Stat_x,
          Avg_Stat_y = avg_Stat_y
        ))
      }
    }
  }

  # Remove the first row which is all NA
  #avg_results = avg_results[-1, ]

  # Save average results to a CSV file
  write.table(avg_results, file = "C:/Users/User/Desktop/avgData.csv", row.names = FALSE, col.names = TRUE, sep = ",", na = "")

  print("Analysis Complete")

  # return
  return(final_output)
}
