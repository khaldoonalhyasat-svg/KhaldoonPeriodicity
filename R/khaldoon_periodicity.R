#' Khaldoon Periodicity Test
#'
#' @description
#' A novel spectral-based statistical test for detecting periodic patterns in
#' time series data. The test combines Fast Fourier Transform (FFT) analysis
#' with bootstrap hypothesis testing.
#'
#' @param x A numeric vector or time series object containing the data
#' @param alpha Significance level for the test (default: 0.05)
#' @param B Number of bootstrap replicates (default: 300)
#' @param method Method for handling the data. Options: "auto", "detrend", "none" (default: "auto")
#' @param window Window function to apply. Options: "hanning", "hamming", "rectangular" (default: "hanning")
#' @param plot Logical. Should diagnostic plots be generated? (default: FALSE)
#'
#' @return An object of class "khaldoon_test" containing:
#' \item{statistic}{The value of the Khaldoon test statistic}
#' \item{p.value}{The bootstrap p-value}
#' \item{method}{Description of the test}
#' \item{data.name}{Name of the data}
#' \item{alpha}{Significance level used}
#' \item{decision}{Test decision (Reject/Fail to reject H0)}
#' \item{peak_to_mean}{Peak-to-mean ratio}
#' \item{significant_peaks}{Number of significant peaks}
#' \item{energy_concentration}{Energy concentration ratio}
#' \item{spectrum}{Power spectrum data}
#' \item{frequencies}{Frequency values}
#' \item{dominant_frequency}{Dominant frequency detected}
#' \item{bootstrap_stats}{Vector of bootstrap statistics}
#'
#' @details
#' The Khaldoon test examines the null hypothesis that the time series contains
#' no systematic periodic patterns. It uses a composite test statistic:
#'
#' \deqn{T_n = R \times \log(1 + S) \times C}
#'
#' where:
#' \itemize{
#'   \item R is the peak-to-mean ratio
#'   \item S is the number of significant peaks
#'   \item C is the energy concentration ratio
#' }
#'
#' @references
#' Alhyasat, K. (2025). Khaldoon Periodicity Test: A Novel Spectral-Based
#' Method for Detecting Periodic Patterns in Time Series Data.
#' \emph{In preparation}.
#'
#' @examples
#' # Simple sine wave (periodic)
#' set.seed(123)
#' t <- seq(0, 10, length.out = 200)
#' x <- sin(2 * pi * t) + rnorm(200, sd = 0.3)
#' result <- khaldoon_periodicity(x)
#' print(result)
#'
#' # White noise (non-periodic)
#' y <- rnorm(200)
#' result2 <- khaldoon_periodicity(y)
#' print(result2)
#'
#' # With plotting
#' result3 <- khaldoon_periodicity(x, plot = TRUE)
#'
#' @export
#' @importFrom stats fft sd rnorm lm residuals
#' @importFrom graphics par plot lines abline legend
khaldoon_periodicity <- function(x,
                                  alpha = 0.05,
                                  B = 300,
                                  method = "auto",
                                  window = "hanning",
                                  plot = FALSE) {
  
  # ========================================
  # Input Validation
  # ========================================
  
  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector")
  }
  
  if (any(is.na(x))) {
    stop("Input 'x' contains missing values (NA)")
  }
  
  n <- length(x)
  if (n < 20) {
    stop("Input 'x' must have at least 20 observations")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("Significance level 'alpha' must be between 0 and 1")
  }
  
  if (B < 100) {
    warning("Bootstrap replications 'B' is less than 100. Consider increasing for more reliable results.")
  }
  
  data_name <- deparse(substitute(x))
  
  # ========================================
  # Data Preprocessing
  # ========================================
  
  # Remove mean
  x_centered <- x - mean(x)
  
  # Apply detrending if needed
  if (method == "auto") {
    # Simple linear trend detection
    time_index <- seq_along(x_centered)
    trend_model <- lm(x_centered ~ time_index)
    if (!is.null(summary(trend_model)$r.squared) && !is.na(summary(trend_model)$r.squared) && summary(trend_model)$r.squared > 0.3) {
      x_centered <- residuals(trend_model)
      message("Linear trend detected and removed")
    }
  } else if (method == "detrend") {
    time_index <- seq_along(x_centered)
    trend_model <- lm(x_centered ~ time_index)
    x_centered <- residuals(trend_model)
  }
  
  # Apply window function
  if (window == "hanning") {
    w <- 0.5 - 0.5 * cos(2 * pi * seq_along(x_centered) / (n + 1))
  } else if (window == "hamming") {
    w <- 0.54 - 0.46 * cos(2 * pi * seq_along(x_centered) / (n + 1))
  } else {
    w <- rep(1, n)
  }
  
  x_windowed <- x_centered * w
  
  # ========================================
  # Compute Test Statistic
  # ========================================
  
  test_stat_result <- compute_khaldoon_statistic(x_windowed)
  observed_statistic <- test_stat_result$statistic
  
  # ========================================
  # Bootstrap Procedure
  # ========================================
  
  bootstrap_stats <- numeric(B)
  
  for (b in 1:B) {
    # Generate random permutation
    x_perm <- sample(x_windowed)
    boot_result <- compute_khaldoon_statistic(x_perm)
    bootstrap_stats[b] <- boot_result$statistic
  }
  
  # Calculate p-value
  p_value <- mean(bootstrap_stats >= observed_statistic)
  
  # ========================================
  # Test Decision
  # ========================================
  
  decision <- ifelse(p_value < alpha,
                     "Reject H0: Evidence of periodicity",
                     "Fail to reject H0: No significant periodicity")
  
  # ========================================
  # Prepare Output
  # ========================================
  
  result <- list(
    statistic = observed_statistic,
    p.value = p_value,
    method = "Khaldoon Periodicity Test",
    data.name = data_name,
    alpha = alpha,
    decision = decision,
    peak_to_mean = test_stat_result$peak_to_mean,
    significant_peaks = test_stat_result$significant_peaks,
    energy_concentration = test_stat_result$energy_concentration,
    spectrum = test_stat_result$spectrum,
    frequencies = test_stat_result$frequencies,
    dominant_frequency = test_stat_result$dominant_frequency,
    bootstrap_stats = bootstrap_stats,
    n = n,
    B = B
  )
  
  class(result) <- "khaldoon_test"
  
  # ========================================
  # Plotting
  # ========================================
  
  if (plot) {
    plot(result)
  }
  
  return(result)
}

#' Compute Khaldoon Test Statistic
#'
#' @description Internal function to compute the composite test statistic
#'
#' @param x Preprocessed time series data
#'
#' @return A list containing the statistic and its components
#'
#' @keywords internal
#' @noRd
compute_khaldoon_statistic <- function(x) {
  
  n <- length(x)
  
  # Check for constant data (zero variance)
  if (stats::var(x) == 0) {
    return(list(
      statistic = 0,
      peak_to_mean = 0,
      significant_peaks = 0,
      energy_concentration = 0,
      spectrum = numeric(0),
      frequencies = numeric(0),
      dominant_frequency = NA
    ))
  }
  
  
  # Compute FFT
  fft_result <- fft(x)
  
  # Compute power spectrum (only positive frequencies)
  spectrum <- Mod(fft_result[1:(n %/% 2)])^2
  frequencies <- (0:(n %/% 2 - 1)) / n
  
  # Remove DC component
  spectrum <- spectrum[-1]
  frequencies <- frequencies[-1]
  
  if (length(spectrum) == 0) {
    return(list(
      statistic = 0,
      peak_to_mean = 0,
      significant_peaks = 0,
      energy_concentration = 0,
      spectrum = numeric(0),
      frequencies = numeric(0),
      dominant_frequency = NA
    ))
  }
  
  # Component 1: Peak-to-Mean Ratio (R)
  peak_to_mean <- max(spectrum) / mean(spectrum)
  
  # Component 2: Number of Significant Peaks (S)
  threshold <- mean(spectrum) + 2 * sd(spectrum)
  significant_peaks <- sum(spectrum > threshold)
  
  # Component 3: Energy Concentration (C)
  top_5_indices <- order(spectrum, decreasing = TRUE)[1:min(5, length(spectrum))]
  energy_concentration <- sum(spectrum[top_5_indices]) / sum(spectrum)
  
  # Composite Statistic
  statistic <- peak_to_mean * log(1 + significant_peaks) * energy_concentration
  
  # Dominant frequency
  dominant_freq_idx <- which.max(spectrum)
  dominant_frequency <- frequencies[dominant_freq_idx]
  
  return(list(
    statistic = statistic,
    peak_to_mean = peak_to_mean,
    significant_peaks = significant_peaks,
    energy_concentration = energy_concentration,
    spectrum = spectrum,
    frequencies = frequencies,
    dominant_frequency = dominant_frequency
  ))
}

#' Print Method for Khaldoon Test
#'
#' @param x An object of class "khaldoon_test"
#' @param ... Additional arguments (not used)
#'
#' @export
print.khaldoon_test <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("         Khaldoon Periodicity Test Results         \n")
  cat("═══════════════════════════════════════════════════\n\n")
  
  cat("Data:                ", x$data.name, "\n")
  cat("Sample size:         ", x$n, "\n")
  cat("Bootstrap replications:", x$B, "\n\n")
  
  cat("───────────────────────────────────────────────────\n")
  cat("Test Statistic:      ", sprintf("%.4f", x$statistic), "\n")
  cat("P-value:             ", sprintf("%.4f", x$p.value), "\n")
  cat("Significance level:  ", x$alpha, "\n")
  cat("───────────────────────────────────────────────────\n\n")
  
  cat("Decision: ", x$decision, "\n\n")
  
  cat("Components:\n")
  cat("  • Peak-to-Mean Ratio:      ", sprintf("%.4f", x$peak_to_mean), "\n")
  cat("  • Significant Peaks:       ", x$significant_peaks, "\n")
  cat("  • Energy Concentration:    ", sprintf("%.4f", x$energy_concentration), "\n")
  
  if (!is.na(x$dominant_frequency)) {
    cat("  • Dominant Frequency:      ", sprintf("%.4f", x$dominant_frequency), "\n")
    if (x$dominant_frequency > 0) {
      period <- 1 / x$dominant_frequency
      cat("  • Estimated Period:        ", sprintf("%.2f observations", period), "\n")
    }
  }
  
  cat("\n═══════════════════════════════════════════════════\n\n")
  
  invisible(x)
}
