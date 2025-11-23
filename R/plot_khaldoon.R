#' Plot Method for Khaldoon Test Results
#'
#' @description
#' Produces diagnostic plots for Khaldoon periodicity test results including
#' power spectrum, bootstrap distribution, and QQ plot.
#'
#' @param x An object of class "khaldoon_test"
#' @param which Which plots to display. Can be 1, 2, 3, or "all" (default: "all")
#'        1 = Power Spectrum
#'        2 = Bootstrap Distribution
#'        3 = Cumulative Spectrum
#' @param ... Additional graphical parameters
#'
#' @examples
#' set.seed(123)
#' x <- sin(2 * pi * seq(0, 10, length.out = 200)) + rnorm(200, sd = 0.3)
#' result <- khaldoon_periodicity(x)
#' plot(result)
#' plot(result, which = 1)  # Only power spectrum
#'
#' @export
#' @importFrom graphics par plot lines abline legend hist
plot.khaldoon_test <- function(x, which = "all", ...) {
  
  if (which == "all") {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    on.exit(par(old_par))
    
    plot_spectrum(x, ...)
    plot_bootstrap(x, ...)
    plot_cumulative_spectrum(x, ...)
    plot_periodogram(x, ...)
    
  } else if (which == 1) {
    plot_spectrum(x, ...)
  } else if (which == 2) {
    plot_bootstrap(x, ...)
  } else if (which == 3) {
    plot_cumulative_spectrum(x, ...)
  } else {
    stop("Invalid 'which' argument. Use 1, 2, 3, or 'all'")
  }
  
  invisible(x)
}


#' Plot Power Spectrum
#'
#' @param x Khaldoon test result object
#' @param ... Additional graphical parameters
#'
#' @keywords internal
#' @noRd
plot_spectrum <- function(x, ...) {
  
  if (length(x$spectrum) == 0) {
    plot.new()
    text(0.5, 0.5, "Insufficient data for spectrum plot", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Compute threshold
  threshold <- mean(x$spectrum) + 2 * sd(x$spectrum)
  
  plot(x$frequencies, x$spectrum,
       type = "h",
       lwd = 2,
       col = "steelblue",
       xlab = "Frequency (cycles per observation)",
       ylab = "Power",
       main = "Power Spectrum",
       ...)
  
  # Add threshold line
  abline(h = threshold, col = "red", lty = 2, lwd = 2)
  
  # Mark dominant frequency
  if (!is.na(x$dominant_frequency)) {
    dominant_idx <- which.max(x$spectrum)
    points(x$frequencies[dominant_idx], x$spectrum[dominant_idx],
           col = "red", pch = 19, cex = 1.5)
  }
  
  # Add legend
  legend("topright",
         legend = c("Spectrum", "Significance Threshold", "Dominant Peak"),
         col = c("steelblue", "red", "red"),
         lty = c(1, 2, NA),
         pch = c(NA, NA, 19),
         lwd = 2,
         bty = "n")
}


#' Plot Bootstrap Distribution
#'
#' @param x Khaldoon test result object
#' @param ... Additional graphical parameters
#'
#' @keywords internal
#' @noRd
plot_bootstrap <- function(x, ...) {
  
  hist(x$bootstrap_stats,
       breaks = 30,
       col = "lightblue",
       border = "white",
       xlab = "Test Statistic",
       main = "Bootstrap Distribution",
       probability = TRUE,
       ...)
  
  # Add observed statistic
  abline(v = x$statistic, col = "red", lwd = 3, lty = 2)
  
  # Add density curve
  dens <- density(x$bootstrap_stats)
  lines(dens, col = "darkblue", lwd = 2)
  
  # Add text
  legend("topright",
         legend = c(
           paste0("Observed: ", sprintf("%.3f", x$statistic)),
           paste0("P-value: ", sprintf("%.4f", x$p.value))
         ),
         col = c("red", NA),
         lty = c(2, NA),
         lwd = c(3, NA),
         bty = "n")
}


#' Plot Cumulative Spectrum
#'
#' @param x Khaldoon test result object
#' @param ... Additional graphical parameters
#'
#' @keywords internal
#' @noRd
plot_cumulative_spectrum <- function(x, ...) {
  
  if (length(x$spectrum) == 0) {
    plot.new()
    text(0.5, 0.5, "Insufficient data for cumulative plot", cex = 1.2)
    return(invisible(NULL))
  }
  
  cumulative_power <- cumsum(x$spectrum) / sum(x$spectrum)
  
  plot(x$frequencies, cumulative_power,
       type = "l",
       lwd = 2,
       col = "darkgreen",
       xlab = "Frequency (cycles per observation)",
       ylab = "Cumulative Power Proportion",
       main = "Cumulative Power Spectrum",
       ylim = c(0, 1),
       ...)
  
  # Add reference lines
  abline(h = c(0.5, 0.8, 0.95), col = "gray", lty = 2)
  
  # Add grid
  grid()
  
  # Add text for energy concentration
  text(max(x$frequencies) * 0.7, 0.1,
       paste0("Energy Concentration:\n", sprintf("%.3f", x$energy_concentration)),
       cex = 0.9,
       adj = 0)
}


#' Plot Periodogram
#'
#' @param x Khaldoon test result object
#' @param ... Additional graphical parameters
#'
#' @keywords internal
#' @noRd
plot_periodogram <- function(x, ...) {
  
  if (length(x$spectrum) == 0) {
    plot.new()
    text(0.5, 0.5, "Insufficient data", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Convert to period scale
  periods <- 1 / x$frequencies
  periods <- periods[is.finite(periods) & periods > 0]
  spectrum_subset <- x$spectrum[is.finite(periods) & periods > 0]
  
  if (length(periods) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid periods", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Limit to reasonable periods
  max_period <- min(length(x$spectrum), 50)
  valid_idx <- periods <= max_period
  
  plot(periods[valid_idx], spectrum_subset[valid_idx],
       type = "h",
       lwd = 2,
       col = "purple",
       xlab = "Period (observations)",
       ylab = "Power",
       main = "Periodogram",
       ...)
  
  # Mark dominant period
  if (!is.na(x$dominant_frequency) && x$dominant_frequency > 0) {
    dominant_period <- 1 / x$dominant_frequency
    if (dominant_period <= max_period) {
      abline(v = dominant_period, col = "red", lwd = 2, lty = 2)
      text(dominant_period, max(spectrum_subset[valid_idx]) * 0.9,
           paste0("Period: ", sprintf("%.1f", dominant_period)),
           pos = 4, col = "red", cex = 0.9)
    }
  }
  
  grid()
}


#' Summary Method for Khaldoon Test
#'
#' @param object An object of class "khaldoon_test"
#' @param ... Additional arguments (not used)
#'
#' @export
summary.khaldoon_test <- function(object, ...) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("     Khaldoon Periodicity Test - Summary          \n")
  cat("═══════════════════════════════════════════════════\n\n")
  
  cat("Call: Khaldoon Periodicity Test\n")
  cat("Data:", object$data.name, "\n\n")
  
  cat("Sample Information:\n")
  cat("  Number of observations:     ", object$n, "\n")
  cat("  Bootstrap replications:     ", object$B, "\n\n")
  
  cat("Test Results:\n")
  cat("  Test statistic:             ", sprintf("%.6f", object$statistic), "\n")
  cat("  P-value:                    ", sprintf("%.6f", object$p.value), "\n")
  cat("  Significance level:         ", object$alpha, "\n")
  cat("  Decision:                   ", object$decision, "\n\n")
  
  cat("Statistic Components:\n")
  cat("  Peak-to-Mean Ratio (R):     ", sprintf("%.4f", object$peak_to_mean), "\n")
  cat("  Significant Peaks (S):      ", object$significant_peaks, "\n")
  cat("  Energy Concentration (C):   ", sprintf("%.4f", object$energy_concentration), "\n\n")
  
  if (!is.na(object$dominant_frequency)) {
    cat("Dominant Frequency Information:\n")
    cat("  Frequency:                  ", sprintf("%.6f", object$dominant_frequency), 
        "cycles/observation\n")
    
    if (object$dominant_frequency > 0) {
      period <- 1 / object$dominant_frequency
      cat("  Estimated Period:           ", sprintf("%.2f", period), "observations\n")
    }
  }
  
  cat("\n")
  cat("Bootstrap Statistics:\n")
  cat("  Min:                        ", sprintf("%.4f", min(object$bootstrap_stats)), "\n")
  cat("  1st Quartile:               ", sprintf("%.4f", quantile(object$bootstrap_stats, 0.25)), "\n")
  cat("  Median:                     ", sprintf("%.4f", median(object$bootstrap_stats)), "\n")
  cat("  Mean:                       ", sprintf("%.4f", mean(object$bootstrap_stats)), "\n")
  cat("  3rd Quartile:               ", sprintf("%.4f", quantile(object$bootstrap_stats, 0.75)), "\n")
  cat("  Max:                        ", sprintf("%.4f", max(object$bootstrap_stats)), "\n")
  
  cat("\n═══════════════════════════════════════════════════\n\n")
  
  invisible(object)
}
