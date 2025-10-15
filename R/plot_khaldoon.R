#' Plot Method for Khaldoon Test Results
#'
#' Visualize periodogram, component statistics, bootstrap distribution, and test results
#'
#' @param x An object of class 'khaldoon_test'
#' @param type Type of plot: "summary" (default), "periodogram", or "bootstrap"
#' @param cex_main Main title size
#' @param cex_lab Axis label size
#' @param cex_axis Axis text size
#' @param cex_text Legend and text size
#' @param ... Additional plotting parameters
#'
#' @return No return value, called for side effects
#'
#' @examples
#' \donttest{
#' data <- sin(seq(0, 4*pi, length.out = 100)) + rnorm(100, 0, 0.3)
#' result <- khaldoon_test(data, n_bootstrap = 50)
#' plot(result)
#' plot(result, type = "periodogram")
#' plot(result, type = "bootstrap")
#' }
#'
#' @export
plot.khaldoon_test <- function(x, type = "summary",
                              cex_main = 1.1,
                              cex_lab = 0.9,
                              cex_axis = 0.85,
                              cex_text = 0.85, ...) {

  if (type == "summary") {
    # Set up 2x2 plot layout with balanced margins
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(
      mfrow = c(2, 2),
      mar = c(3.5, 3.5, 2.5, 1.5),
      mgp = c(2, 0.6, 0),
      cex = 0.85
    )
    on.exit(graphics::par(old_par))

    # Plot 1: Periodogram
    plot(x$periodogram$frequencies,
         x$periodogram$periodogram,
         type = "l",
         main = "Periodogram",
         xlab = "Frequency",
         ylab = "Power",
         cex.main = cex_main,
         cex.lab = cex_lab,
         cex.axis = cex_axis,
         col = "darkblue",
         lwd = 2)
    grid()

    # Plot 2: Component Statistics
    comp_data <- c(x$components$peak_to_mean,
                   x$components$energy_concentration)
    names(comp_data) <- c("Peak-to-Mean", "Energy Conc.")
    barplot(comp_data,
            main = "Component Statistics",
            col = c("skyblue", "lightcoral"),
            ylab = "Value",
            cex.main = cex_main,
            cex.lab = cex_lab,
            cex.names = cex_axis,
            cex.axis = cex_axis,
            ylim = c(0, max(comp_data) * 1.1))

    # Plot 3: Bootstrap Distribution
    hist(x$bootstrap_distribution,
         main = "Bootstrap Distribution",
         xlab = "Test Statistic",
         col = "lightblue",
         breaks = 20,
         probability = TRUE,
         cex.main = cex_main,
         cex.lab = cex_lab,
         cex.axis = cex_axis)
    lines(density(x$bootstrap_distribution), col = "darkblue", lwd = 2)
    abline(v = x$test_statistic, col = "red", lwd = 2.5, lty = 2)
    legend("topright",
           legend = c("Observed", "Bootstrap"),
           col = c("red", "lightblue"),
           lwd = c(2, 8),
           lty = c(2, 1),
           bty = "n",
           cex = cex_text)

    # Plot 4: Test Results
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "",
         main = "Test Results",
         cex.main = cex_main)

    significance_text <- ifelse(x$is_significant, "SIGNIFICANT", "NOT SIGNIFICANT")
    color <- ifelse(x$is_significant, "red", "blue")

    p_value_text <- ifelse(x$p_value < 0.001,
                          format.pval(x$p_value, digits = 2),
                          round(x$p_value, 4))

    text(0.5, 0.7, paste("p-value:", p_value_text), cex = cex_text, font = 2)
    text(0.5, 0.5, significance_text, cex = cex_text * 1.3, col = color, font = 2)
    text(0.5, 0.3, paste("Statistic:", round(x$test_statistic, 4)), cex = cex_text)

  } else if (type == "periodogram") {
    # Periodogram only plot
    plot(x$periodogram$frequencies,
         x$periodogram$periodogram,
         type = "l",
         main = "Periodogram",
         xlab = "Frequency",
         ylab = "Power",
         cex.main = cex_main * 1.2,
         cex.lab = cex_lab,
         cex.axis = cex_axis,
         col = "darkblue",
         lwd = 2,
         ...)
    grid()

  } else if (type == "bootstrap") {
    # Bootstrap only plot
    hist(x$bootstrap_distribution,
         main = "Bootstrap Distribution",
         xlab = "Test Statistic",
         col = "lightblue",
         breaks = 25,
         probability = TRUE,
         cex.main = cex_main * 1.2,
         cex.lab = cex_lab,
         cex.axis = cex_axis,
         ...)
    lines(density(x$bootstrap_distribution), col = "darkblue", lwd = 2)
    abline(v = x$test_statistic, col = "red", lwd = 2.5, lty = 2)
    legend("topright",
           legend = c("Observed", "Bootstrap"),
           col = c("red", "lightblue"),
           lwd = c(2, 8),
           lty = c(2, 1),
           bty = "n",
           cex = cex_text)
  } else {
    stop("Unknown plot type. Use 'summary', 'periodogram', or 'bootstrap'")
  }
}

#' Print Method for Khaldoon Test Results
#'
#' @param x An object of class 'khaldoon_test'
#' @param ... Additional parameters
#'
#' @return No return value, called for side effects
#' @export
print.khaldoon_test <- function(x, ...) {
  cat("Khaldoon Periodicity Test\n")
  cat("=========================\n\n")
  cat("Test Statistic:", round(x$test_statistic, 4), "\n")
  cat("P-value:", round(x$p_value, 4), "\n")
  cat("Significant at", x$significance_threshold, "level:",
      ifelse(x$is_significant, "YES", "NO"), "\n\n")

  cat("Component Statistics:\n")
  cat("  Peak-to-Mean Ratio (R):", round(x$components$peak_to_mean, 4), "\n")
  cat("  Significant Peaks (S):", x$components$significant_peaks, "\n")
  cat("  Energy Concentration (C):", round(x$components$energy_concentration, 4), "\n\n")

  cat("Sample Size:", x$n_observations, "\n")
  cat("Bootstrap Replications:", x$n_bootstrap, "\n")
}

#' Summary Method for Khaldoon Test Results
#'
#' @param object An object of class 'khaldoon_test'
#' @param ... Additional parameters
#'
#' @return No return value, called for side effects
#' @export
summary.khaldoon_test <- function(object, ...) {
  print(object)
}
