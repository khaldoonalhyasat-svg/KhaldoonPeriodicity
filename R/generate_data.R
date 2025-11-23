#' Generate Periodic Time Series Data
#'
#' @description
#' Generate synthetic time series data with various periodic patterns for
#' testing and demonstration purposes.
#'
#' @param n Length of the time series (default: 200)
#' @param type Type of data to generate. Options: "sine", "multi", "seasonal",
#'        "trend_periodic", "white_noise", "random_walk" (default: "sine")
#' @param noise_sd Standard deviation of additive Gaussian noise (default: 0.3)
#' @param frequency Fundamental frequency for periodic signals (default: 0.1)
#' @param trend_coef Trend coefficient for "trend_periodic" type (default: 0.05)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return A numeric vector containing the generated time series
#'
#' @details
#' Available data types:
#' \itemize{
#'   \item \code{"sine"}: Simple sinusoidal pattern
#'   \item \code{"multi"}: Multiple frequency components
#'   \item \code{"seasonal"}: Seasonal pattern with multiple harmonics
#'   \item \code{"trend_periodic"}: Periodic signal with linear trend
#'   \item \code{"white_noise"}: Pure random noise (no periodicity)
#'   \item \code{"random_walk"}: Random walk (no periodicity)
#' }
#'
#' @examples
#' # Simple sine wave
#' x1 <- generate_periodic_data(n = 200, type = "sine")
#' plot(x1, type = "l", main = "Simple Sine Wave")
#'
#' # Multiple frequencies
#' x2 <- generate_periodic_data(n = 200, type = "multi")
#' plot(x2, type = "l", main = "Multiple Frequencies")
#'
#' # White noise (non-periodic)
#' x3 <- generate_periodic_data(n = 200, type = "white_noise")
#' plot(x3, type = "l", main = "White Noise")
#'
#' @export
#' @importFrom stats rnorm
generate_periodic_data <- function(n = 200,
                                    type = "sine",
                                    noise_sd = 0.3,
                                    frequency = 0.1,
                                    trend_coef = 0.05,
                                    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  t <- seq(0, 1, length.out = n)
  
  data <- switch(type,
    "sine" = {
      # Simple sinusoidal pattern
      A <- 2
      sin(2 * pi * frequency * n * t) * A + rnorm(n, sd = noise_sd)
    },
    
    "multi" = {
      # Multiple frequency components
      A1 <- 2
      A2 <- 1.5
      A3 <- 1
      f1 <- frequency
      f2 <- frequency * 2.5
      f3 <- frequency * 4
      A1 * sin(2 * pi * f1 * n * t) +
        A2 * sin(2 * pi * f2 * n * t) +
        A3 * sin(2 * pi * f3 * n * t) +
        rnorm(n, sd = noise_sd)
    },
    
    "seasonal" = {
      # Seasonal pattern (e.g., monthly data with yearly cycle)
      # Assume n represents 3 years of monthly data
      monthly_t <- seq(0, 3, length.out = n)
      A_yearly <- 2
      A_biannual <- 0.8
      A_yearly * sin(2 * pi * monthly_t) +
        A_biannual * sin(2 * pi * 2 * monthly_t) +
        rnorm(n, sd = noise_sd)
    },
    
    "trend_periodic" = {
      # Periodic signal with linear trend
      trend <- trend_coef * n * t
      periodic <- 2 * sin(2 * pi * frequency * n * t)
      trend + periodic + rnorm(n, sd = noise_sd)
    },
    
    "white_noise" = {
      # Pure white noise (no periodicity)
      rnorm(n, mean = 0, sd = 1)
    },
    
    "random_walk" = {
      # Random walk (no periodicity)
      cumsum(rnorm(n, mean = 0, sd = 0.5))
    },
    
    stop("Invalid type. Choose from: 'sine', 'multi', 'seasonal', 'trend_periodic', 'white_noise', 'random_walk'")
  )
  
  return(data)
}


#' Generate Multiple Test Datasets
#'
#' @description
#' Generate a list of multiple synthetic datasets with different characteristics
#' for comprehensive testing.
#'
#' @param n_each Number of observations in each dataset (default: 200)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A named list of numeric vectors, each containing a different type of
#' time series data
#'
#' @examples
#' datasets <- generate_test_datasets()
#' names(datasets)
#'
#' # Test on all datasets
#' results <- lapply(datasets, khaldoon_periodicity)
#'
#' # Print all p-values
#' sapply(results, function(x) x$p.value)
#'
#' @export
generate_test_datasets <- function(n_each = 200, seed = 123) {
  
  set.seed(seed)
  
  datasets <- list(
    simple_sine = generate_periodic_data(n_each, "sine", seed = seed),
    multi_frequency = generate_periodic_data(n_each, "multi", seed = seed + 1),
    seasonal = generate_periodic_data(n_each, "seasonal", seed = seed + 2),
    trend_periodic = generate_periodic_data(n_each, "trend_periodic", seed = seed + 3),
    white_noise = generate_periodic_data(n_each, "white_noise", seed = seed + 4),
    random_walk = generate_periodic_data(n_each, "random_walk", seed = seed + 5)
  )
  
  return(datasets)
}


#' Simulate Power Analysis
#'
#' @description
#' Perform a simulation study to assess the statistical power of the Khaldoon
#' periodicity test under various conditions.
#'
#' @param n Sample size (default: 200)
#' @param n_sim Number of simulations (default: 1000)
#' @param alpha Significance level (default: 0.05)
#' @param signal_strengths Vector of signal-to-noise ratios to test (default: seq(0.5, 3, by = 0.5))
#' @param frequency Frequency of the periodic signal (default: 0.1)
#' @param B Number of bootstrap replicates per test (default: 300)
#' @param seed Random seed (default: NULL)
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{signal_strength}: Signal-to-noise ratio
#'   \item \code{power}: Proportion of rejections (statistical power)
#'   \item \code{n_sim}: Number of simulations
#' }
#'
#' @examples
#' \donttest{
#' # Quick power analysis (use more simulations in practice)
#' power_results <- simulate_power(n = 100, n_sim = 100)
#' print(power_results)
#'
#' # Plot power curve
#' plot(power_results$signal_strength, power_results$power,
#'      type = "b", pch = 19,
#'      xlab = "Signal-to-Noise Ratio",
#'      ylab = "Statistical Power",
#'      main = "Power Analysis for Khaldoon Test",
#'      ylim = c(0, 1))
#' abline(h = 0.8, lty = 2, col = "red")
#' grid()
#' }
#'
#' @export
simulate_power <- function(n = 200,
                           n_sim = 1000,
                           alpha = 0.05,
                           signal_strengths = seq(0.5, 3, by = 0.5),
                           frequency = 0.1,
                           B = 300,
                           seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  results <- data.frame(
    signal_strength = numeric(),
    power = numeric(),
    n_sim = integer()
  )
  
  for (signal_strength in signal_strengths) {
    
    cat("Testing signal strength:", signal_strength, "\n")
    
    rejections <- 0
    
    for (i in 1:n_sim) {
      
      # Generate data with specified signal-to-noise ratio
      t <- seq(0, 1, length.out = n)
      signal <- signal_strength * sin(2 * pi * frequency * n * t)
      noise <- rnorm(n, sd = 1)
      x <- signal + noise
      
      # Perform test
      test_result <- khaldoon_periodicity(x, alpha = alpha, B = B, plot = FALSE)
      
      # Count rejections
      if (test_result$p.value < alpha) {
        rejections <- rejections + 1
      }
    }
    
    power <- rejections / n_sim
    
    results <- rbind(results, data.frame(
      signal_strength = signal_strength,
      power = power,
      n_sim = n_sim
    ))
  }
  
  return(results)
}


#' Simulate Type I Error Rate
#'
#' @description
#' Assess the Type I error rate (false positive rate) of the Khaldoon test
#' under the null hypothesis of no periodicity.
#'
#' @param n Sample size (default: 200)
#' @param n_sim Number of simulations (default: 1000)
#' @param alpha Nominal significance level (default: 0.05)
#' @param B Number of bootstrap replicates (default: 300)
#' @param null_type Type of null data. Options: "white_noise", "random_walk", "linear_trend" (default: "white_noise")
#' @param seed Random seed (default: NULL)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{type_I_error}: Observed Type I error rate
#'   \item \code{n_rejections}: Number of rejections
#'   \item \code{n_sim}: Number of simulations
#'   \item \code{alpha}: Nominal significance level
#'   \item \code{ci_lower}: Lower bound of 95% confidence interval
#'   \item \code{ci_upper}: Upper bound of 95% confidence interval
#' }
#'
#' @examples
#' \donttest{
#' # Assess Type I error rate (use more simulations in practice)
#' type1_results <- simulate_type1_error(n = 200, n_sim = 100)
#' print(type1_results)
#' }
#'
#' @export
simulate_type1_error <- function(n = 200,
                                  n_sim = 1000,
                                  alpha = 0.05,
                                  B = 300,
                                  null_type = "white_noise",
                                  seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  rejections <- 0
  
  cat("Simulating Type I error rate...\n")
  cat("Null hypothesis type:", null_type, "\n")
  cat("Simulations:", n_sim, "\n\n")
  
  for (i in 1:n_sim) {
    
    if (i %% 100 == 0) {
      cat("Progress:", i, "/", n_sim, "\n")
    }
    
    # Generate null data
    x <- generate_periodic_data(n = n, type = null_type)
    
    # Perform test
    test_result <- khaldoon_periodicity(x, alpha = alpha, B = B, plot = FALSE)
    
    # Count rejections
    if (test_result$p.value < alpha) {
      rejections <- rejections + 1
    }
  }
  
  type1_error <- rejections / n_sim
  
  # Calculate 95% confidence interval for Type I error
  se <- sqrt(type1_error * (1 - type1_error) / n_sim)
  ci_lower <- max(0, type1_error - 1.96 * se)
  ci_upper <- min(1, type1_error + 1.96 * se)
  
  result <- list(
    type_I_error = type1_error,
    n_rejections = rejections,
    n_sim = n_sim,
    alpha = alpha,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    null_type = null_type
  )
  
  cat("\n")
  cat("═══════════════════════════════════════\n")
  cat("Type I Error Analysis Results\n")
  cat("═══════════════════════════════════════\n")
  cat("Observed Type I error:  ", sprintf("%.4f", type1_error), "\n")
  cat("Expected (alpha):       ", alpha, "\n")
  cat("95% CI:                 [", sprintf("%.4f", ci_lower), ", ", 
      sprintf("%.4f", ci_upper), "]\n", sep = "")
  cat("Number of rejections:   ", rejections, "/", n_sim, "\n")
  cat("═══════════════════════════════════════\n\n")
  
  return(result)
}
