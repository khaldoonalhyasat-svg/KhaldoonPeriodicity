# Internal Utility Functions
.calculate_periodogram <- function(x) {
  n <- length(x)
  fft_result <- stats::fft(x)
  periodogram <- (Mod(fft_result[2:floor(n/2)])^2) / n
  frequencies <- (1:length(periodogram)) / n
  return(list(periodogram = periodogram, frequencies = frequencies))
}
.calculate_components <- function(periodogram, n) {
  pgram <- periodogram$periodogram
  
  # Peak-to-mean ratio
  peak_to_mean <- max(pgram) / mean(pgram)
  
  # Number of significant peaks (above 95th percentile)
  threshold <- stats::quantile(pgram, 0.95)
  significant_peaks <- sum(pgram > threshold)
  
  # Energy concentration (top 3 peaks)
  sorted_pgram <- sort(pgram, decreasing = TRUE)
  energy_concentration <- sum(sorted_pgram[1:3]) / sum(pgram)
  
  return(list(
    peak_to_mean = peak_to_mean,
    significant_peaks = significant_peaks,
    energy_concentration = energy_concentration
  ))
}
.calculate_test_statistic <- function(components) {
  R <- components$peak_to_mean
  S <- components$significant_peaks
  C <- components$energy_concentration
  
  T_n <- R * log(1 + S) * C
  return(T_n)
}
.perform_bootstrap <- function(x, test_statistic, n_bootstrap, parallel, n_cores) {
  
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    bootstrap_stats <- foreach::foreach(i = 1:n_bootstrap, .combine = c) %dopar% {
      permuted_x <- sample(x)
      permuted_periodogram <- .calculate_periodogram(permuted_x)
      permuted_components <- .calculate_components(permuted_periodogram, length(permuted_x))
      .calculate_test_statistic(permuted_components)
    }
    
    parallel::stopCluster(cl)
  } else {
    bootstrap_stats <- numeric(n_bootstrap)
    for (i in 1:n_bootstrap) {
      permuted_x <- sample(x)
      permuted_periodogram <- .calculate_periodogram(permuted_x)
      permuted_components <- .calculate_components(permuted_periodogram, length(permuted_x))
      bootstrap_stats[i] <- .calculate_test_statistic(permuted_components)
    }
  }
  
  p_value <- mean(bootstrap_stats >= test_statistic)
  
  return(list(
    p_value = p_value,
    bootstrap_distribution = bootstrap_stats
  ))
}




#' Sample Random Data
#'
#' Generate random white noise data for testing purposes
#'
#' @param n Number of observations (default: 100)
#' @param mean Mean of the noise (default: 0)
#' @param sd Standard deviation of the noise (default: 1)
#' @return Numeric vector of white noise
#' @importFrom stats rnorm
#' @export
sample_random_data <- function(n = 100, mean = 0, sd = 1) {
  stats::rnorm(n, mean = mean, sd = sd)
}


