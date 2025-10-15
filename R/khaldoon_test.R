#' Khaldoon Periodicity Test
#'
#' A robust test for detecting periodic patterns in time series data using
#' composite spectral statistics and bootstrap inference.
#'
#' @param x A numeric vector containing the time series data
#' @param n_bootstrap Number of bootstrap replications (default: 1000)
#' @param significance_threshold Significance level for the test (default: 0.05)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: NULL = auto-detect)
#'
#' @return An object of class 'khaldoon_test' containing test results and components
#'
#' @examples
#' # Generate synthetic periodic data
#' set.seed(123)
#' t <- 1:100
#' periodic_data <- sin(2 * pi * t / 20) + rnorm(100, 0, 0.5)
#'
#' # Perform Khaldoon test
#' result <- khaldoon_test(periodic_data, n_bootstrap = 100)
#' print(result)
#' plot(result)
#'
#' @export
#' @importFrom stats na.omit quantile density
#' @importFrom graphics plot grid barplot hist lines abline legend text
#' @importFrom foreach %dopar%
khaldoon_test <- function(x, n_bootstrap = 1000, significance_threshold = 0.05,
                          parallel = FALSE, n_cores = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector")
  }

  if (length(x) < 10) {
    warning("Time series is very short - results may be unreliable")
  }

  # Remove NA values
  x_clean <- stats::na.omit(x)
  n <- length(x_clean)

  # Calculate periodogram
  periodogram <- .calculate_periodogram(x_clean)

  # Calculate component statistics
  components <- .calculate_components(periodogram, n)

  # Calculate test statistic
  test_statistic <- .calculate_test_statistic(components)

  # Bootstrap inference
  bootstrap_results <- .perform_bootstrap(x_clean, test_statistic, n_bootstrap,
                                          parallel, n_cores)

  # Create result object
  result <- list(
    test_statistic = test_statistic,
    p_value = bootstrap_results$p_value,
    components = components,
    periodogram = periodogram,
    bootstrap_distribution = bootstrap_results$bootstrap_distribution,
    significance_threshold = significance_threshold,
    is_significant = bootstrap_results$p_value < significance_threshold,
    n_observations = n,
    n_bootstrap = n_bootstrap
  )

  class(result) <- "khaldoon_test"
  return(result)
}