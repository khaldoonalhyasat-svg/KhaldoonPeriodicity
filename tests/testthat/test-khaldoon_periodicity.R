# Test file for khaldoon_periodicity function

test_that("khaldoon_periodicity works with simple sine wave", {
  set.seed(123)
  t <- seq(0, 10, length.out = 200)
  x <- sin(2 * pi * t) + rnorm(200, sd = 0.3)
  
  result <- khaldoon_periodicity(x, B = 100)
  
  expect_s3_class(result, "khaldoon_test")
  expect_true(result$p.value < 0.05)
  expect_true(result$peak_to_mean > 1)
  expect_type(result$statistic, "double")
})

test_that("khaldoon_periodicity fails to reject for white noise", {
  set.seed(456)
  x <- rnorm(200)
  
  result <- khaldoon_periodicity(x, B = 100)
  
  expect_s3_class(result, "khaldoon_test")
  expect_true(result$p.value >= 0.05)
})

test_that("input validation works", {
  expect_error(khaldoon_periodicity("not numeric"),
               "Input 'x' must be a numeric vector")
  
  expect_error(khaldoon_periodicity(c(1, 2, NA, 4)),
               "Input 'x' contains missing values")
  
  expect_error(khaldoon_periodicity(1:10),
               "Input 'x' must have at least 20 observations")
  
  expect_error(khaldoon_periodicity(rnorm(100), alpha = 1.5),
               "Significance level 'alpha' must be between 0 and 1")
})

test_that("different window functions work", {
  set.seed(789)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  
  result_hanning <- khaldoon_periodicity(x, window = "hanning", B = 50)
  result_hamming <- khaldoon_periodicity(x, window = "hamming", B = 50)
  result_rect <- khaldoon_periodicity(x, window = "rectangular", B = 50)
  
  expect_s3_class(result_hanning, "khaldoon_test")
  expect_s3_class(result_hamming, "khaldoon_test")
  expect_s3_class(result_rect, "khaldoon_test")
})

test_that("detrending method works", {
  set.seed(321)
  t <- seq(0, 10, length.out = 200)
  x <- 0.5 * t + sin(2 * pi * t) + rnorm(200, sd = 0.3)
  
  result_auto <- khaldoon_periodicity(x, method = "auto", B = 100)
  result_detrend <- khaldoon_periodicity(x, method = "detrend", B = 100)
  result_none <- khaldoon_periodicity(x, method = "none", B = 100)
  
  expect_s3_class(result_auto, "khaldoon_test")
  expect_s3_class(result_detrend, "khaldoon_test")
  expect_s3_class(result_none, "khaldoon_test")
})

test_that("print method works", {
  set.seed(111)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  result <- khaldoon_periodicity(x, B = 50)
  
  expect_output(print(result), "Khaldoon Periodicity Test Results")
  expect_output(print(result), "Test Statistic")
  expect_output(print(result), "P-value")
})

test_that("summary method works", {
  set.seed(222)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  result <- khaldoon_periodicity(x, B = 50)
  
  expect_output(summary(result), "Khaldoon Periodicity Test - Summary")
  expect_output(summary(result), "Bootstrap Statistics")
})

test_that("plot method works without errors", {
  set.seed(333)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  result <- khaldoon_periodicity(x, B = 50)
  
  expect_silent(plot(result, which = 1))
  expect_silent(plot(result, which = 2))
  expect_silent(plot(result, which = 3))
})

test_that("bootstrap statistics are reasonable", {
  set.seed(444)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  result <- khaldoon_periodicity(x, B = 100)
  
  expect_length(result$bootstrap_stats, 100)
  expect_true(all(is.finite(result$bootstrap_stats)))
  expect_true(all(result$bootstrap_stats >= 0))
})

test_that("dominant frequency detection works", {
  set.seed(555)
  freq <- 0.1
  t <- seq(0, 1, length.out = 200)
  x <- sin(2 * pi * freq * 200 * t) + rnorm(200, sd = 0.2)
  
  result <- khaldoon_periodicity(x, B = 100)
  
  expect_true(!is.na(result$dominant_frequency))
  expect_true(result$dominant_frequency > 0)
  # Should be close to 0.1 (within tolerance)
  expect_true(abs(result$dominant_frequency - freq) < 0.05)
})

test_that("multiple periodic components are detected", {
  set.seed(666)
  t <- seq(0, 1, length.out = 200)
  x <- sin(2 * pi * 0.1 * 200 * t) +
       0.5 * sin(2 * pi * 0.25 * 200 * t) +
       rnorm(200, sd = 0.2)
  
  result <- khaldoon_periodicity(x, B = 100)
  
  expect_true(result$significant_peaks >= 2)
  expect_true(result$p.value < 0.05)
})

test_that("generate_periodic_data works for all types", {
  types <- c("sine", "multi", "seasonal", "trend_periodic",
             "white_noise", "random_walk")
  
  for (type in types) {
    data <- generate_periodic_data(n = 100, type = type, seed = 123)
    expect_length(data, 100)
    expect_true(all(is.finite(data)))
  }
})

test_that("generate_test_datasets works", {
  datasets <- generate_test_datasets(n_each = 100, seed = 123)
  
  expect_type(datasets, "list")
  expect_length(datasets, 6)
  expect_true(all(sapply(datasets, length) == 100))
})

test_that("edge case: very short series", {
  expect_error(khaldoon_periodicity(rnorm(15)),
               "Input 'x' must have at least 20 observations")
})

test_that("edge case: constant series", {
  x <- rep(1, 100)
  result <- khaldoon_periodicity(x, B = 50)
  
  expect_s3_class(result, "khaldoon_test")
  expect_true(result$p.value >= 0.05)
})

test_that("energy concentration is in valid range", {
  set.seed(777)
  x <- sin(2 * pi * seq(0, 10, length.out = 100)) + rnorm(100, sd = 0.2)
  result <- khaldoon_periodicity(x, B = 50)
  
  expect_true(result$energy_concentration >= 0)
  expect_true(result$energy_concentration <= 1)
})
