test_that("khaldoon_test detects periodicity", {
  set.seed(123)
  t <- 1:200
  periodic <- sin(2*pi*t/24) + rnorm(200, sd = 0.3)
  
  result <- khaldoon_test(periodic, n_bootstrap = 500)
  
  expect_lt(result$p_value, 0.05)
  expect_true(result$is_significant)
})

test_that("khaldoon_test handles white noise", {
  set.seed(456)
  noise <- rnorm(200)
  
  result <- khaldoon_test(noise, n_bootstrap = 500)
  
  expect_true(is.numeric(result$p_value))
  expect_gte(result$p_value, 0)
  expect_lte(result$p_value, 1)
})

test_that("components are valid", {
  set.seed(789)
  x <- sin(2*pi*(1:100)/20) + rnorm(100, 0, 0.3)
  result <- khaldoon_test(x, n_bootstrap = 100)
  
  expect_gte(result$components$peak_to_mean, 0)
  expect_gte(result$components$significant_peaks, 0)
  expect_gte(result$components$energy_concentration, 0)
  expect_lte(result$components$energy_concentration, 1)
})

