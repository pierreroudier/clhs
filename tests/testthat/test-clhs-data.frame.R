context("clhs-data.frame")

test_that("basic clhs on a data.frame works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE),
    stringsAsFactors = TRUE
  )
  
  res <- clhs(df, size = 5, iter = 100, use.cpp = F, progress = FALSE)
  
  expect_equal(res, c(188, 657, 140, 301, 817))
})

test_that("basic clhs using C++ on a data.frame works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE),
    stringsAsFactors = TRUE
  )
  
  res <- clhs(df, size = 5)
  
  expect_equal(res, c(219, 876, 510, 39, 944))
})
