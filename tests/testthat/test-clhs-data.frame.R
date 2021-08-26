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
  
  res <- clhs(df, size = 5, iter = 100, use.cpp = FALSE, progress = FALSE)
  
  # expect_equal(res, c(188, 657, 140, 301, 817))
  expect_equal(res, c(28, 466, 419, 700, 536))
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
  
  expect_equal(res, c(632, 968, 914, 503, 224))
})

test_that("results using C++ and R implementations match", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE),
    stringsAsFactors = TRUE
  )
  
  res_r <- clhs(df, size = 5, use.cpp = FALSE, simple = FALSE, iter = 7500, progress = FALSE)
  res_cpp <- clhs(df, size = 5, use.cpp = TRUE, simple = FALSE, iter = 7500, progress = FALSE)
  
  # Trying to comapre R and C++ results based on objective function values
  # 
  
  # range of values observed
  rg <- range(
    c(range(res_r$obj), range(res_cpp$obj))
  )
  # Compute proportional difference
  obj_diff <- abs(res_r$obj[length(res_r$obj)] - res_cpp$obj[length(res_cpp$obj)])
  prop_obj_diff <- obj_diff / diff(rg)
  
  # Tolerance of a 5% difference
  max_diff <- 0.05
    
  expect_true(prop_obj_diff <= max_diff)
})

