context("clhs-data.frame_method_DLHS")

test_that("DLHS method on a data.frame works (R)", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE),
    stringsAsFactors = TRUE
  )
  
  # Weigth of sampling from 1 for the middle strata to 3 for the edges
  eta <- matrix(c(3, 2, 1, 2, 3), ncol = 2, nrow = 5)
  res <- clhs(df, size = 5, iter = 100, use.cpp = F, progress = FALSE, eta = eta)
  
  expect_equal(res, c(919, 23, 349, 411, 308))
})

test_that("DLHS method on a data.frame works (C++)", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE),
    stringsAsFactors = TRUE
  )
  
  # Weight of sampling from 1 for the middle strata to 3 for the edges
  eta <- matrix(c(3, 2, 1, 2, 3), ncol = 2, nrow = 5)
  res <- clhs(df, size = 5, iter = 5000, progress = FALSE, eta = eta)
  
  expect_equal(res, c(243, 270, 742, 220, 165))
})
