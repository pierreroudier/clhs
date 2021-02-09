context("costs-constrained implemenatation")

test_that("Cost-constrained implementation (R) works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  res <- clhs(mtcars, size = 3, use.cpp = F, cost = "mpg", iter = 250, simple = FALSE)
  
  expect_equal(
    res$index_samples, c(31, 32, 15)
  )
  
  expect_equal(
    round(min(res$obj), digits = 2), 30.92
  )
  
  expect_equal(
    min(res$cost), 46.8
  )
  
})

test_that("Cost-constrained implementation (C++) works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(3.145)
  res <- clhs(mtcars, size = 3, use.cpp = T, cost = "mpg", simple = FALSE)
  
  expect_equal(
    res$index_samples, c(26, 31, 15)
  )
  
  expect_equal(
    round(min(res$obj), digits = 2), 29.21
  )
  
  expect_equal(
    min(res$cost), 52.7
  )
  
})
