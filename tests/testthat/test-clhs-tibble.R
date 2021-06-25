context("tibbles")

test_that("CLHS (R) works for tibbles", {
  
  data(diamonds, package = "ggplot2")
  
  expect_type(clhs(diamonds[1:30,], size = 3, iter = 250, simple = FALSE), "list")
  
})

test_that("CLHS (C++) works for tibbles", {
  
  data(diamonds, package = "ggplot2")
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  res <- clhs(diamonds, size = 3, use.cpp = T, iter = 5000, simple = FALSE)
  
  expect_equal(
    res$index_samples, c(15338, 6921, 30750)
  )
  
})
