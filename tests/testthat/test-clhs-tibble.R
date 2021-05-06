context("tibbles")

test_that("CLHS works for tibbles", {
  
  data(diamonds, package = "ggplot2")
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  res <- clhs(diamonds[1:30,], size = 3, iter = 250, simple = FALSE)
  
  expect_equal(
    res$index_samples, c(2, 9, 15)
  )
  
})
