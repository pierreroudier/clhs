context("tibbles")

test_that("CLHS works for tibbles", {
  
  data(diamonds, package = "ggplot2")
  
  set.seed(1)
  res <- clhs(diamonds[1:30,], size = 3, iter = 250, simple = FALSE)
  
  expect_equal(
    res$index_samples, c(6, 17, 29)
  )
  
})
