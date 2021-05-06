context("tibbles")

test_that("CLHS works for tibbles", {
  
  data(diamonds, package = "ggplot2")
  
  expect_type(clhs(diamonds[1:30,], size = 3, iter = 250, simple = FALSE), "list")
  
})
