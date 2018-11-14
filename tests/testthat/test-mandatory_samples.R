context("mandatory_samples")

test_that("Mandatory samples are actually selected", {
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = sample(LETTERS[1:5], size = 1000, replace = TRUE)
  )
  
  mandatory_idx <- sample(1:nrow(df), size = 3)
  
  # no error
  res <- clhs(df, size = 10, iter = 100, include = mandatory_idx, progress = FALSE)
  expect_true(all(mandatory_idx %in% res))
  
  # error : size <= length(include)
  expect_error(clhs(df, size = 3, iter = 100, include = mandatory_idx, progress = FALSE))
  
})
