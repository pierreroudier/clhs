context("clhs-sf")

test_that("clhs on a sf works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3)),
    x = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3)),
    y = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3))
  )
  
  sf <- sf::st_as_sf(
    df,
    coords = c("x", "y"),
    crs = 4326
  )
  
  res1 <- clhs(sf, size = 5, iter = 100, progress = FALSE, simple = TRUE)
  res2 <- clhs(sf, size = 5, iter = 100, progress = FALSE, simple = TRUE, use.coords = TRUE)
  
  expect_equal(res1, c(809, 559, 264, 898, 891))
  expect_equal(res2, c(618, 547, 123, 230, 233))
})
