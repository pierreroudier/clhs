context("clhs-sp")

test_that("clhs on a SpatialPointsDataFrame works", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  
  df <- data.frame(
    a = runif(1000), 
    b = rnorm(1000), 
    c = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3)),
    x = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3)),
    y = c(rnorm(n = 500, mean = 20, sd = 3), rnorm(n = 500, mean = -20, sd = 3))
  )
  
  spdf <- sp::SpatialPointsDataFrame(
    coords = df[, c("x", "y")],
    data = df[, c("a", "b", "c")],
    proj4string = sp::CRS("+init=epsg:4326")
  )
  
  res1 <- clhs(spdf, size = 5, iter = 100, progress = FALSE, simple = TRUE)
  res2 <- clhs(spdf, size = 5, iter = 100, progress = FALSE, simple = TRUE, use.coords = TRUE)
  
  expect_equal(res1, c(809, 559, 264, 898, 891))
  expect_equal(res2, c(618, 547, 123, 230, 233))
})
