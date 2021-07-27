context("clhs-raster")

test_that("clhs works a raster", {
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  df <- data.frame(
    a = runif(10000), 
    b = rnorm(10000)
  )
  
  coords <- expand.grid(x = 1:100, y = 1:100)
  
  rdf <- data.frame(coords, df)
  r <- raster::rasterFromXYZ(rdf)
  
  set.seed(1)
  res <- clhs(r, size = 5, use.cpp = F, iter = 100, progress = FALSE, simple = FALSE)
  
  expect_equal(
    res$sampled_data@coords, 
    structure(
      c(64, 78, 73, 47, 74, 51, 98, 78, 17, 60), 
      .Dim = c(5L, 2L), 
      .Dimnames = list(NULL, c("x", "y"))
    )
  )
  
})
