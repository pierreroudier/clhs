context("clhs-gower")

test_that("Gower similarity works", {
  
  library(raster)
  library(sp)

  slogo <- stack(system.file("external/rlogo.grd", package = "raster"))
  spdf <- SpatialPointsDataFrame(data.frame(x = 50, y = 50), data = data.frame(ID = 1))

  gw <- similarity_buffer(covs = slogo, pts = spdf, buffer = 1, fac = NA)

  expect_equal(as.numeric(na.omit(values(gw))), c(0.00000000, 0.08841009, 0.70043470, 1.00000000))
})
