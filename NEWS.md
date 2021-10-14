<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

- Fixed slow example in the plot docs
- fixed the way RcppAramadillo dependecny is handled so to pass CRAN tests
- added Rcpp and RcppArmadillo dependencies
- Added C++ implementation (much faster!)


# clhs 0.8.1.9000

- Same as previous version.


# clhs 0.8.1

- the behaviour of the `simple` option is now more consistent for spatial objects
- added `use.coords` option to optionally use spatial coordinates of sf, sp, or raster objects in clhs
- added support for sf objects
- Added test for clhs.Raster

# clhs 0.8-0

- Added a `NEWS.md` file to track changes to the package.
- Implemented an improved logic to drop worse sample (thanks to David Clifford)
- Various minor bug fixes to support R > 4.0

# clhs_0.7-2 
- Maintenance version fixing new RNG sampler call.

# clhs_0.7-0 
- New DLHS method contributed by Benjamin Louis
- Included possibility to include compulsory/existing samples in the set (thanks to Benjamin Louis)
- Added more tests

# clhs_0.6-0
- Added Gower similarity tool (thanks to Colby Brugnard)
- Switched documentation to Roxygen
- Implemented the first few tests

# clhs_0.5-7
- fixed sneaky bug that produced duplicated sampling points (thanks Ankur Gupta for the pointers)

# clhs_0.5-6
- fixed bug when passing a single raster layer (Github issue #1, thanks Github user GreatEmerald)
- Raster* clhs method now returns SpatialPointsDataFrame

# clhs_0.5-3
- corrected dependencies - now depend on R >= 2.14 (due to ggplot2)

# clhs_0.5-2
- version ready for release when ggplot2_0.9.2.1 is rolled out

# clhs_0.5-1
- minor bugfix version fixing compatibility issues with ggplot2 >= 0.9.2
- added ggplot2 in Suggests temporary to fix bug in ggplot2_0.9.2

# clhs_0.5-0 
- introduced cost and cost tracking modes
- general code cleaning
- improved plot function:
  - general buxfixes
  - now using dotplot for factors
  - support cost function
  - new boxplot mode
- cleaner NAMESPACE

# clhs_0.4-3
- Dummy version increment to solve CRAN upload problems

# clhs_0.4-2
- Added the choice between histogram and density plots in
the plot method

# clhs_0.4-1 
- Improved doc
- Complete plot.cLHS_result method using ggplot2
- new reshape2 dependency
- various buxfixes, esp. for spatial classes

# clhs_0.4-0
- Introduced cLHS_result S3 class with associated plot method
- Introduced simple = ... option to the clhs() method. If set to true, returns only the indices of the selected samples, if set to FALSE, returns a cLHS_result object (eg if you want to plot the objective function behaviour).
- added a plot() method. For the moment, it just plots the objective function.

# clhs_0.3-2
- slight improvement on the Raster method using rasterToPoints(...)
- added plotting option for the objective function

# clhs_0.3-1
- Corrected bug on Raster* methods. Now returns a SpatialPointsDataFrame object.

# clhs_0.3-0
- Switch to S4 methods
- Introduced methods for Raster* and SpatialP*DataFrame classes

# clhs_0.2-1
- changed examples

# clhs_0.2-0
- Fixed important bugs in the annealing process.

# clhs_0.1-0
- Initial release of the package.
