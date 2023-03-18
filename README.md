# clhs development

This repo is a fork of the main clhs package, used for development (especially of the C++ functions). The original C++ version is on CRAN, but there are multiple updates to the development version, including bug fixes and specification of minimum distance between points. 

## Installation

The CRAN package currently contains a bug if the C++ version is used along with `must.include`.

You can install the development version using the `devtools` package:

`devtools::install_github("kdaust/clhs")`
