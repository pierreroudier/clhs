#' @importFrom methods isGeneric setGeneric
if (!isGeneric("clhs"))
  setGeneric("clhs", function(x, ...)
    standardGeneric("clhs"))
