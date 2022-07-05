#' Title
#'
#' @param x
#' @param group
#'
#' @return
#'
#' @examples
fastSumID <- function(x,group){
  as.vector(x = rowsum.default(x,group,reorder = FALSE), mode = "numeric")
}
