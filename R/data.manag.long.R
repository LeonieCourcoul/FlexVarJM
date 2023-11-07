#' Management of longitudinal data
#'
#' @param formGroup A formula which indicates the group variable
#' @param formFixed A formula which indicates the fixed effects for the longitudinal submodel
#' @param formRandom A formula which indicates the random effects for the longitudinal submodel
#' @param data.long1 A dataframe with the longitudinal data
#'
#' @return A list with the following components :
#' \describe{
#' \item{\code{data_long}}{a clean dataframe for the longitudinal data}
#' \item{\code{y.new.prog}}{the vector of responses variable}
#' \item{\code{X}}{a matrix with the fixed effects}
#' \item{\code{U}}{a matrix with the random effects}
#' \item{\code{id}}{a vector with the identification of individuals}
#' \item{\code{offset}}{a vector with the number of measurements for each individual}
#' \item{\code{I}}{an integer, the number of individuals}
#' }
#' @importFrom stats model.frame model.matrix

data.manag.long <- function(formGroup, formFixed, formRandom, data.long1){

  data_long <- data.long1[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
  #y.new.prog <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  y.new.prog <- mfX[,1]
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  if(!("id" %in% colnames(data_long))) #To have a column named "id"
    data_long <- cbind(data_long, id = id)

  list.long <- list("data_long"= data_long, "y.new.prog" = y.new.prog, "X" = X, "U" = U,
                    "id" = id, "offset"=offset, "I" = I)

  return(list.long)
}
