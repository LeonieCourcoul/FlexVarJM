#' Management of data for longitudinal submodel
#'
#' @param data.id A dataframe
#' @param Time A vector of Time of events
#' @param formFixed A formula for the fixed effects of the longitudinal submodel
#' @param formRandom  A formula for the random effects of the longitudinal submodel
#' @param timeVar The name of the column of time in data.id. This variable must appears in data.id
#'
#' @return A list with the following components
#' \describe{
#' \item{\code{Xtime}}{a matrix of fixed effects at each time of measure}
#' \item{\code{Utime}}{a matrix of random effects at each time of measure}
#' }
#' @importFrom stats model.frame model.matrix
#' @export
#'
#' @examples
#'
#'
#'
#'
#'
data.time <- function(data.id, Time, formFixed, formRandom, timeVar){

  if (!timeVar %in% names(data.id))
    stop("\n'timeVar' does not correspond to one of the columns in formulas")
  data.id[[timeVar]] <- Time
  mfX.id <- model.frame(formFixed, data = data.id)
  Xtime <- model.matrix(formFixed, mfX.id)
  mfU.id <- model.frame(formRandom, data = data.id)
  Utime <- model.matrix(formRandom, mfU.id)

  list("Xtime" = Xtime, "Utime" = Utime)
}
