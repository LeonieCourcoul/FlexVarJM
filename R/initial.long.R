#' Initialisation of Longitudinal Submodel
#'
#' @param formFixed A formula which indicates the fixed effects for the longitudinal submodel
#' @param formRandom A formula which indicates the random effects for the longitudinal submodel
#' @param idVar A character, indicates the name of the group variable
#' @param data.long1 A dataframe with the longitudinal data
#' @param ncX An integer, the number of columns of matrix X, ie, the number of fixed effects
#' @param nproc An integer, the number of cores for parallel computation
#'
#' @return A list with the following components :
#' \describe{
#' \item{\code{long_model}}{the result of the hlme function}
#' \item{\code{priorMean.beta}}{the estimated parameters for fixed effects in the linear mixed effects model}
#' \item{\code{sigma}}{the estimated sigma of the model}
#' }
#' @importFrom lcmm hlme
#'

initial.long <- function(formFixed, formRandom, idVar, data.long1, ncX, nproc = nproc){

  long_model <- hlme(fixed = formFixed,
                     random= formRandom,
                     subject = idVar,
                     data=data.long1,
                     nproc = nproc,
                     verbose = FALSE)
  priorMean.beta <- long_model$best[1:ncX]
  #priorTau.beta <- diag(rep(precision,length(priorMean.beta)))
  sigma <- long_model$best["stderr"]

  list.init.long <- list("long_model" = long_model, "priorMean.beta" = priorMean.beta,
                         "sigma" = sigma)

  return(list.init.long)

}
