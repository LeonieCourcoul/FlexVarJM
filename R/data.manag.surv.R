#' Management of survival data
#'
#' @param formGroup A formula which indicates the group variable
#' @param formSurv A formula which indicates the variables used in the survival submodel
#' @param data.long1 Database
#' @param formSurv_CompRisk A formula which indicates the variables used in the competing survival submodel
#'
#' @return A list with the following components :
#' \describe{
#' \item{\code{tmp}}{the final database for survival analysis}
#' \item{\code{Time}}{a vector of observed times}
#' \item{\code{event1}}{a vector of first event indicator}
#' \item{\code{nTime}}{length of Time vector}
#' \item{\code{Z}}{matrix of covariables of first survival submodel}
#' \item{\code{event2}}{a vector of second event indicator}
#' \item{\code{Z_CR}}{matrix of covariables of second survival submodel}
#' }
#' @import stats
#' @import survival
#' @export
#'
#' @examples
data.manag.surv <- function(formGroup, formSurv, data.long1,formSurv_CompRisk){
  tmp <- data.long1[unique(c(all.vars(formGroup),all.vars(formSurv),all.vars(formSurv_CompRisk)))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(formSurv)][, 2]  # vector of event indicator (delta)
  event1 <- ifelse(event == 1, 1,0)
  nTime <- length(Time)                   # number of subject having Time
  zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- model.frame(formSurv, data = tmp)
  Z <- model.matrix(formSurv, mfZ)
  event2 = NULL
  Z_CR = NULL
  if(!is.null(formSurv_CompRisk)){
    event2 <- ifelse(event == 2, 1,0)
    mfZ_CR <- model.frame(formSurv_CompRisk, data = tmp)
    Z_CR <- model.matrix(formSurv_CompRisk, mfZ_CR)
  }
  list("tmp" = tmp, "Time" = Time , "event1"= event1, "nTime" = nTime,
       "Z" = Z, "event2" = event2, "Z_CR" = Z_CR)

}
