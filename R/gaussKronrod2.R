#' Initialisation of Survival Data at Gauss Kronrod time points 2
#'
#' @param data.id A database with covariates of interest and 1 line per subject
#' @param k The number of Gauss Kronrod points, by default k = 15
#' @param a First born
#' @param b Second born
#'
#' @return A list with the following components :
#' \describe{
#' \item{\code{K}}{an integer, the number of points}
#' \item{\code{P}}{a vector, of value Time/2}
#' \item{\code{st}}{a matrix with nrow = number of subjects and ncol = k. The new time to compute the survival function}
#' \item{\code{wk}}{a vector of weights}
#' \item{\code{data.id2}}{a database with K lines per subjects}
#' \item{\code{id.GK}}{the vector of IDs}
#' }
#'

data.GaussKronrod2 <- function(data.id, a, b, k = 15){
  #print("new")
  wk <- gaussKronrod()$wk
  sk <- gaussKronrod()$sk
  K <- length(sk)
  P <- (b-a)/2
  st <- outer(P, sk+1)+a
  id.GK <- rep(seq_along(data.id$id), each = K)
  data.id2 <- data.id[id.GK, ]

  list(K = K, P = P, st = st, wk = wk, data.id2 = data.id2,
       id.GK = id.GK)

}
