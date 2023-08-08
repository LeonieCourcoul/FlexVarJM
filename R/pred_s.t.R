#' Predictions computation
#'
#' @param newdata 
#' @param x 
#' @param s 
#' @param window 
#' @param event 
#'
#' @return
#' @export
#'
#' @examples

pred_s.t <- function(object, newdata, s, window, event = 1, IC = NULL, draws=NULL){
  x <- object
  if(!inherits(x, "lsjm")) stop("use only \"lsjm\" objects")
  if(IC<=0 || IC>=1) stop("IC must be between 0 and 1")
  if(!is.null(IC) && (is.null(draws) || draws <=0)) stop("draw must be higher 1")
  param <- x$table.res$Estimation
  pred.true <- pred_s.t.aux(x, newdata, s, window, event, param)
  result <- cbind(unique(newdata[,all.vars(x$control$formGroup)]),pred.true)
  result <- as.data.frame(result)
  colnames(result) <- c("ID", "Prediction")
  if(!is.null(IC)){
    pred.IC <- matrix(ncol = draws, nrow = length(unique(newdata[,all.vars(x$control$formGroup)])))
    Hess <- matrix(rep(0,length(x$result$grad)**2),nrow=length(x$result$grad),ncol=length(x$result$grad))
    Hess[upper.tri(Hess, diag=T)] <- x$result$v
    Hess2 = Hess + t(Hess)
    diag(Hess2) <- diag(Hess2) - diag(Hess)
    for(la in 1:draws){
      tirage <- mvtnorm::rmvnorm(1, mean = x$table.res$Estimation, sigma = Hess2)
      pred.la <- pred_s.t.aux(x, newdata, s, window, event, tirage)
      pred.IC[,la] <- pred.la
    }
    qinf <- (1-IC)/2
    qsup <- (1+IC)/2
    Cinf <- apply(pred.IC, 1, function(x) quantile(x,qinf))
    Csup <- apply(pred.IC, 1, function(x) quantile(x,qsup))
    Med <-  apply(pred.IC, 1, median)
    sdemp <- apply(pred.IC, 1, sd)
    result$Median <- Med
    result$IC_inf <- Cinf
    result$IC_sup <- Csup
    result$sd <- sdemp
  }
  rownames(result) <- c()
  result
}


