#' Function to predict the random effects
#'
#' @param param
#' @param nb.e.a
#' @param X_base_i
#' @param beta
#' @param U_i
#' @param y_i
#' @param Sigma.b
#' @param mu.log.sigma
#' @param tau.log.sigma
#' @param variability_hetero
#' @param alpha.sigma
#' @param competing_risk
#' @param alpha.sigma.CR
#' @param sharedtype
#' @param sharedtype_CR
#' @param Xtime_i
#' @param Utime_i
#' @param Xs_i
#' @param Us_i
#' @param alpha.current
#' @param alpha.current.CR
#' @param hazard_baseline
#' @param wk
#' @param shape
#' @param Time_i
#' @param st_i
#' @param gamma
#' @param B_i
#' @param Bs_i
#' @param Z_i
#' @param alpha
#' @param P_i
#' @param hazard_baseline_CR
#' @param shape.CR
#' @param gamma.CR
#' @param B.CR_i
#' @param Bs.CR_i
#' @param Z.CR_i
#' @param alpha.CR
#' @param event1_i
#' @param event2_i
#'
#' @return
#' @export
#'
#' @examples
#'


predict_ea <- function(param, nb.e.a = 0, X_base_i = NULL, beta = NULL, U_i = NULL, y_i = NULL,
                       Sigma.b = NULL, mu.log.sigma = NULL, tau.log.sigma = NULL,variability_hetero = TRUE,
                       alpha.sigma = NULL, competing_risk = FALSE, alpha.sigma.CR = NULL,sharedtype = NULL, sharedtype_CR = NULL,
                       Xtime_i = NULL, Utime_i = NULL, Xs_i = NULL, Us_i = NULL,
                       alpha.current = NULL, alpha.current.CR = NULL,hazard_baseline = NULL,wk = NULL,
                       shape = 0,Time_i = NULL,st_i = NULL,gamma = NULL, B_i = NULL, Bs_i = NULL,
                       Z_i = NULL, alpha = NULL, P_i = NULL, hazard_baseline_CR = NULL, shape.CR = NULL,
                       gamma.CR = NULL, B.CR_i = NULL, Bs.CR_i = NULL,Z.CR_i = NULL,alpha.CR = NULL,event1_i = NULL, event2_i = NULL
){
  b_ea <- param[1:(nb.e.a)]
  sigma_ea <- exp(param[nb.e.a+1])
  if(is.null(nrow(X_base_i))){
    CV <- (beta%*%X_base_i)[1,1] + b_ea%*%U_i
    f_Y_b_sigma <- dnorm(x=y_i, mean = CV, sd = sigma_ea)#(1/(sqrt(2*pi)*sigma_al))*exp((-1/2)*((y_i-CV)/sigma_al)**2) #
  }
  else{
    f_Y_b_sigma <- 1
    for(k in 1:nrow(X_base_i)){
      CV <- (beta%*%X_base_i[k,])[1,1] + b_ea%*%U_i[k,]
      f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y_i[k], mean = CV, sd = sigma_ea)# exp((-1/2)*((y_i[k]-CV)/sigma_al)**2) #

    }
  }

  ### f(b_i)
  f_b <- mvtnorm::dmvnorm(x = b_ea, mean = rep(0,length(b_ea)), sigma = Sigma.b)

  ### f(sigma_i)
  f_sigma <- (1/sigma_ea)*dnorm(log(sigma_ea), mean = mu.log.sigma, sd = tau.log.sigma)


  ### f(T_i|...)
  ###### lambda(T_i|...)
  h <- 1
  etaBaseline <- 0
  survLong <- 0
  etaBaseline.0 <- 0
  survLong.0 <- 0
  if(variability_hetero){
    h <- h*exp(alpha.sigma*sigma_ea)
    etaBaseline <- etaBaseline + alpha.sigma*sigma_ea
  }
  if(competing_risk){
    h_CR <- 1
    etaBaseline_CR <- 0
    survLong_CR <- 0
    etaBaseline.0_CR <- 0
    survLong.0_CR <- 0
    if(variability_hetero){
      h_CR <- h_CR*exp(alpha.sigma.CR*sigma_ea)
      etaBaseline_CR <- etaBaseline_CR + alpha.sigma.CR*sigma_ea
    }
  }
  if(sharedtype %in% c("CV","CVS") || (competing_risk && sharedtype_CR %in% c("CV","CVS")) ){
    CV <- (beta%*%Xtime_i)[1,1]+b_ea%*%Utime_i
    current.GK <- beta%*%t(Xs_i) + b_ea%*%t(Us_i)
    if(sharedtype %in% c("CV","CVS")){
      h <- h*exp(alpha.current*CV)
      survLong <- survLong + alpha.current*current.GK
    }
    if(competing_risk && sharedtype_CR %in% c("CV","CVS")){
      h_CR <- h_CR*exp(alpha.current.CR*CV)
      survLong_CR <- survLong_CR + alpha.current.CR*current.GK
    }
  }

  ###h0
  if(hazard_baseline == "Exponential"){
    h_0 <- 1
    h_0.GK <- wk
  }

  if(hazard_baseline == "Weibull"){
    h_0 <-  shape*(Time_i**(shape-1))
    h_0.GK <- shape*(st_i**(shape-1))*wk
  }

  if(hazard_baseline == "Splines"){
    h_0 <- exp((gamma%*%B_i)[1,1])
    mat_h0s <- matrix(gamma,ncol=1)
    h_0.GK <- (wk*exp(Bs_i%*%mat_h0s))
  }


  ###hazard function
  if(length(Z_i)==0){
    pred_surv <- 0
  }
  else{
    pred_surv <- (alpha%*%Z_i)[1,1]
  }
  h <- h_0*exp(pred_surv)*h
  etaBaseline <- etaBaseline + pred_surv
  ###GK integration
  survLong <- exp(survLong)


  survLong <- survLong%*%h_0.GK

  Surv <- exp((-exp(etaBaseline)*P_i*survLong))


  if(competing_risk){
    ###h0
    if(hazard_baseline_CR == "Exponential"){
      h_0.CR <- 1
      h_0.GK.CR <- wk
    }

    if(hazard_baseline_CR == "Weibull"){
      h_0.CR <- shape.CR*(Time_i**(shape.CR-1))
      h_0.GK.CR <- shape.CR*(st_i**(shape.CR-1))*wk
    }

    if(hazard_baseline_CR == "Splines"){
      h_0.CR <- exp((gamma.CR%*%B.CR_i)[1,1])
      mat_h0s.CR <- matrix(gamma.CR,ncol=1)
      h_0.GK.CR <- (wk*exp(Bs.CR_i%*%mat_h0s.CR))
    }



    ###hazard function
    if(length(Z.CR_i)==0){
      pred_surv.CR <- 0
    }
    else{
      pred_surv.CR <- (alpha.CR%*%Z.CR_i)[1,1]
    }
    h_CR <- h_0.CR*exp(pred_surv.CR)*h_CR
    etaBaseline_CR <- etaBaseline_CR + pred_surv.CR
    ###GK integration
    survLong_CR <- exp(survLong_CR)
    survLong_CR <- survLong_CR%*%h_0.GK.CR
    Surv.CR <- exp((-exp(etaBaseline_CR)*P_i*survLong_CR))
  }
  if(competing_risk){
    fct <- f_Y_b_sigma*f_sigma*f_b*(h**event1_i)*Surv*(h_CR**event2_i)*Surv.CR #f_sigma*
  }
  else{
    fct <- f_Y_b_sigma*f_b*f_sigma*(h**event1_i)*Surv
  }
  fct
}
