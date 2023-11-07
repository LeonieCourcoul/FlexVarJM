pred.re <- function(param, nb.e.a = NULL, variability_hetero = NULL, nb.e.a.sigma = NULL,
                       Sigma.re = NULL, X_base_i = NULL, U_i = NULL, beta = NULL, omega = NULL, O_base_i = NULL,
                       W_base_i = NULL, y_i = NULL, sigma.epsilon = NULL, Otime_i = NULL, Wtime_i = NULL,
                       Os_i = NULL, Ws_i = NULL, S = NULL, alpha.sigma = NULL, competing_risk = NULL, alpha.sigma.CR = NULL,
                       sharedtype = NULL, sharedtype_CR = NULL, alpha.current = NULL, alpha.current.CR = NULL,
                       alpha.slope = NULL, alpha.slope.CR = NULL, Xtime_i = NULL, Utime_i = NULL, Xs_i = NULL, Us_i = NULL,
                       indices_beta_slope = NULL, hazard_baseline = NULL, wk = NULL, st_i = NULL, gamma = NULL, B_i = NULL, Bs_i = NULL,
                       Z_i = NULL, alpha = NULL, shape = NULL, Time_i = NULL, P_i = NULL, hazard_baseline_CR = NULL, gamma.CR = NULL, B.CR_i = NULL,
                       Bs.CR_i = NULL, Z.CR_i = NULL, alpha.CR = NULL, shape.CR = NULL, event1_i = NULL, event2_i = NULL, Xs.slope_i = NULL, Us.slope_i = NULL,
                       Xslope_i = NULL, Uslope_i = NULL
){
  
  b_re <- param[1:nb.e.a]
  if(variability_hetero){
    tau_re <- param[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]
    #### f(b_i,tau_i)
    f_b_tau <- mvtnorm::dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
    #### f(Y_i|b_i,tau_i)
    if(is.null(nrow(X_base_i))){
      CV <- (beta%*%X_base_i)[1,1] + b_re%*%U_i
      sigma.long <- exp((omega%*%O_base_i)[1,1] + tau_re%*%W_base_i)
      f_Y_b_sigma <- dnorm(x=y_i, mean = CV, sd = sigma.long)
    }
    else{
      f_Y_b_sigma <- 1
      for(k in 1:nrow(X_base_i)){
        sigma.long <- exp((omega%*%O_base_i[k,])[1,1] + tau_re%*%W_base_i[k,])
        CV <- (beta%*%X_base_i[k,])[1,1] + b_re%*%U_i[k,]
        f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y_i[k], mean = CV, sd = sigma.long)
      }
    }
  }
  else{
    #### f(b_i)
    f_b_tau <- mvtnorm::dmvnorm(x = b_re, mean = rep(0,length(b_re)), sigma = Sigma.re)
    #### f(Y_i|b_i)
    if(is.null(nrow(X_base_i))){
      CV <- (beta%*%X_base_i)[1,1] + b_re%*%U_i
      f_Y_b_sigma <- dnorm(x=y_i, mean = CV, sd = sigma.epsilon)
    }
    else{
      f_Y_b_sigma <- 1
      for(k in 1:nrow(X_base_i)){
        CV <- (beta%*%X_base_i[k,])[1,1] + b_re%*%U_i[k,]
        f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y_i[k], mean = CV, sd = sigma.epsilon)
      }
    }
  }
  
  #### f(T_i|...)
  h <- 1
  etaBaseline <- 0
  survLong <- 0
  etaBaseline.0 <- 0
  survLong.0 <- 0
  # variability
  if("variability" %in% sharedtype){
    Sigma.CV <- exp((omega%*%Otime_i)[1,1] + tau_re%*%Wtime_i)
    Sigma.current.GK <- exp(omega%*%t(Os_i) + tau_re%*%t(Ws_i))
    h <- h*exp(alpha.sigma*Sigma.CV)
    survLong <- survLong + alpha.sigma*Sigma.current.GK
  }
  if(competing_risk){
    h_CR <- 1
    etaBaseline_CR <- 0
    survLong_CR <- 0
    etaBaseline.0_CR <- 0
    survLong.0_CR <- 0
    if("variability" %in% sharedtype_CR){
      Sigma.CV <- exp((omega%*%Otime_i)[1,1] + tau_re%*%Wtime_i)
      Sigma.current.GK <- exp(omega%*%t(Os_i) + tau_re%*%t(Ws_i))
      h_CR <- h_CR*exp(alpha.sigma.CR*Sigma.CV)
      survLong_CR <- survLong_CR + alpha.sigma.CR*Sigma.current.GK
    }
  }
  # current value
  if( ("current value" %in% sharedtype) || (competing_risk && "current value" %in% sharedtype_CR)){
    CV <- (beta%*%Xtime_i)[1,1]+b_re%*%Utime_i
    current.GK <- beta%*%t(Xs_i) + b_re%*%t(Us_i)
    if("current value" %in% sharedtype){
      h <- h*exp(alpha.current*CV)
      survLong <- survLong + alpha.current*current.GK
    }
    if(competing_risk && "current value" %in% sharedtype_CR){
      h_CR <- h_CR*exp(alpha.current.CR*CV)
      survLong_CR <- survLong_CR + alpha.current.CR*current.GK
    }
  }
  # slope
  if( ("slope" %in% sharedtype) || (competing_risk && "slope" %in% sharedtype_CR)){
    slope.GK <- beta[indices_beta_slope]%*%t(Xs.slope_i) + b_re[-1]%*%t(Us.slope_i)
    if(length(indices_beta_slope) == 1){
      slope <- (beta[indices_beta_slope]%*%Xslope_i)[1,1]+b_re[-1]*Uslope_i
    }
    else{
      slope <- (beta[indices_beta_slope]%*%Xslope_i)[1,1]+b_re[-1]%*%Uslope_i
    }
    if("slope" %in% sharedtype){
      h <- h*exp(alpha.slope*slope)
      survLong <- survLong + alpha.slope*slope.GK
    }
    if("slope" %in% sharedtype_CR){
      h_CR <- h_CR*exp(alpha.slope.CR*CV)
      survLong_CR <- survLong_CR + alpha.slope.CR*slope.GK
    }
  }
  #h0
  if(hazard_baseline == "Exponential"){
    h_0 <- 1
    h_0.GK <- wk
  }
  if(hazard_baseline == "Weibull"){
    h_0 <- shape*(Time_i**(shape-1))
    h_0.GK <- shape*(st_i**(shape-1))*wk
  }
  if(hazard_baseline == "Splines"){
    h_0 <- exp((gamma%*%B_i)[1,1])
    mat_h0s <- matrix(gamma,ncol=1)
    h_0.GK <- (wk*exp(Bs_i%*%mat_h0s))
  }
  if(length(Z_i)==0){
    pred_surv <- 0
  }
  else{
    pred_surv <- (alpha%*%Z_i)[1,1]
  }
  h <- h_0*exp(pred_surv)*h
  etaBaseline <- etaBaseline + pred_surv
  survLong <- exp(survLong)%*%h_0.GK
  Surv <- (-exp(etaBaseline)*P_i*survLong)
  
  if(competing_risk){
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
    if(length(Z.CR_i)==0){
      pred_surv.CR <- 0
    }
    else{
      pred_surv.CR <- (alpha.CR%*%Z.CR_i)[1,1]
    }
    h_CR <- h_0.CR*exp(pred_surv.CR)*h_CR
    etaBaseline_CR <- etaBaseline_CR + pred_surv.CR
    survLong_CR <- exp(survLong_CR)%*%h_0.GK.CR
    Surv.CR <- (-exp(etaBaseline_CR)*P_i*survLong_CR)
  }
  
  if(competing_risk){
    fct <- log(f_Y_b_sigma*f_b_tau) + log(h**event1_i) + Surv + log(h_CR**event2_i) + Surv.CR
  }
  else{
    fct <- log(f_Y_b_sigma*f_b_tau) + log(h**event1_i) + Surv
  }
  fct
  
}
