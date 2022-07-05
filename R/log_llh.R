#' Log-likelihood computation
#'
#' @param param a vector : paramaters to be estimated
#' @param nb.e.a integer : number of RE
#' @param nb.priorMean.beta integer : number of fixed effects
#' @param nb.alpha integer : number of covariates in survival model
#' @param competing_risk boolean : allow competing risk or not, FALSE by default
#' @param nb.alpha.CR integer : number of covariates in survival model for competing risks
#' @param variability_hetero boolean : allow the heterogeneous variability or not
#' @param S integer : the number of QMC points
#' @param Zq vector : sobol points
#' @param sharedtype vector : dependence structure for survival model : "RE" (random effects) or "CV" (current value) or "CVS" (current value and slope) or "S" (slope)
#' @param sharedtype_CR vector : dependence structure for competing risk survival model : "RE" (random effects) or "CV" (current value) or "CVS" (current value and slope) or "S" (slope)
#' @param hazard_baseline char : baseline hazard function : "Exponential" or "Weibull" or "Splines"
#' @param hazard_baseline_CR char : baseline hazard function, competing risk : "Exponential" or "Weibull" or "Splines"
#' @param ord.splines integer : the order of splines function for baseline hazard function
#' @param Xtime matrix : fixed effects at event time
#' @param Utime matrix : RE at event time
#' @param nb_pointsGK integer : number of points for Gauss-Kronrod approximation, 7 or 15 (default)
#' @param Xs matrix : fixed effects at Gauss-Kronrod times
#' @param Us matrix : RE at Gauss-Kronrod times
#' @param Xslope matrix : fixed effects of slope at event times
#' @param Uslope matrix : RE of slope at event times
#' @param Xs.slope matrix : fixed effects of slope at Gauss-Kronrod times
#' @param Us.slope matrix : RE of slope at Gauss-Kronrod times
#' @param indices_beta_slope vector : position of beta which will be used in the slope computation
#' @param Time vector : observed event times
#' @param st_calc matrix : Gauss-Kronrod times
#' @param B matrix : splines for baseline hazard function of event 1
#' @param Bs matrix : splines for baseline survival function of event 1
#' @param wk vector : Gauss-Kronrod weights
#' @param Z matrix : covariates for survival function of event 1
#' @param P vector : Time/2
#' @param left_trunc boolean : left truncation indicator
#' @param Z_CR matrix : covariates for survival function of event 2
#' @param X_base matrix : fixed effects for longitudinal submodel
#' @param offset vector : number of lines per subjects
#' @param U matrix : RE for longitudinal submodel
#' @param y.new.prog vector : y measures for longitudinal submodel
#' @param event1 vector : event 1 indicator
#' @param event2 vector : event 2 indicator
#' @param Ind integer : number of subjects
#' @param Xs.0 same for left truncation
#' @param Us.0 same for left truncation
#' @param Xs.slope.0 same for left truncation
#' @param Us.slope.0 same for left truncation
#' @param P.0 same for left truncation
#' @param st.0 same for left truncation
#' @param Bs.0 same for left truncation
#' @param B.CR same for left truncation
#' @param Bs.CR same for left truncation
#' @param Bs.0.CR same for left truncation
#'
#' @return The value of the log-likelihood
#' @export
#'
#' @examples
log_llh <- function(param, nb.e.a, nb.priorMean.beta, nb.alpha, competing_risk,
                    nb.alpha.CR, variability_hetero, S,Zq, sharedtype, sharedtype_CR,
                    hazard_baseline, hazard_baseline_CR, ord.splines, Xtime, Utime, nb_pointsGK,
                    Xs,Us, Xslope, Uslope, Xs.slope, Us.slope, indices_beta_slope, Time,
                    st_calc, B, Bs, wk, Z, P, left_trunc, Z_CR, X_base, offset, U, y.new.prog, event1, event2, Ind,
                    Xs.0, Us.0, Xs.slope.0, Us.slope.0, P.0, st.0,Bs.0,B.CR, Bs.CR, Bs.0.CR
){
  #Manage parameters

  borne1 <- choose(n = nb.e.a, k = 2) + nb.e.a
  C <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
  C[lower.tri(C, diag=T)] <- param[1:borne1]
  borne2 <- borne1 + nb.priorMean.beta
  beta <- param[(borne1+1):borne2]
  borne3 <- borne2 + nb.alpha
  if(borne3 != borne2){
    alpha <- param[(borne2+1):borne3]
  }
  else{
    alpha <- 0
  }
  if(competing_risk){
    borne4 <- borne3 + nb.alpha.CR
    if(borne4 != borne3){
      alpha.CR <- param[(borne3+1):borne4]
    }
    else{
      alpha.CR <- 0
    }
  }
  else{
    borne4 <- borne3
  }
  if(variability_hetero){
    mu.log.sigma <- param[borne4+1]
    tau.log.sigma <- abs(param[borne4+2])
    #sigma2.log.sigma <- tau.log.sigma**2
    alpha.sigma <- param[borne4+3]
    borne5 <- borne4+3
    log.sigma_al <- matrix(rep(mu.log.sigma,S), nrow = S, byrow = T) + Zq[,nb.e.a+1]*tau.log.sigma
    sigma_al <- exp(log.sigma_al)
    if(competing_risk){
      alpha.sigma.CR <- param[borne4+4]
      borne5 <- borne4+4
    }
  }
  else{
    sigma.epsilon <- abs(param[borne4+1])
    borne5 <- borne4+1
  }
  curseur <- borne5
  if(sharedtype %in% c("RE")){
    stop("Not implemented yet")
  }
  if(competing_risk && sharedtype_CR %in% c("RE")){
    stop("Not implemented yet")
  }
  if(sharedtype %in% c("CV","CVS")){
    alpha.current <- param[curseur+1]
    curseur <- curseur+1
  }
  if(competing_risk && sharedtype_CR %in% c("CV","CVS")){
    alpha.current.CR <- param[curseur+1]
    curseur <- curseur +1
  }
  if(sharedtype %in% c("CVS","S")){
    alpha.slope <- param[curseur+1]
    curseur <- curseur +1
  }
  if(competing_risk && sharedtype_CR %in% c("CVS","S")){
    alpha.slope.CR <- param[curseur+1]
    curseur <- curseur +1
  }
  #if(hazard_baseline == "Exponential"){
  #  lambda_0 <- alpha[1]
  #}
  if(hazard_baseline == "Weibull"){
    #lambda_0 <- alpha[1]
    shape <- param[curseur+1]**2
    curseur <- curseur +1
  }
  if(hazard_baseline == "Splines"){
    gamma <- param[(curseur+1):(curseur+ord.splines+2)]
    curseur <- curseur + ord.splines + 2
  }
  # if(competing_risk && hazard_baseline_CR == "exponential"){
  #   lambda0.CR <- alpha.CR[1]
  #
  # }
  if(competing_risk && hazard_baseline_CR == "Weibull"){
    #lambda0.CR <- alpha.CR[1]
    shape.CR <- param[curseur+1]**2
    curseur <- curseur +1
  }
  if(competing_risk && hazard_baseline_CR == "Splines"){
    gamma.CR <- param[(curseur+1):(curseur+ord.splines+2)]
    curseur <- curseur + ord.splines + 2
  }

  #Manage random effects
  b_al <- Zq[,1:(nb.e.a)]%*%t(C)

  ll_glob <- 0

  for(i in 1:Ind){#Computation of contribution to the log_lokelihood
    h <- 1
    etaBaseline <- 0
    survLong <- 0
    etaBaseline.0 <- 0
    survLong.0 <- 0
    if(variability_hetero){
      h <- h*exp(alpha.sigma*sigma_al)
      etaBaseline <- etaBaseline + alpha.sigma*sigma_al
      if(left_trunc){
        etaBaseline.0 <- etaBaseline.0 + alpha.sigma*sigma_al
      }
    }
    if(competing_risk){
      h_CR <- 1
      etaBaseline_CR <- 0
      survLong_CR <- 0
      etaBaseline.0_CR <- 0
      survLong.0_CR <- 0
      if(variability_hetero){
        h_CR <- h_CR*exp(alpha.sigma.CR*sigma_al)
        etaBaseline_CR <- etaBaseline_CR + alpha.sigma.CR*sigma_al
        if(left_trunc){
          etaBaseline.0_CR <- etaBaseline.0_CR + alpha.sigma.CR*sigma_al
        }
      }
    }
    if(sharedtype %in% c("CV","CVS") || (competing_risk && sharedtype_CR %in% c("CV","CVS")) ){
      Xtime_i <- Xtime[i,]
      Utime_i <- Utime[i,]
      Xs_i <- Xs[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      Us_i <- Us[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      CV <- (beta%*%Xtime_i)[1,1]+b_al%*%Utime_i
      current.GK <- matrix(rep(beta%*%t(Xs_i),S),nrow=S,byrow = T) + b_al%*%t(Us_i)
      if(left_trunc){
        Xs.0_i <- Xs.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Us.0_i <- Us.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        current.GK.0 <- matrix(rep(beta%*%t(Xs.0_i),S),nrow=S,byrow = T) + b_al%*%t(Us.0_i)
      }
      if(sharedtype %in% c("CV","CVS")){
        h <- h*exp(alpha.current*CV)
        survLong <- survLong + alpha.current*current.GK
        if(left_trunc){
          survLong.0 <- survLong.0 + alpha.current*current.GK.0
        }
      }
      if(competing_risk && sharedtype_CR %in% c("CV","CVS")){
        h_CR <- h_CR*exp(alpha.current.CR*CV)
        survLong_CR <- survLong_CR + alpha.current.CR*current.GK
        if(left_trunc){
          survLong.0_CR <- survLong.0_CR + alpha.current.CR*current.GK.0
        }
      }
    }
    if(sharedtype %in% c("CVS","S") || (competing_risk && sharedtype_CR %in% c("CVS","S") )){
      #print("on passe ici")
      Xslope_i <- Xslope[i,]
      Uslope_i <- Uslope[i,]
      Xs.slope_i <- Xs.slope[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      Us.slope_i <- Us.slope[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      slope.GK <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope_i),S),nrow=S,byrow = T) + b_al[,-1]%*%t(Us.slope_i)
      if(length(indices_beta_slope) == 1){
        slope <- (beta[indices_beta_slope]%*%Xslope_i)[1,1]+b_al[,-1]*Uslope_i
      }
      else{
        slope <- (beta[indices_beta_slope]%*%Xslope_i)[1,1]+b_al[,-1]%*%Uslope_i
      }
      if(left_trunc){
        Xs.slope.0_i <- Xs.slope.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Us.slope.0_i <- Us.slope.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        slope.GK.0 <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope.0_i),S),nrow=S,byrow = T) + b_al[,-1]%*%t(Us.slope.0_i)
      }
      if(sharedtype %in% c("CVS","S")){
        h <- h*exp(alpha.slope*slope)
        survLong <- survLong + alpha.slope*slope.GK
        if(left_trunc){
          survLong.0 <- survLong.0 + alpha.slope*slope.GK.0
        }
      }
      if(competing_risk && sharedtype_CR %in% c("CVS","S") ){
        h_CR <- h_CR*exp(alpha.slope.CR*slope)
        survLong_CR <- survLong_CR + alpha.slope.CR*slope.GK
        if(left_trunc){
          survLong.0_CR <- survLong.0_CR + alpha.slope.CR*slope.GK.0
        }
      }
    }

    ###h0
    if(hazard_baseline == "Exponential"){
      h_0 <- 1
      h_0.GK <- wk
      if(left_trunc){
        h_0.GK.0 <- wk
      }
    }
    if(hazard_baseline == "Weibull"){
      Time_i <- Time[i]
      st_i <- st_calc[i,]
      h_0 <- shape*(Time_i**(shape-1))
      h_0.GK <- shape*(st_i**(shape-1))*wk
      if(left_trunc){
        st.0_i <- st.0[i,]
        h_0.GK.0 <-shape*(st.0_i**(shape-1))*wk
      }
    }
    if(hazard_baseline == "Splines"){
      B_i <- B[i,]
      Bs_i <- Bs[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      h_0 <- exp((gamma%*%B_i)[1,1])
      mat_h0s <- matrix(gamma,ncol=1)
      h_0.GK <- (wk*exp(Bs_i%*%mat_h0s))
      if(left_trunc){
        Bs.0_i <- Bs.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        h_0.GK.0 <- (wk*exp(Bs.0_i%*%mat_h0s))
      }
    }

    ###hazard function
    Z_i <- Z[i,]
    if(length(Z_i)==0){
      pred_surv <- 0
    }
    else{
      pred_surv <- (alpha%*%Z_i)[1,1]
    }
    h <- h_0*exp(pred_surv)*h
    etaBaseline <- etaBaseline + pred_surv
    if(left_trunc){
      etaBaseline.0 <- etaBaseline.0 + pred_surv
    }

    ###GK integration
    survLong <- exp(survLong)
    survLong <- survLong%*%h_0.GK

    P_i <- P[i]
    Surv <- (-exp(etaBaseline)*P_i*survLong)

    if(left_trunc){###Computation of S(T0i)
      #stop("Not implemented yet.")
      survLong.0 <- exp(survLong.0)
      survLong.0 <- survLong.0%*%h_0.GK.0
      P.0_i <- P.0[i]
      Surv.0 <- exp((-exp(etaBaseline.0)*P.0_i*survLong.0))
    }

    if(competing_risk){
      ###h0
      if(hazard_baseline_CR == "Exponential"){
        h_0.CR <- 1
        h_0.GK.CR <- wk
        if(left_trunc){
          h_0.GK.0_CR <- wk
        }
      }
      if(hazard_baseline_CR == "Weibull"){
        Time_i <- Time[i]
        st_i <- st_calc[i,]
        h_0.CR <- shape.CR*(Time_i**(shape.CR-1))
        h_0.GK.CR <- shape.CR*(st_i**(shape.CR-1))*wk
        if(left_trunc){
          st.0_i <- st.0[i,]
          h_0.GK.0_CR <- shape.CR*(st.0_i**(shape.CR-1))*wk
        }
      }
      if(hazard_baseline_CR == "Splines"){
        B.CR_i <- B.CR[i,]
        Bs.CR_i <- Bs.CR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        h_0.CR <- exp((gamma.CR%*%B.CR_i)[1,1])
        mat_h0s.CR <- matrix(gamma.CR,ncol=1)
        h_0.GK.CR <- (wk*exp(Bs.CR_i%*%mat_h0s.CR))
        if(left_trunc){
          Bs.0.CR_i <- Bs.0.CR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
          h_0.GK.0_CR <- (wk*exp(Bs.0.CR_i%*%mat_h0s.CR))
        }
      }

      ###hazard function
      Z.CR_i <- Z_CR[i,]
      if(length(Z.CR_i)==0){
        pred_surv.CR <- 0
      }
      else{
        pred_surv.CR <- (alpha.CR%*%Z.CR_i)[1,1]
      }
      h_CR <- h_0.CR*exp(pred_surv.CR)*h_CR
      etaBaseline_CR <- etaBaseline_CR + pred_surv.CR
      if(left_trunc){
        etaBaseline.0_CR <- etaBaseline.0_CR + pred_surv.CR
      }

      ###GK integration
      survLong_CR <- exp(survLong_CR)
      survLong_CR <- survLong_CR%*%h_0.GK.CR
      P_i <- P[i]
      Surv.CR <- (-exp(etaBaseline_CR)*P_i*survLong_CR)

      if(left_trunc){###Computation of S(T0i)
        #stop("Not implemented yet.")
        survLong.0_CR <- exp(survLong.0_CR)
        survLong.0_CR <- survLong.0_CR%*%h_0.GK.0_CR
        P.0_i <- P.0[i]
        Surv.0.CR <- exp((-exp(etaBaseline.0_CR)*P.0_i*survLong.0_CR))
      }
    }

    #Longitudinal part
    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    U_i <- U[offset[i]:(offset[i+1]-1),]
    y_i <- y.new.prog[offset[i]:(offset[i+1]-1)]
    if(variability_hetero  == T){
      sigma.long <- sigma_al
    }
    else{
      sigma.long <- sigma.epsilon
    }
    if(is.null(nrow(X_base_i))){
      CV <- (beta%*%X_base_i)[1,1] + b_al%*%U_i
      f_Y_b_sigma <- log((1/(sqrt(2*pi)*sigma.long))) + (-1/2)*((y_i-CV)/sigma.long)**2 #dnorm(x=y_i, mean = CV, sd = sigma.long)
    }
    else{
      f_Y_b_sigma <- rep(0,S)
      for(k in 1:nrow(X_base_i)){
        CV <- (beta%*%X_base_i[k,])[1,1] + b_al%*%U_i[k,]
        f_Y_b_sigma <- f_Y_b_sigma + (-1/2)*((y_i[k]-CV)/sigma.long)**2 #*dnorm(x = y_i[k], mean = CV, sd = sigma.long)
      }
      add <- nrow(X_base_i)*log((1/(sqrt(2*pi)*sigma.long)))
      f_Y_b_sigma <- f_Y_b_sigma+add
    }
    event1_i <- event1[i]
    if(competing_risk){
      event2_i <- event2[i]
      #log_dens <- log(sum((h**event1_i)*(h_CR**event2_i)*Surv*Surv.CR*f_Y_b_sigma)) - log(S)
      log_dens_int <- f_Y_b_sigma + log(h**event1_i)+log(h_CR**event2_i)+Surv+Surv.CR
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp


      log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)

    }
    else{
      log_dens_int <- f_Y_b_sigma + log(h**event1_i)+Surv
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp
      log_dens <- Clogexp  +log(sum(exp(log_dens_int))) - log(S)

    }

    if(left_trunc){
      #stop("Not yet implemented.")
      if(competing_risk){
        den <- log(sum(Surv.0*Surv.0.CR))-log(S)
      }
      else{
        den <- log(sum(Surv.0))-log(S)
        #print(den)
      }
      log_dens <- log_dens - den
      #print(log_dens)
    }
    ll_glob <- ll_glob + log_dens
    #print(log_dens)
  }

  #print(ll_glob)
  #if(is.na(ll_glob)){
  # ll_glob <- -10E8
  #}
  ll_glob
}
