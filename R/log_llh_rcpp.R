#' Log-likelihood computation in RCPP
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
#' @param Os matrix : fixed effects of variability at Gauss-Kronrod times
#' @param Ws matrix : random effects of variability at Gauss-Kronrod times
#' @param Otime matrix : fixed effects of variability at event time
#' @param Wtime matrix : RE of variability at event time
#' @param O_base matrix : fixed effects for variability
#' @param W_base matrix : fixed effects for variability
#' @param Os.0 matrix : same for left truncation
#' @param Ws.0 matrix : same for left truncation
#' @param nb.e.a.sigma integer : number of RE for variability
#' @param nb.omega integer : number of fixed effects for variability
#' @param correlated_re boolean : indicator to allow all the random effects to be correlated
#'
#' @import Rcpp
#' @return The value of the log-likelihood
#'


log_llh_rcpp <- function(param, nb.e.a, nb.priorMean.beta, nb.alpha, competing_risk,
                         nb.alpha.CR, variability_hetero, S,Zq, sharedtype, sharedtype_CR,
                         hazard_baseline, hazard_baseline_CR, ord.splines, Xtime, Utime, nb_pointsGK,
                         Xs,Us, Xslope, Uslope, Xs.slope, Us.slope, indices_beta_slope, Time,
                         st_calc, B, Bs, wk, Z, P, left_trunc, Z_CR, X_base, offset, U, y.new.prog, event1, event2, Ind,
                         Xs.0, Us.0, Xs.slope.0, Us.slope.0, P.0, st.0,Bs.0,B.CR, Bs.CR, Bs.0.CR,
                         nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega, Otime = Otime, Wtime = Wtime,
                         Os = Os, Ws = Ws, O_base = O_base, W_base=W_base, correlated_re = correlated_re, Os.0 = Os.0, Ws.0 = Ws.0
){

  #initialisation des paramètres
  #browser()
  Otime_i <- c(1); Wtime_i <- c(1); Os_i <- as.matrix(1); Ws_i <- as.matrix(1); omega <- c(1)
  b_om <- as.matrix(1); alpha.sigma <- 0; Os.0_i <- as.matrix(1); Ws.0_i <- as.matrix(1); Xtime_i <- c(1); Utime_i <- c(1);
  alpha.sigma.CR <- 0;  Xs_i <- as.matrix(1); Us_i <- as.matrix(1); Xs.0_i <- as.matrix(1); Us.0_i<-as.matrix(1)
  alpha.current <- 0; alpha.current.CR <- 0; alpha.sigma.CR <- 0; beta_slope <- c(1); Xslope_i <- c(1); b_al_slope <- as.matrix(1)
  Uslope_i <- c(1); Xs.slope_i <- as.matrix(1); Us.slope_i <- as.matrix(1); Xs.slope.0_i <- as.matrix(1)
  Us.slope.0_i <- as.matrix(1); alpha.slope <- 0; alpha.slope.CR <- 0; Time_i <- 0; st_i <- c(1); st.0_i <- c(1)
  shape <- 0; gamma <- c(1); B_i <- c(1); Bs_i <- as.matrix(1); Bs.0_i <- as.matrix(1); Z_i <- c(1); alpha <- c(1)
  P.0_i <- 0; shape.CR <- 0; gamma.CR <- c(1); B.CR_i <- c(1); Bs.CR_i <- as.matrix(1); Bs.0.CR_i <- as.matrix(1)
  Z.CR_i <- c(1); alpha.CR <- c(1); O_base_i <- as.matrix(1); W_base_i <- as.matrix(1); sigma.epsilon <- 0
  
  #Manage parameters
  curseur <- 1
  #Evenement 1 :
  ## Risque de base :
  if(hazard_baseline == "Weibull"){
    #if(scaleWeibull == "square"){
    #  alpha_weib <- param[curseur]**2
    #  curseur <- curseur + 1
    #  shape <- param[curseur]**2
    #  curseur <- curseur + 1
    #}
    #else{
    #  alpha_weib <- exp(param[curseur])
    #  curseur <- curseur + 1
    #  shape <- exp(param[curseur])
    #  curseur <- curseur + 1
    #}
    shape <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline == "Splines"){
    gamma <- param[(curseur):(curseur+ord.splines+1)]
    curseur <- curseur + ord.splines + 2
  }
  ## Covariables :
  if(nb.alpha >=1){
    alpha <- param[(curseur):(curseur+nb.alpha-1)]
    curseur <- curseur+nb.alpha
  }
  ## Association :
  if("current value" %in% sharedtype){
    alpha.current <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    alpha.current <- 0
  }
  if("slope" %in% sharedtype){
    alpha.slope <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    alpha.slope <- 0
  }
  if("variability" %in% sharedtype){
    alpha.sigma <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    alpha.sigma <- 0
  }
  #if(sharedtype %in% c("RE")){
  #  stop("Not implemented yet")
  #}
  #if(sharedtype %in% c("CV","CVS")){
  #  alpha.current <- param[curseur]
  #  curseur <- curseur + 1
  #}
  #if(sharedtype %in%  c("CVS","S")){
  #  alpha.slope <- param[curseur]
  #  curseur <- curseur + 1
  #}
  #if(variability_hetero){
  #  alpha.sigma <- param[curseur]
  #  curseur <- curseur + 1
  #}
  # Evenement 2
  if(competing_risk){
    ## Risque de base :
    if(hazard_baseline_CR == "Weibull"){
      #if(scaleWeibull == "square"){
      #  alpha_weib.CR <- param[curseur]**2
      #  curseur <- curseur + 1
      #  shape.CR <- param[curseur]**2
      #  curseur <- curseur + 1
      #}
      #else{
      #  alpha_weib.CR <- exp(param[curseur])
      #  curseur <- curseur + 1
      #  shape.CR <- exp(param[curseur])
      #  curseur <- curseur + 1
      #}
      shape.CR <- param[curseur]**2
      curseur <- curseur + 1
    }
    if(hazard_baseline_CR == "Splines"){
      gamma.CR <- param[(curseur):(curseur+ord.splines+1)]
      curseur <- curseur + ord.splines + 2
    }
    ## Covariables :
    if(nb.alpha.CR >=1){
      alpha.CR <- param[(curseur):(curseur+nb.alpha.CR-1)]
      curseur <- curseur+nb.alpha.CR
    }
    ## Association :
    if("current value" %in% sharedtype_CR){
      alpha.current.CR <- param[curseur]
      curseur <- curseur + 1
    }
    else{
      alpha.current.CR <- 0
    }
    if("slope" %in% sharedtype_CR){
      alpha.slope.CR <- param[curseur]
      curseur <- curseur + 1
    }
    else{
      alpha.slope <- 0
    }
    if("variability" %in% sharedtype_CR){
      alpha.sigma.CR <- param[curseur]
      curseur <- curseur + 1
    }
    else{
      alpha.sigma.CR <- 0
    }
    #if(sharedtype_CR %in% c("RE")){
    #  stop("Not implemented yet")
    #}
    #if(sharedtype_CR %in% c("CV","CVS")){
    #  alpha.current.CR <- param[curseur]
    #  curseur <- curseur + 1
    #}
    #if(sharedtype_CR %in%  c("CVS","S")){
    #  alpha.slope.CR <- param[curseur]
    #  curseur <- curseur + 1
    #}
    #if(variability_hetero){
    #  alpha.sigma.CR <- param[curseur]
    #  curseur <- curseur + 1
    #}
  }
  # Marqueur :
  ## Effets fixes trend :
  beta <- param[curseur:(curseur+nb.priorMean.beta-1)]
  if( "slope" %in% sharedtype || "slope" %in% sharedtype_CR){
    beta_slope <- beta[indices_beta_slope]
  }
  curseur <- curseur+nb.priorMean.beta
  ## Effets fixes var :
  if(variability_hetero){
    omega <- param[curseur:(curseur+nb.omega-1)]
    curseur <- curseur + nb.omega
  }
  else{
    sigma.epsilon <- abs(param[curseur])
    curseur <- curseur + 1
  }
  ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
  if(variability_hetero){
    if(correlated_re){
      C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      MatCov <- C1
      MatCov <- as.matrix(MatCov)
      random.effects <- Zq%*%t(MatCov)
      b_al <- random.effects[,1:nb.e.a]
      b_al <- matrix(b_al, ncol = nb.e.a)
      b_om <- random.effects[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]
      b_om <- matrix(b_om, ncol = nb.e.a.sigma)
    }
    else{
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
      C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
      MatCovb <- as.matrix(C1)
      MatCovSig <- as.matrix(C3)
      b_al <- Zq[,1:nb.e.a]%*%t(MatCovb)
      b_al <- matrix(b_al, ncol = nb.e.a)
      b_om <- Zq[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]%*%t(MatCovSig)
      b_om <- matrix(b_om, ncol = nb.e.a.sigma)
    }
  }
  else{
    borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    MatCov <- as.matrix(C1)
    b_al <- Zq%*%t(MatCov)
    b_al <- matrix(b_al, ncol = nb.e.a)
  }
  #Manage random effects
  
  if("slope" %in% sharedtype || (competing_risk && "slope" %in% sharedtype_CR)){
    b_al_slope <- as.matrix(b_al[,-1])
  }

  ll_glob <- 0
  sht <- list(sharedtype, sharedtype_CR)
  HB <- list(hazard_baseline, hazard_baseline_CR)
  list_nb_points_int <- list(S , nb_pointsGK)
  sharedtype_bool <- c("current value" %in% sharedtype, "slope" %in% sharedtype, "variability" %in% sharedtype)
  sharedtype_CR_bool <- c("current value" %in% sharedtype_CR, "slope" %in% sharedtype_CR, "variability" %in% sharedtype_CR)
  for(i in 1:Ind){
    if(variability_hetero){
      Otime_i <- as.matrix(Otime[i,])
      Wtime_i <- as.matrix(Wtime[i,])
      Os_i <- as.matrix(Os[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      Ws_i <- as.matrix(Ws[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      O_base_i <- O_base[offset[i]:(offset[i+1]-1),]
      O_base_i <- matrix(O_base_i, nrow = offset[i+1]-offset[i])
      W_base_i <- W_base[offset[i]:(offset[i+1]-1),]
      W_base_i <- matrix(W_base_i, nrow = offset[i+1]-offset[i])
      if(left_trunc){
        Os.0_i <- as.matrix(Os.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
        Ws.0_i <- as.matrix(Ws.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      }
    }
    if("current value" %in% sharedtype ||(competing_risk && "current value" %in% sharedtype_CR)){
      Xtime_i <- as.matrix(Xtime[i,])
      Utime_i <- as.matrix(Utime[i,])
      Xs_i <- as.matrix(Xs[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      Us_i <- as.matrix(Us[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      if(left_trunc){
        Xs.0_i <- as.matrix(Xs.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
        Us.0_i <- as.matrix(Us.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      }
    }
    if("slope" %in% sharedtype ||(competing_risk && "slope" %in% sharedtype_CR)){
      Xslope_i <- as.matrix(Xslope[i,])
      Uslope_i <- as.matrix(Uslope[i,])
      Xs.slope_i <- as.matrix(Xs.slope[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      Us.slope_i <- as.matrix(Us.slope[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      if(left_trunc){
        Xs.slope.0_i <- as.matrix(Xs.slope.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
        Us.slope.0_i <- as.matrix(Us.slope.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),])
      }
    }
    if(hazard_baseline == "Weibull"){
      Time_i <- Time[i]
      st_i <- st_calc[i,]
      if(left_trunc){
        st.0_i <- st.0[i,]
      }
    }
    if(hazard_baseline == "Splines"){
      B_i <- B[i,]
      Bs_i <- Bs[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      if(left_trunc){
        Bs.0_i <- Bs.0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      }
    }
    ###hazard function
    Z_i <- Z[i,]
    P_i <- P[i]
    if(left_trunc){
      P.0_i <- P.0[i]
    }
    
    if(competing_risk){
      if(hazard_baseline_CR == "Weibull"){
        Time_i <- Time[i]
        st_i <- st_calc[i,]
        if(left_trunc){
          st.0_i <- st.0[i,]
        }
      }
      if(hazard_baseline_CR == "Splines"){
        B.CR_i <- B.CR[i,]
        Bs.CR_i <- Bs.CR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Bs.0.CR_i <- Bs.0.CR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      ###hazard function
      Z.CR_i <- Z_CR[i,]
    }
    
    #Longitudinal part
    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
    n_row_X <- nrow(X_base_i)
    U_i <- U[offset[i]:(offset[i+1]-1),]
    U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
    y_i <- y.new.prog[offset[i]:(offset[i+1]-1)]
    event1_i <- event1[i]
    event2_i <- 0
    if(competing_risk){
      event2_i <- event2[i]
    }
    list.event <- list(event1_i,event2_i)
    log_ind <- log_llh_ind(variability_hetero, Otime_i, Wtime_i, Os_i, Ws_i, omega, b_om, list_nb_points_int,
                           alpha.sigma, left_trunc, Os.0_i, Ws.0_i, competing_risk, alpha.sigma.CR, beta,
                           Xtime_i, Utime_i, b_al,Xs_i, Us_i, Xs.0_i, Us.0_i, alpha.current, sharedtype_bool,sharedtype_CR_bool,
                           alpha.current.CR, beta_slope, Xslope_i, b_al_slope, Uslope_i, Xs.slope_i, 
                           Us.slope_i,Xs.slope.0_i, Us.slope.0_i, alpha.slope, alpha.slope.CR, HB, 
                           wk, Time_i, st_i, st.0_i, shape, gamma, B_i, Bs_i, Bs.0_i, Z_i, alpha, P_i, P.0_i, 
                           shape.CR, gamma.CR, B.CR_i, Bs.CR_i, Bs.0.CR_i,
                           Z.CR_i, alpha.CR, X_base_i, O_base_i, W_base_i, U_i, y_i, n_row_X, sigma.epsilon, list.event)
    
    ll_glob <- ll_glob + log_ind
  }
  if(is.na(ll_glob)){
    ll_glob <- -10E8
  }
  
  print(ll_glob)
  ll_glob
}