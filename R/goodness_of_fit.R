#' Predictions for the goodness of fit, of the random effects, the current value for each individuals and the cumulative hazard function for both events
#'
#' @param object an object of class lsjm
#' @param graph a boolean to indicate to print graphics, False by default
#' @param break.times a vector of times for the time points of longitudinal graphic
#'
#' @return A list which contains the following elements :
#' \describe{
#' \item{\code{tables}}{A list with the table of the predicted random effect, the table of the predicted current value, table(s) of predictive cumulative hazard function(s)}
#' \item{\code{graphs}}{A list with 2 or 3 graphs : one for the longitudinal adjustment and one for each risk function}
#'
#' }
#' @export
#'
#' @examples
#' 
#' \donttest{
#' 
#'
#' #Fit a joint model with competing risks and subject-specific variability
#' example <- lsjm(formFixed = y~visit,
#' formRandom = ~ visit,
#' formGroup = ~ID,
#' formSurv = Surv(time, event ==1 ) ~ 1,
#' timeVar = "visit",
#' data.long = Data_toy,
#' variability_hetero = TRUE,
#' formFixedVar =~visit,
#' formRandomVar =~visit,
#' correlated_re = TRUE,
#' sharedtype = c("current value", "variability"),
#' hazard_baseline = "Weibull",
#' formSlopeFixed =~1,
#' formSlopeRandom = ~1,
#' indices_beta_slope = c(2), 
#' competing_risk = TRUE,
#' formSurv_CR = Surv(time, event ==2 ) ~ 1,
#' hazard_baseline_CR = "Weibull",
#' sharedtype_CR = c("current value", "variability"),
#' S1 = 100,
#' S2 = 1000,
#' nproc = 1,
#' maxiter = 100,
#' Comp.Rcpp = TRUE
#' )
#' 
#' #Assesment of the goodness of fit:
#' gof <- goodness_of_fit(example, graph = TRUE)
#' gof$tables
#' gof$graphs
#' }
#'
goodness_of_fit <- function(object, graph = FALSE, break.times = NULL){
  x <- object
  if(!inherits(x, "lsjm")) stop("use only \"lsjm\" objects")
  if(x$result$istop != 1) stop("The estimation didn't reach convergence \n")
  message("Computation of predictions")
  Xtime <- NULL
  Utime <- NULL
  Xs <- NULL
  Us <- NULL
  Xslope <- NULL
  Uslope <- NULL
  Xs.slope <- NULL
  Us.slope <- NULL
  wk <- NULL
  P <- NULL
  st_calc <- NULL
  B <- NULL
  Bs <- NULL
  #LT
  Xs.0 <- NULL
  Us.0 <- NULL
  Xs.slope.0 <- NULL
  Us.slope.0 <- NULL
  st.0 <- NULL
  Bs.0 <- NULL
  P.0 <- NULL
  #CR
  event2 <- NULL
  Z_CR <- NULL
  B.CR <- NULL
  Bs.CR <- NULL
  Bs.0.CR <- NULL
  st.0.CR <- NULL
  Bs.0.CR <- NULL
  gamma.CR <- NULL
  rr.CR <- NULL
  O_base <- NULL
  nb.e.a.sigma <- NULL
  nb.omega <- NULL
  Otime <- NULL
  Wtime <- NULL
  Os <- NULL
  Ws <- NULL
  W_base <- NULL
  shape <- NULL
  shape.CR <- NULL
  
  #data management
  data.long <- x$control$data.long
  id <- as.integer(data.long[all.vars(x$control$formGroup)][,1])
  if(!("id" %in% colnames(x$control$data.long))){ #To have a column named "id"
    data.long <- cbind(x$control$data.long, id = id)
  }
  idVar = "id"
  
  #Longitudinal part
  list.long <- data.manag.long(x$control$formGroup,x$control$formFixed, x$control$formRandom,x$control$data.long)
  X_base <- list.long$X
  U <- list.long$U
  nb.e.a <- ncol(U)
  y.new.prog <- list.long$y.new.prog
  list.var <- data.manag.sigma(x$control$formGroup,x$control$formFixedVar, x$control$formRandomVar,x$control$data.long)
  O_base <- list.var$X
  W_base <- list.var$U
  nb.omega <- ncol(O_base)
  #print(nb.omega)
  nb.e.a.sigma <- ncol(W_base)
  #data.long <- cbind(data.long,y.new.prog)
  data.long <- as.data.frame(data.long)
  offset <- list.long$offset
  Ind <- list.long$I
  
  # Survival part
  list.surv <- data.manag.surv(x$control$formGroup, x$control$formSurv, data.long, formSurv_CompRisk = x$control$formSurv_CR)
  event1 <- list.surv$event1
  event2 <- list.surv$event2
  Time <- list.surv$Time
  
  #Dependence
  data.id <- data.long[!duplicated(id),]
  data.id <- cbind(data.id,event1)
  if(c("random effects") %in% x$control$sharedtype){
    stop("Not implemented yet")
  }
  else{
    list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = x$control$nb_pointsGK)
    wk <- list.GaussKronrod$wk
    st_calc <- list.GaussKronrod$st
    P <- list.GaussKronrod$P
    id.GK <- list.GaussKronrod$id.GK
    if(x$control$left_trunc){
      list.GaussKronrod.0 <- data.GaussKronrod(data.id, x$control$Time.0, k = x$control$nb_pointsGK)
      st.0 <- list.GaussKronrod.0$st
      P.0 <- list.GaussKronrod.0$P
    }
    if(c("current value") %in% x$control$sharedtype){
      list.data.current.time <- data.time(data.id, list.surv$Time, x$control$formFixed, x$control$formRandom,x$control$timeVar)
      list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        x$control$formFixed, x$control$formRandom,x$control$timeVar)
      Xtime <- list.data.current.time$Xtime
      Utime <- list.data.current.time$Utime
      Xs <- list.data.GK.current$Xtime
      Us <- list.data.GK.current$Utime
      if(x$control$left_trunc){
        list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            x$control$formFixed, x$control$formRandom,x$control$timeVar)
        Xs.0 <- list.data.GK.current.0$Xtime
        Us.0 <- list.data.GK.current.0$Utime
      }
    }
    if(c("slope") %in% x$control$sharedtype){
      list.data.slope.time <- data.time(data.id, list.surv$Time, x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
      list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                      x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
      Xslope <- list.data.slope.time$Xtime
      Uslope <- list.data.slope.time$Utime
      Xs.slope <- list.data.GK.slope$Xtime
      Us.slope <- list.data.GK.slope$Utime
      if(x$control$left_trunc){
        list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                          x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
        Xs.slope.0 <- list.data.GK.slope.0$Xtime
        Us.slope.0 <- list.data.GK.slope.0$Utime
      }
    }
  }
  
  
  if(x$control$hazard_baseline == "Exponential"){
    Z <- list.surv$Z
  }
  else{
    if(x$control$hazard_baseline == "Weibull"){
      Z <- list.surv$Z
    }
    else{
      if(x$control$hazard_baseline == "Splines"){
        Z <- list.surv$Z[,-1]
        pp <- seq(0,1, length.out = x$control$ord.splines)
        pp <- utils::tail(utils::head(pp,-1),-1)
        tt1 <- as.data.frame(cbind(Time,event1))
        tt <- tt1$Time[which(tt1$event1 == 1)]
        kn <- quantile(tt, pp, names = FALSE)
        kn <- kn[kn<max(Time)]
        rr <- sort(c(rep(range(Time,0), 4L), kn))
        B <- splines::splineDesign(rr, Time, ord = 4L)
        Bs <- splines::splineDesign(rr, c(t(st_calc)), ord = 4L)
        if(x$control$left_trunc){
          Bs.0 <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
        }
      }
      else{
        stop("This type of base survival function is not implemented.")
      }
    }
    
  }
  nb.alpha.CR <- 0
  if(x$control$competing_risk){
    data.id <- cbind(data.id,event2)
    if(c("random effects") %in% x$control$sharedtype_CR){
      stop("Not implemented yet")
    }
    else{
      list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = x$control$nb_pointsGK)
      wk <- list.GaussKronrod$wk
      st_calc <- list.GaussKronrod$st
      P <- list.GaussKronrod$P
      id.GK <- list.GaussKronrod$id.GK
      if(x$control$left_trunc){
        list.GaussKronrod.0 <- data.GaussKronrod(data.id, x$control$Time.0, k = x$control$nb_pointsGK)
        st.0 <- list.GaussKronrod.0$st
        P.0 <- list.GaussKronrod.0$P
      }
      if(c("current value") %in% x$control$sharedtype_CR){
        list.data.current.time <- data.time(data.id, list.surv$Time, x$control$formFixed, x$control$formRandom,x$control$timeVar)
        list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                          x$control$formFixed, x$control$formRandom,x$control$timeVar)
        Xtime <- list.data.current.time$Xtime
        Utime <- list.data.current.time$Utime
        Xs <- list.data.GK.current$Xtime
        Us <- list.data.GK.current$Utime
        if(x$control$left_trunc){
          list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                              x$control$formFixed, x$control$formRandom,x$control$timeVar)
          Xs.0 <- list.data.GK.current.0$Xtime
          Us.0 <- list.data.GK.current.0$Utime
        }
      }
      if(c("slope") %in% x$control$sharedtype_CR){
        list.data.slope.time <- data.time(data.id, list.surv$Time, x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
        list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
        Xslope <- list.data.slope.time$Xtime
        Uslope <- list.data.slope.time$Utime
        Xs.slope <- list.data.GK.slope$Xtime
        Us.slope <- list.data.GK.slope$Utime
        if(x$control$left_trunc){
          list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
          Xs.slope.0 <- list.data.GK.slope.0$Xtime
          Us.slope.0 <- list.data.GK.slope.0$Utime
        }
        
      }
      
      if(x$control$hazard_baseline_CR == "Exponential"){
        Z_CR <- list.surv$Z_CR
      }
      else{
        if(x$control$hazard_baseline_CR == "Weibull"){
          Z_CR <- list.surv$Z_CR
        }
        else{
          if(x$control$hazard_baseline_CR == "Splines"){
            Z_CR <- list.surv$Z_CR[,-1]
            pp <- seq(0,1, length.out = x$control$ord.splines)
            pp <- utils::tail(utils::head(pp,-1),-1)
            tt2 <- as.data.frame(cbind(Time,event2))
            tt <- tt2$Time[which(tt2$event2 == 1)]
            kn <- quantile(tt, pp, names = FALSE)
            kn <- kn[kn<max(Time)]
            rr <- sort(c(rep(range(Time,0), 4L), kn))
            B.CR <- splines::splineDesign(rr, Time, ord = 4L)
            Bs.CR <- splines::splineDesign(rr, c(t(st_calc)), ord = 4L)
            if(x$control$left_trunc){
              Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
            }
          }
          else{
            stop("This type of base survival function is not implemented.")
          }
        }
        
      }
    }
  }
  if(c("variability") %in% x$control$sharedtype){
    list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = x$control$nb_pointsGK)
    wk <- list.GaussKronrod$wk
    st_calc <- list.GaussKronrod$st
    P <- list.GaussKronrod$P
    id.GK <- list.GaussKronrod$id.GK
    
    list.data.current.sigma.time <- data.time(data.id, list.surv$Time, x$control$formFixedVar, x$control$formRandomVar,x$control$timeVar)
    list.data.GK.current.sigma <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                            x$control$formFixedVar, x$control$formRandomVar,x$control$timeVar)
    Otime <- list.data.current.sigma.time$Xtime
    Wtime <- list.data.current.sigma.time$Utime
    Os <- list.data.GK.current.sigma$Xtime
    Ws <- list.data.GK.current.sigma$Utime
  }
  #Manage parameters
  estim_param <- x$table.res$Estimation
  curseur <- 1
  #Evenement 1 :
  ## Risque de base :
  if(x$control$hazard_baseline == "Weibull"){
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
    shape <- estim_param[curseur]**2
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline == "Splines"){
    gamma <- estim_param[(curseur):(curseur+x$control$ord.splines+1)]
    curseur <- curseur + x$control$ord.splines + 2
  }
  ## Covariables :
  if(x$control$nb.alpha >=1){
    alpha <- estim_param[(curseur):(curseur+x$control$nb.alpha-1)]
    curseur <- curseur+x$control$nb.alpha
  }
  ## Association :
  if(c("random effects") %in% x$control$sharedtype){
    stop("Not implemented yet")
  }
  if(c("current value") %in% x$control$sharedtype){
    alpha.current <- estim_param[curseur]
    curseur <- curseur + 1
  }
  if(c("slope") %in% x$control$sharedtype){
    alpha.slope <- estim_param[curseur]
    curseur <- curseur + 1
  }
  if(c("variability") %in% x$control$sharedtype){
    alpha.sigma <- estim_param[curseur]
    curseur <- curseur + 1
  }
  # Evenement 2
  if(x$control$competing_risk){
    ## Risque de base :
    if(x$control$hazard_baseline_CR == "Weibull"){
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
      shape.CR <- estim_param[curseur]**2
      curseur <- curseur + 1
    }
    if(x$control$hazard_baseline_CR == "Splines"){
      gamma.CR <- estim_param[(curseur):(curseur+x$control$ord.splines+1)]
      curseur <- curseur + x$control$ord.splines + 2
    }
    ## Covariables :
    if(x$control$nb.alpha.CR >=1){
      alpha.CR <- estim_param[(curseur):(curseur+x$control$nb.alpha.CR-1)]
      curseur <- curseur+x$control$nb.alpha.CR
    }
    ## Association :
    if(c("random effects") %in% x$control$sharedtype_CR){
      stop("Not implemented yet")
    }
    if(c("current value") %in% x$control$sharedtype_CR){
      alpha.current.CR <- estim_param[curseur]
      curseur <- curseur + 1
    }
    if(c("slope") %in% x$control$sharedtype_CR){
      alpha.slope.CR <- estim_param[curseur]
      curseur <- curseur + 1
    }
    if(c("variability") %in% x$control$sharedtype_CR){
      alpha.sigma.CR <- estim_param[curseur]
      curseur <- curseur + 1
    }
  }
  # Marqueur :
  ## Effets fixes trend :
  beta <- estim_param[curseur:(curseur+x$control$nb.priorMean.beta-1)]
  curseur <- curseur+x$control$nb.priorMean.beta
  ## Effets fixes var :
  if(x$control$variability_hetero){
    omega <- estim_param[curseur:(curseur+nb.omega-1)]
    curseur <- curseur + nb.omega
  }
  else{
    sigma.epsilon <- estim_param[curseur]
    curseur <- curseur + 1
  }
  ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
  if(x$control$variability_hetero){
    if(x$control$correlated_re){
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- estim_param[curseur:borne1]
      C2 <- matrix(estim_param[(borne1+1):(borne1+nb.e.a.sigma*x$control$nb.e.a)],nrow=nb.e.a.sigma,ncol=x$control$nb.e.a, byrow = TRUE)
      borne2 <- borne1+nb.e.a.sigma*x$control$nb.e.a + 1
      borne3 <- borne2 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma - 1
      C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- estim_param[borne2:borne3]
      C4 <- matrix(rep(0,x$control$nb.e.a*nb.e.a.sigma),nrow=x$control$nb.e.a,ncol=nb.e.a.sigma)
      MatCov <- rbind(cbind(C1,C4),cbind(C2,C3))
      MatCov <- as.matrix(MatCov)
      diag(MatCov) <- abs(diag(MatCov))
    }
    else{
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- estim_param[curseur:borne1]
      borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
      C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- estim_param[(borne1+1):borne3]
      C4 <- matrix(rep(0,x$control$nb.e.a*nb.e.a.sigma),nrow=x$control$nb.e.a,ncol=nb.e.a.sigma)
      C2.bis <- matrix(rep(0,x$control$nb.e.a*nb.e.a.sigma),nrow=nb.e.a.sigma,ncol=x$control$nb.e.a)
      MatCov <- rbind(cbind(C1,C4),cbind(C2.bis,C3))
      MatCov <- as.matrix(MatCov)
      diag(MatCov) <- abs(diag(MatCov))
    }
  }
  else{
    borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- estim_param[curseur:borne1]
    MatCov <- C1
  }
  Cum_risk2 <- c()
  Cum_risk1 <- c()
  Time.sort.unique <- unique(sort(Time))
  data.GaussKronrod.sort.unique <- data.GaussKronrod(data.id = data.id, Time = Time.sort.unique, k = x$control$nb_pointsGK)
  st_calc.sort.unique <- data.GaussKronrod.sort.unique$st
  P.sort.unique <- data.GaussKronrod.sort.unique$P
  pred.r.e.table <- c()
  pred.CV <- c()
  pb <- utils::txtProgressBar(min = 0,
                              max = Ind,
                              initial = 0,
                              char = "*",
                              style = 3)
  for(i in 1:Ind){
    Cum_risk_2i <- c()
    Cum_risk_1i <- c()
    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    U_i <- U[offset[i]:(offset[i+1]-1),]
    y_i <- y.new.prog[offset[i]:(offset[i+1]-1)]
    Xtime_i <- Xtime[i,]
    Utime_i <- Utime[i,]
    Xs_i <- Xs[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Us_i <- Us[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Xslope_i <- Xslope[i,]
    Uslope_i <- Uslope[i,]
    Xs.slope_i <- Xs.slope[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Us.slope_i <- Us.slope[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Time_i <- Time[i]
    st_i <- st_calc[i,]
    B_i <- B[i,]
    Bs_i <- Bs[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Z_i <- Z[i,]
    P_i <- P[i]
    event1_i <- event1[i]
    event2_i <- event2[i]
    B.CR_i <- B.CR[i,]
    Bs.CR_i <- Bs.CR[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Z.CR_i <- Z_CR[i,]
    if(is.null(dim(U_i))){
      U_i <- matrix(U_i, nrow= 1)
      X_base_i <- matrix(X_base_i, nrow= 1)
    }
    O_base_i <- O_base[offset[i]:(offset[i+1]-1),]
    W_base_i <- W_base[offset[i]:(offset[i+1]-1),]
    if(is.null(dim(W_base_i))){
      W_base_i <- matrix(W_base_i, nrow= 1)
      O_base_i <- matrix(O_base_i, nrow= 1)
    }
    Otime_i <- Otime[i,]
    Wtime_i <- Wtime[i,]
    Os_i <- Os[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    Ws_i <- Ws[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),]
    
    Sigma.re <- MatCov%*%t(MatCov)
    if(x$control$variability_hetero){
      binit <- mvtnorm::rmvnorm(1, mean = rep(0, x$control$nb.e.a+nb.e.a.sigma), Sigma.re)
    }
    else{
      binit <- mvtnorm::rmvnorm(1, mean = rep(0, x$control$nb.e.a), Sigma.re)
    }
    
    #Longitudinal prediction
    
    pred.r.e <- marqLevAlg(binit, fn = pred.re, minimize = FALSE,
                           nb.e.a = x$control$nb.e.a, variability_hetero = x$control$variability_hetero, nb.e.a.sigma = nb.e.a.sigma,
                           Sigma.re = Sigma.re, X_base_i = X_base_i, U_i = U_i, beta = beta, omega = omega, O_base_i = O_base_i,
                           W_base_i = W_base_i, y_i = y_i, sigma.epsilon = sigma.epsilon, Otime_i = Otime_i, Wtime_i = Wtime_i,
                           Os_i = Os_i, Ws_i = Ws_i, S = x$control$S2, alpha.sigma = alpha.sigma, competing_risk = x$control$competing_risk, alpha.sigma.CR = alpha.sigma.CR,
                           sharedtype = x$control$sharedtype, sharedtype_CR = x$control$sharedtype_CR, alpha.current = alpha.current, alpha.current.CR = alpha.current.CR,
                           alpha.slope = alpha.slope, alpha.slope.CR = alpha.slope.CR, Xtime_i = Xtime_i, Utime_i = Utime_i, Xs_i = Xs_i, Us_i = Us_i,
                           indices_beta_slope = x$control$indices_beta_slope, hazard_baseline = x$control$hazard_baseline, wk = wk, st_i = st_i, gamma = gamma, B_i = B_i, Bs_i = Bs_i,
                           Z_i = Z_i, alpha = alpha, shape = shape, Time_i = Time_i, P_i = P_i, hazard_baseline_CR = x$control$hazard_baseline_CR, gamma.CR = gamma.CR, B.CR_i = B.CR_i,
                           Bs.CR_i = Bs.CR_i, Z.CR_i = Z.CR_i, alpha.CR = alpha.CR, shape.CR = shape.CR, event1_i = event1_i, event2_i = event2_i,Xs.slope_i = Xs.slope_i, Us.slope_i = Us.slope_i,
                           Xslope_i = Xslope_i, Uslope_i = Uslope_i, nproc = 1,
                           clustertype = x$control$clustertype,
                           maxiter = x$control$maxiter, print.info = F,blinding = TRUE, epsa = x$control$epsa,
                           epsb = x$control$epsb, epsd = x$control$epsd)
    
    
    pred.r.e.table <- rbind(pred.r.e.table,c(data.id$ID[i], pred.r.e$b))
    CV <- X_base_i%*%beta + U_i%*%pred.r.e$b[1:(x$control$nb.e.a)]
    if(x$control$variability_hetero){
      Var <- O_base_i%*%omega + W_base_i%*%pred.r.e$b[(x$control$nb.e.a+1):(x$control$nb.e.a+nb.e.a.sigma)]
    }
    if(x$control$variability_hetero){
      pred.CV <- rbind(pred.CV,cbind(rep(data.id$ID[i],length(CV)), CV,Var,X_base_i[,2]))
    }
    else{
      pred.CV <- rbind(pred.CV,cbind(rep(data.id$ID[i],length(CV)), CV,X_base_i[,2]))
    }
    
    ### Survival goodness-of-fit
    #print("Survival part")
    
    for(j in 1:nrow(st_calc.sort.unique)){
      pred_haz <- 0
      if((c("variability") %in% x$control$sharedtype)|| (x$control$competing_risk && c("variability") %in% x$control$sharedtype_CR) ){
        list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),], st_calc.sort.unique[j,],
                                                            x$control$formFixedVar, x$control$formRandomVar,x$control$timeVar)
        Os.j <- list.data.GK.current.sigma.sort.unique$Xtime
        Ws.j <- list.data.GK.current.sigma.sort.unique$Utime
        Sigma.current.GK <- exp(omega%*%t(Os_i) + pred.r.e$b[(x$control$nb.e.a+1):(x$control$nb.e.a+nb.e.a.sigma)]%*%t(Ws_i))
        if(c("variability") %in% x$control$sharedtype){
          pred_haz <- pred_haz +  alpha.sigma*Sigma.current.GK
        }
      }
      
      if((c("current value") %in% x$control$sharedtype) || (x$control$competing_risk && c("current value") %in% x$control$sharedtype_CR) ){
        
        list.data.GK.current.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),], st_calc.sort.unique[j,],
                                                      x$control$formFixed, x$control$formRandom,x$control$timeVar)
        
        Xs.j <- list.data.GK.current.sort.unique$Xtime
        Us.j <- list.data.GK.current.sort.unique$Utime
        
        current.GK <-beta%*%t(Xs.j) + pred.r.e$b[1:(x$control$nb.e.a)]%*%t(Us.j)
        
        if(c("current value") %in% x$control$sharedtype){
          pred_haz <- pred_haz +  alpha.current*current.GK
        }
      }
      
      if((c("slope") %in% x$control$sharedtype)|| (x$control$competing_risk && c("slope") %in% x$control$sharedtype_CR)){
        list.data.GK.slope.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(i-1)+1):(x$control$nb_pointsGK*i),], st_calc.sort.unique[j,],
                                                    x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$timeVar)
        
        Xs.slope.j <- list.data.GK.slope.sort.unique$Xtime
        Us.slope.j <- list.data.GK.slope.sort.unique$Utime
        
        slope.GK <- beta[x$control$indices_beta_slope]%*%t(Xs.slope.j) +  pred.r.e$b[2:(x$control$nb.e.a)]%*%t(Us.slope.j)
        
        if(c("slope") %in% x$control$sharedtype){
          pred_haz <- pred_haz +  alpha.slope*slope.GK
        }
      }
      
      if(x$control$hazard_baseline == "Exponential"){
        h_0 <- 1
        h_0.GK <- wk
      }
      if(x$control$hazard_baseline == "Weibull"){
        st_j <- st_calc.sort.unique[j,]
        h_0.GK <- shape*(st_j**(shape-1))*wk
      }
      if(x$control$hazard_baseline == "Splines"){
        st_j <- st_calc.sort.unique[j,]
        Bs_j <- splines::splineDesign(x$control$knots.hazard_baseline.splines, st_j, ord = 4L)
        #Bs_j <- Bs[(x$control$nb_pointsGK*(j-1)+1):(x$control$nb_pointsGK*j),]
        mat_h0s <- matrix(gamma,ncol=1)
        h_0.GK <- (wk*exp(Bs_j%*%mat_h0s))
      }
      
      ###hazard function
      Z_i <- Z[i,]
      if(length(Z_i)==0){
        pred_surv <- 0
      }
      else{
        pred_surv <- (alpha%*%Z_i)[1,1]
      }
      
      pred_haz <- pred_haz + pred_surv
      
      Cum_risk_1i <- c(Cum_risk_1i, P.sort.unique[j]*sum(exp(pred_haz)%*%h_0.GK))
      
      if(x$control$competing_risk){
        if(c("variability") %in% x$control$sharedtype_CR){
          pred_haz.CR <- alpha.sigma.CR*Sigma.current.GK
        }
        else{
          pred_haz.CR <- 0
        }
        
        if(c("current value") %in% x$control$sharedtype_CR){
          pred_haz.CR <- pred_haz.CR + alpha.current.CR*current.GK
        }
        
        if(c("slope") %in% x$control$sharedtype_CR){
          pred_haz.CR <- pred_haz.CR + alpha.slope.CR*slope.GK
        }
        
        
        if(x$control$hazard_baseline_CR == "Exponential"){
          h_0.GK.CR <- wk
        }
        if(x$control$hazard_baseline_CR == "Weibull"){
          st_j <- st_calc.sort.unique[j,]
          h_0.GK.CR <- shape.CR*(st_j**(shape.CR-1))*wk
        }
        if(x$control$hazard_baseline_CR == "Splines"){
          st_j <- st_calc.sort.unique[j,]
          Bs_j <- splines::splineDesign(x$control$knots.hazard_baseline.splines.CR, st_j, ord = 4L)
          #Bs_j <- Bs.CR[(x$control$nb_pointsGK*(j-1)+1):(x$control$nb_pointsGK*j),]
          mat_h0s <- matrix(gamma.CR,ncol=1)
          h_0.GK.CR <- (wk*exp(Bs_j%*%mat_h0s))
        }
        
        ###hazard function
        Z.CR_i <- Z_CR[i,]
        if(length(Z.CR_i)==0){
          pred_surv <- 0
        }
        else{
          pred_surv.CR <- (alpha.CR%*%Z.CR_i)[1,1]
        }
        
        pred_haz.CR <- pred_haz.CR + pred_surv.CR
        
        Cum_risk_2i <- c(Cum_risk_2i, P.sort.unique[j]*sum(exp(pred_haz.CR)%*%h_0.GK.CR))
        
      }
    }
    
    Cum_risk1 <- rbind(Cum_risk1,Cum_risk_1i)
    Cum_risk2 <- rbind(Cum_risk2,Cum_risk_2i)
    
    utils::setTxtProgressBar(pb,i)
  }
  result <- list(pred.r.e.table = pred.r.e.table,
                 pred.CV = pred.CV, Cum_risk1 = Cum_risk1, Cum_risk2 = Cum_risk2)
  graphs <- NULL
  if(graph){
    timeInterv <- range(data.long[,x$control$timeVar])
    if(is.null(break.times)) break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
    graphs <- plot_goodnessoffit(data.long,data.id,pred.CV,break.times, formFixed = x$control$formFixed, formSurv = x$control$formSurv,
                                 timeVar = x$control$timeVar,Cum_risk1, competing_risk = x$control$competing_risk, formSurv_CR = x$control$formSurv_CR, Cum_risk2)
  }
  
  result.final <- list(tables = result,
                       graphs = graphs)
  result.final
}


