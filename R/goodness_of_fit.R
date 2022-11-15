#' Predictions for the goodness of fit, of the random effects, the current value for each individuals and the cumulative hazard function for both events
#'
#' @param FlexVarJM an object of class FlexVarJM
#'
#' @return
#' @export
#'
#' @examples
#'
goodness_of_fit <- function(FlexVarJM, graph = F, break.times =NULL ){
  cat("Computation of prediction \n")
  time.prog1 <- Sys.time()

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
  #data management

  data.long <- FlexVarJM$control$data.long
  id <- as.integer(data.long[all.vars(FlexVarJM$control$formGroup)][,1])
  #print(id)
  if(!("id" %in% colnames(FlexVarJM$control$data.long))) #To have a column named "id"
    data.long <- cbind(FlexVarJM$control$data.long, id = id)
  idVar = "id"

  ##longitudinal part
  list.long <- data.manag.long(FlexVarJM$control$formGroup,FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,data.long)
  X_base <- list.long$X
  U <- list.long$U
  y.new.prog <- list.long$y.new.prog
  data.long <- cbind(data.long,y.new.prog)
  data.long <- as.data.frame(data.long)
  offset <- list.long$offset
  Ind <- list.long$I

  ## survival part
  list.surv <- data.manag.surv(FlexVarJM$control$formGroup, FlexVarJM$control$formSurv, data.long, formSurv_CompRisk = FlexVarJM$control$formSurv_CR)
  event1 <- list.surv$event1
  event2 <- list.surv$event2
  Time <- list.surv$Time

  ##dependance
  data.id <- data.long[!duplicated(id),]
  data.id <- cbind(data.id,event1)
  lag = 0
  if(FlexVarJM$control$sharedtype %in% c("random effects")){
    stop("Not implemented yet")
  }
  else{
    list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = FlexVarJM$control$nb_pointsGK)
    wk <- list.GaussKronrod$wk
    st_calc <- list.GaussKronrod$st
    P <- list.GaussKronrod$P
    id.GK <- list.GaussKronrod$id.GK
    if(FlexVarJM$control$left_trunc){
      list.GaussKronrod.0 <- data.GaussKronrod(data.id, FlexVarJM$control$Time.0, k = FlexVarJM$control$nb_pointsGK)
      st.0 <- list.GaussKronrod.0$st
      P.0 <- list.GaussKronrod.0$P
    }
    if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
      list.data.current.time <- data.time(data.id, list.surv$Time, FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      Xtime <- list.data.current.time$Xtime
      Utime <- list.data.current.time$Utime
      Xs <- list.data.GK.current$Xtime
      Us <- list.data.GK.current$Utime
      if(FlexVarJM$control$left_trunc){
        list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xs.0 <- list.data.GK.current.0$Xtime
        Us.0 <- list.data.GK.current.0$Utime
      }
    }
    if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
      list.data.slope.time <- data.time(data.id, list.surv$Time, FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                      FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      Xslope <- list.data.slope.time$Xtime
      Uslope <- list.data.slope.time$Utime
      Xs.slope <- list.data.GK.slope$Xtime
      Us.slope <- list.data.GK.slope$Utime
      if(FlexVarJM$control$left_trunc){
        list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                          FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        Xs.slope.0 <- list.data.GK.slope.0$Xtime
        Us.slope.0 <- list.data.GK.slope.0$Utime
      }
    }
  }


  if(FlexVarJM$control$hazard_baseline == "Exponential"){
    Z <- list.surv$Z
  }
  else{
    if(FlexVarJM$control$hazard_baseline == "Weibull"){
      Z <- list.surv$Z
    }
    else{
      if(FlexVarJM$control$hazard_baseline == "Splines"){
        Z <- list.surv$Z[,-1]
        pp <- seq(0,1, length.out = FlexVarJM$control$ord.splines)
        pp <- tail(head(pp,-1),-1)
        tt1 <- as.data.frame(cbind(Time,event1))
        tt <- tt1$Time[which(tt1$event1 == 1)]
        kn <- quantile(tt, pp, names = FALSE)
        kn <- kn[kn<max(Time)]
        rr <- sort(c(rep(range(Time,0), 4L), kn))
        B <- splineDesign(rr, Time, ord = 4L)
        Bs <- splineDesign(rr, c(t(st_calc)), ord = 4L)
        if(FlexVarJM$control$left_trunc){
          Bs.0 <- splineDesign(rr, c(t(st.0)), ord = 4L)
        }
      }
      else{
        stop("This type of base survival function is not implemented.")
      }
    }

  }
  nb.alpha.CR <- 0
  if(FlexVarJM$control$competing_risk){
    data.id <- cbind(data.id,event2)
    if(FlexVarJM$control$sharedtype_CR %in% c("random effects")){
      stop("Not implemented yet")
    }
    else{
      list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = FlexVarJM$control$nb_pointsGK)
      wk <- list.GaussKronrod$wk
      st_calc <- list.GaussKronrod$st
      P <- list.GaussKronrod$P
      id.GK <- list.GaussKronrod$id.GK
      if(FlexVarJM$control$left_trunc){
        list.GaussKronrod.0 <- data.GaussKronrod(data.id, FlexVarJM$control$Time.0, k = FlexVarJM$control$nb_pointsGK)
        st.0 <- list.GaussKronrod.0$st
        P.0 <- list.GaussKronrod.0$P
      }
      if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
        list.data.current.time <- data.time(data.id, list.surv$Time, FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                          FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xtime <- list.data.current.time$Xtime
        Utime <- list.data.current.time$Utime
        Xs <- list.data.GK.current$Xtime
        Us <- list.data.GK.current$Utime
        if(FlexVarJM$control$left_trunc){
          list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                              FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
          Xs.0 <- list.data.GK.current.0$Xtime
          Us.0 <- list.data.GK.current.0$Utime
        }
      }
      if(FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
        list.data.slope.time <- data.time(data.id, list.surv$Time, FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        Xslope <- list.data.slope.time$Xtime
        Uslope <- list.data.slope.time$Utime
        Xs.slope <- list.data.GK.slope$Xtime
        Us.slope <- list.data.GK.slope$Utime
        if(FlexVarJM$control$left_trunc){
          list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
          Xs.slope.0 <- list.data.GK.slope.0$Xtime
          Us.slope.0 <- list.data.GK.slope.0$Utime
        }

      }

      if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
        Z_CR <- list.surv$Z_CR
      }
      else{
        if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
          Z_CR <- list.surv$Z_CR
        }
        else{
          if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
            Z_CR <- list.surv$Z_CR[,-1]
            pp <- seq(0,1, length.out = FlexVarJM$control$ord.splines)
            pp <- tail(head(pp,-1),-1)
            tt2 <- as.data.frame(cbind(Time,event2))
            tt <- tt2$Time[which(tt2$event2 == 1)]
            kn <- quantile(tt, pp, names = FALSE)
            kn <- kn[kn<max(Time)]
            rr <- sort(c(rep(range(Time,0), 4L), kn))
            B.CR <- splineDesign(rr, Time, ord = 4L)
            Bs.CR <- splineDesign(rr, c(t(st_calc)), ord = 4L)
            if(FlexVarJM$control$left_trunc){
              Bs.0.CR <- splineDesign(rr, c(t(st.0)), ord = 4L)
            }
          }
          else{
            stop("This type of base survival function is not implemented.")
          }
        }

      }
    }
  }

  #cat("Prediction \n")
  shape <- 0
  shape.CR <- 0
  estim_param <- FlexVarJM$table.res$Estimation
  pred.e.a.table <- c()
  pred.CV <- c()

  borne1 <- choose(n = FlexVarJM$control$nb.e.a, k = 2) + FlexVarJM$control$nb.e.a
  C <- matrix(rep(0,(FlexVarJM$control$nb.e.a)**2),nrow=FlexVarJM$control$nb.e.a,ncol=FlexVarJM$control$nb.e.a)
  C[lower.tri(C, diag=T)] <- estim_param[1:borne1]
  Sigma.b <- C%*%t(C)
  borne2 <- borne1 + FlexVarJM$control$nb.priorMean.beta
  beta <- estim_param[(borne1+1):borne2]
  borne3 <- borne2 + FlexVarJM$control$nb.alpha
  if(borne3 != borne2){
    alpha <- estim_param[(borne2+1):borne3]
  }
  else{
    alpha <- 0
  }
  alpha.CR <- 0
  if(FlexVarJM$control$competing_risk){
    borne4 <- borne3 + FlexVarJM$control$nb.alpha.CR
    if(borne4 != borne3){
      alpha.CR <- estim_param[(borne3+1):borne4]
    }
    else{
      alpha.CR <- 0
    }
  }
  else{
    borne4 <- borne3
  }
  mu.log.sigma <- 0
  tau.log.sigma <- 0
  alpha.sigma <- 0
  alpha.sigma.CR <- 0
  sigma.epsilon <- 0
  if(FlexVarJM$control$variability_hetero){
    mu.log.sigma <- estim_param[borne4+1]
    tau.log.sigma <- abs(estim_param[borne4+2])
    alpha.sigma <- estim_param[borne4+3]
    borne5 <- borne4+3
    if(FlexVarJM$control$competing_risk){
      alpha.sigma.CR <- estim_param[borne4+4]
      borne5 <- borne4+4
    }
  }
  else{
    sigma.epsilon <- abs(estim_param[borne4+1])
    borne5 <- borne4+1
  }
  curseur <- borne5
  alpha.current <- 0
  alpha.current.CR <- 0
  alpha.slope <- 0
  alpha.slope.CR <- 0
  shape <- 0
  gamma <- 0
  shape.CR <- 0
  gamma.CR <- 0
  if(FlexVarJM$control$sharedtype %in% c("random effects")){
    stop("Not implemented yet")
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("random effects")){
    stop("Not implemented yet")
  }
  if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
    alpha.current <- estim_param[curseur+1]
    curseur <- curseur+1
  }
  if(FlexVarJM$control$competing_risk &&FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
    alpha.current.CR <- estim_param[curseur+1]
    curseur <- curseur +1
  }
  if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
    alpha.slope <- estim_param[curseur+1]
    curseur <- curseur +1
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
    alpha.slope.CR <- estim_param[curseur+1]
    curseur <- curseur +1
  }
  if(FlexVarJM$control$hazard_baseline == "Weibull"){
    shape <- estim_param[curseur+1]**2
    curseur <- curseur +1
  }
  if(FlexVarJM$control$hazard_baseline == "Splines"){
    gamma <- estim_param[(curseur+1):(curseur+FlexVarJM$control$ord.splines+2)]
    curseur <- curseur + FlexVarJM$control$ord.splines + 2
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Weibull"){
    shape.CR <- estim_param[curseur+1]**2
    curseur <- curseur +1
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Splines"){
    gamma.CR <- estim_param[(curseur+1):(curseur+FlexVarJM$control$ord.splines+2)]
    curseur <- curseur + FlexVarJM$control$ord.splines + 2
  }

  #print(head(X_base))
  Cum_risk2 <- c()
  Cum_risk1 <- c()
  Time.sort.unique <- unique(sort(Time))
  data.GaussKronrod.sort.unique <- data.GaussKronrod(data.id = data.id, Time = Time.sort.unique, k = FlexVarJM$control$nb_pointsGK)
  st_calc.sort.unique <- data.GaussKronrod.sort.unique$st
  P.sort.unique <- data.GaussKronrod.sort.unique$P
  for(i in 1:Ind){
    #print(i)
    Cum_risk_2i <- c()
    Cum_risk_1i <- c()
    #print(i)
    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    U_i <- U[offset[i]:(offset[i+1]-1),]
    y_i <- y.new.prog[offset[i]:(offset[i+1]-1)]
    Xtime_i <- Xtime[i,]
    Utime_i <- Utime[i,]
    Xs_i <- Xs[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Us_i <- Us[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Xslope_i <- Xslope[i,]
    Uslope_i <- Uslope[i,]
    Xs.slope_i <- Xs.slope[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Us.slope_i <- Us.slope[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Time_i <- Time[i]
    st_i <- st_calc[i,]
    B_i <- B[i,]
    Bs_i <- Bs[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Z_i <- Z[i,]
    P_i <- P[i]
    event1_i <- event1[i]
    event2_i <- event2[i]
    B.CR_i <- B.CR[i,]
    Bs.CR_i <- Bs.CR[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),]
    Z.CR_i <- Z_CR[i,]
    sigma_i <- mu.log.sigma
    if(is.null(dim(U_i))){
      U_i <- matrix(U_i, nrow= 1)
      X_base_i <- matrix(X_base_i, nrow= 1)
    }
    bi <- Sigma.b%*%t(U_i)%*%solve(U_i%*%Sigma.b%*%t(U_i) + diag(rep(exp(sigma_i),nrow(U_i)), nrow = nrow(U_i), ncol = nrow(U_i)))%*%(y_i - X_base_i%*%beta)

    ### Longitudinal prediction
    #print("Longitudinal part")
    binit <- c(bi,sigma_i)

    pred.e.a <- marqLevAlg(binit, fn = predict_re, minimize = FALSE,
                           nb.e.a = FlexVarJM$control$nb.e.a, X_base_i = X_base_i, beta = beta, U_i =U_i, y_i=y_i,
                           Sigma.b = Sigma.b, mu.log.sigma = mu.log.sigma, tau.log.sigma = tau.log.sigma,
                           variability_hetero = FlexVarJM$control$variability_hetero,
                           alpha.sigma = alpha.sigma, competing_risk = FlexVarJM$control$competing_risk, alpha.sigma.CR = alpha.sigma.CR,
                           sharedtype = FlexVarJM$control$sharedtype, sharedtype_CR = FlexVarJM$control$sharedtype_CR,
                           Xtime_i = Xtime_i, Utime_i = Utime_i, Xs_i = Xs_i, Us_i = Us_i,
                           alpha.current = alpha.current, alpha.current.CR = alpha.current.CR,
                           Xslope_i = Xslope_i, Uslope_i = Uslope_i, Xs.slope_i = Xs.slope_i, Us.slope_i = Us.slope_i,
                           alpha.slope = alpha.slope, alpha.slope.CR = alpha.slope.CR, indices_beta_slope = FlexVarJM$control$indices_beta_slope,
                           hazard_baseline = FlexVarJM$control$hazard_baseline,  wk = wk, shape = shape,
                           Time_i = Time_i,  st_i = st_i,gamma = gamma, B_i = B_i, Bs_i = Bs_i, Z_i = Z_i, alpha = alpha, P_i = P_i,
                           hazard_baseline_CR = FlexVarJM$control$hazard_baseline_CR, shape.CR = shape.CR, gamma.CR = gamma.CR, B.CR_i = B.CR_i, Bs.CR_i = Bs.CR_i,
                           Z.CR_i = Z.CR_i,alpha.CR = alpha.CR,event1_i = event1_i, event2_i=event2_i,
                           nproc = FlexVarJM$control$nproc,
                           clustertype = FlexVarJM$control$clustertype,
                           maxiter = FlexVarJM$control$maxiter, print.info = F,blinding = TRUE, epsa = FlexVarJM$control$epsa,
                           epsb = FlexVarJM$control$epsb, epsd = FlexVarJM$control$epsd)



    #print(pred.e.a$b)
    pred.e.a.table <- rbind(pred.e.a.table,c(data.id$ID[i], pred.e.a$b))
    CV <- X_base_i%*%beta + U_i%*%pred.e.a$b[1:(FlexVarJM$control$nb.e.a)]
    pred.CV <- rbind(pred.CV,cbind(rep(data.id$ID[i],length(CV)), CV))

    ### Survival goodness-of-fit
    #print("Survival part")

    for(j in 1:nrow(st_calc.sort.unique)){

      if(FlexVarJM$control$variability_hetero){
        pred_haz <- alpha.sigma*exp(pred.e.a$b[FlexVarJM$control$nb.e.a+1])
      }
      else{
        pred_haz <- 0
      }

      if(FlexVarJM$control$sharedtype %in% c("CV","CVS") || (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")) ){

        list.data.GK.current.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),], st_calc.sort.unique[j,],
                                                      FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)

        Xs.j <- list.data.GK.current.sort.unique$Xtime
        Us.j <- list.data.GK.current.sort.unique$Utime

        current.GK <-beta%*%t(Xs.j) + pred.e.a$b[1:(FlexVarJM$control$nb.e.a)]%*%t(Us.j)

        if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
          pred_haz <- pred_haz +  alpha.current*current.GK
        }
      }

      if(FlexVarJM$control$sharedtype %in% c("CVS","S") || (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS","S") )){
        list.data.GK.slope.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(FlexVarJM$control$nb_pointsGK*(i-1)+1):(FlexVarJM$control$nb_pointsGK*i),], st_calc.sort.unique[j,],
                                                    FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)

        Xs.slope.j <- list.data.GK.slope.sort.unique$Xtime
        Us.slope.j <- list.data.GK.slope.sort.unique$Utime

        slope.GK <- beta[indices_beta_slope]%*%t(Xs.slope.j) +  pred.e.a$b[2:(FlexVarJM$control$nb.e.a)]%*%t(Us.slope.j)

        if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
          pred_haz <- pred_haz +  alpha.slope*slope.GK
        }
      }

      if(FlexVarJM$control$hazard_baseline == "Exponential"){
        h_0 <- 1
        h_0.GK <- wk
      }
      if(FlexVarJM$control$hazard_baseline == "Weibull"){
        st_j <- st_calc.sort.unique[j,]
        h_0.GK <- shape*(st_j**(shape-1))*wk
      }
      if(FlexVarJM$control$hazard_baseline == "Splines"){
        st_j <- st_calc.sort.unique[j,]
        Bs_j <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, st_j, ord = 4L)
        #Bs_j <- Bs[(FlexVarJM$control$nb_pointsGK*(j-1)+1):(FlexVarJM$control$nb_pointsGK*j),]
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

      if(FlexVarJM$control$competing_risk){
        if(FlexVarJM$control$variability_hetero){
          pred_haz.CR <- alpha.sigma.CR*exp(pred.e.a$b[FlexVarJM$control$nb.e.a+1])
        }
        else{
          pred_haz.CR <- 0
        }

        if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
          pred_haz.CR <- pred_haz.CR + alpha.current.CR*current.GK
        }

        if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
          pred_haz.CR <- pred_haz.CR + alpha.slope.CR*slope.GK
        }


        if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
          h_0.GK.CR <- wk
        }
        if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
          st_j <- st_calc.sort.unique[j,]
          h_0.GK.CR <- shape.CR*(st_j**(shape.CR-1))*wk
        }
        if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
          st_j <- st_calc.sort.unique[j,]
          Bs_j <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, st_j, ord = 4L)
          #Bs_j <- Bs.CR[(FlexVarJM$control$nb_pointsGK*(j-1)+1):(FlexVarJM$control$nb_pointsGK*j),]
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
  }
  result <- list(pred.e.a.table = pred.e.a.table,
                 pred.CV = pred.CV, Cum_risk1 = Cum_risk1, Cum_risk2 = Cum_risk2)
  graphs <- NULL
  if(graph){
    timeInterv <- range(data.long[,FlexVarJM$control$timeVar])
    if(is.null(break.times)) break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
    graphs <- plot.goodnessoffit(data.long,data.id,pred.CV,break.times, formFixed = FlexVarJM$control$formFixed, formSurv = FlexVarJM$control$formSurv,
                                 timeVar = FlexVarJM$control$timeVar,Cum_risk1, competing_risk = FlexVarJM$control$competing_risk, formSurv_CR = FlexVarJM$control$formSurv_CR, Cum_risk2)
  }

  result.final <- list(tables = result,
                       graphs = graphs)
  result.final
}
