#' Title
#'
#' @param newdata
#' @param FlexVarJM
#' @param s
#' @param t
#' @param event
#' @param L
#'
#' @return
#' @export
#'
#' @examples
pred_s.t <- function(newdata,FlexVarJM, s, t, event, L){#L: nb of bootsrap
  result <- c()
  newdata <- as.data.frame(newdata)
  newdata$id <- 1

  data.long.until.time.s <-newdata[which(newdata[,FlexVarJM$control$timeVar] <=s),]
  name.time.event <- all.vars(FlexVarJM$control$formSurv)[1]
  name.event.event <- all.vars(FlexVarJM$control$formSurv)[2]
  data.long.until.time.s[which(data.long.until.time.s[,name.time.event]>=s),name.event.event] <- 0
  data.long.until.time.s[which(data.long.until.time.s[,name.time.event]>=s),name.time.event] <- max(data.long.until.time.s[,FlexVarJM$control$timeVar])

  data.long.until.time.s.id <- data.long.until.time.s[1,]

  # Longitudinal data
  list.long <- data.manag.long(FlexVarJM$control$formGroup,FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,data.long.until.time.s)
  X_base <- list.long$X
  U <- list.long$U
  y.new.prog <- list.long$y.new.prog

  # Survival data
  ### Between s and s+t
  data.GaussKronrod.1 <- data.GaussKronrod2(data.long.until.time.s.id,a=s,b=s+t,k = FlexVarJM$control$nb_pointsGK)
  P.1 <- data.GaussKronrod.1$P
  st.1 <- data.GaussKronrod.1$st
  wk.1 <- data.GaussKronrod.1$wk
  data.id.1 <- data.GaussKronrod.1$data.id2

  ##########Computing little lambda############
  ### Matrix for current value and slope
  if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                      FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
    Xs <- list.data.GK.current$Xtime
    Us <- list.data.GK.current$Utime
    #if(FlexVarJM$control$left_trunc){
    #  list.data.GK.current.0 <-  data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
    #                                      FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
    #  Xs.0 <- list.data.GK.current.0$Xtime
    #  Us.0 <- list.data.GK.current.0$Utime
    #}
  }
  if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
    list.data.GK.slope <-  data.time(data.id.1, c(t(st.1)),
                                    FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
    Xs.slope <- list.data.GK.slope$Xtime
    Us.slope <- list.data.GK.slope$Utime
    #if(FlexVarJM$control$left_trunc){
    #  list.data.GK.slope.0 <-  data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
    #                                    FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
    #  Xs.slope.0 <- list.data.GK.slope.0$Xtime
    #  Us.slope.0 <- list.data.GK.slope.0$Utime
    #}
  }

  #### lambda0
  if(FlexVarJM$control$hazard_baseline == "Exponential"){
    mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
    Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
  }else{
    if(FlexVarJM$control$hazard_baseline == "Weibull"){
      mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
      Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
    }else{
      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
        Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
        Z <- Z[,-1]
        #print((st.1))
        Bs <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, c(t(st.1)), ord = 4L)
        #print("ok")
        if(FlexVarJM$control$left_trunc){
          Bs.0 <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
        }
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }

  }

  ####Same for competing risks
  if(FlexVarJM$control$competing_risk){
    if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
      list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                        FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      Xs <- list.data.GK.current$Xtime
      Us <- list.data.GK.current$Utime
      #if(FlexVarJM$control$left_trunc){
      #  list.data.GK.current.0 <-  data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
      #                                      FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      #  Xs.0 <- list.data.GK.current.0$Xtime
      #  Us.0 <- list.data.GK.current.0$Utime
      #}
    }
    if(FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
      list.data.GK.slope <-  data.time(data.id.1, c(t(st.1)),
                                      FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      Xs.slope <- list.data.GK.slope$Xtime
      Us.slope <- list.data.GK.slope$Utime
      #if(FlexVarJM$control$left_trunc){
      #  list.data.GK.slope.0 <-  data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
      #                                    FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      #  Xs.slope.0 <- list.data.GK.slope.0$Xtime
      #  Us.slope.0 <- list.data.GK.slope.0$Utime
      #}

    }

    if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
      mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
      Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
    }else{
      if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
        mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
        Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
      }else{
        if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
          mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
          Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
          Z_CR <- Z_CR[,-1]
          Bs.CR <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, c(t(st.1)), ord = 4L)
          if(FlexVarJM$control$left_trunc){
            Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
          }
        }else{
          stop("This type of base survival function is not implemented.")
        }
      }

    }
  }

  ### At s
  if(FlexVarJM$control$sharedtype %in% c("random effects")){
    stop("Not implemented yet")
  }else{
    list.GaussKronrod <- data.GaussKronrod(data.long.until.time.s.id, s, k = FlexVarJM$control$nb_pointsGK)
    wk.den <- list.GaussKronrod$wk
    st_calc.den <- list.GaussKronrod$st
    P.den <- list.GaussKronrod$P
    id.GK.den <- list.GaussKronrod$id.GK
    if(FlexVarJM$control$sharedtype %in% c("CV","CVS") ){
      list.data.current.time <-  data.time(data.long.until.time.s.id,s, FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      list.data.GK.current <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                        FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      Xtime.den <- list.data.current.time$Xtime
      Utime.den <- list.data.current.time$Utime
      Xs.den <- list.data.GK.current$Xtime
      Us.den <- list.data.GK.current$Utime
      #if(left_trunc){
      #  list.data.GK.current.0 <-  data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
      #                                      formFixed, formRandom,timeVar)
      #  Xs.0 <- list.data.GK.current.0$Xtime
      #  Us.0 <- list.data.GK.current.0$Utime
      #}
    }
    if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
      list.data.slope.time <-  data.time(data.long.until.time.s.id, s, FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      list.data.GK.slope <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                       FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      Xslope.den <- list.data.slope.time$Xtime
      Uslope.den <- list.data.slope.time$Utime
      Xs.slope.den <- list.data.GK.slope$Xtime
      Us.slope.den <- list.data.GK.slope$Utime
    }
  }
     if(FlexVarJM$control$hazard_baseline == "Splines"){
       Bs.den <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, c(t(st_calc.den)), ord = 4L)
      }
  if(FlexVarJM$control$competing_risk){
    if(FlexVarJM$control$sharedtype_CR %in% c("RE")){
      stop("Not implemented yet")
    }else{
      list.GaussKronrod <- data.GaussKronrod(data.long.until.time.s.id, s, k = nb_pointsGK)
      wk.den <- list.GaussKronrod$wk
      st_calc.den <- list.GaussKronrod$st
      P.den <- list.GaussKronrod$P
      id.GK.den <- list.GaussKronrod$id.GK
      if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
        #list.data.current.time <-  data.time(data.long.until.time.s.id, s, FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        list.data.GK.current <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                          FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xs.den <- list.data.GK.current$Xtime
        Us.den <- list.data.GK.current$Utime
      }
      if(FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
        #list.data.slope.time <-  data.time(data.until.time.s.id, list.surv$Time, formSlopeFixed, formSlopeRandom,timeVar)
        list.data.GK.slope <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                         FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        #Xslope <- list.data.slope.time$Xtime
        #Uslope <- list.data.slope.time$Utime
        Xs.slope.den <- list.data.GK.slope$Xtime
        Us.slope.den <- list.data.GK.slope$Utime
      }
    }
        if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
          Bs.CR.den <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, c(t(st_calc.den)), ord = 4L)
        }
  }
  if(FlexVarJM$control$variability_hetero){
    Zq <- sobol(FlexVarJM$control$S2, FlexVarJM$control$nb.e.a+1, normal = TRUE, scrambling = 1)
  }else{
    Zq <- sobol(FlexVarJM$control$S2, FlexVarJM$control$nb.e.a, normal = TRUE, scrambling = 1)
  }
  Ncpus <- FlexVarJM$control$nproc
  cl <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(cl)
  Hess <- matrix(rep(0,length(FlexVarJM$result$grad)**2),nrow=length(FlexVarJM$result$grad),ncol=length(FlexVarJM$result$grad))
  Hess[upper.tri(Hess, diag=T)] <- FlexVarJM$result$v
  Hess2 = Hess + t(Hess)
  diag(Hess2) <- diag(Hess2) - diag(Hess)
  result <- c()
  res <- foreach::foreach(l=1:L, .combine='c',.packages=c("survival","splines","FlexVarJM","mvtnorm")) %dopar%{
  #for(l in 1:L){

    tirage <- rmvnorm(1, mean = FlexVarJM$table.res$Estimation, sigma = Hess2)
    borne1 <- choose(n = FlexVarJM$control$nb.e.a, k = 2) + FlexVarJM$control$nb.e.a
    C <- matrix(rep(0,(FlexVarJM$control$nb.e.a)**2),nrow=FlexVarJM$control$nb.e.a,ncol=FlexVarJM$control$nb.e.a)
    C[lower.tri(C, diag=T)] <- tirage[1:borne1]
    Sigma.b <- C%*%t(C)
    borne2 <- borne1 + FlexVarJM$control$nb.priorMean.beta
    beta <- tirage[(borne1+1):borne2]
    borne3 <- borne2 + FlexVarJM$control$nb.alpha
    if(borne3 != borne2){
      alpha <- tirage[(borne2+1):borne3]
    }else{
      alpha <- 0
    }
    alpha.CR <- 0
    if(FlexVarJM$control$competing_risk){
      borne4 <- borne3 + FlexVarJM$control$nb.alpha.CR
      if(borne4 != borne3){
        alpha.CR <- tirage[(borne3+1):borne4]
      }else{
        alpha.CR <- 0
      }
    }else{
      borne4 <- borne3
    }
    mu.log.sigma <- 0
    tau.log.sigma <- 0
    alpha.sigma <- 0
    alpha.sigma.CR <- 0
    sigma.epsilon <- 0
    if(FlexVarJM$control$variability_hetero){
      mu.log.sigma <- tirage[borne4+1]
      tau.log.sigma <- abs(tirage[borne4+2])
      alpha.sigma <- tirage[borne4+3]
      borne5 <- borne4+3
      if(FlexVarJM$control$competing_risk){
        alpha.sigma.CR <- tirage[borne4+4]
        borne5 <- borne4+4
      }
    }else{
      sigma.epsilon <- abs(tirage[borne4+1])
      borne5 <- borne4+1
    }
    curseur <- borne5
    if(FlexVarJM$control$sharedtype %in% c("random effects")){
      stop("Not implemented yet")
    }
    if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("random effects")){
      stop("Not implemented yet")
    }
    if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
      alpha.current <- tirage[curseur+1]
      curseur <- curseur+1
    }
    if(FlexVarJM$control$competing_risk &&FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
      alpha.current.CR <- tirage[curseur+1]
      curseur <- curseur +1
    }
    if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
      alpha.slope <- tirage[curseur+1]
      curseur <- curseur +1
    }
    if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
      alpha.slope.CR <- tirage[curseur+1]
      curseur <- curseur +1
    }
    if(FlexVarJM$control$hazard_baseline == "Weibull"){
      shape <- tirage[curseur+1]**2
      curseur <- curseur +1
    }
    if(FlexVarJM$control$hazard_baseline == "Splines"){
      gamma <- tirage[(curseur+1):(curseur+FlexVarJM$control$ord.splines+2)]
      curseur <- curseur + FlexVarJM$control$ord.splines + 2
    }
    if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Weibull"){
      shape.CR <- tirage[curseur+1]**2
      curseur <- curseur +1
    }
    if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Splines"){
      gamma.CR <- tirage[(curseur+1):(curseur+FlexVarJM$control$ord.splines+2)]
      curseur <- curseur + FlexVarJM$control$ord.splines + 2
    }


    b_al <- Zq[,1:(FlexVarJM$control$nb.e.a)]%*%t(C)
    log.sigma_al <- matrix(rep(mu.log.sigma,FlexVarJM$control$S2), nrow = FlexVarJM$control$S2, byrow = T) + Zq[,FlexVarJM$control$nb.e.a+1]*tau.log.sigma
    sigma_al <- exp(log.sigma_al)


    ###########################################################################
    ###########################################################################
    ############################ NUMERATEUR ##################################
    ###########################################################################
    ###########################################################################
    #####compute lambda
    h <- 1
    etaBaseline <- 0
    survLong <- 0
    etaBaseline.0 <- 0
    survLong.0 <- 0
    if(event==1){
      if(FlexVarJM$control$variability_hetero){
        h <- h*exp(alpha.sigma*sigma_al)
      }
      if(FlexVarJM$control$sharedtype %in% c("CV") ){
        current.GK <- matrix(rep(beta%*%t(Xs),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al%*%t(Us)
        h <- matrix(rep(h,ncol(current.GK)),ncol=ncol(current.GK))*exp(alpha.current*current.GK)
      }
      if(FlexVarJM$control$sharedtype %in% c("CVS") ){
        current.GK <- matrix(rep(beta%*%t(Xs),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al%*%t(Us)
        h <- matrix(rep(h,ncol(current.GK)),ncol=ncol(current.GK))*exp(alpha.current*current.GK)
        h <- h*exp(alpha.slope*slope.GK)
      }
      if(FlexVarJM$control$sharedtype %in% c("S") ){
        slope.GK <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope)
        h <- matrix(rep(h,ncol(slope.GK)),ncol=ncol(slope.GK))*exp(alpha.slope*slope.GK)
      }
      ###h0
      if(FlexVarJM$control$hazard_baseline == "Exponential"){
        h_0.GK <- wk.1
      }

      if(FlexVarJM$control$hazard_baseline == "Weibull"){
        h_0.GK <- shape*(st.1**(shape-1))*wk.1
      }

      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mat_h0s <- matrix(gamma,ncol=1)
        h_0.GK <- (wk.1*exp(Bs%*%mat_h0s))
      }

      ###hazard function
      if(length(Z)==0){
        pred_surv <- 0
      }else{
        pred_surv <- (alpha%*%Z)[1,1]
      }
      h <- h*exp(pred_surv)
      #h <- as.vector(h)
      #h_0.GK <- as.vector(h_0.GK)
      h <- matrix(rep(h_0.GK,nrow(h)),nrow = nrow(h),byrow = T)*h
    }else{
      if(FlexVarJM$control$variability_hetero){
        h <- h*exp(alpha.sigma.CR*sigma_al)
      }
      if(FlexVarJM$control$sharedtype %in% c("CV") ){
        current.GK <- matrix(rep(beta%*%t(Xs),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al%*%t(Us)
        h <- matrix(rep(h,ncol(current.GK)),ncol=ncol(current.GK))*exp(alpha.current.CR*current.GK)
      }
      if(FlexVarJM$control$sharedtype %in% c("CVS") ){
        current.GK <- matrix(rep(beta%*%t(Xs),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al%*%t(Us)
        h <- matrix(rep(h,ncol(current.GK)),ncol=ncol(current.GK))*exp(alpha.current.CR*current.GK)
        h <- h*exp(alpha.slope.CR*slope.GK)
      }
      if(FlexVarJM$control$sharedtype %in% c("S") ){
        slope.GK <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope)
        h <- matrix(rep(h,ncol(slope.GK)),ncol=ncol(slope.GK))*exp(alpha.slope.CR*slope.GK)
      }
      ###h0
      if(FlexVarJM$control$hazard_baseline == "Exponential"){
        h_0.GK <- wk.1
      }

      if(FlexVarJM$control$hazard_baseline == "Weibull"){
        h_0.GK <- shape.CR*(st.1**(shape.CR-1))*wk.1
      }

      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mat_h0s <- matrix(gamma,ncol=1)
        h_0.GK <- (wk.1*exp(Bs.CR%*%mat_h0s))
      }

      ###hazard function
      if(length(Z_CR)==0){
        pred_surv <- 0
      }else{
        pred_surv <- (alpha.CR%*%Z_CR)[1,1]
      }
      h <- h*exp(pred_surv)
      #h <- as.vector(h)
      #h_0.GK <- as.vector(h_0.GK)
      h <- matrix(rep(h_0.GK,nrow(h)),nrow = nrow(h),byrow = T)*h
    }

    #########################################
    ##### Computing LAMBDA1 and LAMBDA2 #####
    #########################################
    Gamma1 <- c()
    Gamma2 <- c()
    for(t2 in st.1){

      data.GaussKronrod.2 <-  data.GaussKronrod(data.long.until.time.s.id,t2,k = FlexVarJM$control$nb_pointsGK)
      P.2 <- data.GaussKronrod.2$P
      st.2 <- data.GaussKronrod.2$st
      wk.2 <- data.GaussKronrod.2$wk
      data.id.2 <- data.GaussKronrod.2$data.id2

      if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
        list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                            FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xs.2 <- list.data.GK.current.2$Xtime
        Us.2 <- list.data.GK.current.2$Utime
      }
      if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
        list.data.GK.slope.2 <-  data.time(data.id.2, c(t(st.2)),
                                          FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        Xs.slope.2 <- list.data.GK.slope.2$Xtime
        Us.slope.2 <- list.data.GK.slope.2$Utime
      }

      #### lambda0
      if(FlexVarJM$control$hazard_baseline == "Exponential"){
        mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
        Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
      }else{
        if(FlexVarJM$control$hazard_baseline == "Weibull"){
          mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
          Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
        }else{
          if(FlexVarJM$control$hazard_baseline == "Splines"){
            mfZ <- model.frame(FlexVarJM$control$formSurv, data = data.long.until.time.s.id)
            Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
            Z <- Z[,-1]
            Bs.2 <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, c(t(st.2)), ord = 4L)
            if(FlexVarJM$control$left_trunc){
              Bs.0 <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
            }
          }else{
            stop("This type of base survival function is not implemented.")
          }
        }

      }

      ####Same for competing risks
      if(FlexVarJM$control$competing_risk){
        if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
          list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                              FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
          Xs.2 <- list.data.GK.current.2$Xtime
          Us.2 <- list.data.GK.current.2$Utime
        }
        if(FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
          list.data.GK.slope.2 <-  data.time(data.id.2, c(t(st.2)),
                                            FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
          Xs.slope.2 <- list.data.GK.slope.2$Xtime
          Us.slope.2 <- list.data.GK.slope.2$Utime
        }

        if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
          mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
          Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
        }else{
          if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
            mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
            Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
          }else{
            if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
              mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = data.long.until.time.s.id)
              Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
              Z_CR <- Z_CR[,-1]
              Bs.CR.2 <- splines::splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, c(t(st.2)), ord = 4L)
              if(FlexVarJM$control$left_trunc){
                Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
              }
            }else{
              stop("This type of base survival function is not implemented.")
            }
          }

        }
      }


      h.2.1 <- 1
      h.2.2 <- 1
      if(FlexVarJM$control$variability_hetero){
        h.2.1 <- h.2.1*exp(alpha.sigma*sigma_al)
        if(FlexVarJM$control$competing_risk){
          h.2.2 <- h.2.2*exp(alpha.sigma.CR*sigma_al)
        }
      }
      if(FlexVarJM$control$sharedtype %in% c("CV")|| (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV"))){
        current.GK.2 <- matrix(rep(beta%*%t(Xs.2),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T)+ b_al%*%t(Us.2)
        if(FlexVarJM$control$sharedtype %in% c("CV")){
          h.2.1 <- matrix(rep(h.2.1,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current*current.GK.2)
        }
        if((FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV"))){
          h.2.2 <- matrix(rep(h.2.2,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current.CR*current.GK.2)
        }
      }
      if(FlexVarJM$control$sharedtype %in% c("CVS")|| (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS"))){
        current.GK.2 <- matrix(rep(beta%*%t(Xs.2),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T)+ b_al%*%t(Us.2)
        slope.GK.2 <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope.2),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T)+ b_al[,-1]%*%t(Us.slope.2)
        if(FlexVarJM$control$sharedtype %in% c("CVS")){
          h.2.1 <- matrix(rep(h.2.1,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current*current.GK.2)
          h.2.1 <- h.2.1*exp(alpha.slope*slope.GK.2)
        }
        if((FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS"))){
          h.2.2 <- matrix(rep(h.2.2,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current.CR*current.GK.2)
          h.2.2 <- h.2.2*exp(alpha.slope.CR*slope.GK.2)
        }
      }
      if(FlexVarJM$control$sharedtype %in% c("S")|| (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("S"))){
        slope.GK.2 <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope.2),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T)+ b_al[,-1]%*%t(Us.slope.2)
        if(FlexVarJM$control$sharedtype %in% c("S")){
          h.2.1 <- matrix(rep(h.2.1,ncol(slope.GK.2)),ncol=ncol(slope.GK.2))*exp(alpha.slope*slope.GK.2)

        }
        if((FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("S"))){
          h.2.2 <- matrix(rep(h.2.2,ncol(slope.GK.2)),ncol=ncol(slope.GK.2))*exp(alpha.slope.CR*slope.GK.2)
        }
      }

      ###h0
      if(FlexVarJM$control$hazard_baseline == "Exponential"){
        h_0.GK.2 <- wk.2
      }

      if(FlexVarJM$control$hazard_baseline == "Weibull"){
        h_0.GK.2 <- shape*(st.2**(shape-1))*wk.2
      }

      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mat_h0s <- matrix(gamma,ncol=1)
        h_0.GK.2 <- (wk.2*exp(Bs.2%*%mat_h0s))
      }

      ###hazard function
      if(length(Z)==0){
        pred_surv <- 0
      }else{
        pred_surv <- (alpha%*%Z)[1,1]
      }
      h.2.1 <- h.2.1*exp(pred_surv)
      h.2.1 <- matrix(rep(h_0.GK.2,nrow(h.2.1)),nrow = nrow(h.2.1),byrow = T)*h.2.1

      if(FlexVarJM$control$competing_risk){
        if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
          h_0.GK.2.CR <- wk.2
        }

        if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
          h_0.GK.2.CR <- shape.CR*(st.2**(shape.CR-1))*wk.2
        }

        if(FlexVarJM$control$hazard_baseline == "Splines"){
          mat_h0s <- matrix(gamma.CR,ncol=1)
          h_0.GK.2.CR <- (wk.2*exp(Bs.CR.2%*%mat_h0s))
        }

        ###hazard function
        if(length(Z_CR)==0){
          pred_surv.CR <- 0
        }else{
          pred_surv.CR <- (alpha.CR%*%Z_CR)[1,1]
        }
        h.2.2 <- h.2.2*exp(pred_surv.CR)
        h.2.2 <- matrix(rep(h_0.GK.2.CR,nrow(h.2.2)),nrow = nrow(h.2.2),byrow = T)*h.2.2
      }
      Gamma1 <- cbind(Gamma1, (t2/2)*rowSums(h.2.1))
      #matrix(rep(wk.2,nrow(h.2.1)), nrow = nrow(h.2.1), byrow = T)*
      if(FlexVarJM$control$competing_risk){
        Gamma2 <- cbind(Gamma2,(t2/2)*rowSums(h.2.2))
        #matrix(rep(wk.2,nrow(h.2.2)), nrow = nrow(h.2.2), byrow = T)*
      }

    }

    int <- exp(-Gamma1-Gamma2)*h

    surv.num <- P.1*rowSums(matrix(rep(wk.1,nrow(int)), nrow = nrow(int), byrow = T)*int)


    #Longitudinal part
    if(FlexVarJM$control$variability_hetero  == T){
      sigma.long <- sigma_al
    }else{
      sigma.long <- sigma.epsilon
    }
    if(is.null(nrow(X_base))){
      CV <- (beta%*%X_base)[1,1] + b_al%*%U
      f_Y_b_sigma <- dnorm(x=y.new.prog, mean = CV, sd = sigma.long)
    }else{
      f_Y_b_sigma <- rep(1,FlexVarJM$control$S2)
      for(k in 1:nrow(X_base)){
        CV <- (beta%*%X_base[k,])[1,1] + b_al%*%U[k,]
        f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y.new.prog[k], mean = CV, sd = sigma.long)
      }
    }

    numerateur <- surv.num*f_Y_b_sigma
    numerateur <- mean(numerateur)

    ###########################################################################
    ###########################################################################
    ############################ DENOMINATEUR##################################
    ###########################################################################
    ###########################################################################

    etaBaseline <- 0
    survLong <- 0
    #etaBaseline.0 <- 0
    #survLong.0 <- 0
    if(FlexVarJM$control$variability_hetero){
      etaBaseline <- etaBaseline + alpha.sigma*sigma_al
    }
    if(FlexVarJM$control$competing_risk){
      etaBaseline_CR <- 0
      survLong_CR <- 0
      if(FlexVarJM$control$variability_hetero){
        etaBaseline_CR <- etaBaseline_CR + alpha.sigma.CR*sigma_al
      }
    }
    if(FlexVarJM$control$sharedtype %in% c("CV","CVS") || (FlexVarJM$control$competing_risk && FlexVarJM$sharedtype_CR %in% c("CV","CVS")) ){
      current.GK <- matrix(rep(beta%*%t(Xs.den),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al%*%t(Us.den)
      if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
        survLong <- survLong + alpha.current*current.GK
      }
      if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
        survLong_CR <- survLong_CR + alpha.current.CR*current.GK
      }
    }
    if(FlexVarJM$control$sharedtype %in% c("CVS","S") || (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS","S")) ){
      slope.GK <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope.den),FlexVarJM$control$S2),nrow=FlexVarJM$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope.den)
      if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
        survLong <- survLong + alpha.slope*slope.GK
      }
      if(FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
        survLong_CR <- survLong_CR + alpha.slope.CR*slope.GK
      }
    }

    ###h0
    if(FlexVarJM$control$hazard_baseline == "Exponential"){
      h_0 <- 1
      h_0.GK <- wk.den
      if(left_trunc){
        h_0.GK.0 <- wk.den
      }
    }
    if(FlexVarJM$control$hazard_baseline == "Weibull"){
      h_0 <- shape*(s**(shape-1))
      h_0.GK <- shape*(st_calc.den**(shape-1))*wk.den
    }
    if(FlexVarJM$control$hazard_baseline == "Splines"){
      mat_h0s <- matrix(gamma,ncol=1)
      h_0.GK <- (wk.den*exp(Bs.den%*%mat_h0s))
    }

    ###hazard function
    if(length(Z)==0){
      pred_surv <- 0
    }else{
      pred_surv <- (alpha%*%Z)[1,1]
    }
    etaBaseline <- etaBaseline + pred_surv
    #if(left_trunc){
    #  etaBaseline.0 <- etaBaseline.0 + pred_surv
    #}

    ###GK integration
    survLong <- exp(survLong)
    h_0.GK <- as.vector(h_0.GK)
    survLong <- survLong%*%h_0.GK


    Surv <- (-exp(etaBaseline)*P.den*survLong)

    #if(left_trunc){###Computation of S(T0i)
    #  #stop("Not implemented yet.")
    #  survLong.0 <- exp(survLong.0)
    #  survLong.0 <- survLong.0%*%h_0.GK.0
    #  P.0_i <- P.0[i]
    #  Surv.0 <- exp((-exp(etaBaseline.0)*P.0_i*survLong.0))
    #}

    if(FlexVarJM$control$competing_risk){
      ###h0
      if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
        h_0.CR <- 1
        h_0.GK.CR <- wk.den
        if(left_trunc){
          h_0.GK.0_CR <- wk.den
        }
      }
      if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
        h_0.GK.CR <- shape.CR*(st_calc.den**(shape.CR-1))*wk.den
        #if(left_trunc){
        #  st.0_i <- st.0[i,]
        #  h_0.GK.0_CR <- shape.CR*(st.0_i**(shape.CR-1))*wk
        #}
      }
      if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
        mat_h0s.CR <- matrix(gamma.CR,ncol=1)
        h_0.GK.CR <- (wk.den*exp(Bs.CR.den%*%mat_h0s.CR))
        #if(left_trunc){
        #  Bs.0.CR_i <- Bs.0.CR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        #  h_0.GK.0_CR <- (wk*exp(Bs.0.CR_i%*%mat_h0s.CR))
        #}
      }

      ###hazard function
      if(length(Z_CR)==0){
        pred_surv.CR <- 0
      }else{
        pred_surv.CR <- (alpha.CR%*%Z_CR)[1,1]
      }
      etaBaseline_CR <- etaBaseline_CR + pred_surv.CR
      #if(left_trunc){
      #  etaBaseline.0_CR <- etaBaseline.0_CR + pred_surv.CR
      #}

      ###GK integration
      survLong_CR <- exp(survLong_CR)
      h_0.GK.CR <- as.vector(h_0.GK.CR)
      survLong_CR <- survLong_CR%*%h_0.GK.CR
      Surv.CR <- (-exp(etaBaseline_CR)*P.den*survLong_CR)

      #if(left_trunc){###Computation of S(T0i)
      #  #stop("Not implemented yet.")
      #  survLong.0_CR <- exp(survLong.0_CR)
      #  survLong.0_CR <- survLong.0_CR%*%h_0.GK.0_CR
      #  P.0_i <- P.0[i]
      #  Surv.0.CR <- exp((-exp(etaBaseline.0_CR)*P.0_i*survLong.0_CR))
      #}
    }

    denominateur <- exp(Surv+Surv.CR)*f_Y_b_sigma
    denominateur <- mean(denominateur)
    numerateur/denominateur
    result <- c(result, numerateur/denominateur)

  }
  parallel::stopCluster(cl)
  res
}












