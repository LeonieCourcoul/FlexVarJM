pred_surv <- function(FlexVarJM, s, t, re.b, re.sigma, newdata.id, event){

  newdata.id$id <- 1
  data.GaussKronrod.1 <- data.GaussKronrod2(newdata.id,a=s,b=s+t,k = FlexVarJM$control$nb_pointsGK)
  P.1 <- data.GaussKronrod.1$P
  st.1 <- data.GaussKronrod.1$st
  wk.1 <- data.GaussKronrod.1$wk
  data.id.1 <- data.GaussKronrod.1$data.id2

  ##########Computing little lambda############
  ### Matrix for current value and slope
  if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
    list.data.GK.current <- data.time(data.id.1, c(t(st.1)),
                                      FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
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
    list.data.GK.slope <- data.time(data.id.1, c(t(st.1)),
                                    FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
    Xs.slope <- list.data.GK.slope$Xtime
    Us.slope <- list.data.GK.slope$Utime
    if(FlexVarJM$control$left_trunc){
      list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                        FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      Xs.slope.0 <- list.data.GK.slope.0$Xtime
      Us.slope.0 <- list.data.GK.slope.0$Utime
    }
  }

  #### lambda0
  if(FlexVarJM$control$hazard_baseline == "Exponential"){
    mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
    Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
  }
  else{
    if(FlexVarJM$control$hazard_baseline == "Weibull"){
      mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
      Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
    }
    else{
      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
        Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
        Z <- Z[,-1]
        #print((st.1))
        Bs <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, c(t(st.1)), ord = 4L)
        #print("ok")
        if(FlexVarJM$control$left_trunc){
          Bs.0 <- splineDesign(rr, c(t(st.0)), ord = 4L)
        }
      }
      else{
        stop("This type of base survival function is not implemented.")
      }
    }

  }

  ####Same for competing risks
  nb.alpha.CR <- 0
  if(FlexVarJM$control$competing_risk){
      if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
        list.data.GK.current <- data.time(data.id.1, c(t(st.1)),
                                          FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
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
        list.data.GK.slope <- data.time(data.id.1, c(t(st.1)),
                                        FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
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
        mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
        Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
      }
      else{
        if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
          mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
          Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
        }
        else{
          if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
            mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
            Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
            Z_CR <- Z_CR[,-1]
            Bs.CR <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, c(t(st.1)), ord = 4L)
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

  shape <- 0
  shape.CR <- 0
  estim_param <- FlexVarJM$table.res$Estimation
  pred.e.a.table <- c()
  pred.CV <- c()

  borne1 <- choose(n = FlexVarJM$control$nb.e.a, k = 2) + FlexVarJM$control$nb.e.a
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
    gamma <- estim_param[(curseur+1):(curseur+FlexVarJM$control$nb.knots.splines+2+2)]
    #print(gamma)
    curseur <- curseur + FlexVarJM$control$nb.knots.splines+2 + 2
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Weibull"){
    shape.CR <- estim_param[curseur+1]**2
    curseur <- curseur +1
  }
  if(FlexVarJM$control$competing_risk && FlexVarJM$control$hazard_baseline_CR == "Splines"){
    gamma.CR <- estim_param[(curseur+1):(curseur+FlexVarJM$control$nb.knots.splines+2+2)]
    curseur <- curseur + FlexVarJM$control$nb.knots.splines +2 + 2
  }

  bi <- re.b
  sigmai <- re.sigma

  #####compute lambda
  h <- 1
  etaBaseline <- 0
  survLong <- 0
  etaBaseline.0 <- 0
  survLong.0 <- 0
  if(event==1){
    if(FlexVarJM$control$variability_hetero){
      #print(alpha.sigma)
      #print(sigmai)
      h <- h*exp(alpha.sigma*sigmai)
    }
    if(FlexVarJM$control$sharedtype %in% c("CV","CVS") ){
      current.GK <- beta%*%t(Xs) + bi%*%t(Us)
      h <- h*exp(alpha.current*current.GK)
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
    }
    else{
      pred_surv <- (alpha%*%Z)[1,1]
    }
    h <- h*exp(pred_surv)
    h <- as.vector(h)
    h_0.GK <- as.vector(h_0.GK)
    h <- h_0.GK*h
  }
  else{
    if(FlexVarJM$control$variability_hetero){
      h <- h*exp(alpha.sigma.CR*sigmai)
    }
    if((FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")) ){
      current.GK <- beta%*%t(Xs) + bi%*%t(Us)
      h <- h*exp(alpha.current.CR*current.GK)
    }
    ###h0
    if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
      h_0.GK <- wk.1
    }

    if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
      h_0.GK <- shape.CR*(st.1**(shape.CR-1))*wk.1
    }

    if(FlexVarJM$control$hazard_baseline == "Splines"){
      mat_h0s <- matrix(gamma.CR,ncol=1)
      h_0.GK <- (wk.1*exp(Bs.CR%*%mat_h0s))
    }

    ###hazard function
    if(length(Z_CR)==0){
      pred_surv <- 0
    }
    else{
      pred_surv <- (alpha.CR%*%Z_CR)[1,1]
    }
    h <- h*exp(pred_surv)
    h <- as.vector(h)
    h_0.GK <- as.vector(h_0.GK)
    h <- h_0.GK*h
  }

  #########################################
  ##### Computing LAMBDA1 and LAMBDA2 #####
  #########################################
  Gamma1 <- c()
  Gamma2 <- c()
  for(t2 in st.1){

    data.GaussKronrod.2 <- data.GaussKronrod(newdata.id,t2,k = FlexVarJM$control$nb_pointsGK)
    P.2 <- data.GaussKronrod.2$P
    st.2 <- data.GaussKronrod.2$st
    wk.2 <- data.GaussKronrod.2$wk
    data.id.2 <- data.GaussKronrod.2$data.id2

    if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
      list.data.GK.current.2 <- data.time(data.id.2, c(t(st.2)),
                                        FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
      Xs.2 <- list.data.GK.current.2$Xtime
      Us.2 <- list.data.GK.current.2$Utime
      if(FlexVarJM$control$left_trunc){
        list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xs.0 <- list.data.GK.current.0$Xtime
        Us.0 <- list.data.GK.current.0$Utime
      }
    }
    if(FlexVarJM$control$sharedtype %in% c("CVS","S")){
      list.data.GK.slope.2 <- data.time(data.id.2, c(t(st.2)),
                                      FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
      Xs.slope.2 <- list.data.GK.slope.2$Xtime
      Us.slope.2 <- list.data.GK.slope.2$Utime
      if(FlexVarJM$control$left_trunc){
        list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                          FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        Xs.slope.0 <- list.data.GK.slope.0$Xtime
        Us.slope.0 <- list.data.GK.slope.0$Utime
      }
    }

    #### lambda0
    if(FlexVarJM$control$hazard_baseline == "Exponential"){
      mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
      Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
    }
    else{
      if(FlexVarJM$control$hazard_baseline == "Weibull"){
        mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
        Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
      }
      else{
        if(FlexVarJM$control$hazard_baseline == "Splines"){
          mfZ <- model.frame(FlexVarJM$control$formSurv, data = newdata.id)
          Z <- model.matrix(FlexVarJM$control$formSurv, mfZ)
          Z <- Z[,-1]
          Bs.2 <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines, c(t(st.2)), ord = 4L)
          if(FlexVarJM$control$left_trunc){
            Bs.0 <- splineDesign(rr, c(t(st.0)), ord = 4L)
          }
        }
        else{
          stop("This type of base survival function is not implemented.")
        }
      }

    }

    ####Same for competing risks
    nb.alpha.CR <- 0
    if(FlexVarJM$control$competing_risk){
      if(FlexVarJM$control$sharedtype_CR %in% c("CV","CVS")){
        list.data.GK.current.2 <- data.time(data.id.2, c(t(st.2)),
                                          FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
        Xs.2 <- list.data.GK.current.2$Xtime
        Us.2 <- list.data.GK.current.2$Utime
        if(FlexVarJM$control$left_trunc){
          list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                              FlexVarJM$control$formFixed, FlexVarJM$control$formRandom,FlexVarJM$control$timeVar)
          Xs.0 <- list.data.GK.current.0$Xtime
          Us.0 <- list.data.GK.current.0$Utime
        }
      }
      if(FlexVarJM$control$sharedtype_CR %in% c("CVS","S")){
        list.data.GK.slope.2 <- data.time(data.id.2, c(t(st.2)),
                                        FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
        Xs.slope.2 <- list.data.GK.slope.2$Xtime
        Us.slope.2 <- list.data.GK.slope.2$Utime
        if(FlexVarJM$control$left_trunc){
          list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            FlexVarJM$control$formSlopeFixed, FlexVarJM$control$formSlopeRandom,FlexVarJM$control$timeVar)
          Xs.slope.0 <- list.data.GK.slope.0$Xtime
          Us.slope.0 <- list.data.GK.slope.0$Utime
        }

      }

      if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
        mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
        Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
      }
      else{
        if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
          mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
          Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
        }
        else{
          if(FlexVarJM$control$hazard_baseline_CR == "Splines"){
            mfZ.CR <- model.frame(FlexVarJM$control$formSurv_CR, data = newdata.id)
            Z_CR <- model.matrix(FlexVarJM$control$formSurv_CR, mfZ.CR)
            Z_CR <- Z_CR[,-1]
            Bs.CR.2 <- splineDesign(FlexVarJM$control$knots.hazard_baseline.splines.CR, c(t(st.2)), ord = 4L)
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


    h.2.1 <- 1
    h.2.2 <- 1
      if(FlexVarJM$control$variability_hetero){
        h.2.1 <- h.2.1*exp(alpha.sigma*sigmai)
        if(FlexVarJM$control$competing_risk){
          h.2.2 <- h.2.2*exp(alpha.sigma.CR*sigmai)
        }
      }
      if(FlexVarJM$control$sharedtype %in% c("CV","CVS")|| (FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV","CVS"))){
        current.GK.2 <- beta%*%t(Xs.2) + bi%*%t(Us.2)
        if(FlexVarJM$control$sharedtype %in% c("CV","CVS")){
          h.2.1 <- h.2.1*exp(alpha.current*current.GK.2)
        }
        if((FlexVarJM$control$competing_risk && FlexVarJM$control$sharedtype_CR %in% c("CV","CVS"))){
          h.2.2 <- h.2.2*exp(alpha.current.CR*current.GK.2)
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
    }
    else{
      pred_surv <- (alpha%*%Z)[1,1]
    }
    h.2.1 <- h.2.1*exp(pred_surv)
    h.2.1 <- as.vector(h.2.1)
    h_0.GK.2 <- as.vector(h_0.GK.2)
    h.2.1 <- h_0.GK.2*h.2.1

    if(FlexVarJM$control$competing_risk){
      if(FlexVarJM$control$hazard_baseline_CR == "Exponential"){
        h_0.GK.2.CR <- wk.2
      }

      if(FlexVarJM$control$hazard_baseline_CR == "Weibull"){
        h_0.GK.2.CR <- shape.CR*(st.2**(shape.CR-1))*wk.2
      }

      if(FlexVarJM$control$hazard_baseline == "Splines"){
        mat_h0s <- matrix(gamma,ncol=1)
        h_0.GK.2.CR <- (wk.2*exp(Bs.CR.2%*%mat_h0s))
      }

      ###hazard function
      if(length(Z_CR)==0){
        pred_surv.CR <- 0
      }
      else{
        pred_surv.CR <- (alpha.CR%*%Z)[1,1]
      }
      h.2.2 <- h.2.2*exp(pred_surv.CR)
      h.2.2 <- as.vector(h.2.2)
      h_0.GK.2.CR <- as.vector(h_0.GK.2.CR)
      h.2.2 <- h_0.GK.2.CR*h.2.2
    }
    Gamma1 <- c(Gamma1, (t2/2)*sum(wk.2*h.2.1))
    if(FlexVarJM$control$competing_risk){
      Gamma2 <- c(Gamma2,(t2/2)*sum(wk.2*h.2.2))
    }

  }

  int <- exp(-Gamma1-Gamma2)*h

  res <- (t/2)*sum(wk.1*int)
  res
}
