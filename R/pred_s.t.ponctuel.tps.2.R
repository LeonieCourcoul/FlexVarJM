#' Predictions computation
#'
#' @param newdata data frame : collected data for a new individual
#' @param object lsjm object : estimation of the model
#' @param s numeric : the time to begin prediction
#' @param window numeric : the side of the prediction window
#' @param event integer (0, 1 or 2) : the event of interest for the prediction
#'
#'

pred_s.t.ponctuel.tps.2 <- function(newdata,object, s, window, event = 1){
  
  newdata <- as.data.frame(newdata)
  table.predyn.ponct <- c()
  #Ncpus <- 40
  
  #cl <- parallel::makeCluster(Ncpus)
  #doParallel::registerDoParallel(cl)
  #id.pred.to <- unique(newdata[,all.vars(object$control$formGroup)])
  #res.pred <- foreach(id.pred.new=1:length(id.pred.to ), .combine='c', .packages = c("survival")) %dopar% {
  #print(id.pred.new)
  #########
  #.packages plus tard
  
  ##################################################################
  for(id.pred.new in unique(newdata[,all.vars(object$control$formGroup)])){
    #print(id.pred.new)
    newdata.id <- subset(newdata, get(all.vars(object$control$formGroup)) == id.pred.new)
    newdata.id$id <- id.pred.new
    newdata.id <- as.data.frame(newdata.id)
    data.long.until.time.s <-subset(newdata.id, get(object$control$timeVar)<=s)
    name.time.event <- all.vars(object$control$formSurv)[1]
    name.event.event <- all.vars(object$control$formSurv)[2]
    data.long.until.time.s[which(data.long.until.time.s[,name.time.event]>=s),name.event.event] <- 0
    data.long.until.time.s[which(data.long.until.time.s[,name.time.event]>=s),name.time.event] <- max(data.long.until.time.s[,object$control$timeVar])
    
    data.long.until.time.s.id <- data.long.until.time.s[1,]
    
    ################
    ###Parameters###
    ################
    curseur <- 1
    param <- object$result_step2$result_step2$b
    #Manage parameters
    #Evenement 1 :
    ## Risque de base :
    if(object$control$hazard_baseline == "Weibull"){
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
    if(object$control$hazard_baseline == "Splines"){
      gamma <- param[(curseur):(curseur+object$control$ord.splines+1)]
      curseur <- curseur + object$control$ord.splines + 2
    }
    ## Covariables :
    if(object$control$nb.alpha >=1){
      alpha <- param[(curseur):(curseur+object$control$nb.alpha-1)]
      curseur <- curseur+object$control$nb.alpha
    }
    ## Association :
    if("current value" %in% object$control$sharedtype){
      alpha.current <- param[curseur]
      curseur <- curseur + 1
    }
    else{
      alpha.current <- 0
    }
    if("slope" %in% object$control$sharedtype){
      alpha.slope <- param[curseur]
      curseur <- curseur + 1
    }
    else{
      alpha.slope <- 0
    }
    if("variability" %in% object$control$sharedtype){
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
    if(object$control$competing_risk){
      ## Risque de base :
      if(object$control$hazard_baseline_CR == "Weibull"){
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
      if(object$control$hazard_baseline_CR == "Splines"){
        gamma.CR <- param[(curseur):(curseur+object$control$ord.splines+1)]
        curseur <- curseur + object$control$ord.splines + 2
      }
      ## Covariables :
      if(object$control$nb.alpha.CR >=1){
        alpha.CR <- param[(curseur):(curseur+object$control$nb.alpha.CR-1)]
        curseur <- curseur+object$control$nb.alpha.CR
      }
      ## Association :
      if("current value" %in% object$control$sharedtype_CR){
        alpha.current.CR <- param[curseur]
        curseur <- curseur + 1
      }
      else{
        alpha.current.CR <- 0
      }
      if("slope" %in% object$control$sharedtype_CR){
        alpha.slope.CR <- param[curseur]
        curseur <- curseur + 1
      }
      else{
        alpha.slope <- 0
      }
      if("variability" %in% object$control$sharedtype_CR){
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
    beta <- param[curseur:(curseur+object$control$nb.priorMean.beta-1)]
    if( "slope" %in% object$control$sharedtype || "slope" %in% object$control$sharedtype_CR){
      beta_slope <- beta[object$control$indices_beta_slope]
    }
    curseur <- curseur+object$control$nb.priorMean.beta
    ## Effets fixes var :
    if(object$control$variability_hetero){
      omega <- param[curseur:(curseur+object$control$nb.omega-1)]
      curseur <- curseur + object$control$nb.omega
    }
    else{
      sigma.epsilon <- param[curseur]
      curseur <- curseur + 1
    }
    ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
    if(object$control$variability_hetero){
      Zq <- randtoolbox::sobol(object$control$S2, object$control$nb.e.a+object$control$nb.e.a.sigma, normal = TRUE, scrambling = 1)
    }else{
      Zq <- randtoolbox::sobol(object$control$S2, object$control$nb.e.a, normal = TRUE, scrambling = 1)
    }
    if(object$control$variability_hetero){
      if(object$control$correlated_re){
        C1 <- matrix(rep(0,(object$control$nb.e.a+object$control$nb.e.a.sigma)**2),nrow=object$control$nb.e.a+object$control$nb.e.a.sigma,ncol=object$control$nb.e.a+object$control$nb.e.a.sigma)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        MatCov <- C1
        MatCov <- as.matrix(MatCov)
        random.effects <- Zq%*%t(MatCov)
        b_al <- random.effects[,1:object$control$nb.e.a]
        b_al <- matrix(b_al, ncol = object$control$nb.e.a)
        b_om <- random.effects[,(object$control$nb.e.a+1):(object$control$nb.e.a+object$control$nb.e.a.sigma)]
        b_om <- matrix(b_om, ncol = object$control$nb.e.a.sigma)
      }
      else{
        borne1 <- curseur + choose(n = object$control$nb.e.a, k = 2) + object$control$nb.e.a - 1
        C1 <- matrix(rep(0,(object$control$nb.e.a)**2),nrow=object$control$nb.e.a,ncol=object$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        borne3 <- borne1 + choose(n = object$control$nb.e.a.sigma, k = 2) + object$control$nb.e.a.sigma
        C3 <- matrix(rep(0,(object$control$nb.e.a.sigma)**2),nrow=object$control$nb.e.a.sigma,ncol=object$control$nb.e.a.sigma)
        C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
        MatCovb <- as.matrix(C1)
        MatCovSig <- as.matrix(C3)
        b_al <- Zq[,1:object$control$nb.e.a]%*%t(MatCovb)
        b_al <- matrix(b_al, ncol = object$control$nb.e.a)
        b_om <- Zq[,(object$control$nb.e.a+1):(object$control$nb.e.a+object$control$nb.e.a.sigma)]%*%t(MatCovSig)
        b_om <- matrix(b_om, ncol = object$control$nb.e.a.sigma)
      }
    }
    else{
      borne1 <- curseur + choose(n = object$control$nb.e.a, k = 2) + object$control$nb.e.a - 1
      C1 <- matrix(rep(0,(object$control$nb.e.a)**2),nrow=object$control$nb.e.a,ncol=object$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      MatCov <- as.matrix(C1)
      b_al <- Zq%*%t(MatCov)
      b_al <- matrix(b_al, ncol = object$control$nb.e.a)
    }
    
    
    #####################
    # Longitudinal part #
    #####################
    list.long <- data.manag.long(object$control$formGroup,object$control$formFixed, object$control$formRandom,data.long.until.time.s)
    X_base <- list.long$X
    U <- list.long$U
    y.new.prog <- list.long$y.new.prog
    if(object$control$variability_hetero){
      list.var <- data.manag.sigma(object$control$formGroup,object$control$formFixedVar, object$control$formRandomVar,data.long.until.time.s)
      O_base <- list.var$X
      W_base <- list.var$U
    }
    if(is.null(nrow(X_base))){
      if(object$control$variability_hetero){
        sigma.long <- exp((omega%*%O_base)[1,1] + b_om%*%W_base)
      }
      else{
        sigma.long <- sigma.epsilon
      }
      CV <- (beta%*%X_base)[1,1] + b_al%*%U
      f_Y_b_sigma <- dnorm(x=y.new.prog, mean = CV, sd = sigma.long)
    }else{
      f_Y_b_sigma <- rep(1,object$control$S2)
      for(k in 1:nrow(X_base)){
        if(object$control$variability_hetero){
          sigma.long <- exp((omega%*%O_base[k,])[1,1] + b_om%*%W_base[k,])
        }
        else{
          sigma.long <- sigma.epsilon
        }
        CV <- (beta%*%X_base[k,])[1,1] + b_al%*%U[k,]
        f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y.new.prog[k], mean = CV, sd = sigma.long)
      }
    }
    
    #####################
    # Survival part #
    #####################
    
    ######## Table creation
    ### Between s and s+t and between 0 and s
    #browser()
    data.GaussKronrod.1 <- data.GaussKronrod2(data.long.until.time.s.id,a=s,b=s+window,k = object$control$nb_pointsGK)
    P.1 <- data.GaussKronrod.1$P
    st.1 <- data.GaussKronrod.1$st
    wk.1 <- data.GaussKronrod.1$wk
    data.id.1 <- data.GaussKronrod.1$data.id2
    data.GaussKronrod.den <- data.GaussKronrod(data.long.until.time.s.id,Time=s,k = object$control$nb_pointsGK)
    P.den <- data.GaussKronrod.den$P
    st.den <- data.GaussKronrod.den$st
    wk.den <- data.GaussKronrod.den$wk
    data.id.den <- data.GaussKronrod.den$data.id2
    
    ### Matrix for current value and slope
    if((c("variability") %in% object$control$sharedtype )|| (object$control$competing_risk && c("variability") %in% object$control$sharedtype_CR )){
      list.data.GK.current.sigma <- data.time(data.id.1, c(t(st.1)),
                                              object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
      Os <- list.data.GK.current.sigma$Xtime
      Ws <- list.data.GK.current.sigma$Utime
      #Sigma.current.GK <- exp(matrix(rep(omega%*%t(Os),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(Ws))
      
      list.data.GK.current.sigma.den <- data.time(data.id.den, c(t(st.den)),
                                              object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
      Os.den <- list.data.GK.current.sigma.den$Xtime
      Ws.den <- list.data.GK.current.sigma.den$Utime
      #Sigma.current.den <- exp(matrix(rep(omega%*%t(Os.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(Ws.den))
      
      #Sigma.current.0_u <- matrix(rep(beta%*%t(X_0_LT_i),S),nrow=S,byrow = T) + b_y%*%t(U_0_LT_i)
    }
    if((c("current value") %in% object$control$sharedtype )|| (object$control$competing_risk && c("current value") %in% object$control$sharedtype_CR )){
      list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                         object$control$formFixed, object$control$formRandom,object$control$timeVar)
      Xs <- list.data.GK.current$Xtime
      Us <- list.data.GK.current$Utime
      
      list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                         object$control$formFixed, object$control$formRandom,object$control$timeVar)
      Xs.den <- list.data.GK.current.den$Xtime
      Us.den <- list.data.GK.current.den$Utime
    }
    if((c("slope") %in% object$control$sharedtype )|| (object$control$competing_risk && c("slope") %in% object$control$sharedtype_CR )){
      list.data.GK.slope <-  data.time(data.id.1, c(t(st.1)),
                                       object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
      Xs.slope <- list.data.GK.slope$Xtime
      Us.slope <- list.data.GK.slope$Utime
      
      list.data.GK.slope.den <-  data.time(data.id.den, c(t(st.den)),
                                       object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
      Xs.slope.den <- list.data.GK.slope.den$Xtime
      Us.slope.den <- list.data.GK.slope.den$Utime
      
      #slope.GK <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope)
    }
    #### lambda0
    if(object$control$hazard_baseline == "Exponential"){
      mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
      Z <- model.matrix(object$control$formSurv, mfZ)
    }else{
      if(object$control$hazard_baseline == "Weibull" || object$control$hazard_baseline == "Gompertz"){
        mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
        Z <- model.matrix(object$control$formSurv, mfZ)
      }else{
        if(object$control$hazard_baseline == "Splines"){
          mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
          Z <- model.matrix(object$control$formSurv, mfZ)
          Z <- Z[,-1]
          Bs <- splines::splineDesign(object$control$knots.hazard_baseline.splines, c(t(st.1)), ord = 4L)
          Bs.den <- splines::splineDesign(object$control$knots.hazard_baseline.splines, c(t(st.den)), ord = 4L)
        }else{
          stop("This type of base survival function is not implemented.")
        }
      }
    }
    
    ### Same for competing risk
    if(object$control$competing_risk){
      
      if(object$control$hazard_baseline_CR == "Exponential"){
        mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
        Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
      }else{
        if(object$control$hazard_baseline_CR == "Weibull" || object$control$hazard_baseline_CR == "Gompertz"){
          mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
          Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
        }else{
          if(object$control$hazard_baseline_CR == "Splines"){
            mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
            Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
            Z_CR <- Z_CR[,-1]
            Bs.CR <- splines::splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st.1)), ord = 4L)
            Bs.CR.den <- splines::splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st.den)), ord = 4L)
            # if(object$control$left_trunc){
            #   Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
            # }
          }else{
            stop("This type of base survival function is not implemented.")
          }
        }
      }
    }
    
    ### Between 0 and u (double integral)
    st_0_u <- c(); X_0_u <- c(); U_0_u <- c(); Xslope_0_u <- c(); Uslope_0_u <- c(); Bs_0_u <- c(); Bs_CR_0_u <- c(); O_0_u <- c(); W_0_u <- c()
    data.id.integrale <- data.long.until.time.s.id
    for(st.integrale in st.1){
      list.GK_0_st.2 <- data.GaussKronrod(data.id.integrale, Time = st.integrale, k = nb_pointsGK)
      st.2 <- list.GK_0_st.2$st
      st_0_u <- rbind(st_0_u, st.2)
      if(("variability" %in% object$control$sharedtype) || (object$control$competing_risk && "variability" %in% object$control$sharedtype_CR)){
        list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
        O_0_st_u <- list.data.GK_0_u$Xtime; W_0_st_u <- list.data.GK_0_u$Utime
        O_0_u <- rbind(O_0_u,O_0_st_u); W_0_u <- rbind(W_0_u,W_0_st_u)
      }
      if(("current value" %in% object$control$sharedtype) || (object$control$competing_risk && "current value" %in% object$control$sharedtype_CR)){
        list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),object$control$formFixed, object$control$formRandom,object$control$timeVar)
        X_0_st_u <- list.data.GK_0_u$Xtime; U_0_st_u <- list.data.GK_0_u$Utime
        X_0_u <- rbind(X_0_u,X_0_st_u); U_0_u <- rbind(U_0_u,U_0_st_u)
      }
      if(("slope" %in% object$control$sharedtype) || (object$control$competing_risk && "slope" %in% object$control$sharedtype_CR)){
        list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
        Xslope_0_st_u <- list.data.GK_0_u$Xtime; Uslope_0_st_u <- list.data.GK_0_u$Utime
        Xslope_0_u <- rbind(Xslope_0_u,Xslope_0_st_u); Uslope_0_u <- rbind(Uslope_0_u,Uslope_0_st_u)
      }
      if(object$control$hazard_baseline == "Splines"){
        Bs_0_u <- rbind(Bs_0_u,splineDesign(object$control$knots.hazard_baseline.splines, c(t(st.2)), ord = 4L))
      }
      if(object$control$competing_risk && object$control$hazard_baseline == "Splines"){
        Bs_CR_0_u <- rbind(Bs_CR_0_u,splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st.2)), ord = 4L))
      }
    }
    
    ### Computation
    etaBaseline_s_t <- 0; survLong_s_t <- 0; etaBaseline_0_u <- 0; survLong_0_u <- 0;
    etaBaseline_0_s <- 0; survLong_0_s <- 0; etaBaseline_0_u.CR <- 0; survLong_0_u.CR <- 0;
    etaBaseline_0_s.CR <- 0; survLong_0_s.CR <- 0

    if((c("current value") %in% object$control$sharedtype )|| (object$control$competing_risk && c("current value") %in% object$control$sharedtype_CR )){
      current.GK <- matrix(rep(beta%*%t(Xs),object$control$S2),nrow=object$control$S2,byrow = T) + b_al%*%t(Us)
      current.GK.den <- matrix(rep(beta%*%t(Xs.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_al%*%t(Us.den)
      current.GK.0_u <- matrix(rep(beta%*%t(X_0_u),object$control$S2),nrow=object$control$S2,byrow = T) + b_al%*%t(U_0_u)
      survLong_0_s <- survLong_0_s + alpha.current*current.GK.den
      survLong_0_u <- survLong_0_u + alpha.current*current.GK.0_u
      if(event == 1){
        survLong_s_t <- survLong_s_t + alpha.current*current.GK
      }
      else{
        survLong_s_t <- survLong_s_t + alpha.current.CR*current.GK
      }
      if(object$control$competing_risk){
        survLong_0_s.CR <- survLong_0_s.CR + alpha.current.CR*current.GK.den
        survLong_0_u.CR <- survLong_0_u.CR + alpha.current.CR*current.GK.0_u
      }
    }
    if((c("slope") %in% object$control$sharedtype )|| (object$control$competing_risk && c("slope") %in% object$control$sharedtype_CR )){
      slope.GK <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope)
      slope.GK.den <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope.den)
      slope.GK.0_u <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xslope_0_u),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Uslope_0_u)
      survLong_0_s <- survLong_0_s + alpha.slope*slope.GK.den
      survLong_0_u <- survLong_0_u + alpha.slope*slope.GK.0_u
      if(event == 1){
        survLong_s_t <- survLong_s_t + alpha.slope*slope.GK
      }
      else{
        survLong_s_t <- survLong_s_t + alpha.slope.CR*slope.GK
      }
      if(object$control$competing_risk){
        survLong_0_s.CR <- survLong_0_s.CR + alpha.slope.CR*slope.GK.den
        survLong_0_u.CR <- survLong_0_u.CR + alpha.slope.CR*slope.GK.0_u
      }
    }
    if((c("variability") %in% object$control$sharedtype )|| (object$control$competing_risk && c("variability") %in% object$control$sharedtype_CR )){
      variability.GK.s_t <- matrix(rep(omega%*%t(Os),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(Ws)
      variability.GK.den <- matrix(rep(omega%*%t(Os.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(Ws.den)
      variability.GK.0_u <- matrix(rep(omega%*%t(O_0_u),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(W_0_u)
      survLong_0_s <- survLong_0_s + alpha.sigma*variability.GK.den
      survLong_0_u <- survLong_0_u + alpha.sigma*variability.GK.0_u
      if(event == 1){
        survLong_s_t <- survLong_s_t + alpha.sigma*variability.GK.s_t
      }
      else{
        survLong_s_t <- survLong_s_t + alpha.sigma.CR*variability.GK.s_t
      }
      if(object$control$competing_risk){
        survLong_0_s.CR <- survLong_0_s.CR + alpha.sigma.CR*variability.GK.den
        survLong_0_u.CR <- survLong_0_u.CR + alpha.sigma.CR*variability.GK.0_u
      }
    }
    wk <- wk.1
    if(object$control$hazard_baseline == "Exponential"){
      h_0.GK_0_s <- wk.1
      h_0.GK_0_u <- rep(wk.1, length(wk))
      if(event == 1){
        h_0.GK_s_t <- wk
      }
    }
    if(object$control$hazard_baseline == "Weibull"){
      h_0.GK_0_s <- shape*(st.den**(shape-1))*wk
      h_0.GK_0_u <- shape*(st_0_u**(shape-1))*wk
      if(event == 1){
        h_0.GK_s_t <- shape*(st.1**(shape-1))*wk
      }
    }
    else{
      if(object$control$hazard_baseline == "Gompertz"){
        h_0.GK_0_s <- Gompertz.1*exp(Gompertz.2*st.den)*wk
        h_0.GK_0_u <- Gompertz.1*exp(Gompertz.2*st_0_u)*wk
        if(event == 1){
          h_0.GK_s_t <- Gompertz.1*exp(Gompertz.2*st.1)*wk
        }
      }
      else{
        if(object$control$hazard_baseline == "Splines"){
          mat_h0s <- matrix(gamma,ncol=1)
          h_0.GK_0_s <- (wk*exp(Bs.den%*%mat_h0s))
          h_0.GK_0_u <- exp(Bs_0_u%*%mat_h0s)*rep(wk, length(wk))
          if(event == 1){
            h_0.GK_s_t <- (wk*exp(Bs%*%mat_h0s))
          }
        }
      }
    }
    
    if(object$control$competing_risk){
      if(object$control$hazard_baseline_CR == "Exponential"){
        h_0.GK_0_s.CR <- wk
        h_0.GK_0_u.CR <- rep(wk, length(wk))
        if(event == 2){
          h_0.GK_s_t <- wk
        }
      }
      if(object$control$hazard_baseline_CR == "Weibull"){
        h_0.GK_0_s.CR <- shape.CR*(st.den**(shape.CR-1))*wk
        h_0.GK_0_u.CR <- shape.CR*(st_0_u**(shape.CR-1))*wk
        if(event == 2){
          h_0.GK_s_t <- shape.CR*(st.1**(shape.CR-1))*wk
        }
      }
      else{
        if(object$control$hazard_baseline_CR == "Gompertz"){
          h_0.GK_0_s.CR <- Gompertz.1.CR*exp(Gompertz.2.CR*st.den)*wk
          h_0.GK_0_u.CR <- Gompertz.1.CR*exp(Gompertz.2.CR*st_0_u)*wk
          if(event == 2){
            h_0.GK_s_t <- Gompertz.1.CR*exp(Gompertz.2.CR*st.1)*wk
          }
        }
        else{
          if(object$control$hazard_baseline_CR == "Splines"){
            mat_h0s <- matrix(gamma.CR,ncol=1)
            h_0.GK_0_s.CR <- (wk*exp(Bs.CR.den%*%mat_h0s))
            h_0.GK_0_u.CR <- exp(Bs_CR_0_u%*%mat_h0s)*rep(wk, length(wk))
            if(event == 2){
              h_0.GK_s_t <- (wk*exp(Bs.CR%*%mat_h0s))
            }
          }
        }
      }
      
      
    }
    
    if(length(Z)==0){
      pred_surv <- 0
    }
    else{
      pred_surv <- (alpha%*%Z)[1,1]
    }
    
    if(object$control$competing_risk){
      if(length(Z_CR)==0){
        pred_surv_CR <- 0
      }
      else{
        pred_surv_CR <- (alpha.CR%*%Z_CR)[1,1]
      }
    }
    
    if(event == 1){
      etaBaseline_s_t <- etaBaseline_s_t + pred_surv
    }
    else{etaBaseline_s_t <- etaBaseline_s_t + pred_surv_CR}
    
    survLong_s_t <- exp(survLong_s_t)*matrix(rep(t(h_0.GK_s_t), object$control$S2), nrow = object$control$S2, byrow = TRUE)
    h <- exp(etaBaseline_s_t)*survLong_s_t
    
    etaBaseline_0_s <- etaBaseline_0_s + pred_surv
    survLong_0_s <- exp(survLong_0_s)%*%h_0.GK_0_s
    A1_0_s <- exp(etaBaseline_0_s)*P.den*survLong_0_s
    
    etaBaseline_0_u <- etaBaseline_0_u + pred_surv
    survLong_0_u <- exp(survLong_0_u)*rep(c(h_0.GK_0_u),each = object$control$S2)
    ck_15 <- rep(st.1, each = object$control$nb_pointsGK)
    A1_0_u <- exp(etaBaseline_0_u)*survLong_0_u*matrix(rep(ck_15, object$control$S2), nrow = object$control$S2, byrow = T)
    A1_0_u_red <- c()
    for(nb.col in 1:object$control$nb_pointsGK){
      A1_0_u_red <- cbind(A1_0_u_red, rowSums(A1_0_u[,(object$control$nb_pointsGK*(nb.col-1)+1):(object$control$nb_pointsGK*nb.col)]))
    }
    A1_0_u_red <- 0.5*A1_0_u_red
    
    if(object$control$competing_risk){
      etaBaseline_0_s.CR <- etaBaseline_0_s.CR + pred_surv_CR
      survLong_0_s.CR <- exp(survLong_0_s.CR)%*%h_0.GK_0_s.CR
      A1_0_s.CR <- exp(etaBaseline_0_s.CR)*P.den*survLong_0_s.CR
      
      etaBaseline_0_u.CR <- etaBaseline_0_u.CR + pred_surv_CR
      survLong_0_u.CR <- exp(survLong_0_u.CR)*rep(c(h_0.GK_0_u.CR),each = object$control$S2)
      ck_15 <- rep(st.1, each = object$control$nb_pointsGK)
      A1_0_u.CR <- exp(etaBaseline_0_u.CR)*survLong_0_u.CR*rep(ck_15, each = object$control$S2)
      A1_0_u_red.CR <- c()
      for(nb.col in 1:object$control$nb_pointsGK){
        A1_0_u_red.CR <- cbind(A1_0_u_red.CR, rowSums(A1_0_u.CR[,(object$control$nb_pointsGK*(nb.col-1)+1):(object$control$nb_pointsGK*nb.col)]))
      }
      A1_0_u_red.CR <- 0.5*A1_0_u_red.CR
    }
    if(object$control$competing_risk){
      Surv.den <- exp(-A1_0_s-A1_0_s.CR)
      Surv.num <- P.1*rowSums(h*exp(-A1_0_u_red - A1_0_u_red.CR))
    }
    else{
      Surv.den <- exp(-A1_0_s)
      Surv.num <- P.1*rowSums(h*exp(-A1_0_u_red))
    }
    
    numerateur <- Surv.num*f_Y_b_sigma
    denominateur <- Surv.den*f_Y_b_sigma
    pred.current <- mean(numerateur)/mean(denominateur)
    table.predyn.ponct <- rbind(table.predyn.ponct,pred.current)
    print(table.predyn.ponct)
  }
  table.predyn.ponct
}

  