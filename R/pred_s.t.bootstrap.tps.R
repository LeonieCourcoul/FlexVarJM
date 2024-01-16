#' Predictions computation
#'
#' @param newdata data frame : collected data for a new individual
#' @param object lsjm object : estimation of the model
#' @param s numeric : the time to begin prediction
#' @param window numeric : the side of the prediction window
#' @param event integer (0, 1 or 2) : the event of interest for the prediction
#' @param nb.draws integer : the number of draws to compute the IC 
#'

pred_s.t.bootstrap.tps <- function(newdata,object, s, window, event = 1, nb.draws ){
  
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
    #Ncpus <- object$control$nproc
    #cl <- parallel::makeCluster(Ncpus)
    #doParallel::registerDoParallel(cl)
    Hess <- matrix(rep(0,length(object$result_step2$result_step2$grad)**2),nrow=length(object$result_step2$result_step2$grad),ncol=length(object$result_step2$result_step2$grad))
    Hess[upper.tri(Hess, diag=T)] <- object$result_step2$result_step2$v
    Hess2 = Hess + t(Hess)
    diag(Hess2) <- diag(Hess2) - diag(Hess)
    result <- c()
    #res <- foreach::foreach(l=1:nb.draws, .combine='c',.packages=c("survival","splines","FlexVarJM","mvtnorm")) %dopar%{
    for(l in 1:nb.draws){
      param <- mvtnorm::rmvnorm(1, mean = object$result_step2$result_step2$b, sigma = Hess2)
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
          C1 <- matrix(rep(0,(nb.e.a)**2),nrow=object$control$nb.e.a,ncol=object$control$nb.e.a)
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
      
      # Survival data
      ### Between s and s+t
      data.GaussKronrod.1 <- data.GaussKronrod2(data.long.until.time.s.id,a=s,b=s+window,k = object$control$nb_pointsGK)
      P.1 <- data.GaussKronrod.1$P
      st.1 <- data.GaussKronrod.1$st
      wk.1 <- data.GaussKronrod.1$wk
      data.id.1 <- data.GaussKronrod.1$data.id2
      
      ##########Computing little lambda############
      ### Matrix for current value and slope
      if((c("variability") %in% object$control$sharedtype )|| (c("variability") %in% object$control$sharedtype_CR )){
        list.data.GK.current.sigma <- data.time(data.id.1, c(t(st.1)),
                                                object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
        Os <- list.data.GK.current.sigma$Xtime
        Ws <- list.data.GK.current.sigma$Utime
        Sigma.current.GK <- exp(matrix(rep(omega%*%t(Os),object$control$S2),nrow=object$control$S2,byrow = T) + b_om%*%t(Ws))
      }
      if((c("current value") %in% object$control$sharedtype )|| (c("current value") %in% object$control$sharedtype_CR )){
        list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                           object$control$formFixed, object$control$formRandom,object$control$timeVar)
        Xs <- list.data.GK.current$Xtime
        Us <- list.data.GK.current$Utime
        current.GK <- matrix(rep(beta%*%t(Xs),object$control$S2),nrow=object$control$S2,byrow = T) + b_al%*%t(Us)
      }
      if((c("slope") %in% object$control$sharedtype )|| (c("slope") %in% object$control$sharedtype_CR )){
        list.data.GK.slope <-  data.time(data.id.1, c(t(st.1)),
                                         object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
        Xs.slope <- list.data.GK.slope$Xtime
        Us.slope <- list.data.GK.slope$Utime
        slope.GK <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope)
      }
      #### lambda0
      if(object$control$hazard_baseline == "Exponential"){
        mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
        Z <- model.matrix(object$control$formSurv, mfZ)
      }else{
        if(object$control$hazard_baseline == "Weibull"){
          mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
          Z <- model.matrix(object$control$formSurv, mfZ)
        }else{
          if(object$control$hazard_baseline == "Splines"){
            mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
            Z <- model.matrix(object$control$formSurv, mfZ)
            Z <- Z[,-1]
            Bs <- splines::splineDesign(object$control$knots.hazard_baseline.splines, c(t(st.1)), ord = 4L)
           # if(object$control$left_trunc){
           #   Bs.0 <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
           # }
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
          if(object$control$hazard_baseline_CR == "Weibull"){
            mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
            Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
          }else{
            if(object$control$hazard_baseline_CR == "Splines"){
              mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
              Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
              Z_CR <- Z_CR[,-1]
              Bs.CR <- splines::splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st.1)), ord = 4L)
              #if(object$control$left_trunc){
              #  Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
              #}
            }else{
              stop("This type of base survival function is not implemented.")
            }
          }
        }
      }
      
      h <- 1
      etaBaseline <- 0
      survLong <- 0
      etaBaseline.0 <- 0
      survLong.0 <- 0
      if(event==1){
        if(c("variability") %in% object$control$sharedtype){
          h <- h*exp(alpha.sigma*Sigma.current.GK)
        }
        if(c("current value") %in% object$control$sharedtype){          
          h <- h*exp(alpha.current*current.GK)
        }
        
        if(c("slope") %in% object$control$sharedtype){
          h <- h*exp(alpha.slope*slope.GK)
        }
        ###h0
        if(object$control$hazard_baseline == "Exponential"){
          h_0.GK <- wk.1
        }
        
        if(object$control$hazard_baseline == "Weibull"){
          h_0.GK <- shape*(st.1**(shape-1))*wk.1
        }
        
        if(object$control$hazard_baseline == "Splines"){
          mat_h0s <- matrix(gamma,ncol=1)
          h_0.GK <- (wk.1*exp(Bs%*%mat_h0s))
        }
        
        ###hazard function
        if(length(Z)==0){
          pred_surv <- 0
        }else{
          pred_surv <- (alpha%*%Z)[1,1]
        }
      }else{
        if(c("variability") %in% object$control$sharedtype_CR){
          h <- h*exp(alpha.sigma.CR*Sigma.current.GK)
        }
        if(c("current value") %in% object$control$sharedtype_CR ){
          h <- h*exp(alpha.current.CR*current.GK)
        }
        if(c("slope") %in% object$control$sharedtype_CR){
          h <- h*exp(alpha.slope.CR*slope.GK)
        }
        ###h0
        if(object$control$hazard_baseline_CR == "Exponential"){
          h_0.GK <- wk.1
        }
        
        if(object$control$hazard_baseline_CR == "Weibull"){
          h_0.GK <- shape.CR*(st.1**(shape.CR-1))*wk.1
        }
        
        if(object$control$hazard_baseline_CR == "Splines"){
          mat_h0s <- matrix(gamma.CR,ncol=1)
          h_0.GK <- (wk.1*exp(Bs.CR%*%mat_h0s))
        }
        
        ###hazard function
        if(length(Z_CR)==0){
          pred_surv <- 0
        }else{
          pred_surv <- (alpha.CR%*%Z_CR)[1,1]
        }
      }
      
      h <- h*exp(pred_surv)
      h <- matrix(rep(h_0.GK,nrow(h)),nrow = nrow(h),byrow = T)*h
      
      Gamma1 <- c()
      Gamma2 <- c()
      for(t2 in st.1){
        
        data.GaussKronrod.2 <-  data.GaussKronrod(data.long.until.time.s.id,t2,k = object$control$nb_pointsGK)
        P.2 <- data.GaussKronrod.2$P
        st.2 <- data.GaussKronrod.2$st
        wk.2 <- data.GaussKronrod.2$wk
        data.id.2 <- data.GaussKronrod.2$data.id2
        
        if(c("variability") %in% object$control$sharedtype){
          list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                               object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
          Os.2 <- list.data.GK.current.2$Xtime
          Ws.2 <- list.data.GK.current.2$Utime
        }
        
        if(c("current value") %in% object$control$sharedtype){
          list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                               object$control$formFixed, object$control$formRandom,object$control$timeVar)
          Xs.2 <- list.data.GK.current.2$Xtime
          Us.2 <- list.data.GK.current.2$Utime
        }
        if(c("slope") %in% object$control$sharedtype){
          list.data.GK.slope.2 <-  data.time(data.id.2, c(t(st.2)),
                                             object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
          Xs.slope.2 <- list.data.GK.slope.2$Xtime
          Us.slope.2 <- list.data.GK.slope.2$Utime
        }
        
        #### lambda0
        if(object$control$hazard_baseline == "Exponential"){
          mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
          Z <- model.matrix(object$control$formSurv, mfZ)
        }else{
          if(object$control$hazard_baseline == "Weibull"){
            mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
            Z <- model.matrix(object$control$formSurv, mfZ)
          }else{
            if(object$control$hazard_baseline == "Splines"){
              mfZ <- model.frame(object$control$formSurv, data = data.long.until.time.s.id)
              Z <- model.matrix(object$control$formSurv, mfZ)
              Z <- Z[,-1]
              Bs.2 <- splines::splineDesign(object$control$knots.hazard_baseline.splines, c(t(st.2)), ord = 4L)
             # if(object$control$left_trunc){
             #   Bs.0 <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
             # }
            }else{
              stop("This type of base survival function is not implemented.")
            }
          }
        }
        ####Same for competing risks
        if(object$control$competing_risk){
          if(c("variability") %in% object$control$sharedtype_CR){
            list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                                 object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
            Os.2 <- list.data.GK.current.2$Xtime
            Ws.2 <- list.data.GK.current.2$Utime
          }
          if(c("current value") %in% object$control$sharedtype_CR){
            list.data.GK.current.2 <-  data.time(data.id.2, c(t(st.2)),
                                                 object$control$formFixed, object$control$formRandom,object$control$timeVar)
            Xs.2 <- list.data.GK.current.2$Xtime
            Us.2 <- list.data.GK.current.2$Utime
          }
          if(c("slope") %in% object$control$sharedtype_CR){
            list.data.GK.slope.2 <-  data.time(data.id.2, c(t(st.2)),
                                               object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
            Xs.slope.2 <- list.data.GK.slope.2$Xtime
            Us.slope.2 <- list.data.GK.slope.2$Utime
          }
          
          if(object$control$hazard_baseline_CR == "Exponential"){
            mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
            Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
          }else{
            if(object$control$hazard_baseline_CR == "Weibull"){
              mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
              Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
            }else{
              if(object$control$hazard_baseline_CR == "Splines"){
                mfZ.CR <- model.frame(object$control$formSurv_CR, data = data.long.until.time.s.id)
                Z_CR <- model.matrix(object$control$formSurv_CR, mfZ.CR)
                Z_CR <- Z_CR[,-1]
                Bs.CR.2 <- splines::splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st.2)), ord = 4L)
               # if(object$control$left_trunc){
               #   Bs.0.CR <- splines::splineDesign(rr, c(t(st.0)), ord = 4L)
               # }
              }else{
                stop("This type of base survival function is not implemented.")
              }
            }
          }
        }
        
        h.2.1 <- 1
        h.2.2 <- 1
        if(c("variability") %in% object$control$sharedtype){
          sigma.GK.2 <- exp(matrix(rep(omega%*%t(Os.2),object$control$S2),nrow=object$control$S2,byrow = T)+ b_om%*%t(Ws.2))
          h.2.1 <- h.2.1*exp(alpha.sigma*sigma.GK.2)
          if(object$control$competing_risk && (c("variability") %in% object$control$sharedtype_CR)){
            h.2.2 <- h.2.2*exp(alpha.sigma.CR*sigma.GK.2)
          }
        }
        if((c("current value") %in% object$control$sharedtype)|| (c("current value") %in% object$control$sharedtype_CR)){
          current.GK.2 <- matrix(rep(beta%*%t(Xs.2),object$control$S2),nrow=object$control$S2,byrow = T)+ b_al%*%t(Us.2)
          if((c("current value") %in% object$control$sharedtype)){
            h.2.1 <- h.2.1*exp(alpha.current*current.GK.2)
            #h.2.1 <- matrix(rep(h.2.1,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current*current.GK.2)
          }
          if((object$control$competing_risk && (c("current value") %in% object$control$sharedtype_CR))){
            h.2.2 <- h.2.2*exp(alpha.current.CR*current.GK.2)
            #h.2.2 <- matrix(rep(h.2.2,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current.CR*current.GK.2)
          }
        }
        if((c("slope") %in% object$control$sharedtype)||(c("slope") %in% object$control$sharedtype_CR)){
          #current.GK.2 <- matrix(rep(beta%*%t(Xs.2),object$control$S2),nrow=object$control$S2,byrow = T)+ b_al%*%t(Us.2)
          slope.GK.2 <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope.2),object$control$S2),nrow=object$control$S2,byrow = T)+ b_al[,-1]%*%t(Us.slope.2)
          if((c("slope") %in% object$control$sharedtype)){
            #h.2.1 <- h.2.1*exp(alpha.current*current.GK.2)
            #h.2.1 <- matrix(rep(h.2.1,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current*current.GK.2)
            h.2.1 <- h.2.1*exp(alpha.slope*slope.GK.2)
          }
          if((object$control$competing_risk && (c("slope") %in% object$control$sharedtype_CR))){
            #h.2.2 <- h.2.2*exp(alpha.current.CR*current.GK.2)
            #h.2.2 <- matrix(rep(h.2.2,ncol(current.GK.2)),ncol=ncol(current.GK.2))*exp(alpha.current.CR*current.GK.2)
            h.2.2 <- h.2.2*exp(alpha.slope.CR*slope.GK.2)
          }
        }
        
        ###h0
        if(object$control$hazard_baseline == "Exponential"){
          h_0.GK.2 <- wk.2
        }
        if(object$control$hazard_baseline == "Weibull"){
          h_0.GK.2 <- shape*(st.2**(shape-1))*wk.2
        }
        if(object$control$hazard_baseline == "Splines"){
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
        
        if(object$control$competing_risk){
          if(object$control$hazard_baseline_CR == "Exponential"){
            h_0.GK.2.CR <- wk.2
          }
          if(object$control$hazard_baseline_CR == "Weibull"){
            h_0.GK.2.CR <- shape.CR*(st.2**(shape.CR-1))*wk.2
          }
          if(object$control$hazard_baseline_CR == "Splines"){
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
        if(object$control$competing_risk){
          Gamma2 <- cbind(Gamma2,(t2/2)*rowSums(h.2.2))
        }
      }
      
      if(object$control$competing_risk){
        int <- exp(-Gamma1-Gamma2)*h
      }
      else{
        int <- exp(-Gamma1)*h
      }
      
      surv.num <- P.1*rowSums(int)
      
      numerateur <- surv.num*f_Y_b_sigma
      
      numerateur <- mean(numerateur)
      
      ###### Denominateur #######
      ### At s
      if((c("random effect") %in% object$control$sharedtype)){
        stop("Not implemented yet")
      }else{
        list.GaussKronrod <- data.GaussKronrod(data.long.until.time.s.id, s, k = object$control$nb_pointsGK)
        wk.den <- list.GaussKronrod$wk
        st_calc.den <- list.GaussKronrod$st
        P.den <- list.GaussKronrod$P
        id.GK.den <- list.GaussKronrod$id.GK
        if((c("variability") %in% object$control$sharedtype)){
          list.data.current.time <-  data.time(data.long.until.time.s.id,s, object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
          list.data.GK.current <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                             object$control$formFixedVar, object$control$formRandomVar,object$control$timeVar)
          Otime.den <- list.data.current.time$Xtime
          Wtime.den <- list.data.current.time$Utime
          Os.den <- list.data.GK.current$Xtime
          Ws.den <- list.data.GK.current$Utime
        }
        if((c("current value") %in% object$control$sharedtype)){
          list.data.current.time <-  data.time(data.long.until.time.s.id,s, object$control$formFixed, object$control$formRandom,object$control$timeVar)
          list.data.GK.current <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                             object$control$formFixed, object$control$formRandom,object$control$timeVar)
          Xtime.den <- list.data.current.time$Xtime
          Utime.den <- list.data.current.time$Utime
          Xs.den <- list.data.GK.current$Xtime
          Us.den <- list.data.GK.current$Utime
        }
        if((c("slope") %in% object$control$sharedtype)){
          list.data.slope.time <-  data.time(data.long.until.time.s.id, s, object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
          list.data.GK.slope <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                           object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
          Xslope.den <- list.data.slope.time$Xtime
          Uslope.den <- list.data.slope.time$Utime
          Xs.slope.den <- list.data.GK.slope$Xtime
          Us.slope.den <- list.data.GK.slope$Utime
        }
      }
      if(object$control$hazard_baseline == "Splines"){
        Bs.den <- splines::splineDesign(object$control$knots.hazard_baseline.splines, c(t(st_calc.den)), ord = 4L)
      }
      if(object$control$competing_risk){
        if((c("random effect") %in% object$control$sharedtype_CR)){
          stop("Not implemented yet")
        }else{
          list.GaussKronrod <- data.GaussKronrod(data.long.until.time.s.id, s, k = object$control$nb_pointsGK)
          wk.den <- list.GaussKronrod$wk
          st_calc.den <- list.GaussKronrod$st
          P.den <- list.GaussKronrod$P
          id.GK.den <- list.GaussKronrod$id.GK
          if((c("current value") %in% object$control$sharedtype_CR)){
            #list.data.current.time <-  data.time(data.long.until.time.s.id, s, object$control$formFixed, object$control$formRandom,object$control$timeVar)
            list.data.GK.current <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                               object$control$formFixed, object$control$formRandom,object$control$timeVar)
            Xs.den <- list.data.GK.current$Xtime
            Us.den <- list.data.GK.current$Utime
          }
          if((c("slope") %in% object$control$sharedtype_CR)){
            #list.data.slope.time <-  data.time(data.until.time.s.id, list.surv$Time, formSlopeFixed, formSlopeRandom,timeVar)
            list.data.GK.slope <-  data.time(list.GaussKronrod$data.id2, c(t(st_calc.den)),
                                             object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$timeVar)
            #Xslope <- list.data.slope.time$Xtime
            #Uslope <- list.data.slope.time$Utime
            Xs.slope.den <- list.data.GK.slope$Xtime
            Us.slope.den <- list.data.GK.slope$Utime
          }
        }
        if(object$control$hazard_baseline_CR == "Splines"){
          Bs.CR.den <- splines::splineDesign(object$control$knots.hazard_baseline.splines.CR, c(t(st_calc.den)), ord = 4L)
        }
      }
      
      etaBaseline <- 0
      survLong <- 0
      
      if((c("variability") %in% object$control$sharedtype)){
        sigma.GK.den <- exp(matrix(rep(omega%*%t(Os.den),object$control$S2),nrow=object$control$S2,byrow = T)+ b_om%*%t(Ws.den))
        survLong <- survLong + alpha.sigma*sigma.GK.den
      }
      if(object$control$competing_risk){
        etaBaseline_CR <- 0
        survLong_CR <- 0
        if((c("variability") %in% object$control$sharedtype_CR)){
          survLong_CR <- survLong_CR + alpha.sigma.CR*sigma.GK.den
        }
      }
      if((c("current value") %in% object$control$sharedtype) || (c("current value") %in% object$control$sharedtype_CR)){
        current.GK <- matrix(rep(beta%*%t(Xs.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_al%*%t(Us.den)
        if((c("current value") %in% object$control$sharedtype)){
          survLong <- survLong + alpha.current*current.GK
        }
        if(object$control$competing_risk && (c("current value") %in% object$control$sharedtype_CR)){
          survLong_CR <- survLong_CR + alpha.current.CR*current.GK
        }
      }
      if((c("slope") %in% object$control$sharedtype)|| (c("slope") %in% object$control$sharedtype_CR)){
        slope.GK <- matrix(rep(beta[object$control$indices_beta_slope]%*%t(Xs.slope.den),object$control$S2),nrow=object$control$S2,byrow = T) + b_al[,-1]%*%t(Us.slope.den)
        if((c("slope") %in% object$control$sharedtype)){
          survLong <- survLong + alpha.slope*slope.GK
        }
        if(object$control$competing_risk && (c("slope") %in% object$control$sharedtype_CR)){
          survLong_CR <- survLong_CR + alpha.slope.CR*slope.GK
        }
      }
      if(object$control$hazard_baseline == "Exponential"){
        h_0 <- 1
        h_0.GK <- wk.den
      }
      if(object$control$hazard_baseline == "Weibull"){
        h_0 <- shape*(s**(shape-1))
        h_0.GK <- shape*(st_calc.den**(shape-1))*wk.den
      }
      if(object$control$hazard_baseline == "Splines"){
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
      
      ###GK integration
      survLong <- exp(survLong)
      h_0.GK <- as.vector(h_0.GK)
      survLong <- survLong%*%h_0.GK
      
      Surv <- (-exp(etaBaseline)*P.den*survLong)
      
      
      if(object$control$competing_risk){
        ###h0
        if(object$control$hazard_baseline_CR == "Exponential"){
          h_0.CR <- 1
          h_0.GK.CR <- wk.den
          if(object$control$left_trunc){
            h_0.GK.0_CR <- wk.den
          }
        }
        if(object$control$hazard_baseline_CR == "Weibull"){
          h_0.GK.CR <- shape.CR*(st_calc.den**(shape.CR-1))*wk.den
        }
        if(object$control$hazard_baseline_CR == "Splines"){
          mat_h0s.CR <- matrix(gamma.CR,ncol=1)
          h_0.GK.CR <- (wk.den*exp(Bs.CR.den%*%mat_h0s.CR))
        }
        
        ###hazard function
        if(length(Z_CR)==0){
          pred_surv.CR <- 0
        }else{
          pred_surv.CR <- (alpha.CR%*%Z_CR)[1,1]
        }
        etaBaseline_CR <- etaBaseline_CR + pred_surv.CR
        
        ###GK integration
        survLong_CR <- exp(survLong_CR)
        h_0.GK.CR <- as.vector(h_0.GK.CR)
        survLong_CR <- survLong_CR%*%h_0.GK.CR
        Surv.CR <- (-exp(etaBaseline_CR)*P.den*survLong_CR)
        
      }
      if(object$control$competing_risk){
        denominateur <- exp(Surv+Surv.CR)*f_Y_b_sigma
      }
      else{
        denominateur <- exp(Surv)*f_Y_b_sigma
      }
      
      denominateur <- mean(denominateur)
      #print(denominateur)
      pred.current <- numerateur/denominateur
      result <- c(result, numerateur/denominateur)
  } 
  #parallel::stopCluster(cl)
  }
  result
}

  