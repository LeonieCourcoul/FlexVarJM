#' @export

summary.lsjm <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsjm")) stop("use only \"lsjm\" objects")

  cat("Joint model for quantitative outcome and competing risks", "\n")
  cat("with heterogenous variability and fitted by maximum likelihood method", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Number of subjects:", x$control$Ind),"\n")
  cat(paste("    Number of observations:", nrow(x$control$data.long)),"\n")
  #cat(paste("    Number of events 1:", ),"\n")
  #cat(paste("    Number of events 2:", ),"\n")

  cat("\n")
  cat("Iteration process:", "\n")

  if(x$control$conv==1) cat("    Convergence criteria satisfied")
  if(x$control$conv==2) cat("    Maximum number of iteration reached without convergence")
  if(x$control$conv==4) {cat("    The program stopped abnormally. No results can be displayed. \n")
  }
  else{
    cat("\n")
    cat(paste("     Number of iterations: ",x$control$niter), "\n")
    cat(paste("     Convergence criteria: parameters =" ,signif(x$control$convcrit[1],3)), "\n")
    cat(paste("                         : likelihood =" ,signif(x$control$convcrit[2],3)), "\n")
    cat(paste("                         : second derivatives =" ,signif(x$control$convcrit[3],3)), "\n")
    cat(paste("     Time of computation :" ,format(x$time.compute)))
  }

  cat("\n")
  cat("\n")
  cat("Goodness-of-fit statistics:")
  cat("\n")
  cat(paste("    Likelihood: ", x$control$likelihood_value),"\n")
  cat(paste("    AIC: ", 2*nrow(x$table.res) - 2* x$control$likelihood_value),"\n")

  cat("\n")
  cat("Maximum Likelihood Estimates:")
  
  #Manage parameters
  curseur <- 1
  param <- x$table.res$Estimation
  param.se <- x$table.res$SE
  param.names <- rownames(x$table.res)
  #browser()
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
    shape <- param[curseur]
    shape.se <- param.se[curseur]
    shape.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline == "Splines"){
    gamma <- param[(curseur):(curseur+x$control$ord.splines+1)]
    gamma.se <- param.se[(curseur):(curseur+x$control$ord.splines+1)]
    gamma.name <- param.names[(curseur):(curseur+x$control$ord.splines+1)]
    curseur <- curseur + x$control$ord.splines + 2
  }
  ## Covariables :
  if(x$control$nb.alpha >=1){
    alpha <- param[(curseur):(curseur+x$control$nb.alpha-1)]
    alpha.se <- param.se[(curseur):(curseur+x$control$nb.alpha-1)]
    alpha.name <- param.names[(curseur):(curseur+x$control$nb.alpha-1)]
    curseur <- curseur+x$control$nb.alpha
  }
  ## Association :
  if(c("random effect") %in% x$control$sharedtype ){
    stop("Not implemented yet")
  }
  if(c("current value") %in% x$control$sharedtype){
    alpha.current <- param[curseur]
    alpha.current.se <- param.se[curseur]
    alpha.current.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(c("current slope") %in% x$control$sharedtype){
    alpha.slope <- param[curseur]
    alpha.slope.se <- param.se[curseur]
    alpha.slope.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(c("variability") %in% x$control$sharedtype){
    alpha.sigma <- param[curseur]
    alpha.sigma.se <- param.se[curseur]
    alpha.sigma.name <- param.names[curseur]
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
      shape.CR <- param[curseur]
      shape.CR.se <- param.se[curseur]
      shape.CR.name <- param.names[curseur]
      curseur <- curseur + 1
    }
    if(x$control$hazard_baseline_CR == "Splines"){
      gamma.CR <- param[(curseur):(curseur+x$control$ord.splines+1)]
      gamma.CR.se <- param.se[(curseur):(curseur+x$control$ord.splines+1)]
      gamma.CR.name <- param.names[(curseur):(curseur+x$control$ord.splines+1)]
      curseur <- curseur + x$control$ord.splines + 2
    }
    ## Covariables :
    if(x$control$nb.alpha.CR >=1){
      alpha.CR <- param[(curseur):(curseur+x$control$nb.alpha.CR-1)]
      alpha.CR.se <- param.se[(curseur):(curseur+x$control$nb.alpha.CR-1)]
      alpha.CR.name <- param.names[(curseur):(curseur+x$control$nb.alpha.CR-1)]
      curseur <- curseur+x$control$nb.alpha.CR
    }
    ## Association :
    if(c("random effect") %in% x$control$sharedtype_CR){
      stop("Not implemented yet")
    }
    if(c("current value") %in% x$control$sharedtype_CR){
      alpha.current.CR <- param[curseur]
      alpha.current.CR.se <- param.se[curseur]
      alpha.current.CR.name <- param.names[curseur]
      curseur <- curseur + 1
    }
    if(c("current slope") %in% x$control$sharedtype_CR){
      alpha.slope.CR <- param[curseur]
      alpha.slope.CR.se <- param.se[curseur]
      alpha.slope.CR.name <- param.names[curseur]
      curseur <- curseur + 1
    }
    if(c("variability") %in% x$control$sharedtype){
      alpha.sigma.CR <- param[curseur]
      alpha.sigma.CR.se <- param.se[curseur]
      alpha.sigma.CR.name <- param.names[curseur]
      curseur <- curseur + 1
    }
  }
  # Marqueur :
  ## Effets fixes trend :
  beta <- param[curseur:(curseur+x$control$nb.priorMean.beta-1)]
  beta.se <- param.se[curseur:(curseur+x$control$nb.priorMean.beta-1)]
  beta.name <- param.names[curseur:(curseur+x$control$nb.priorMean.beta-1)]
  curseur <- curseur+x$control$nb.priorMean.beta
  ## Effets fixes var :
  if(x$control$variability_hetero){
    omega <- param[curseur:(curseur+x$control$nb.omega-1)]
    omega.se <- param.se[curseur:(curseur+x$control$nb.omega-1)]
    omega.name <- param.names[curseur:(curseur+x$control$nb.omega-1)]
    curseur <- curseur + x$control$nb.omega
  }
  else{
    sigma.epsilon <- param[curseur]
    sigma.epsilon.se <- param.se[curseur]
    sigma.epsilon.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
  if(x$control$variability_hetero){
    if(x$control$correlated_re){
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <- matrix(param[(borne1+1):(borne1+x$control$nb.e.a.sigma*x$control$nb.e.a)],nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a, byrow = TRUE)
      borne2 <- borne1+x$control$nb.e.a.sigma*x$control$nb.e.a + 1
      borne3 <- borne2 + choose(n = x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a.sigma - 1
      C3 <- matrix(rep(0,(x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- param[borne2:borne3]
      C4 <- matrix(rep(0,x$control$nb.e.a*x$control$nb.e.a.sigma),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a.sigma)
      MatCov <- rbind(cbind(C1,C4),cbind(C2,C3))
      MatCov <- as.matrix(MatCov)
    }
    else{
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      borne3 <- borne1 + choose(n = x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a.sigma
      C3 <- matrix(rep(0,(x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
      C4 <- matrix(rep(0,x$control$nb.e.a*x$control$nb.e.a.sigma),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a.sigma)
      C2.bis <- matrix(rep(0,x$control$nb.e.a*x$control$nb.e.a.sigma),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a)
      MatCov <- rbind(cbind(C1,C4),cbind(C2.bis,C3))
      MatCov <- as.matrix(MatCov)
    }
  }
  else{
    borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    MatCov <- C1
  }
  
  
  
  cat("\n")
  cat("Longitudinal model:")
  cat("\n")
  
  cat("      Fixed effects:")

 # betas_tab <- x$table.res[grep("^beta", rownames(x$table.res)),]
  #print(rownames(betas_tab))
  #r.name.betas <- strsplit(rownames(betas_tab), "#§#_")
  #r.name.betas2 <- c()
  #for(l in r.name.betas){
  #  r.name.betas2 <- c(r.name.betas2, l[-1])
  #}
  #print(r.name.betas2)
  betas_tab <- matrix(nrow = length(beta), ncol = 4)
  betas_tab[,1] <- beta
  betas_tab[,2] <- beta.se
  betas_tab[,3] <- betas_tab[,1]/betas_tab[,2]
  betas_tab[,4] <- 1 - pchisq(betas_tab[,3]**2,1)
  betas_tab <- as.data.frame(betas_tab)
  rownames(betas_tab) <- beta.name
  colnames(betas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  cat("\n")
  print(betas_tab)

  if(x$control$variability_hetero){
    cat("\n")
    cat("     Fixed effects of the linear predictor associated with variability:")
    var_tab <- matrix(nrow = length(omega), ncol = 4)
    var_tab[,1] <- omega
    var_tab[,2] <- omega.se
    var_tab[,3] <- var_tab[,1]/var_tab[,2]
    var_tab[,4] <- 1 - pchisq(var_tab[,3]**2,1)
    var_tab <- as.data.frame(var_tab)
    rownames(var_tab) <- omega.name
    colnames(var_tab) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_tab)
  }
  else{
    cat("\n")
    var_tab <- matrix(nrow = length(sigma.epsilon), ncol = 4)
    var_tab[,1] <- sigma.epsilon
    var_tab[,2] <- sigma.epsilon.se
    var_tab[,3] <- var_tab[,1]/var_tab[,2]
    var_tab[,4] <- 1 - pchisq(var_tab[,3]**2,1)
    var_tab <- as.data.frame(var_tab)
    rownames(var_tab) <- sigma.epsilon.name
    colnames(var_tab) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_tab)
  }
  

  cat("\n")
  
  cat("     Covariance matrix of the random effects:")
  cat("\n")
  
  #chol_mat <- x$table.res$Estimation[grep("^chol", rownames(x$table.res))]
  #chol_mat2 <- matrix(0,nrow = quad(1,1,-2*length(chol_mat)), ncol =  quad(1,1,-2*length(chol_mat)) )
  #chol_mat2[lower.tri(chol_mat2, diag=T)] <- chol_mat
  #print(chol_mat2,quote=FALSE,na.print="")
  print(MatCov%*%t(MatCov),quote=FALSE,na.print="")
  
  cat("\n")

  cat("Survival model(s):")
  cat("\n")
  cat("    First event:")
  e1_var_tab <- NULL
  e1_share_current_tab <- NULL
  e1_share_slope_tab <- NULL
  e1_alpha_tab <- NULL
  e1_names_tab <- c()
  #browser()
  if(c("variability") %in% x$control$sharedtype){
    e1_var_tab <- matrix(nrow = 1, ncol = 4)
    e1_var_tab[,1] <- alpha.sigma
    e1_var_tab[,2] <- alpha.sigma.se
    e1_var_tab[,3] <- e1_var_tab[,1]/e1_var_tab[,2]
    e1_var_tab[,4] <- 1 - pchisq(e1_var_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.sigma.name)
  }
  if(c("random effect") %in% x$control$sharedtype){
    print("Not implemented yet")
  }
  if(c("current value") %in% x$control$sharedtype){
    e1_share_current_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_current_tab[,1] <- alpha.current
    e1_share_current_tab[,2] <- alpha.current.se
    e1_share_current_tab[,3] <- e1_share_current_tab[,1]/e1_share_current_tab[,2]
    e1_share_current_tab[,4] <- 1 - pchisq(e1_share_current_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.current.name)
  }
  #leonie la plus belle <3 <3 <3
  if(c("slope") %in% x$control$sharedtype){
    e1_share_slope_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_slope_tab[,1] <- c(alpha.slope)
    e1_share_slope_tab[,2] <- c(alpha.slope.se)
    e1_share_slope_tab[,3] <- e1_share_slope_tab[,1]/e1_share_slope_tab[,2]
    e1_share_slope_tab[,4] <- 1 - pchisq(e1_share_slope_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.slope.name)
  }
  if(x$control$nb.alpha >=1){
    e1_alpha_tab <- matrix(nrow = length(alpha), ncol = 4)
    e1_alpha_tab[,1] <- alpha
    e1_alpha_tab[,2] <- alpha.se
    e1_alpha_tab[,3] <- e1_alpha_tab[,1]/e1_alpha_tab[,2]
    e1_alpha_tab[,4] <- 1 - pchisq(e1_alpha_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.name)
  }
  
  
  
  e1_bas_tab <- NULL
  if(x$control$hazard_baseline == "Exponential"){
    e1_bas_tab <- matrix(nrow = 1, ncol = 4)
    e1_bas_tab[1,] <- e1_alpha_tab[1,]
    e1_alpha_tab <- e1_alpha_tab[-1,]
    e1_names_tab <- c(e1_names_tab, alpha.name[-1])
    rownames(e1_bas_tab) <- c("intercept")
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline == "Weibull"){
    e1_bas_tab <- matrix(nrow = 2, ncol = 4)
    e1_bas_tab[1,] <- e1_alpha_tab[1,]
    e1_alpha_tab <- e1_alpha_tab[-1,]
    e1_names_tab <- c(e1_names_tab, alpha.name[-1])
    e1_bas_tab[2,1] <- shape
    e1_bas_tab[2,2] <- shape.se
    e1_bas_tab[2,3] <- e1_bas_tab[2,1]/e1_bas_tab[2,2]
    e1_bas_tab[2,4] <- 1 - pchisq(e1_bas_tab[2,3]**2,1)
    rownames(e1_bas_tab) <- c("intercept",shape.name)
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline == "Splines"){
    e1_bas_tab <- matrix(nrow = length(gamma), ncol = 4)
    e1_bas_tab[,1] <- gamma
    e1_bas_tab[,2] <- gamma.se
    e1_bas_tab[,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    rownames(e1_bas_tab) <- gamma.name
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  e1_surv_tab <- rbind(e1_var_tab, e1_share_current_tab, e1_share_slope_tab, e1_alpha_tab)
  rownames(e1_surv_tab) <- e1_names_tab
  colnames(e1_surv_tab) <- c("Coeff", "SE", "Wald", "P-value")
  
  if(nrow(e1_bas_tab)!=0){
    cat("\n")
    cat("       Regression:")
    cat("\n")
    print(e1_surv_tab)
  }
  cat("\n")
  cat(paste("     Baseline: ",x$control$hazard_baseline), "\n")
  cat("\n")
  print(e1_bas_tab)

  cat("\n")

  if(x$control$competing_risk){
    cat("    Second event:")
    e2_var_tab <- NULL
    e2_share_current_tab <- NULL
    e2_share_slope_tab <- NULL
    e2_alpha_tab <- NULL
    e2_names_tab <- c()
    if(c("variability") %in% x$control$sharedtype_CR){
      e2_var_tab <- matrix(nrow = 1, ncol = 4)
      e2_var_tab[,1] <- alpha.sigma.CR
      e2_var_tab[,2] <- alpha.sigma.CR.se
      e2_var_tab[,3] <- e2_var_tab[,1]/e2_var_tab[,2]
      e2_var_tab[,4] <- 1 - pchisq(e2_var_tab[,3]**2,1)
      e2_names_tab <- c(e2_names_tab, alpha.sigma.CR.name)
    }
    if(c("random effect") %in% x$control$sharedtype_CR){
      print("Not implemented yet")
    }
    if(c("current value") %in% x$control$sharedtype_CR){
      e2_share_current_tab <- matrix(nrow = 1, ncol = 4)
      e2_share_current_tab[,1] <- alpha.current.CR
      e2_share_current_tab[,2] <- alpha.current.CR.se
      e2_share_current_tab[,3] <- e2_share_current_tab[,1]/e2_share_current_tab[,2]
      e2_share_current_tab[,4] <- 1 - pchisq(e2_share_current_tab[,3]**2,1)
      e2_names_tab <- c(e2_names_tab, alpha.current.CR.name)
    }
    if(c("slope") %in% x$control$sharedtype_CR){
      e2_share_slope_tab <- matrix(nrow = 1, ncol = 4)
      e2_share_slope_tab[,1] <- alpha.slope.CR
      e2_share_slope_tab[,2] <- alpha.slope.CR.se
      e2_share_slope_tab[,3] <- e2_share_slope_tab[,1]/e2_share_slope_tab[,2]
      e2_share_slope_tab[,4] <- 1 - pchisq(e2_share_slope_tab[,3]**2,1)
      e2_names_tab <- c(e2_names_tab, alpha.slope.CR.name)
    }
    if(x$control$nb.alpha.CR >=1){
      e2_alpha_tab <- matrix(nrow = length(alpha.CR), ncol = 4)
      e2_alpha_tab[,1] <- alpha.CR
      e2_alpha_tab[,2] <- alpha.CR.se
      e2_alpha_tab[,3] <- e2_alpha_tab[,1]/e2_alpha_tab[,2]
      e2_alpha_tab[,4] <- 1 - pchisq(e2_alpha_tab[,3]**2,1)
    }
    
    
    
    e2_bas_tab <- NULL
    if(x$control$hazard_baseline_CR == "Exponential"){
      e2_bas_tab <- matrix(nrow = 1, ncol = 4)
      e2_bas_tab[1,] <- e2_alpha_tab[1,]
      e2_alpha_tab <- e2_alpha_tab[-1,]
      e2_names_tab <- c(e2_names_tab, alpha.CR.name[-1])
      rownames(e2_bas_tab) <- c("intercept")
      colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
    }
    if(x$control$hazard_baseline_CR == "Weibull"){
      e2_bas_tab <- matrix(nrow = 2, ncol = 4)
      e2_bas_tab[1,] <- e2_alpha_tab[1,]
      e2_alpha_tab <- e2_alpha_tab[-1,]
      e2_names_tab <- c(e2_names_tab, alpha.CR.name[-1])
      e2_bas_tab[2,1] <- shape.CR
      e2_bas_tab[2,2] <- shape.CR.se
      e2_bas_tab[2,3] <- e2_bas_tab[2,1]/e2_bas_tab[2,2]
      e2_bas_tab[2,4] <- 1 - pchisq(e2_bas_tab[2,3]**2,1)
      rownames(e2_bas_tab) <- c("intercept",shape.CR.name)
      colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
    }
    if(x$control$hazard_baseline_CR == "Splines"){
      e2_bas_tab <- matrix(nrow = length(gamma.CR), ncol = 4)
      e2_bas_tab[,1] <- gamma.CR
      e2_bas_tab[,2] <- gamma.CR.se
      e2_bas_tab[,3] <- e2_bas_tab[,1]/e2_bas_tab[,2]
      e2_bas_tab[,4] <- 1 - pchisq(e2_bas_tab[,3]**2,1)
      rownames(e2_bas_tab) <- gamma.CR.name
      colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
    }
    
    e2_surv_tab <- rbind(e2_var_tab, e2_share_current_tab, e2_share_slope_tab, e2_alpha_tab)
    rownames(e2_surv_tab) <- e2_names_tab
    colnames(e2_surv_tab) <- c("Coeff", "SE", "Wald", "P-value")
    
    if(nrow(e2_bas_tab)!=0){
      cat("\n")
      cat("       Regression:")
      cat("\n")
      print(e2_surv_tab)
    }
    cat("\n")
    cat(paste("     Baseline: ",x$control$hazard_baseline_CR), "\n")
    cat("\n")
    print(e2_bas_tab)
    
    cat("\n")
    
    
  }



}
