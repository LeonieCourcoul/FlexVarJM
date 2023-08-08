#' @export

summary.lsmm <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm")) stop("use only \"lsmm\" objects")
  
  cat("Linear mixed-effect model for quantitative outcome", "\n")
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
  
}
