#' @export

summary.FlexVarJM <- function(object,...)
{
  x <- object
  if(!inherits(x, "FlexVarJM")) stop("use only \"FlexVarJM\" objects")

  cat("Joint model for quantitative outcome and competing risks", "\n")
  cat("with heterogenous variability and fitted by maximum likelihood method", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Nuumber of subjects:", x$control$Ind),"\n")
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
    cat(paste("     Time of computation :" ,round(x$time.compute,3)))
  }

  cat("\n")
  cat("\n")
  cat("Goodness-of-fit statistics:")
  cat("\n")
  cat(paste("    Likelihood: ", x$control$likelihood_value),"\n")

  cat("\n")
  cat("Maximum Likelihood Estimates:")
  cat("\n")
  cat("Longitudinal model:")
  cat("\n")
  cat("     Cholesky matrix of the random effects:")
  cat("\n")

  chol_mat <- x$table.res$Estimation[grep("^chol", rownames(x$table.res))]
  chol_mat2 <- matrix(0,nrow = quad(1,1,-2*length(chol_mat)), ncol =  quad(1,1,-2*length(chol_mat)) )
  chol_mat2[lower.tri(chol_mat2, diag=T)] <- chol_mat
  print(chol_mat2,quote=FALSE,na.print="")

  cat("\n")
  cat("      Fixed effects:")

  betas_tab <- x$table.res[grep("^beta", rownames(x$table.res)),]
  #print(rownames(betas_tab))
  r.name.betas <- strsplit(rownames(betas_tab), "#§#_")
  r.name.betas2 <- c()
  for(l in r.name.betas){
    r.name.betas2 <- c(r.name.betas2, l[-1])
  }
  print(r.name.betas2)
  rownames(betas_tab) <- r.name.betas2
  betas_tab$Wald <- betas_tab$Estimation/betas_tab$SE
  betas_tab$pvalue <- 1 - pchisq(betas_tab$Wald**2,1)
  colnames(betas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  cat("\n")
  print(betas_tab)

  cat("\n")
  cat("     Variability:")
  var_tab <- x$table.res[grep("log.sigma", rownames(x$table.res)),]
  var_tab$Wald <- var_tab$Estimation/var_tab$SE
  var_tab$pvalue <- 1 - pchisq(var_tab$Wald**2,1)
  colnames(var_tab) <- c("Coeff", "SE", "Wald", "P-value")
  cat("\n")
  print(var_tab)

  cat("\n")

  cat("Survival model(s):")
  cat("\n")
  cat("    First event:")
  e1_surv_tab <- x$table.res[grep(".e1", rownames(x$table.res)),]
  e1_surv_tab$Wald <- e1_surv_tab$Estimation/e1_surv_tab$SE
  e1_surv_tab$pvalue <- 1 - pchisq(e1_surv_tab$Wald**2,1)

  e1_reg <- e1_surv_tab[grep("alpha.e1_", rownames(e1_surv_tab)),]
  e1_reg_other <- rbind(e1_surv_tab[grep("alpha.sigma", rownames(e1_surv_tab)),],
                        e1_surv_tab[grep("alpha.current", rownames(e1_surv_tab)),],
                        e1_surv_tab[grep("alpha.slope", rownames(e1_surv_tab)),])
  r.name.e1_reg <- strsplit(rownames(e1_reg), "#§#_")
  r.name.e1_reg2 <- c()
  for(l in r.name.e1_reg){
    r.name.e1_reg2 <- c(r.name.e1_reg2, l[-1])
  }
  rownames(e1_reg) <- r.name.e1_reg2
  e1_reg <- rbind(e1_reg,e1_reg_other)

  e1_baz <- rbind(e1_surv_tab[grep("weibull.", rownames(e1_surv_tab)),],
                  e1_surv_tab[grep("sp", rownames(e1_surv_tab)),])
  colnames(e1_baz) <- c("Coeff", "SE", "Wald", "P-value")
  colnames(e1_reg) <- c("Coeff", "SE", "Wald", "P-value")

  if(nrow(e1_baz)!=0){
    cat("\n")
    cat("       Baseline:")
    print(e1_baz)
  }
  cat("\n")
  cat("        Regression:")
  cat("\n")
  print(e1_reg)

  cat("\n")

  if(x$control$competing_risk){
    cat("    Second event:")
    e2_surv_tab <- x$table.res[grep(".e2", rownames(x$table.res)),]
    e2_surv_tab$Wald <- e2_surv_tab$Estimation/e2_surv_tab$SE
    e2_surv_tab$pvalue <- 1 - pchisq(e2_surv_tab$Wald**2,1)

    e2_reg <- e2_surv_tab[grep("alpha.e2_", rownames(e2_surv_tab)),]
    e2_reg_other <- rbind(e2_surv_tab[grep("alpha.sigma", rownames(e2_surv_tab)),],
                          e2_surv_tab[grep("alpha.current", rownames(e2_surv_tab)),],
                          e2_surv_tab[grep("alpha.slope", rownames(e2_surv_tab)),])

    r.name.e2_reg <- strsplit(rownames(e2_reg), "_")
    r.name.e2_reg2 <- c()
    for(l in r.name.e2_reg){
      r.name.e2_reg2 <- c(r.name.e2_reg2, l[-1])
    }
    rownames(e2_reg) <- r.name.e2_reg2
    e2_reg <- rbind(e2_reg,e2_reg_other)

    e2_baz <- rbind(e2_surv_tab[grep("weibull.", rownames(e2_surv_tab)),],
                    e2_surv_tab[grep("sp", rownames(e2_surv_tab)),])
    colnames(e2_baz) <- c("Coeff", "SE", "Wald", "P-value")
    colnames(e2_reg) <- c("Coeff", "SE", "Wald", "P-value")

    if(nrow(e2_baz)!=0){
      cat("\n")
      cat("       Baseline:")
      print(e2_baz)
    }
    cat("\n")
    cat("        Regression:")
    cat("\n")
    print(e2_reg)

    cat("\n")

  }



}
