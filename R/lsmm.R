lsmm <- function(formFixed, formRandom, formGroup, timeVar, data.long,
                 variability_hetero = TRUE,  formFixedVar, formRandomVar, correlated_re = FALSE,S1 = 1000, S2= 5000, nproc = 1, clustertype = "SOCK", maxiter = 100,
                 print.info = FALSE, file = NULL, epsa = 1e-03, epsb = 1e-03, epsd = 1e-03, binit = NULL
                       
){
  time.prog1 <- Sys.time()
  precision = 0.01
  #Check enter parameters
  if(missing(formFixed)) stop("The argument formFixed must be specified")
  if(missing(formRandom)) stop("The argument formRandom must be specified")
  if(!inherits(class(formFixed),"formula")) stop("The argument formFixed must be a formula")
  if(!inherits(class(formRandom),"formula")) stop("The argument formRandom must be a formula")
  if(missing(formGroup)) stop("The argument formGroup must be specified")
  if(!inherits(class(formGroup),"formula")) stop("The argument formGroup must be a formula")
  if(missing(timeVar)) stop("The argument timeVar must be specified")
  if(!inherits(class(timeVar),"character")) stop("The argument timeVar must be a character")
  if(length(timeVar) != 1) stop("The argument timeVar must be of length 1")
  if(missing(data.long)) stop("The argument data.long must be specified")
  if(!inherits(class(data.long) ,"data.frame") )stop("The argument data.long must be a data frame")
  if(nrow(data.long) == 0) stop("Data should not be empty")
  if(!(timeVar %in% colnames(data.long))) stop("Unable to find variable 'timeVar' in 'data.long'")
  if(!inherits(class(variability_hetero) ,"logical") )stop("The argument 'varability_hetero' must be a logical")
  if(!inherits(class(precision) , "numeric") )stop("The argument precision must be a numeric")
  if(!inherits(class(S1),"numeric")) stop("The argument S1 must be a numeric")
  if(!inherits(class(S2),"numeric")) stop("The argument S2 must be a numeric")
  #if(!(all.vars(formFixed) %in% colnames(data.long))) stop("All variables used in the argument formFixed must be in data.long")
  #if(!(all.vars(formRandom) %in% colnames(data.long))) stop("All variables used in the argument formRandom must be in data.long")
  #if(!(all.vars(formGroup) %in% colnames(data.long))) stop("All variables used in the argument formGroup must be in data.long")
  #if(!(all.vars(formSurv) %in% colnames(data.long))) stop("All variables used in the argument formSurv must be in data.long")
  #if(competing_risk && !(all.vars(formSurv_CR) %in% colnames(data.long))) stop("All variables used in the argument formSurv_CR must be in data.long")
  
  
  
  
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
  rr <- NULL
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
  #data management
  id <- as.integer(data.long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data.long))) #To have a column named "id"
    data.long <- cbind(data.long, id = id)
  else{
    data.long$id <- as.integer(data.long$id)
  }
  idVar = "id"
  
  ##longitudinal part
  cat("Longitudinal initialisation \n")
  list.long <- data.manag.long(formGroup,formFixed, formRandom,data.long)
  X_base <- list.long$X
  U <- list.long$U
  nb.e.a <- ncol(U)
  y.new.prog <- list.long$y.new.prog
  list.var <- data.manag.sigma(formGroup,formFixedVar, formRandomVar,data.long)
  O_base <- list.var$X
  W_base <- list.var$U
  nb.omega <- ncol(O_base)
  nb.e.a.sigma <- ncol(W_base)
  data.long <- cbind(data.long,y.new.prog)
  data.long <- data.long %>% group_by(id) %>% dplyr::mutate(sd.emp = sd(y.new.prog),
                                                            VC.emp = mean(y.new.prog) )%>% ungroup()
  data.long <- as.data.frame(data.long)
  offset <- list.long$offset
  Ind <- list.long$I
  if(is.null(binit)){
    list.init.long <- initial.long(formFixed, formRandom, idVar, data.long,
                                   ncol(list.long$X), nproc = nproc)
    sigma_epsilon <- list.init.long$sigma
    mu.log.sigma <- log(sigma_epsilon)
    tau.log.sigma <- precision
    cholesky_b <- list.init.long$long_model$cholesky
    priorMean.beta <- list.init.long$priorMean.beta
  }
  if(is.null(binit)){
    # Marqueur :
    ## Effets fixes trend :
    binit <- c(binit, priorMean.beta)
    names_param <- c(names_param, paste(colnames(X_base),"Y",sep = "_"))
    ## Effets fixes var :
    if(variability_hetero){
      binit <- c(binit,mu.log.sigma,rep(0,nb.omega-1))
      names_param <- c(names_param, paste(colnames(O_base),"Var",sep = "_"))
    }
    else{
      binit <- c(binit, sigma_epsilon)
      names_param <- c(names_param, sigma)
    }
    ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
    if(variability_hetero){
      if(correlated_re){
        binit <- c(binit,
                   cholesky_b,
                   rep(0, nb.e.a*nb.e.a.sigma),
                   rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))
        for(i in 1:length(c(cholesky_b,
                            rep(0, nb.e.a*nb.e.a.sigma),
                            rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma)))){
          names_param <- c(names_param, paste("chol", i, sep = "_"))
        }
      }
      else{
        binit <- c(binit,
                   cholesky_b,
                   rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))
        for(i in 1:length(cholesky_b)){
          names_param <- c(names_param, paste("chol_b", i, sep = "_"))
        }
        for(i in 1:length(rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))){
          names_param <- c(names_param, paste("chol_tau", i, sep = "_"))
        }
      }
      
    }
    else{
      binit <- c(binit,
                 cholesky_b)
      for(i in 1:length(cholesky_b)){
        names_param <- c(names_param, paste("chol_b", i, sep = "_"))
      }
    }
    nb.priorMean.beta = length(priorMean.beta)
  }
  if(variability_hetero){
    Zq <- sobol(S1, nb.e.a+nb.e.a.sigma, normal = TRUE, scrambling = 1)
  }
  else{
    Zq <- sobol(S1, nb.e.a, normal = TRUE, scrambling = 1)
  }
  cat("First estimation  \n")
  estimation <- marqLevAlg(binit, fn = log_llh_lsmm, minimize = FALSE,
                           nb.e.a = nb.e.a, nb.priorMean.beta = nb.priorMean.beta,variability_hetero = variability_hetero, S = S1,Zq = Zq, X_base = X_base, offset = offset, 
                           U = U, y.new.prog = y.new.prog, Ind = Ind,
                           nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega, 
                           O_base = O_base, W_base=W_base,
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = TRUE, epsa = epsa, epsb = epsb, epsd = epsd)
  #S2 <- 10000
  if(variability_hetero){
    Zq <- sobol(S2, nb.e.a+1, normal = TRUE, scrambling = 1)
  }
  else{
    Zq <- sobol(S2, nb.e.a, normal = TRUE, scrambling = 1)
  }
  cat("Second estimation  \n")
  estimation2 <- marqLevAlg(estimation$b, fn = log_llh_lsmm, minimize = FALSE,
                            nb.e.a = nb.e.a, nb.priorMean.beta = nb.priorMean.beta,variability_hetero = variability_hetero, S = S1,Zq = Zq, X_base = X_base, offset = offset, 
                            U = U, y.new.prog = y.new.prog, Ind = Ind,
                            nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega, 
                            O_base = O_base, W_base=W_base,
                            nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info, file = file,
                            blinding = TRUE, epsa = epsa, epsb = epsb, epsd = epsd)
  var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
  var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
  sd.param <- sqrt(diag(var_trans))
  param_est <-  estimation2$b
  table.res <- cbind(param_est, sd.param)
  table.res <- as.data.frame(table.res)
  colnames(table.res) <- c("Estimation", "SE")
  rownames(table.res) <- names_param
  
  time.prog2 <- Sys.time()
  time.prog.fin <- difftime(time.prog2, time.prog1)
  
  final_object <- list( result = estimation2,
                        table.res = table.res,
                        time.compute = time.prog.fin,
                        control = list( formFixed = formFixed,
                                        formRandom = formRandom,
                                        formGroup = formGroup,
                                        timeVar = timeVar,
                                        variability_hetero = variability_hetero,
                                        data.long = data.long,
                                        S1 = S1, S2 = S2, nb.e.a = nb.e.a,
                                          nproc = nproc, clustertype = clustertype,
                                        maxiter = maxiter, print.info = print.info,
                                        file = file, epsa = epsa, epsb = epsb, epsd = epsd,
                                        nb.priorMean.beta = nb.priorMean.beta, 
                                        Ind = Ind, conv = estimation2$istop, niter = estimation2$ni,
                                        convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm),
                                        names_long = colnames(X_base),
                                        likelihood_value = estimation2$fn.value)
                        
  )
  class(final_object) <- c("lsmm")
  return(final_object)
  
}
