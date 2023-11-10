#' lsmm : Estimation of location scale mixed model
#' 
#' This function fits complex mixed effects model with a time and covariate dependent variance.
#' We suppose that the variance of the residual error is time-dependent and subject-specific.
#' Parameters are estimated simultaneously through a maximum likelihood method, using a Marquardt-Levenberg algorithm.
#'
#' @details
#' 
#' The model is defined by :
#' #' \eqn{\quad\left\{\begin{array}{ll}
#' Y_{ij} = Y_{i}(t_{ij}) = \widetilde{Y}_i(t_{ij}) + \epsilon_{ij} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij}, \\
#' \epsilon_{ij}(t_{ij}) \sim \mathcal{N}(0,\sigma_i^2(t_{ij})) \hspace{3mm} \text{with} \hspace{3mm} \log(\sigma_i(t_{ij}))  = O_{ij}^{\top} \mu+M_{ij}^{\top} \tau_{i}
#' \end{array}
#' \right.}
#' 
#' with \eqn{X_{ij}}, \eqn{O_{ij}}, \eqn{Z_{ij}} and \eqn{M_{ij}} four vectors of explanatory variables for subject \eqn{i} at visit \eqn{j}, 
#' respectively associated with the fixed-effect vectors \eqn{\beta} and \eqn{\mu}, and the subject-specific random-effect vector \eqn{b_i} and \eqn{\tau_i}, such as 
#' \eqn{\quad\left(\begin{array}{c}
#'              b_{i} \\
#'              \tau_i
#'              \end{array}\right) \sim N\left(\left(\begin{array}{c}
#'                                                   0 \\
#'                                                   0
#'                                                   \end{array}\right),\left(\begin{array}{cc}
#'                                                                            \Sigma_{b} & \Sigma_{\tau b} \\
#'                                                                            \Sigma_{\tau b}' & \Sigma_{\tau}
#' \end{array}\right)\right)}
#' -------------------------------------------------------------------------------------------------------------------------------------------------------------
#' \eqn{Y_{i}(t_{ij}) = \tilde{Y}_i(t_{ij}) + \epsilon_{ij} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij}}
#'
#' with \eqn{X_{ij}} and \eqn{Z_{ij}} two covariate vectors for subject i at visit j,
#' respectively associated with the vector of fixed effects \eqn{\beta} and the vector of
#' subject-specific individual random effects \eqn{b_i}.
#' The vector \eqn{b_i} is assumed to be normally distributed and a specific-subject random effect on the
#' variance of the measure error can be added: \eqn{\epsilon_{ij} \sim \mathcal{N}(0,\sigma_i^2)} and
#'
#' \eqn{\quad\left(\begin{array}{c}
#' b_{i} \\
#' \log \sigma_{i}
#' \end{array}\right) \sim \mathcal{N}\left(\left(\begin{array}{c}
#'                                                0 \\
#'                                                \mu_{\sigma}
#'                                                \end{array}\right),\left(\begin{array}{cc}
#'                                                                         \Sigma_{b} & 0 \\
#'                                                                         0 & \tau_{\sigma}^{2}
#'                                                                         \end{array}\right)\right)}
#'
#'
#'
#'
#' @param formFixed A formula for the fixed effects of the longitudinal submodel
#' @param formRandom A formula for the random effects of the longitudinal submodel
#' @param formGroup A formula which indicates the group variable
#' @param timeVar The name of the column of time in data.long. This variable must appears in data.long
#' @param data.long A dataframe with the longitudinal data
#' @param variability_hetero A logical to indicate if we suppose a subject_specific variability
#' @param formFixedVar A formula for the fixed effects of the variance predictor
#' @param formRandomVar A formula for the random effects of the variance predictor
#' @param correlated_re A logical to indicate if the random effects of the marker and the variance predictors are correlated (By default there are supposed to be independent) 
#' @param S1 An integer : the number of QMC draws for the first step
#' @param S2 An integer : the number of QMC draws for the second step
#' @param nproc An integer : the number of processors for parallel computing
#' @param clustertype one of the supported types from \code{makeCluster} function
#' @param maxiter optional maximum number of iterations for the marqLevAlg iterative algorithm.
#' @param print.info logical indicating if the outputs of each iteration should be written
#' @param file optional character giving the name of the file where the outputs of each iteration should be written (if print.info=TRUE)
#' @param epsa optional threshold for the convergence criterion based on the parameter stability.
#' @param epsb optional threshold for the convergence criterion based on the objective function stability.
#' @param epsd optional threshold for the relative distance to maximum. This criterion has the nice interpretation of estimating the ratio of the approximation error over the statistical error, thus it can be used for stopping the iterative process whatever the problem.
#' @param binit optional initials parameters.
#'
#' @return A FlexVarJoint object which contains the following elements :
#' \describe{
#' \item{\code{result}}{A marqLevAlg object with the results of the estimation.}
#' \item{\code{table.res}}{The table of results : Estimation and SE}
#' \item{\code{time.compute}}{Computation time}
#' \item{\code{control}}{A list of control elements}
#'
#' }
#' @import dplyr
#' @import marqLevAlg
#' @import splines
#' @importFrom randtoolbox sobol
#' @export
#'
#' @examples
#'
#' \donttest{
#'
#' 
#' #fit a joint model with competing risks and subject-specific variability
#' example <- lsmm(formFixed = y~visit,
#' formRandom = ~ visit,
#' formGroup = ~ID,
#' timeVar = "visit",
#' data.long = Data_toy,
#' variability_hetero = TRUE,
#' formFixedVar =~visit,
#' formRandomVar =~visit,
#' correlated_re = TRUE,
#' S1 = 100,
#' S2 = 1000,
#' nproc = 1,
#' maxiter = 100
#' )
#' 
#' summary(example)
#' }
#'
lsmm <- function(formFixed, formRandom, formGroup, timeVar, data.long,
                 variability_hetero = TRUE,  formFixedVar, formRandomVar, correlated_re = FALSE,S1 = 1000, S2= 5000, nproc = 1, clustertype = "SOCK", maxiter = 100,
                 print.info = FALSE, file = NULL, epsa = 1e-03, epsb = 1e-03, epsd = 1e-03, binit = NULL
                       
){
  time.prog1 <- Sys.time()
  precision = 0.01
  #Check enter parameters
  if(missing(formFixed)) stop("The argument formFixed must be specified")
  if(missing(formRandom)) stop("The argument formRandom must be specified")
  if(!inherits((formFixed),"formula")) stop("The argument formFixed must be a formula")
  if(!inherits((formRandom),"formula")) stop("The argument formRandom must be a formula")
  if(missing(formGroup)) stop("The argument formGroup must be specified")
  if(!inherits((formGroup),"formula")) stop("The argument formGroup must be a formula")
  if(missing(timeVar)) stop("The argument timeVar must be specified")
  if(!inherits((timeVar),"character")) stop("The argument timeVar must be a character")
  if(length(timeVar) != 1) stop("The argument timeVar must be of length 1")
  if(missing(data.long)) stop("The argument data.long must be specified")
  if(!inherits((data.long) ,"data.frame") )stop("The argument data.long must be a data frame")
  if(nrow(data.long) == 0) stop("Data should not be empty")
  if(!(timeVar %in% colnames(data.long))) stop("Unable to find variable 'timeVar' in 'data.long'")
  if(!inherits((variability_hetero) ,"logical") )stop("The argument 'varability_hetero' must be a logical")
  if(!inherits((precision) , "numeric") )stop("The argument precision must be a numeric")
  if(!inherits((S1),"numeric")) stop("The argument S1 must be a numeric")
  if(!inherits((S2),"numeric")) stop("The argument S2 must be a numeric")
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
  message("Longitudinal initialisation")
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
    names_param <- c()
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
  message("First estimation")
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
  message("Second estimation")
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
  #rownames(table.res) <- names_param
  
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
