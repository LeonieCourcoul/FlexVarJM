#' FlexVar_JM : Estimation of joint model for longitudinal data with a subject-specific time-dependent variability and time-to-event data.
#'
#' This function fits complex joint models with shared random effects.
#' The longitudinal submodel estimates longitudinal data with a mixed-effects model in which
#' we suppose that the variance of the residual error is time-dependent and subject-specific.
#' The survival submodel handles right-censored and left-truncated time-to-event data and competing risks.
#' The dependence structure between the longitudinal and the survival data can be the random effects from the mixed
#' model or the current value of the marker and/or the slope of the marker. We can also adjust on the current variance of the marker.
#' (See below)
#' Parameters are estimated simultaneously through a maximum likelihood method, using a Marquardt-Levenberg algorithm.
#'
#' @details
#' A. LONGITUDINAL SUBMODEL
#'
#' The longitudinal submodel is defined by a linear mixed effects model with the residual variance which coulb be supposed to be time-dependent and subject-specific :
#' \eqn{\quad\left\{\begin{array}{ll}
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
#' B. SURVIVAL SUBMODEL
#' 
#' The risk function for the event $k \in \{1,2\}$ is defined by:
#' \eqn{\lambda_{ik}(t)=\lambda_{0k}(t) \exp \left(W_{i}^{\top} \gamma_{k}+\alpha_{1k}\tilde{y}_i(t)+\\
#' \alpha_{2k}\tilde{y}'_i(t)+ \alpha_{\sigma k} \sigma_i(t) \right)}
#' 
#' with \eqn{\lambda_{0k}(t)} the baseline risk function, \eqn{W_{i}} a vector of baseline covariates associated with the regression coefficient \eqn{\gamma_k}, 
#' and \eqn{\alpha_{1k}}, \eqn{\alpha_{2k}} and \eqn{\alpha_{\sigma k}} the regression coefficients associated with the current value \eqn{\tilde{y}_i(t)}, 
#' the current slope \eqn{\tilde{y}'_i(t)} and the current variability \eqn{\sigma_i(t)} of the marker, respectively. 
#' Different parametric forms for the baseline risk function can be considered, such as exponential, Weibull, or, for more flexibility, a B-splines base.
#' 
#'
#'
#' @param formFixed A formula for the fixed effects of the longitudinal submodel
#' @param formRandom A formula for the random effects of the longitudinal submodel
#' @param formGroup A formula which indicates the group variable
#' @param formSurv A formula which indicates the variables used in the survival submodel
#' @param timeVar The name of the column of time in data.long. This variable must appears in data.long
#' @param data.long A dataframe with the longitudinal data
#' @param variability_hetero A logical to indicate if we suppose a subject_specific variability
#' @param formFixedVar A formula for the fixed effects of the variance predictor
#' @param formRandomVar A formula for the random effects of the variance predictor
#' @param correlated_re A logical to indicate if the random effects of the marker and the variance predictors are correlated (By default there are supposed to be independent) 
#' @param sharedtype char : dependence structure for survival model : "RE" (random effects) or "CV" (current value) or "CVS" (current value and slope) or "S" (slope)
#' @param hazard_baseline char : baseline hazard function : "Exponential" or "Weibull" or "Splines"
#' @param formSlopeFixed A formula for the fixed effects of the slope of the longitudinal submodel : the derivative of the formFixed
#' @param formSlopeRandom A formula for the random effects of the slope of the longitudinal submodel : the derivative of the formRandom
#' @param indices_beta_slope A vector of index indicating which beta of the formFixed formula is used in the formSlopeFixed formula
#' @param nb_pointsGK the number of points for Gauss-Kronrod approximation : choice between 7 and 15. 15 by default.
#' @param ord.splines  A numeric, the order of splines for the baseline risk function (3 by default)
#' @param competing_risk A logical indicating if the model handles with competing risks
#' @param formSurv_CR In case of competing risk A formula which indicates the variables used in the survival submodel for the second event
#' @param hazard_baseline_CR In case of competing risk : a character for the baseline hazard function of the second event
#' @param sharedtype_CR In case of competing risk ; a character for the dependence structure
#' @param left_trunc A logical indicating if the model handles with left truncated data
#' @param Time.0 In case of left truncation : a vector of entry times
#' @param S1 An integer : the number of QMC draws for the first step
#' @param S2 An integer : the number of QMC draws for the second step
#' @param nproc An integer : the number of processors for parallel computing
#' @param clustertype ...
#' @param maxiter optional maximum number of iterations for the marqLevAlg iterative algorithm.
#' @param print.info logical indicating if the outputs of each iteration should be written
#' @param file optional character giving the name of the file where the outputs of each iteration should be written (if print.info=TRUE)
#' @param epsa optional threshold for the convergence criterion based on the parameter stability.
#' @param epsb optional threshold for the convergence criterion based on the objective function stability.
#' @param epsd optional threshold for the relative distance to maximum. This criterion has the nice interpretation of estimating the ratio of the approximation error over the statistical error, thus it can be used for stopping the iterative process whathever the problem.
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
#' @import survival
#' @import marqLevAlg
#' @importFrom survival Surv
#' @importFrom randtoolbox sobol
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' #fit a joint model with competing risks and subaject-specific variability
#' exemple <- FlexVar_JM(formFixed = y~visit,
#' formRandom = ~ visit,
#' formGroup = ~ID,
#' formSurv = Surv(time, event ==1 ) ~ 1,
#' timeVar = "visit",
#' data.long = Data_toy,
#' variability_hetero = TRUE,
#' formFixedVar =~visit, 
#' formRandomVar =~visit,
#' correlated_re = TRUE,
#' sharedtype = "CV",
#' hazard_baseline = "Weibull",
#' competing_risk = TRUE,
#' formSurv_CR = Surv(time, event ==2 ) ~ 1,
#' hazard_baseline_CR = "Weibull",
#' sharedtype_CR = "CV",
#' S1 = 1000,
#' S2 = 8000,
#' nproc = 5
#' )
#' 
#' summary.FlexVarJM(exemple)
#' 
#' #Predictions for the goodness-of-fit :
#' goodness <- goodness_of_fit(exemple, graph = T)
#'
#' #Predictions for a (new) individual to have event 1 between time
#' s = 3 and time S+t = 4 years given that he did not experience any event
#' before time s, its trajectory of marker until time s ans the set
#' of estimated parameters :
#'
#' newdata <- Data_toy[which(Data_toy$ID == 15),]
#' predictions <- pred_s.t(newdata, exemple, s = 3, window = 1, event = 1, tirage = NULL)
#'
#'
#'
#' }
#'
FlexVar_JM <- function(formFixed, formRandom, formGroup, formSurv, timeVar, data.long,
                       variability_hetero = TRUE, formFixedVar, formRandomVar, correlated_re = FALSE, sharedtype = "CV", hazard_baseline = "Exponential",
                       formSlopeFixed = NULL, formSlopeRandom = NULL, indices_beta_slope = NULL,
                       nb_pointsGK = 15, ord.splines = 3, competing_risk = FALSE, formSurv_CR = NULL,
                       hazard_baseline_CR = "Exponential", sharedtype_CR = "CV", left_trunc = FALSE,
                       Time.0 = NULL, S1 = 1000, S2= 5000, nproc = 1, clustertype = "SOCK", maxiter = 100,
                       print.info = FALSE, file = NULL, epsa = 1e-03, epsb = 1e-03, epsd = 1e-03, binit = NULL
                       
){
  time.prog1 <- Sys.time()
  precision = 0.01
  #Check enter parameters
  if(missing(formFixed)) stop("The argument formFixed must be specified")
  if(missing(formRandom)) stop("The argument formRandom must be specified")
  if(class(formFixed)!="formula") stop("The argument formFixed must be a formula")
  if(class(formRandom)!="formula") stop("The argument formRandom must be a formula")
  if(missing(formGroup)) stop("The argument formGroup must be specified")
  if(class(formGroup)!="formula") stop("The argument formGroup must be a formula")
  if(missing(timeVar)) stop("The argument timeVar must be specified")
  if(class(timeVar) != "character") stop("The argument timeVar must be a character")
  if(length(timeVar) != 1) stop("The argument timeVar must be of length 1")
  if(missing(data.long)) stop("The argument data.long must be specified")
  if(class(data.long) != "data.frame") stop("The argument data.long must be a data frame")
  if(nrow(data.long) == 0) stop("Data should not be empty")
  if(!(timeVar %in% colnames(data.long))) stop("Unable to find variable 'timeVar' in 'data.long'")
  if(length(sharedtype) != 1 || !(sharedtype %in% c("RE", "CV", "CVS", "S"))) stop("The value of argument 'sharedtype' must be of lenght 1 and must be 'RE' or 'CV' or 'CVS' or 'S'")
  if(class(variability_hetero) != "logical") stop("The argument 'varability_hetero' must be a logical")
  if(sharedtype %in% c("CVS", "S") && missing(formSlopeFixed)) stop("The argument formSlopeFixed must be specified when the 'sharedtype' variable has the value CVS or S")
  if(sharedtype %in% c("CVS", "S") && class(formSlopeFixed) != "formula") stop("The argument formSlopeFixed must be a formula")
  if(sharedtype %in% c("CVS", "S") && missing(formSlopeRandom)) stop("The argument formSlopeRandom must be specified when the 'sharedtype' variable has the value CVS or S")
  if(sharedtype %in% c("CVS", "S") && class(formSlopeRandom) != "formula") stop("The argument formSlopeRandom must be a formula")
  if(missing(formSurv)) stop("The argument formSurv must be specified")
  if(class(formSurv)!="formula") stop("The argument formSurv must be a formula")
  if(class(precision) != "numeric") stop("The argument precision must be a numeric")
  if(!(nb_pointsGK %in% c(7,15))) stop("The argument nb_pointsGK must be equal to 7 or 15.")
  if(length(hazard_baseline) != 1 || !(hazard_baseline %in% c("Weibull", "Splines","Exponential"))) stop("The value of argument 'hazard_baseline' must be of lenght 1 and must be 'Exponential' or 'Weibull' or 'Splines'")
  if(class(S1)!="numeric") stop("The argument S1 must be a numeric")
  if(class(S2)!="numeric") stop("The argument S2 must be a numeric")
  #if(missing(nb.e.a)) stop("The argument nb.e.a must be specified : it is the number of random effects + 1 when variability_hetero is TRUE")
  #if(class(nb.e.a)!="numeric") stop("The argument nb.e.a must be a numeric : it is the number of random effects + 1 when variability_hetero is TRUE")
  if(hazard_baseline == "Splines" && class(ord.splines) != "numeric") stop("The argument ord.splines must be a numeric : the order of splines for the baseline hazard function")
  if(class(left_trunc) != "logical") stop("The argument left_trunc (to take into account left truncation/delay entry) must be a logical")
  if(class(competing_risk) != "logical") stop("The argument competing_risk (to take into account two competing events) must be a logical")
  if(competing_risk && missing(formSurv_CR)) stop("The argument formSurv_CR must be specified when the argument competing_risk is TRUE")
  if(competing_risk && class(formSurv_CR)!="formula") stop("The argument formSurv_CR must be a formula when the argument competing_risk is TRUE")
  if(competing_risk && length(sharedtype_CR) != 1 ) stop("The value of argument 'sharedtype' must be of lenght 1 and must be 'RE' or 'CV' or 'CVS' or 'S'")
  if(competing_risk && !(sharedtype %in% c("RE", "CV", "CVS", "S"))) stop("The value of argument 'sharedtype' must be of lenght 1 and must be 'RE' or 'CV' or 'CVS' or 'S'")
  if(competing_risk && length(hazard_baseline_CR) != 1 ) stop("The value of argument 'hazard_baseline_CR' must be of lenght 1 and must be 'Exponential' or 'Weibull' or 'Splines'")
  if(competing_risk && !(hazard_baseline_CR %in% c("Weibull", "Splines","Exponential"))) stop("The value of argument 'hazard_baseline_CR' must be of lenght 1 and must be 'Exponential' or 'Weibull' or 'Splines'")
  if(left_trunc && missing(Time.0)) stop("The argument Time.0 (time of entry into the study) must be specified when left_trunc is TRUE")
  if(left_trunc && class(Time.0) !="numeric") stop("The argument Time.0 (time of entry into the study) must be a numeric when left_trunc is TRUE")
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
  O_base <- NULL
  nb.e.a.sigma <- NULL
  nb.omega <- NULL
  Otime <- NULL
  Wtime <- NULL
  Os <- NULL
  Ws <- NULL
  W_base <- NULL
  
  
  #data management
  id <- as.integer(data.long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data.long))) #To have a column named "id"
    data.long <- cbind(data.long, id = id)
  else{
    data.long$id <- as.integer(data.long$id)
  }
  idVar = "id"
  
  ##longitudinal part
  cat("Longitudinal management \n")
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
                                                            VC.emp = mean(y.new.prog))%>% ungroup()
  data.long <- as.data.frame(data.long)
  offset <- list.long$offset
  Ind <- list.long$I
  if(1==1){
    cat("Longitudinal initialisation  \n")
    list.init.long <- initial.long(formFixed, formRandom, idVar, data.long,
                                   ncol(list.long$X), nproc = nproc)
    sigma_epsilon <- list.init.long$sigma
    mu.log.sigma <- log(sigma_epsilon)
    tau.log.sigma <- precision
    cholesky_b <- list.init.long$long_model$cholesky
    priorMean.beta <- list.init.long$priorMean.beta
  }
  
  
  
  ## survival part
  cat("Survival management  \n")
  list.surv <- data.manag.surv(formGroup, formSurv, data.long, formSurv_CompRisk = formSurv_CR)
  event1 <- list.surv$event1
  event2 <- list.surv$event2
  Time <- list.surv$Time
  
  formSurv_dep <- formSurv
  if(variability_hetero){
    formSurv_dep <- update(formSurv_dep, ~. + sd.emp)
  }
  ##dependance
  data.id <- data.long[!duplicated(id),]
  data.id <- cbind(data.id,event1)
  lag = 0
  if(sharedtype %in% c("random effects")){
    stop("Not implemented yet")
  }
  else{
    list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = nb_pointsGK)
    wk <- list.GaussKronrod$wk
    st_calc <- list.GaussKronrod$st
    P <- list.GaussKronrod$P
    id.GK <- list.GaussKronrod$id.GK
    if(left_trunc){
      list.GaussKronrod.0 <- data.GaussKronrod(data.id, Time.0, k = nb_pointsGK)
      st.0 <- list.GaussKronrod.0$st
      P.0 <- list.GaussKronrod.0$P
    }
    if(sharedtype %in% c("CV","CVS") ){
      formSurv_dep <- update(formSurv_dep, ~. + VC.emp)
      list.data.current.time <- data.time(data.id, list.surv$Time, formFixed, formRandom,timeVar)
      list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        formFixed, formRandom,timeVar)
      Xtime <- list.data.current.time$Xtime
      Utime <- list.data.current.time$Utime
      Xs <- list.data.GK.current$Xtime
      Us <- list.data.GK.current$Utime
      if(left_trunc){
        list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            formFixed, formRandom,timeVar)
        Xs.0 <- list.data.GK.current.0$Xtime
        Us.0 <- list.data.GK.current.0$Utime
      }
    }
    if(sharedtype %in% c("CVS","S")){
      list.data.slope.time <- data.time(data.id, list.surv$Time, formSlopeFixed, formSlopeRandom,timeVar)
      list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                      formSlopeFixed, formSlopeRandom,timeVar)
      Xslope <- list.data.slope.time$Xtime
      Uslope <- list.data.slope.time$Utime
      Xs.slope <- list.data.GK.slope$Xtime
      Us.slope <- list.data.GK.slope$Utime
      if(left_trunc){
        list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                          formSlopeFixed, formSlopeRandom,timeVar)
        Xs.slope.0 <- list.data.GK.slope.0$Xtime
        Us.slope.0 <- list.data.GK.slope.0$Utime
      }
    }
  }
  
  
  if(hazard_baseline == "Exponential"){
    mod_surv <- survreg(formSurv_dep, data = data.id, dist = "exponential")
    #lambda0 <- mod_surv$coefficients[1]
    Z <- list.surv$Z
    alpha <- mod_surv$coefficients[1:ncol(Z)]
    alpha[1] <- -alpha[1]
  }
  else{
    if(hazard_baseline == "Weibull"){
      mod_surv <- survreg(formSurv_dep, data = data.id, dist = "exponential")
      #Z <- list.surv$Z#[,-1]
      ##lambda0 <- mod_surv$coefficients[1]
      #shape_racine <- 1    
      #alpha_weib <- -mod_surv$coefficients[1]
      #if(length(mod_surv$coefficients >=2)){
      #    alpha <- mod_surv$coefficients[2:(ncol(Z)+1)]
      #}
      #else{
      #  alpha <- NULL
      #}
      #alpha[1] <- -alpha[1]     
      shape_racine <- 1
      Z <- list.surv$Z
      alpha <- mod_surv$coefficients[1:ncol(Z)]
      alpha[1] <- -alpha[1]
      
    }
    else{
      if(hazard_baseline == "Splines"){
        Z <- list.surv$Z[,-1]
        pp <- seq(0,1, length.out = ord.splines)
        pp <- tail(head(pp,-1),-1)
        tt1 <- as.data.frame(cbind(Time,event1))
        tt <- tt1$Time[which(tt1$event1 == 1)]
        kn <- quantile(tt, pp, names = FALSE)
        kn <- kn[kn<max(Time)]
        rr <- sort(c(rep(range(Time,0), 4L), kn))
        B <- splineDesign(rr, Time, ord = 4L)
        Bs <- splineDesign(rr, c(t(st_calc)), ord = 4L)
        opt_splines <- optim(rep(0,ncol(B)), fn2,event = event1,W2 = B,P = P,wk = wk,
                             Time = Time,W2s = Bs,id.GK = id.GK, method="BFGS", hessian = T)
        tmp_model <- coxph(formSurv_dep,
                           data = data.id,
                           x = TRUE)
        if(length(tmp_model$coefficients) == 2){
          alpha <- NULL
        }
        else{
          alpha <- tmp_model$coefficients[1:ncol(Z)]
        }
        if(left_trunc){
          Bs.0 <- splineDesign(rr, c(t(st.0)), ord = 4L)
        }
      }
      else{
        stop("This type of base survival function is not implemented.")
      }
    }
    
  }
  nb.alpha.CR <- 0
  if(competing_risk){
    data.id <- cbind(data.id,event2)
    formSurv_dep_CR <- formSurv_CR
    if(variability_hetero){
      formSurv_dep_CR <- update(formSurv_dep_CR, ~. + sd.emp)
    }
    if(sharedtype_CR %in% c("RE")){
      stop("Not implemented yet")
    }
    else{
      list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = nb_pointsGK)
      wk <- list.GaussKronrod$wk
      st_calc <- list.GaussKronrod$st
      P <- list.GaussKronrod$P
      id.GK <- list.GaussKronrod$id.GK
      if(left_trunc){
        list.GaussKronrod.0 <- data.GaussKronrod(data.id, Time.0, k = nb_pointsGK)
        st.0 <- list.GaussKronrod.0$st
        P.0 <- list.GaussKronrod.0$P
      }
      if(sharedtype_CR %in% c("CV","CVS")){
        formSurv_dep_CR <- update(formSurv_dep_CR, ~. + VC.emp)
        list.data.current.time <- data.time(data.id, list.surv$Time, formFixed, formRandom,timeVar)
        list.data.GK.current <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                          formFixed, formRandom,timeVar)
        Xtime <- list.data.current.time$Xtime
        Utime <- list.data.current.time$Utime
        Xs <- list.data.GK.current$Xtime
        Us <- list.data.GK.current$Utime
        if(left_trunc){
          list.data.GK.current.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                              formFixed, formRandom,timeVar)
          Xs.0 <- list.data.GK.current.0$Xtime
          Us.0 <- list.data.GK.current.0$Utime
        }
      }
      if(sharedtype_CR %in% c("CVS","S")){
        list.data.slope.time <- data.time(data.id, list.surv$Time, formSlopeFixed, formSlopeRandom,timeVar)
        list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                        formSlopeFixed, formSlopeRandom,timeVar)
        Xslope <- list.data.slope.time$Xtime
        Uslope <- list.data.slope.time$Utime
        Xs.slope <- list.data.GK.slope$Xtime
        Us.slope <- list.data.GK.slope$Utime
        if(left_trunc){
          list.data.GK.slope.0 <- data.time(list.GaussKronrod.0$data.id2, c(t(list.GaussKronrod.0$st)),
                                            formSlopeFixed, formSlopeRandom,timeVar)
          Xs.slope.0 <- list.data.GK.slope.0$Xtime
          Us.slope.0 <- list.data.GK.slope.0$Utime
        }
      }
    }
    
    if(hazard_baseline_CR == "Exponential"){
      mod_surv <- survreg(formSurv_dep_CR, data = data.id, dist = "exponential")
      #lambda0_CR <- mod_surv$coefficients[1]
      Z_CR <- list.surv$Z_CR
      alpha_CR <- mod_surv$coefficients[1:ncol(Z_CR)]
      alpha_CR[1] <- -alpha_CR[1]
    }
    else{
      if(hazard_baseline_CR == "Weibull"){
        mod_surv <- survreg(formSurv_dep_CR, data = data.id, dist = "exponential")
        #    Z_CR <- list.surv$Z_CR[,-1]
        #    #lambda0_CR <- mod_surv$coefficients[1]
        #    #shape_CR <- 1/mod_surv$scale
        #    shape_CR_racine <- 1
        #    alpha_weib_CR <- -mod_surv$coefficients[1]
        #  if(length(mod_surv$coefficients >=2)){
        #      alpha_CR <- mod_surv$coefficients[2:ncol((Z_CR+1))]
        #  }
        #  else{
        #    alpha_CR <- NULL
        #  }
        shape_CR_racine <- 1
        Z_CR <- list.surv$Z_CR
        alpha_CR <- mod_surv$coefficients[1:ncol(Z_CR)]
        alpha_CR[1] <- -alpha_CR[1]
        #alpha_CR <- mod_surv$coefficients[1:ncol(Z_CR)]
        #alpha_CR[1] <- -alpha_CR[1]
        
      }
      else{
        if(hazard_baseline_CR == "Splines"){
          Z_CR <- list.surv$Z_CR[,-1]
          pp.CR <- seq(0,1, length.out = ord.splines)
          pp.CR <- tail(head(pp.CR,-1),-1)
          tt2.CR <- as.data.frame(cbind(Time,event2))
          tt.CR <- tt2.CR$Time[which(tt2.CR$event2 == 1)]
          kn.CR <- quantile(tt.CR, pp.CR, names = FALSE)
          kn.CR <- kn.CR[kn.CR<max(Time)]
          rr.CR <- sort(c(rep(range(Time,0), 4L), kn.CR))
          B.CR <- splineDesign(rr.CR, Time, ord = 4L)
          Bs.CR <- splineDesign(rr.CR, c(t(st_calc)), ord = 4L)
          opt_splines_CR <- optim(rep(0,ncol(B)), fn2,event = event2,W2 = B.CR,P = P,wk = wk,Time = Time,W2s = Bs.CR,id.GK = id.GK, method="BFGS", hessian = T)
          tmp_model <- coxph(formSurv_CR,
                             data = data.id,
                             x = TRUE)
          if(length(tmp_model$coefficients) == 2){
            alpha_CR <- NULL
          }
          else{
            alpha_CR <- tmp_model$coefficients[1:ncol(Z_CR)]
          }
          if(left_trunc){
            Bs.0.CR <- splineDesign(rr, c(t(st.0)), ord = 4L)
          }
        }
        else{
          stop("This type of base survival function is not implemented.")
        }
      }
      
    }
    nb.alpha.CR <- length(alpha_CR)
  }
  
  if(variability_hetero){
    list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = nb_pointsGK)
    wk <- list.GaussKronrod$wk
    st_calc <- list.GaussKronrod$st
    P <- list.GaussKronrod$P
    id.GK <- list.GaussKronrod$id.GK
    
    list.data.current.sigma.time <- data.time(data.id, list.surv$Time, formFixedVar, formRandomVar,timeVar)
    list.data.GK.current.sigma <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st)),
                                            formFixedVar, formRandomVar,timeVar)
    Otime <- list.data.current.sigma.time$Xtime
    Wtime <- list.data.current.sigma.time$Utime
    Os <- list.data.GK.current.sigma$Xtime
    Ws <- list.data.GK.current.sigma$Utime
  }
  
  #browser()
  
  
  if(is.null(binit)){
    alpha.sigma <- 0
    alpha.current <- 0
    alpha.slope <- 0
    alpha.shared.effects <- rep(0, nb.e.a)
    alpha.sigma.CR <- 0
    alpha.current.CR <- 0
    alpha.slope.CR <- 0
    alpha.shared.effects.CR <- rep(0, nb.e.a)
    
    binit <- c()
    names_param <- c()
    #Evenement 1 :
    ## Risque de base :
    if(hazard_baseline == "Weibull"){
      binit <- c(binit,shape_racine)
      names_param <- c(names_param, "Shape")
    }
    if(hazard_baseline == "Splines"){
      binit <- c(binit,opt_splines$par)
      for(i in 1:length(opt_splines$par)){
        names_param <- c(names_param, paste("splines", i, sep = "_"))
      }
    }
    ## Covariables :
    binit <- c(binit, alpha)
    names_param <- c(names_param, paste(colnames(Z),"E1",sep = "_"))
    ## Association :
    if(sharedtype %in% c("RE")){
      stop("Not implemented yet")
    }
    if(sharedtype %in% c("CV","CVS")){
      binit <- c(binit, alpha.current)
      names_param <- c(names_param, "Current Value")
    }
    if(sharedtype %in%  c("CVS","S")){
      binit <- c(binit, alpha.slope)
      names_param <- c(names_param, "Current Slope")
    }
    if(variability_hetero){
      binit <- c(binit,alpha.sigma)
      names_param <- c(names_param, "Current Variance")
    }
    # Evenement 2
    if(competing_risk){
      ## Risque de base :
      if(hazard_baseline_CR == "Weibull"){
        binit <- c(binit,shape_CR_racine)
        names_param <- c(names_param, "Shape (CR)")
      }
      if(hazard_baseline_CR == "Splines"){
        binit <- c(binit,opt_splines_CR$par)
        for(i in 1:length(opt_splines$par)){
          names_param <- c(names_param, paste("splines (CR)", i, sep = "_"))
        }
      }
      ## Covariables :
      binit <- c(binit, alpha_CR)
      names_param <- c(names_param, paste(colnames(Z_CR),"E2",sep = "_"))
      ## Association :
      if(sharedtype %in% c("RE")){
        stop("Not implemented yet")
      }
      if(sharedtype %in% c("CV","CVS")){
        binit <- c(binit, alpha.current.CR)
        names_param <- c(names_param, "Current Value (CR)")
      }
      if(sharedtype %in%  c("CVS","S")){
        binit <- c(binit, alpha.slope.CR)
        names_param <- c(names_param, "Current Slope (CR)")
      }
      if(variability_hetero){
        binit <- c(binit,alpha.sigma.CR)
        names_param <- c(names_param, "Current Variance (CR)")
      }
    }
    # Marqueur :
    ## Effets fixes trend :
    binit <- c(binit, priorMean.beta)
    names_param <- c(names_param, paste(colnames(X_base),"Y",sep = "_"))
    ## Effets fixes var :
    if(variability_hetero){
      binit <- c(binit,mu.log.sigma,rep(0,nb.omega-1))
      names_param <- c(names_param, paste(colnames(W_base),"Var",sep = "_"))
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
    nb.alpha = length(alpha)
  }
  if(variability_hetero){
    Zq <- randtoolbox::sobol(S1,  nb.e.a+nb.e.a.sigma, normal = TRUE, scrambling = 1)

  }
  else{
    Zq <- randtoolbox::sobol(S1, dim = nb.e.a, mixed = TRUE, normal = TRUE)
  }
  nb.priorMean.beta = length(priorMean.beta)
  nb.alpha = length(alpha)
  cat("First estimation  \n")
  estimation <- marqLevAlg(binit, fn = log_llh, minimize = FALSE,
                           nb.e.a = nb.e.a, nb.priorMean.beta = nb.priorMean.beta,nb.alpha = nb.alpha,
                           competing_risk = competing_risk,
                           nb.alpha.CR = nb.alpha.CR, variability_hetero = variability_hetero, S = S1,Zq = Zq, sharedtype = sharedtype,
                           sharedtype_CR = sharedtype_CR,
                           hazard_baseline = hazard_baseline, hazard_baseline_CR = hazard_baseline_CR, ord.splines = ord.splines, Xtime = Xtime,
                           Utime = Utime, nb_pointsGK = nb_pointsGK,
                           Xs = Xs,Us = Us, Xslope = Xslope, Uslope = Uslope, Xs.slope = Xs.slope, Us.slope = Us.slope,
                           indices_beta_slope = indices_beta_slope, Time =Time,
                           st_calc = st_calc, B = B, Bs = Bs, wk = wk, Z = Z, P = P, left_trunc = left_trunc,
                           Z_CR = Z_CR, X_base = X_base, offset = offset, U = U, y.new.prog = y.new.prog, event1 = event1, event2 = event2, Ind = Ind,
                           Xs.0 = Xs.0, Us.0 = Us.0, Xs.slope.0 = Xs.slope.0, Us.slope.0 = Us.slope.0, P.0 = P.0, st.0 = st.0,Bs.0 = Bs.0,
                           B.CR = B.CR, Bs.CR = Bs.CR, Bs.0.CR = Bs.0.CR,
                           
                           nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega, Otime = Otime, Wtime = Wtime,
                           Os = Os, Ws = Ws, O_base = O_base, W_base=W_base, correlated_re = correlated_re,
                           
                           
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)
  if(variability_hetero){
    Zq <- randtoolbox::sobol(S2,  nb.e.a+nb.e.a.sigma, normal = TRUE,scrambling = 1)
  }
  else{
    Zq <- randtoolbox::sobol(S2, dim = nb.e.a, mixed = TRUE, normal = TRUE)
  }
  cat("Second estimation  \n")
  estimation2 <- marqLevAlg(estimation$b, fn = log_llh, minimize = FALSE,
                            nb.e.a = nb.e.a, nb.priorMean.beta = nb.priorMean.beta,nb.alpha = nb.alpha,
                            competing_risk = competing_risk,
                            nb.alpha.CR = nb.alpha.CR, variability_hetero = variability_hetero, S = S2,Zq = Zq, sharedtype = sharedtype,
                            sharedtype_CR = sharedtype_CR,
                            hazard_baseline = hazard_baseline, hazard_baseline_CR = hazard_baseline_CR, ord.splines = ord.splines, Xtime = Xtime,
                            Utime = Utime, nb_pointsGK = nb_pointsGK,
                            Xs = Xs,Us = Us, Xslope = Xslope, Uslope = Uslope, Xs.slope = Xs.slope, Us.slope = Us.slope,
                            indices_beta_slope = indices_beta_slope, Time =Time,
                            st_calc = st_calc, B = B, Bs = Bs, wk = wk, Z = Z, P = P, left_trunc = left_trunc,
                            Z_CR = Z_CR, X_base = X_base, offset = offset, U = U, y.new.prog = y.new.prog, event1 = event1, event2 = event2, Ind = Ind,
                            Xs.0 = Xs.0, Us.0 = Us.0, Xs.slope.0 = Xs.slope.0, Us.slope.0 = Us.slope.0, P.0 = P.0, st.0 = st.0,Bs.0 = Bs.0,
                            B.CR = B.CR, Bs.CR = Bs.CR, Bs.0.CR = Bs.0.CR,
                            
                            nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega, Otime = Otime, Wtime = Wtime,
                            Os = Os, Ws = Ws, O_base = O_base, W_base=W_base, correlated_re = correlated_re,
                            
                            
                            nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info, file = file,
                            blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)
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
                                        formFixedVar = formFixedVar,
                                        formRandomVar = formRandomVar,
                                        formGroup = formGroup,
                                        timeVar = timeVar,
                                        variability_hetero = variability_hetero,
                                        formSlopeFixed = formSlopeFixed,
                                        formSlopeRandom = formSlopeRandom,
                                        indices_beta_slope = indices_beta_slope,
                                        formSurv = formSurv,
                                        data.long = data.long,
                                        precision = precision,
                                        sharedtype = sharedtype,
                                        nb_pointsGK = nb_pointsGK,
                                        hazard_baseline = hazard_baseline,
                                        S1 = S1, S2 = S2, nb.e.a = nb.e.a,
                                        ord.splines = ord.splines,
                                        left_trunc = left_trunc,
                                        competing_risk = competing_risk,
                                        formSurv_CR = formSurv_CR,
                                        hazard_baseline_CR = hazard_baseline_CR,
                                        sharedtype_CR = sharedtype_CR,
                                        Time.0 = Time.0,
                                        nproc = nproc, clustertype = clustertype,
                                        maxiter = maxiter, print.info = print.info,
                                        file = file, epsa = epsa, epsb = epsb, epsd = epsd,
                                        nb.priorMean.beta = nb.priorMean.beta, nb.alpha = nb.alpha,
                                        nb.alpha.CR = nb.alpha.CR,
                                        knots.hazard_baseline.splines = rr,
                                        knots.hazard_baseline.splines.CR = rr.CR,
                                        nb.omega = nb.omega,
                                        nb.e.a.sigma = nb.e.a.sigma,
                                        correlated_re = correlated_re,
                                        Ind = Ind, conv = estimation2$istop, niter = estimation2$ni,
                                        convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm),
                                        names_long = colnames(X_base), names_surv = colnames(Z),
                                        names_surv2 = colnames(Z_CR),
                                        likelihood_value = estimation2$fn.value,
                                        names_param = names_param)
                        
  )
  class(final_object) <- c("FlexVar_JM")
  return(final_object)
  
}
