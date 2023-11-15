log_llh_lsmm <- function(param, nb.e.a, nb.priorMean.beta, variability_hetero, S,Zq, X_base, offset, U, y.new.prog,Ind,
                         nb.e.a.sigma, nb.omega , 
                         O_base , W_base
                    
){
  #Manage parameters
  
  borne1 <- choose(n = nb.e.a, k = 2) + nb.e.a
  C <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
  C[lower.tri(C, diag=T)] <- param[1:borne1] #1,2,3
  borne2 <- borne1 + nb.priorMean.beta #5
  beta <- param[(borne1+1):borne2] #4,5
  if(variability_hetero){
    borne3 <- borne2+nb.omega #7
    omega <- param[(borne2+1):(borne3)] #6,7
    borne4 <- choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma + borne3 #10
    C.sigma <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
    C.sigma[lower.tri(C.sigma, diag=T)] <- param[(borne3+1):borne4] #8,9,10
  }
  else{
    sigma.epsilon <- abs(param[borne2+1]) #6
  }
  
  #Manage random effects
  b_al <- Zq[,1:(nb.e.a)]%*%t(C)
  b_om <- Zq[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]%*%t(C.sigma)
  
  ll_glob <- 0
  
  for(i in 1:Ind){#Computation of contribution to the log_lokelihood
    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
    O_base_i <- O_base[offset[i]:(offset[i+1]-1),]
    O_base_i <- matrix(O_base_i, nrow = offset[i+1]-offset[i])
    W_base_i <- W_base[offset[i]:(offset[i+1]-1),]
    W_base_i <- matrix(W_base_i, nrow = offset[i+1]-offset[i])
    U_i <- U[offset[i]:(offset[i+1]-1),]
    U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
    y_i <- y.new.prog[offset[i]:(offset[i+1]-1)]
    if(is.null(nrow(X_base_i))){
      if(variability_hetero){
        sigma.long <- exp((omega%*%O_base_i)[1,1] + b_om%*%W_base_i)
      }
      else{
        sigma.long <- sigma.epsilon
      }
      CV <- (beta%*%X_base_i)[1,1] + b_al%*%U_i
      f_Y_b_sigma <- log((1/(sqrt(2*pi)*sigma.long))) + (-1/2)*((y_i-CV)/sigma.long)**2 #dnorm(x=y_i, mean = CV, sd = sigma.long)
    }
    else{
      f_Y_b_sigma <- rep(0,S)
      for(k in 1:nrow(X_base_i)){
        if(variability_hetero){
          sigma.long <- exp((omega%*%O_base_i[k,])[1,1] + b_om%*%W_base_i[k,])
        }
        else{
          sigma.long <- sigma.epsilon
        }
        CV <- (beta%*%X_base_i[k,])[1,1] + b_al%*%U_i[k,]
        f_Y_b_sigma <- log((1/(sqrt(2*pi)*sigma.long)))+f_Y_b_sigma + (-1/2)*((y_i[k]-CV)/sigma.long)**2 #*dnorm(x = y_i[k], mean = CV, sd = sigma.long)
      }
    }
    log_dens_int <- f_Y_b_sigma
    Clogexp <- max(log_dens_int) - 500
    log_dens_int <- log_dens_int - Clogexp
    log_dens <- Clogexp +log(sum(exp(log_dens_int))) - log(S)
    ll_glob <- ll_glob + log_dens
    }
  ll_glob
}
