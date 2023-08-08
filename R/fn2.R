fn2 <- function(Bs.gammas,event,W2,P,wk,Time,W2s,id.GK){
  -sum(event*drop(W2%*%Bs.gammas)-P*fastSumID(rep(wk,length(Time))*exp(drop(W2s%*%Bs.gammas)),id.GK))
}
