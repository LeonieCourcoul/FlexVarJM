data.manag.sigma <- function(formGroup, formFixed, formRandom, data.long1){
  
  data_long <- data.long1[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
  #y.new.prog <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  list.long <- list("X" = X, "U" = U)
  
  return(list.long)
}