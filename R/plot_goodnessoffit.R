plot_goodnessoffit<-function(data.long,data.id,pred.CV,break.times, formFixed, formSurv,timeVar,Cum_risk1, competing_risk, formSurv_CR,Cum_risk2){

  #Longitudinal part
  #if(is.null(break.times)) break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
  value.var <- as.character(formFixed[[2]])
  data.long$window <- cut(data.long[,timeVar], break.times, include.lowest = T)
  mean.obs <- by(data.long[,value.var], data.long$window, mean)
  sd.obs <- by(data.long[,value.var], data.long$window, sd)
  length.obs <- by(data.long[,value.var], data.long$window, length)
  IC.inf <- mean.obs - 1.96*sd.obs/sqrt(length.obs)
  IC.sup <- mean.obs + 1.96*sd.obs/sqrt(length.obs)
  prediction <- cbind(pred.CV, data.long$window)
  mean.pred <- by(prediction[,2], prediction[,ncol(prediction)], mean)
  obstime.mean <- by(data.long[,timeVar], data.long$window, mean)
  df <- cbind(obstime.mean, mean.obs, IC.sup, IC.inf, mean.pred)
  df <- as.data.frame(df)
  oldpar <- graphics::par(no.readonly = TRUE) # code line i
  on.exit(graphics::par(oldpar)) # code line i + 1
  k <- ggplot2::ggplot(df,  ggplot2::aes(obstime.mean, mean.obs, ymin = IC.sup, ymax = IC.inf))
  graph.long <- k +  ggplot2::geom_pointrange( ggplot2::aes(ymin = IC.sup, ymax = IC.inf), shape =1) +
     ggplot2::geom_point(ggplot2::aes(obstime.mean, mean.pred), size = 3, shape = 17) +
     ggplot2::scale_x_continuous(name = "Time") +
     ggplot2::scale_y_continuous(name = "Current Value") +
     ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
    ggplot2::ggtitle("Longitudinal goodness-of-fit")

  #Survival part

  data.id$e1.new <- ifelse(data.id[,all.vars(formSurv)[2]] == 1, 1, 0)
  C1.sort <- data.id[order(data.id[,all.vars(formSurv)[1]]),]
  Cum.pred1 <- apply(Cum_risk1, 2, mean)
  Cum.pred1 <- cbind(Cum.pred1, unique(sort(data.id[,all.vars(formSurv)[1]])))
  Cum.pred1 <- as.data.frame(Cum.pred1)
  Cum.pred1.sort <- Cum.pred1[order(Cum.pred1[,2]),]
  colnames(Cum.pred1.sort) <- c("pred","timeFormSurv")
  timeFormSurv <- Cum.pred1.sort$timeFormSurv
  pred <- Cum.pred1.sort$pred
    #arrange(Cum.pred1, V2)
  Surv.fit1 <- survminer::surv_fit(as.formula(paste(as.character(formSurv)[2],1,sep="~")), data = C1.sort)
  graph.surv.1 <- survminer::ggsurvplot(Surv.fit1, data = C1.sort, fun = "cumhaz", 
                             conf.int.style = "step", legend = "none", xlab = "Time")$plot +
    ggplot2::geom_line(ggplot2::aes(timeFormSurv,pred),
              data = Cum.pred1.sort,
              color = "black",
              linetype = "3313",
              size = 1)
  graph.surv.2 <- NULL
  if(competing_risk){
    data.id$e2.new <- ifelse(data.id[,all.vars(formSurv_CR)[2]] == 2, 1, 0)
    C2.sort <- data.id[order(data.id[,all.vars(formSurv_CR)[1]]),]
    Cum.pred2 <- apply(Cum_risk2, 2, mean)
    Cum.pred2 <- cbind(Cum.pred2, unique(sort(data.id[,all.vars(formSurv)[1]])))
    Cum.pred2 <- as.data.frame(Cum.pred2)
    Cum.pred2.sort <- Cum.pred2[order(Cum.pred2[,2]),]#arrange(Cum.pred2, V2)
    colnames(Cum.pred2.sort) <- c("pred","timeFormSurv")
    timeFormSurv <- Cum.pred2.sort$timeFormSurv
    pred <- Cum.pred2.sort$pred
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1
    Surv.fit2 <- survminer::surv_fit(as.formula(paste(as.character(formSurv_CR)[2],1,sep="~")), data = C2.sort)
    graph.surv.2 <- survminer::ggsurvplot(Surv.fit2, data = C2.sort, fun = "cumhaz",
                               conf.int.style = "step", legend = "none", xlab = "Time")$plot +
      ggplot2::geom_line(ggplot2::aes(timeFormSurv,pred),
                data = Cum.pred2.sort,
                color = "black",
                linetype = "3313",
                size = 1)

  }

  graphs <- list(graph.long = graph.long, graph.surv.1 = graph.surv.1, graph.surv.2 = graph.surv.2)
  graphs
}
