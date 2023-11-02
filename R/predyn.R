#' Dynamic prediction for a new individual
#'
#' @param newdata data frame : collected data for a new individual
#' @param object lsjm object : estimation of the model
#' @param s numeric : the time to begin prediction
#' @param times numeric vector : future times to calculate predictions
#' @param event integer (0, 1 or 2) : the event of interest for the prediction
#' @param IC integer : percentage of confidence for the interval confidence (between 0 and 100), 95 by default, NULL if no IC
#' @param nb.draws integer : the number of simulations to compute the interval confidence (by bootstrap), 500 by default
#' @param graph boolean : indicator to plot the graphs or not
#'
#' @return
#'
#' @examples

predyn <- function(newdata, object, s, times, event = 1, IC = 95, nb.draws = 500, graph = F){
  bootstrap <- c()
  pred.ponct <- c()
  id <- as.integer(newdata[all.vars(object$control$formGroup)][,1])
  newdata$id <- id
  window <- times - s
  if(!is.null(IC)){
    newdata.id <- newdata[which(newdata$id == i),]
    for(i in unique(id)){
      for(t in window){
        pred.t <- pred_s.t.ponctuel.tps(newdata.id,object,s,t,event)
        pred.t.bootstrap.tps <- pred_s.t.bootstrap(newdata.id,object,s,t,event,L)
        pred.ponct <- c(pred.ponct, pred.t)
        bootstrap<- cbind(bootstrap,pred.t.bootstrap)
      }
      table.pred.id <- cbind(i, times, pred.ponct, apply(bootstrap,2, function(x) quantile(x, 0.50)),
                             apply(bootstrap,2, function(x) quantile(x, ((100-IC)/2)/100)),
                             apply(bootstrap,2, function(x) quantile(x, 1-((100-IC)/2)/100)),
                             apply(bootstrap,2, sd)
      )
      colnames(table.pred.id) <- c("ID","Time","Prediction","Median","ICinf","ICsup", "Empirical SD")
      if(graph){
        x.axe <- c(newdata.id[,all.vars(object$control$formFixed)[2]],times)
        #print(x.axe)
        y.axe <- c(newdata.id[,all.vars(object$control$formFixed)[1]], rep(NA,length(window)))
        #print(y.axe)
        y.axe2 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$ICinf)
        #print(y.axe2)
        y.axe3 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$Median)
        #print(y.axe3)
        y.axe4 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$ICsup)
        y.axe5 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$Prediction)
        plot(x = x.axe, y = y.axe,xlim = c(0,s+window),
             xlab = "Time", ylab = "Marker", cex.lab = 1, col = "black",
             main = "Prediction of event", pch = 20, cex = 1, font = 1, font.lab = 1, cex.lab = 1, cex.main = 1)
        graphics::par(new = TRUE, font = 1, cex.lab = 1)
        plot(x.axe, y.axe3, axes = FALSE,col = "black", type = "l", ylim = c(0.000001, max(table.pred.id$ICsup, na.rm = T)),ylab = "", xlab = "", lwd =2, font.lab = 1, cex.lab = 1 )
        graphics::lines(x.axe, y.axe5, col = "red", lty=2, lwd = 2)
        graphics::lines(x.axe, y.axe2, col = "black", lty=2, lwd = 2)
        graphics::lines(x.axe, y.axe4, col = "black", lty=2, lwd = 2)
        graphics::axis(side= 4, cex = 2)
        graphics::abline(v = s, lty = 3)
        graphics::mtext("Probability of event", cex = 1, side = 4, line = 3, font.lab = 1)
        
      }
      table.pred <- rbind(table.pred, table.pred.id)
    }
    table.pred <- as.data.frame(table.pred)
    colnames(table.pred) <- c("ID","Time","Prediction","Median","ICinf","ICsup", "Empirical SD")
    graph.predyn <- NULL
  }
  else{
    for(i in unique(id)){
      newdata.id <- newdata[which(newdata$id == i),]
      for(t in window){
        pred.t <- pred_s.t.ponctuel.tps(newdata.id,object,s,t,event)
        pred.ponct <- c(pred.ponct, pred.t)
      }
      table.pred.id <- cbind(i, times, pred.ponct
      )
      table.pred <- rbind(table.pred, table.pred.id)
    }
    table.pred <- as.data.frame(table.pred)
    colnames(table.pred) <- c("ID","Time","Prediction")
    graph.predyn <- NULL
    if(graph){
      x.axe <- c(newdata.id[,all.vars(object$control$formFixed)[2]],times)
      #print(x.axe)
      y.axe <- c(newdata.id[,all.vars(object$control$formFixed)[1]], rep(NA,length(window)))
      #print(y.axe)
      #y.axe2 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$ICinf)
      #print(y.axe2)
      #y.axe3 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$Median)
      #print(y.axe3)
      #y.axe4 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$ICsup)
      y.axe5 <- c(rep(NA,length(newdata.id[,all.vars(object$control$formFixed)[2]])),0,table.pred.id$Prediction)
      plot(x = x.axe, y = y.axe,xlim = c(0,s+window),
           xlab = "Time", ylab = "Marker", cex.lab = 1, col = "black",
           main = "Prediction of event", pch = 20, cex = 1, font = 1, font.lab = 1, cex.lab = 1, cex.main = 1)
      graphics::par(new = TRUE, font = 1, cex.lab = 1)
      plot(x.axe, y.axe5, axes = FALSE,col = "black", type = "l", ylim = c(0.000001, max(table.pred.id$ICsup, na.rm = T)),ylab = "", xlab = "", lwd =2, font.lab = 1, cex.lab = 1 )
      #graphics::lines(x.axe, y.axe5, col = "red", lty=2, lwd = 2)
      #graphics::lines(x.axe, y.axe2, col = "black", lty=2, lwd = 2)
      #graphics::lines(x.axe, y.axe4, col = "black", lty=2, lwd = 2)
      graphics::axis(side= 4, cex = 2)
      graphics::abline(v = s, lty = 3)
      graphics::mtext("Probability of event", cex = 1, side = 4, line = 3, font.lab = 1)
      
    }
  }
  

  result <- table.pred
  result
}
