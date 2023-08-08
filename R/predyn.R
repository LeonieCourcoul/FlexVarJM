#' Dynamic prediction for a new individual
#'
#' @param newdata
#' @param FlexVarJM
#' @param s
#' @param window
#' @param event
#' @param L
#' @param graph
#'
#' @return
#'
#' @examples

predyn <- function(newdata, FlexVarJM, s, window =1, event = 1, L =500, graph = F){
  bootstrap <- c()
  pred.ponct <- c()
  for( t in seq(0.1,window,0.1)){
    pred.t <- pred_s.t.ponctuel.tps(newdata,FlexVarJM,s,t,event)
    pred.t.bootstrap.tps <- pred_s.t.bootstrap(newdata,FlexVarJM,s,t,event,L)
    pred.ponct <- c(pred.ponct, pred.t)
    bootstrap<- cbind(bootstrap,pred.t.bootstrap)
  }
  table.pred <- cbind(apply(bootstrap,2, function(x) quantile(x, 0.025)),
                      pred.ponct,
                      apply(bootstrap,2, function(x) quantile(x, 0.975)))
  table.pred <- as.data.frame(table.pred)
  colnames(table.pred) <- c("ICinf", "Prediction", "ICsup")
  graph.predyn <- NULL

  if(graph){
    x.axe <- c(newdata[,all.vars(FlexVarJM$control$formFixed)[2]],seq(s,window+s,0.1))
    print(x.axe)
    y.axe <- c(newdata[,all.vars(FlexVarJM$control$formFixed)[1]], rep(NA,length(seq(s,window+s,0.1))))
    print(y.axe)
    y.axe2 <- c(rep(NA,length(newdata[,all.vars(FlexVarJM$control$formFixed)[2]])),0,table.pred$ICinf)
    print(y.axe2)
    y.axe3 <- c(rep(NA,length(newdata[,all.vars(FlexVarJM$control$formFixed)[2]])),0,table.pred$Prediction)
    print(y.axe3)
    y.axe4 <- c(rep(NA,length(newdata[,all.vars(FlexVarJM$control$formFixed)[2]])),0,table.pred$ICsup)
    plot(x = x.axe, y = y.axe,xlim = c(0,s+window),
        xlab = "Time", ylab = "Marker", cex.lab = 1, col = "black",
         main = "Prediction of event", pch = 20, cex = 1, font = 1, font.lab = 1, cex.lab = 1, cex.main = 1)
    graphics::par(new = TRUE, font = 1, cex.lab = 1)
      plot(x.axe, y.axe3, axes = FALSE,col = "black", type = "l", ylim = c(0.000001, max(table.pred$ICsup, na.rm = T)),ylab = "", xlab = "", lwd =2, font.lab = 1, cex.lab = 1 )
      graphics::lines(x.axe, y.axe2, col = "black", lty=2, lwd = 2)
      graphics::lines(x.axe, y.axe4, col = "black", lty=2, lwd = 2)
      graphics::axis(side= 4, cex = 2)
      graphics::abline(v = s, lty = 3)
      graphics::mtext("Probability of event", cex = 1, side = 4, line = 3, font.lab = 1)

  }

  result <- table.pred
  result
}
