##' Predict function for list of \code{\link[DiceKriging]{km}} models.
##' @param model list of \code{\link[DiceKriging]{km}} models
##' @param newdata,type,se.compute,cov.compute,light.return,bias.correct,checkNames,... see \code{\link[DiceKriging]{predict.km}}
##' @export 
##' @details So far only \code{light.return = TRUE} and \code{cov.compute = FALSE} handled.
predict_kms <- function(model, newdata, type, se.compute = TRUE, 
                       cov.compute = FALSE, light.return = FALSE,
                       bias.correct = FALSE, checkNames = TRUE, ...){
  
  pred <- lapply(model, FUN=predict, newdata=newdata, checkNames=FALSE, type=type, cov.compute = FALSE, light.return = TRUE)
  
  return(list(mean =  Reduce(rbind, lapply(pred, function(alist) alist$mean)), sd = Reduce(rbind, lapply(pred, function(alist) alist$sd))))
}

