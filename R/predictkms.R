#' Predict function for list of \code{\link[DiceKriging]{km}} models.
#' @param model list of \code{\link[DiceKriging]{km}} models
#' @param newdata,type,se.compute,cov.compute,light.return,bias.correct,checkNames,... see \code{\link[DiceKriging]{predict.km}}
#' @export 
#' @details So far only \code{light.return = TRUE} handled.
predict_kms <- function(model, newdata, type, se.compute = TRUE, 
                        cov.compute = FALSE, light.return = TRUE,
                        bias.correct = FALSE, checkNames = FALSE, ...){
  
  pred <- lapply(model, FUN=predict, newdata=newdata, checkNames=checkNames, type=type, cov.compute = cov.compute, light.return = light.return)
  
  if(cov.compute)
    return(list(mean =  Reduce(rbind, lapply(pred, function(alist) alist$mean)),
                sd = Reduce(rbind, lapply(pred, function(alist) alist$sd)),
                cov = lapply(pred, function(alist) alist$cov)))
  return(list(mean =  Reduce(rbind, lapply(pred, function(alist) alist$mean)),
              sd = Reduce(rbind, lapply(pred, function(alist) alist$sd))))
}

