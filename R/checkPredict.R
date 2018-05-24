##' Check that the new point is not too close to already known observations to avoid numerical issues.
##' Closeness can be estimated with several distances.
##' @title Prevention of numerical instability for a new observation
##' @param x a vector representing the input to check, alternatively a matrix with one point per row,
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param threshold optional value for the minimal distance to an existing observation, default to \code{1e-4},
##' @param distance selection of the distance between new observations, between "\code{euclidean}", "\code{none}",
##'  "\code{covdist}" (default) and "\code{covratio}", see details,
##' @param type "\code{SK}" or "\code{UK}" (default), depending whether uncertainty related to trend estimation has to be taken into account.
##' @details If the distance between \code{x} and the closest observations in \code{model} is below
##'  \code{threshold}, \code{x} should not be evaluated to avoid numerical instabilities.
##' The distance can simply be the Euclidean distance or the canonical distance associated with the kriging covariance k: 
##'   \deqn{d(x,y) = \sqrt{k(x,x) - 2k(x,y) + k(y,y)}.}{d(x,y) = \sqrt(k(x,x) - 2k(x,y) + k(y,y)).} 
##'   The last solution is the ratio between the prediction variance at \code{x} and the variance of the process.
##'   \code{none} can be used, e.g., if points have been selected already.
##' @return \code{TRUE} if the point should not be tested.
##' @export
checkPredict <- function(x, model, threshold = 1e-4, distance = "covdist", type = "UK"){
  
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1) 
  }
  
  if(is.null(distance))
    distance <- "covdist"
  if(is.null(threshold))
    threshold <- 1e-4

  if(distance == "none") return(rep(FALSE, nrow(x)))
    
  if(distance == "euclidean"){
    # tp1 <- as.numeric(t(model[[1]]@X))
    # tp2 <- matrix(tp1 - as.numeric(x), nrow = model[[1]]@n, byrow = TRUE)^2
    mindist <- distcpp_2(x, model[[1]]@X)
    mindist <- apply(mindist, 1, min)
    mindist <- sqrt(mindist)

  }else{
    # d <- length(which(sapply(model, class) == "km")) # to filter the fastobj (no troubles to update)
    if(distance == "covratio"){
      
      mindist <- Inf
      for(i in 1:length(model)){
        if(class(model[[i]]) != "fastfun"){
          pred.sd <- predict(object = model[[i]], newdata = x, type = type, checkNames = FALSE, light.return = TRUE)
          mindist <- pmin(mindist, pred.sd$sd/sqrt(model[[i]]@covariance@sd2))
        }
      }
      
    }else{
      
      ## covdist
      
      # Note : general expression but it is much faster for stationary kernels

      mindist <- Inf
      for(i in 1:length(model)){
        if(class(model[[i]]) != "fastfun"){
          
          if(class(model[[i]]@covariance) %in% c("covTensorProduct", "covIso")){
            kxy <- covMat1Mat2(model[[i]]@covariance, x, model[[i]]@X)
            kxx <- model[[i]]@covariance@sd2 # k(x,x) = k(y,y) = variance term
            mindist <- pmin(mindist, sqrt(pmax(0, 2 * kxx - 2 * apply(kxy, 1, max))/model[[i]]@covariance@sd2))
            
          }else{
            kxx <- diag(covMatrix(model[[i]]@covariance, x)$C)
            kxy <- covMat1Mat2(model[[i]]@covariance, x, model[[i]]@X)
            kyy <- diag(covMatrix(model[[i]]@covariance, model[[i]]@X)$C)
            mindist <- pmin(mindist,
                            sqrt(pmax(0, kxx - 2 * apply(kxy - matrix(kyy, nrow(x), model[[i]]@n, byrow = T), 1, max))/model[[i]]@covariance@sd2))
          }
          
          # kxx <- drop(covMat1Mat2(model[[i]]@covariance, x, matrix(x, nrow = 1)))
          # kyy <- drop(diag(covMat1Mat2(model[[i]]@covariance, model[[i]]@X, model[[i]]@X)))
          # mindist <- sqrt(pmax(0, min(mindist, min(kxx - 2*kxy + kyy) / model[[i]]@covariance@sd2)))
        }
      } 
    }
  }
  
  return(mindist < threshold)
  # if(mindist < threshold){
  #   return(TRUE)
  # }else{
  #   return(FALSE)
  # }
}