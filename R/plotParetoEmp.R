##' Plot the Pareto front with step functions.
##' 
##' @title Pareto front visualization
##' @param nondominatedPoints points considered to plot the Pareto front with segments, matrix with one point per row,
##' @param add optional boolean indicating whether a new graphic should be drawn,
##' @param max optional boolean indicating whether to display a Pareto front in a maximization context,
##' @param bounds for 3D, optional 2*nobj matrix of boundaries
##' @param alpha for 3D, optional value in [0,1] for transparency
##' @param ... additional values to be passed to the \code{\link[graphics]{lines}} function.
##' @export
##' @importFrom rgl plot3d quads3d
##' @examples
##' #------------------------------------------------------------
##' # Simple example
##' #------------------------------------------------------------
##' 
##' x <- c(0.2, 0.4, 0.6, 0.8)
##' y <- c(0.8, 0.7, 0.5, 0.1)
##' 
##' plot(x, y, col = "green", pch = 20) 
##' 
##' plotParetoEmp(cbind(x, y), col = "green")
##' ## Alternative
##' plotParetoEmp(cbind(x, y), col = "red", add = FALSE)
##' 
##' ## With maximization
##' plotParetoEmp(cbind(x, y), col = "blue", max = TRUE)
##' 
##' ## 3D plots
##' library(rgl)
##' set.seed(5)
##' X <- matrix(runif(60), ncol=3)
##' Xnd <- t(nondominated_points(t(X)))
##' plot3d(X)
##' plot3d(Xnd, col="red", size=8, add=TRUE)
##' plot3d(x=min(Xnd[,1]), y=min(Xnd[,2]), z=min(Xnd[,3]), col="green", size=8, add=TRUE)
##' X.range <- diff(apply(X,2,range))
##' bounds <- rbind(apply(X,2,min)-0.1*X.range,apply(X,2,max)+0.1*X.range)
##' plotParetoEmp(nondominatedPoints = Xnd, add=TRUE, bounds=bounds, alpha=0.5)
##' 

plotParetoEmp <- function(nondominatedPoints, add = TRUE, max = FALSE, bounds=NULL, alpha=0.5, ...){
  if (length(dim(nondominatedPoints)) != 2) {
    cat("The nondominatedPoints argument should be a matrix \n")
  } else {
    
    #---- 2D CASE ---------------
    if (ncol(nondominatedPoints) == 2) {
      if(add == FALSE){
        plot(nondominatedPoints, ...)
      }
      
      temp <- nondominatedPoints[order(nondominatedPoints[,1]), , drop = FALSE]
      
      if(max){
        # Limiting points
        if(add){
          lim_left <- par('usr')[1]
          lim_bottom <- par('usr')[3]
        }else{
          lim_left <- -abs(20*temp[nrow(temp),1])
          lim_bottom <- -abs(20*temp[1,2])
        }
        
        lines(c(lim_left, temp[1, 1]), c(temp[1, 2], temp[1, 2]), ...)
        
        lines(c(temp[dim(temp)[1], 1], temp[dim(temp)[1], 1]), c(temp[dim(temp)[1], 2], lim_bottom), ...)
        
        # Segments in between
        if(nrow(temp) > 1){
          for (i in 1:(nrow(temp) - 1)) {
            lines(c(temp[i, 1], temp[i , 1]), c(temp[i, 2], temp[i+ 1, 2]), ...)
            lines(c(temp[i, 1], temp[i + 1, 1]), c(temp[i+1, 2],  temp[i + 1, 2]), ...)
          }
        }
        
      }else{
        # Limiting points
        if(add){
          lim_right <- par('usr')[2]
          lim_up <- par('usr')[4]
        } else {
          lim_right <- abs(20*temp[nrow(temp),1])
          lim_up <- abs(20*temp[1,2])
        }
        
        lines(c(temp[1,1], temp[1,1]),c(lim_up, temp[1,2]),...) #first segment
        
        lines(c(temp[dim(temp)[1], 1], lim_right),
              c(temp[dim(temp)[1], 2], temp[dim(temp)[1], 2]),...) #last segment
        
        # Segments in between
        if(nrow(temp) > 1){
          for(i in 1:(nrow(temp)-1)){
            lines(c(temp[i,1], temp[i+1,1]),c(temp[i,2], temp[i,2]),...) #horizontal part
            lines(c(temp[i+1,1], temp[i+1,1]),c(temp[i,2], temp[i+1,2]),...) #vertical part
          }
        }
      }
    } else if (ncol(nondominatedPoints) == 3) {
      
      if (max) {
        cat("plotParetoEmp with 3 objectives only works for minimization \n")
      } else {
        #---- 3D CASE ---------------
        paretoFront <- nondominatedPoints
        Front.range <- diff(apply(paretoFront,2,range))
        if (is.null(bounds)) bounds <- rbind(apply(paretoFront,2,min)-0.1*Front.range,apply(paretoFront,2,max)+0.1*Front.range)
        n.pareto <- nrow(paretoFront)
        
        if (add == FALSE){
          plot3d(x=paretoFront[,1], y=paretoFront[,2], z=paretoFront[,3],
                 xlim=bounds[,1],ylim=bounds[,2],zlim=bounds[,3], xlab="f1", ylab="f2",zlab="f3", col="red", size=4)
        }
        
        # First plot (xz) and (yz) facets
        paretoFront <- paretoFront[order(paretoFront[,3], decreasing = TRUE),,drop=FALSE]
        
        for (i in 1:(n.pareto)){
          if (i==1) z2 <- bounds[2,3] else z2 <- paretoFront[i-1,3]
          z1 <- paretoFront[i,3]
          
          # Subset of the Pareto front
          front2D <- t(nondominated_points(t(paretoFront[(i:n.pareto),1:2, drop=FALSE])))
          # Reorder by increasing first objective
          front2D <- front2D[order(front2D[,1]),, drop=FALSE]
          
          if (nrow(front2D)>1) {
            for (j in 1:(nrow(front2D)-1)) {
              quads3d(x=front2D[c(0,1,1,0) + j,1], y=front2D[rep(j,4),2], z=c(z1,z1,z2,z2), col="blue", alpha=alpha)
              quads3d(x=front2D[rep(j+1,4),1], y=front2D[c(0,1,1,0) + j,2], z=c(z1,z1,z2,z2), col="darkblue", alpha=alpha)
            }
          }
          quads3d(x=c(max(front2D[,1]), bounds[2,1], bounds[2,1], max(front2D[,1])),
                  y=rep(min(front2D[,2]),4), z=c(z1,z1,z2,z2), col="blue", alpha=alpha)
          
          quads3d(y=c(max(front2D[,2]), bounds[2,2], bounds[2,2], max(front2D[,2])),
                  x=rep(min(front2D[,1]),4), z=c(z2,z2,z1,z1), col="darkblue", alpha=alpha)
        }
        
        # Now plot (xy) facets
        paretoFront <- paretoFront[order(paretoFront[,1], decreasing = TRUE),,drop=FALSE]
        
        for (i in 1:(n.pareto)){
          if (i==1) x2 <- bounds[2,1] else x2 <- paretoFront[i-1,1]
          x1 <- paretoFront[i,1]
          
          # Subset of the Pareto front
          front2D <- t(nondominated_points(t(paretoFront[(i:n.pareto),2:3, drop=FALSE])))
          # Reorder by increasing first objective
          front2D <- front2D[order(front2D[,1]),, drop=FALSE]
          
          if (nrow(front2D)>1) {
            for (j in 1:(nrow(front2D)-1)) {
              quads3d(y=front2D[c(0,1,1,0) + j,1], z=front2D[rep(j,4),2], x=c(x1,x1,x2,x2), col="blue", alpha=alpha)
            }
          }
          quads3d(y=c(max(front2D[,1]), bounds[2,2], bounds[2,2], max(front2D[,1])),
                  z=rep(min(front2D[,2]),4), x=c(x1,x1,x2,x2), col="blue", alpha=alpha)
        }
      }
    } else {
      cat("plotParetoEmp only works with 2 or 3 objectives \n")
    }
  }
}