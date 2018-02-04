## ' Computes the Expected reduction of the volume of the excursion sets behind a 2D Pareto front.
## ' @title Expected reduction of the volume of excursion
## ' @param phi.x.bar see below,
## ' @param phi.x.tilde see below,
## ' @param  phi.eta.x see below, 
## ' @param  phi.y.bar see below, 
## ' @param  phi.y.tilde see below, 
## ' @param  phi.eta.y see below, 
## ' @param  phi2.x.x see below, 
## ' @param  phi2.x.eta see below, 
## ' @param  phi2.y.y see below, 
## ' @param  phi2.y.eta vectors of size n.pareto*n.integration.points. 'x' refers to the first objective, 'y' to the second.
## ' @details To be called by \code{\link[GPareto]{crit_SUR}}
## ' 
## ' @return list with three elements: 
## ' \code{crit} is the expected reduction
## ' \code{pijold} is the current probability of non-domination
## ' \code{pij} is the expected future probability of non-domination
## ' @references
## ' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
## ' \emph{Statistics and Computing}
## ' @export

EEV.2D.computation <- function(phi.x.bar, phi.x.tilde, phi.eta.x, phi.y.bar, phi.y.tilde, phi.eta.y, 
                               phi2.x.x, phi2.x.eta, phi2.y.y, phi2.y.eta)
{
  ############################################################################
  n.pareto <- length(phi.x.bar)
  n.integration.points <- ncol(phi.x.tilde)
  
  pijold  <- diff  <- matrix(0, n.pareto+1, n.integration.points)
  pijold[1,] <- phi.x.tilde[1,]
  diff[1,] <- (phi2.x.x[1,] - phi2.x.eta[1,])*(1 - phi.eta.y)
  
  if (n.pareto > 1){
    pijold[2:(n.pareto),] <- (phi.x.tilde[2:(n.pareto),] - phi.x.tilde[1:(n.pareto-1),])*phi.y.tilde[1:(n.pareto-1),]
    diff[2:(n.pareto),]   <- (phi2.x.x[2:(n.pareto),] - phi2.x.x[1:(n.pareto-1),] + phi2.x.eta[1:(n.pareto-1),] - phi2.x.eta[2:(n.pareto),] )*
      (phi2.y.y[1:(n.pareto-1),] - phi2.y.eta[1:(n.pareto-1),]) 
  }
  pijold[n.pareto+1,] <- (1 - phi.x.tilde[n.pareto,])*phi.y.tilde[n.pareto,]
  diff[n.pareto+1,]   <- (1 - phi.eta.x + phi2.x.eta[n.pareto,] - phi2.x.x[n.pareto,])*(phi2.y.y[n.pareto,] - phi2.y.eta[n.pareto,] )
  
  pij <- pijold - diff
  
  crit <- sum(diff)
  return(list(crit, pijold, pij))
}