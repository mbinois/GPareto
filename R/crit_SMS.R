##' Computes a slightly modified infill Criterion of the SMS-EGO. 
##' To avoid numerical instabilities, an additional penalty is added to the new point if it is too close to an existing observation.
##' @title Analytical expression of the SMS-EGO criterion with m>1 objectives
##' @param x a vector representing the input for which one wishes to calculate the criterion,
##' @param model a list of objects of class \code{\link[DiceKriging]{km}} (one for each objective),
##' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, or any reference set of observations,  
##' @param critcontrol list with arguments: 
##' \itemize{
##' \item \code{currentHV} current hypervolume;
##' \item \code{refPoint} reference point for hypervolume computations;
##' \item \code{extendper} if no reference point \code{refPoint} is provided,
##'  for each objective it is fixed to the maximum over the Pareto front plus extendper times the range. 
##'  Default value to \code{0.2}, corresponding to \code{1.1} for a scaled objective with a Pareto front in \code{[0,1]^n.obj};
##' \item \code{epsilon} optional value to use in additive epsilon dominance;
##' \item \code{gain} optional gain factor for sigma.
##' }
##'        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend 
##'        estimation has to be taken into account.
##' @references 
##' W. Ponweiser, T. Wagner, D. Biermann, M. Vincze (2008), Multiobjective Optimization on a Limited Budget of Evaluations Using Model-Assisted S-Metric Selection,
##' \emph{Parallel Problem Solving from Nature}, pp. 784-794. Springer, Berlin. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, pp. 718-727. Springer, Berlin.
##' @seealso \code{\link[GPareto]{crit_EHI}}, \code{\link[GPareto]{crit_SUR}}, \code{\link[GPareto]{crit_EMI}}.
##' @return Value of the criterion.
##' @examples
##' #---------------------------------------------------------------------------
##' # SMS-EGO surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' f_name <- "P1" 
##' n.grid <- 26
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' n_appr <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
##' response.grid <- t(apply(design.grid, 1, f_name))
##' PF <- t(nondominated_points(t(response.grid)))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' 
##' model <- list(mf1, mf2)
##' critcontrol <- list(refPoint = c(300, 0), currentHV = dominated_hypervolume(t(PF), c(300, 0)))
##' SMSEGO_grid <- apply(test.grid, 1, crit_SMS, model = model,
##'                      paretoFront = PF, critcontrol = critcontrol)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                matrix(pmax(0, SMSEGO_grid), nrow = n.grid), nlevels = 50,
##'                main = "SMS-EGO criterion (positive part)", xlab = expression(x[1]),
##'                ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1],design.grid[,2], pch = 21, bg = "white")
##'                             }
##'               )
##' @export

crit_SMS <- function(x, model, paretoFront=NULL, critcontrol=NULL, type="UK")
{
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  
  distp <- 0  # penalty if too close in the checkPredict sense
  
  if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    # return(0) Not compatible with penalty with SMS
    distp <- 1 # may be changed
  }#else{
  
  refPoint  <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  epsilon   <- critcontrol$epsilon
  gain      <- critcontrol$gain
  nsteps.remaining <- critcontrol$nsteps.remaining
  
  if(is.null(paretoFront) || is.null(refPoint)) {
    observations <- Reduce(cbind, lapply(model, slot, "y"))
  }
  if(is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)){
    if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
    # refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  n.pareto <- nrow(paretoFront)
  
  if (is.null(currentHV)) currentHV <- dominated_hypervolume(points=t(paretoFront), ref=refPoint)
  if (is.null(nsteps.remaining)) nsteps.remaining <- 1
  if (is.null(gain))  gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
  
  if (is.null(epsilon)) {
    if (n.pareto < 2){
      epsilon <- rep(0, n.obj)
    } else {
      spread <- apply(paretoFront,2,max) - apply(paretoFront,2,min)
      c <- 1 - (1 / (2^n.obj) )
      epsilon <- spread / (n.pareto + c * (nsteps.remaining-1))
    }
  }
  
  # mu    <- rep(NaN, n.obj)
  # sigma <- rep(NaN, n.obj)
  # for (i in 1:n.obj){    
  #   pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  #   mu[i]    <- pred$mean
  #   sigma[i] <- pred$sd
  # }
  pred <- predict_kms(model, newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.numeric(pred$mean)
  sigma <- as.numeric(pred$sd)
  
  potSol <- mu - gain*sigma
  penalty <- distp
  for (j in 1:n.pareto){
    # assign penalty to all epsilon-dominated solutions
    if (min(paretoFront[j,] <= potSol + epsilon)){
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j,], rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0){
    # non epsilon-dominated solution
    potFront <- rbind(paretoFront, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points=t(mypoints), ref=refPoint)
    f    <- currentHV - myhv
  } else{
    f <- penalty
  }
  #}
  
  return(-f)
}