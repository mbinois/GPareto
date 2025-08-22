#' Multi-objective Expected Hypervolume Improvement with respect to the
#' current Pareto front. With two objectives the analytical formula is used, while 
#' Sample Average Approximation (SAA) is used with more objectives.
#' To avoid numerical instabilities, the new point is penalized if it is too close to an existing observation.
#' 
#' @title Expected Hypervolume Improvement with m objectives
#' 
#' @param x a vector representing the input for which one wishes to calculate \code{EHI}, alternatively a matrix with one point per row,
#' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
#' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, or any reference set of observations, 
#' @param critcontrol optional list with arguments:
#' \itemize{
#' \item \code{nb.samp} number of random samples from the posterior distribution (with more than two objectives),
#'        default to \code{50}, increasing gives more reliable results at the cost of longer computation time;
#' \item \code{seed} seed used for the random samples (with more than two objectives);
#' \item \code{refPoint} reference point for Hypervolume Expected Improvement;
#' \item \code{extendper} if no reference point \code{refPoint} is provided,
#'  for each objective it is fixed to the maximum over the Pareto front plus extendper times the range, 
#'  Default value to \code{0.2}, corresponding to \code{1.1} for a scaled objective with a Pareto front in \code{[0,1]^n.obj}.
#' } 
#'   
#'        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
#' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
#'        has to be taken into account. 
#' @return The Expected Hypervolume Improvement at \code{x}.
#' @seealso \code{[DiceOptim]{EI}} from package DiceOptim, \code{\link[GPareto]{crit_EMI}}, \code{\link[GPareto]{crit_SUR}}, \code{\link[GPareto]{crit_SMS}}.
#' 
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom emoa hypervolume_indicator
#' @useDynLib GPareto, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @details
#' The computation of the analytical formula with two objectives is adapted from the Matlab source code by Michael Emmerich and Andre Deutz, LIACS, 
#'          Leiden University, 2010 available here :
#'          \url{http://liacs.leidenuniv.nl/~csmoda/code/HV_based_expected_improvement.zip}.
#' @references 
#' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State University, PhD thesis.  \cr \cr
#' M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
#' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
#' 
#' @examples
#' #---------------------------------------------------------------------------
#' # Expected Hypervolume Improvement surface associated with the "P1" problem at a 15 points design
#' #---------------------------------------------------------------------------
#' set.seed(25468)
#' library(DiceDesign)
#' 
#' n_var <- 2 
#' f_name <- "P1" 
#' n.grid <- 26
#' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
#' n_appr <- 15 
#' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
#' response.grid <- t(apply(design.grid, 1, f_name))
#' Front_Pareto <- t(nondominated_points(t(response.grid)))
#' mf1 <- km(~., design = design.grid, response = response.grid[,1])
#' mf2 <- km(~., design = design.grid, response = response.grid[,2])
#' 
#' EHI_grid <- crit_EHI(x = as.matrix(test.grid), model = list(mf1, mf2),
#'          critcontrol = list(refPoint = c(300,0)))
#' 
#' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
#'                matrix(EHI_grid, n.grid), main = "Expected Hypervolume Improvement",
#'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
#'                plot.axes = {axis(1); axis(2);
#'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
#'                             }
#'               )
crit_EHI <- function(x, model, paretoFront = NULL,
                     critcontrol = list(nb.samp = 50, seed = 42),
                     type = "UK"){
  nobj     <- length(model)
  if(nobj < 2){
    cat("Incorrect Number of objectives \n")
    return(NA)
  }
  
  if(nobj == 2){
    return(EHI_2d(x, model, critcontrol, type, paretoFront))
  } else {
    critcontrol$type <- "hypervolume"
    if (is.null(critcontrol$nb.samp)) critcontrol$nb.samp <- 50
    if (is.null(critcontrol$seed)) critcontrol$seed <- 42
    
    return(SAA_mEI(x = x, model = model, critcontrol = critcontrol, type = type, paretoFront = paretoFront))
  }
}

#' Parallel Multi-objective Expected Hypervolume Improvement with respect to the
#' current Pareto front, via Sample Average Approximation (SAA).
#' 
#' @title Batch Expected Hypervolume Improvement with m objectives
#' 
#' @param x a matrix representing the inputs (one per row) for which one wishes to calculate \code{EHI},
#' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
#' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, or any reference set of observations, 
#' @param critcontrol optional list with arguments:
#' \itemize{
#' \item \code{nb.samp} number of random samples from the posterior distribution (with more than two objectives),
#'        default to \code{50}, increasing gives more reliable results at the cost of longer computation time;
#' \item \code{seed} seed used for the random samples (with more than two objectives);
#' \item \code{refPoint} reference point for Hypervolume Expected Improvement;
#' \item \code{extendper} if no reference point \code{refPoint} is provided,
#'  for each objective it is fixed to the maximum over the Pareto front plus extendper times the range, 
#'  Default value to \code{0.2}, corresponding to \code{1.1} for a scaled objective with a Pareto front in \code{[0,1]^n.obj}.
#' } 
#'   
#' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
#'        has to be taken into account. 
#' @return The Expected Hypervolume Improvement at \code{x}.
#' @seealso \code{[DiceOptim]{EI}} from package DiceOptim, \code{\link[GPareto]{crit_EMI}}, \code{\link[GPareto]{crit_SUR}}, \code{\link[GPareto]{crit_SMS}}.
#' 
#' @export
#' @details
#' 
#' The batch EHI is computed by simulated \code{nb.samp} conditional simulation of the objectives on the candidate points,
#' before averaging the corresponding hypervolume improvement of each set of m simulations.
#' 
#' @references 
#' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State University, PhD thesis.  \cr \cr
#' 
#' @examples
#' #---------------------------------------------------------------------------
#' # Expected Hypervolume Improvement surface associated with the "P1" problem at a 15 points design
#' #---------------------------------------------------------------------------
#' set.seed(25468)
#' library(DiceDesign)
#' 
#' n_var <- 2 
#' f_name <- "P1" 
#' n.grid <- 26
#' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
#' test.grid <- as.matrix(test.grid)
#' n_appr <- 15
#' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
#' response.grid <- t(apply(design.grid, 1, f_name))
#' Front_Pareto <- t(nondominated_points(t(response.grid)))
#' mf1 <- km(~., design = design.grid, response = response.grid[,1])
#' mf2 <- km(~., design = design.grid, response = response.grid[,2])
#' 
#' refPoint <- c(300, 0)
#' EHI_grid <- crit_EHI(x = test.grid, model = list(mf1, mf2),
#'          critcontrol = list(refPoint = refPoint))
#' 
#' ## Create potential batches
#' q <- 3
#' ncb <- 500
#' Xbcands <- array(NA, dim = c(ncb, q, n_var))
#' for(i in 1:ncb) Xbcands[i,,] <- test.grid[sample(1:nrow(test.grid), q, prob = pmax(0,EHI_grid)),]
#' 
#' qEHI_grid <- apply(Xbcands, 1, crit_qEHI, model = list(mf1, mf2),
#'   critcontrol = list(refPoint = refPoint))
#'   
#' Xq <- Xbcands[which.max(qEHI_grid),,]
#' 
#' ## For further optimization (gradient may not be reliable)
#' # sol <- optim(as.vector(Xq), function(x, ...) crit_qEHI(matrix(x, q), ...),
#' #   model = list(mf1, mf2), lower = c(0,0), upper = c(1,1), method = "L-BFGS-B",
#' #   control = list(fnscale = -1), critcontrol = list(refPoint = refPoint, nb.samp = 10000))
#' # sol <- psoptim(as.vector(Xq), function(x, ...) crit_qEHI(matrix(x, q), ...),
#' #  model = list(mf1, mf2), lower = c(0,0), upper = c(1,1), 
#' #  critcontrol = list(refPoint = refPoint, nb.samp = 100),
#' #  control = list(fnscale = -1))
#' 
#' # Plot EHI surface and selected designs for parallel evaluation
#' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
#'                matrix(EHI_grid, n.grid), 
#'                main = "Expected Hypervolume Improvement surface and best 3-EHI designs",
#'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
#'                plot.axes = {axis(1); axis(2);
#'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
#'                             points(Xq, pch = 20, col = 2)
#' #                            points(matrix(sol$par, q), col = 4)
#'                             }
#'               )
crit_qEHI <- function(x, model, paretoFront = NULL,
                      critcontrol = list(nb.samp = 50, seed = 42),
                      type = "UK"){
  nobj     <- length(model)
  if(nobj < 2){
    cat("Incorrect Number of objectives \n")
    return(NA)
  }
  critcontrol$type <- "hypervolume"
  if (is.null(critcontrol$nb.samp)) critcontrol$nb.samp <- 50
  if (is.null(critcontrol$seed)) critcontrol$seed <- 42
  
  if(is.null(dim(x))) x <- matrix(x, ncol = model[[1]]@d)
  
  return(SAA_mEI(x = x, model = model, critcontrol = critcontrol, type = type,
                 paretoFront = paretoFront, batch = TRUE))
  
}


