##' Executes \code{nsteps} iterations of multi-objective EGO methods to objects of class \code{\link[DiceKriging]{km}}.
##' At each step, kriging models are re-estimated (including covariance parameters re-estimation)
##'  based on the initial design points plus the points visited during all previous iterations;
##'  then a new point is obtained by maximizing one of the four multi-objective Expected Improvement criteria available. 
##'  Handles noiseless and noisy objective functions.
##' @title Sequential multi-objective Expected Improvement maximization and model re-estimation,
##'  with a number of iterations fixed in advance by the user
##' @details Extension of the function \code{\link[DiceOptim]{EGO.nsteps}} for multi-objective optimization.\cr
##' Available infill criteria with \code{crit} are: \cr
##' \itemize{
##' \item Expected Hypervolume Improvement (\code{EHI}) \code{\link[GPareto]{crit_EHI}},
##' \item SMS criterion (\code{SMS}) \code{\link[GPareto]{crit_SMS}},
##' \item Expected Maximin Improvement (\code{EMI}) \code{\link[GPareto]{crit_EMI}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[GPareto]{crit_SUR}}.
##' }
##' 
##' Depending on the selected criterion, parameters such as reference point for \code{SMS} and \code{EHI} or arguments for \code{\link[GPareto]{integration_design_optim}} 
##' with \code{SUR} can be given with \code{critcontrol}.
##' Also options for \code{\link[GPareto]{checkPredict}} are available.
##' More precisions are given in the corresponding help pages.\cr
##' 
##' The \code{reinterpolation=TRUE} setting can be used to handle noisy objective functions. It works with all criteria and is the recommended option. 
##' If \code{reinterpolation=FALSE} and \code{noise.var!=NULL}, the criteria are used based on a "denoised" Pareto front.
##' 
##' If \code{noise.var="given_by_fn"}, \code{fn} must return a list of two vectors, the first being the objective functions and the second 
##' the corresponding noise variances (see examples).
##' 
##' Display of results and various post-processings are available with \code{\link[GPareto]{plotGPareto}}.  
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param fn the multi-objective function to be minimized (vectorial output), found by a call to \code{\link[base]{match.fun}},
##' @param cheapfn optional additional fast-to-evaluate objective function (handled next with class \code{\link[GPareto]{fastfun}}), which does not need a kriging model, 
##' handled by a call to \code{\link[base]{match.fun}},
##' @param crit choice of multi-objective improvement function: "\code{SMS}", "\code{EHI}", "\code{EMI}" or "\code{SUR}",
##' see details below,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation has to be taken into account,
##'  see \code{\link[DiceKriging]{km}}
##' @param cov.reestim optional boolean specifying if the kriging hyperparameters should be re-estimated at each iteration,
##' @param noise.var noise variance (of the objective functions). Either \code{NULL} (noiseless objectives), a scalar (constant noise, identical for all objectives), 
##' a vector (constant noise, different for each objective) or 
##' a function (type closure) with vectorial output (variable noise, different for each objective). Alternatively, set \code{noise.var="given_by_fn"}, see details. 
##' If not provided but km models are based on noisy observations, \code{noise.var} is taken as the average of \code{model@noise.var}.
##' @param reinterpolation Boolean: for noisy problems, indicates whether a reinterpolation model is used, see details,
##' @param critcontrol optional list of parameters for criterion \code{crit}, see details,
##' @param optimcontrol an optional list of control parameters for optimization of the selected infill criterion:
##' "\code{method}" can be set to "\code{discrete}", "\code{pso}", "\code{genoud}" or a user defined method name (passed to \code{\link[base]{match.fun}}). For "\code{discrete}", 
##' a matrix \code{candidate.points} must be given.
##' For "\code{pso}" and "\code{genoud}", specific parameters to the chosen method can also be specified  (see \code{\link[rgenoud]{genoud}} and \code{\link[pso]{psoptim}}).
##' A user defined method must have arguments like the default \code{\link[stats]{optim}} method, i.e. \code{par}, \code{fn}, \code{lower}, \code{upper}, 
##' \code{...} and eventually \code{control}.\cr
##'  A trace \code{trace} argument is available, it can be set to \code{0} to suppress all messages, to \code{1} (default) for displaying the optimization progresses,
##'  and \code{>1} for the highest level of details. 
##' 
##' @param ncores number of CPU available (> 1 makes mean parallel \code{TRUE}). Only used with \code{discrete} optimization for now.
##' @param ... additional parameters to be given to the objective \code{fn}.
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: a data frame representing the additional points visited during the algorithm,}
##' \item{\code{values}}{: a data frame representing the response values at the points given in \code{par},}
##' \item{\code{nsteps}}{: an integer representing the desired number of iterations (given in argument),}
##' \item{\code{lastmodel}}{: a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted.}
##' \item{\code{observations.denoised}}{:  if \code{noise.var!=NULL}, a matrix representing the mean values of the \code{\link[DiceKriging]{km}} models at observation points.}
##' If a problem occurs during either model updates or criterion maximization, the last working model and corresponding values are returned.
##' }
##' 
##' @references 
##' M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
##' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
##' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}, 25(6), 1265-1280\cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, 718-727, Springer, Berlin. \cr \cr
##' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State university, PhD thesis. 
##' V. Picheny and D. Ginsbourger (2013), \emph{Noisy kriging-based optimization methods: A unified implementation within the DiceOptim package}, 
##' \emph{Computational Statistics & Data Analysis}, 71: 1035-1053.
##' 
##' @importFrom stats runif pnorm qnorm
##' 
##' @examples
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' ################################################
##' # NOISELESS PROBLEMS
##' ################################################
##' d <- 2 
##' fname <- ZDT3
##' n.grid <- 21
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 15 
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' Front_Pareto <- t(nondominated_points(t(response.grid)))
##' 
##' mf1 <- km(~1, design = design.grid, response = response.grid[, 1], lower=c(.1,.1))
##' mf2 <- km(~., design = design.grid, response = response.grid[, 2], lower=c(.1,.1))
##' model <- list(mf1, mf2)
##' 
##' nsteps <- 2
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' 
##' # Optimization 1: EHI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO1 <- GParetoptim(model = model, fn = fname, crit = "EHI", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO1$par)
##' print(omEGO1$values)
##' 
##' \dontrun{
##' nsteps <- 10
##' # Optimization 2: SMS with discrete search
##' optimcontrol <- list(method = "discrete", candidate.points = test.grid)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO2 <- GParetoptim(model = model, fn = fname, crit = "SMS", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO2$par)
##' print(omEGO2$values)
##' 
##' # Optimization 3: SUR with genoud
##' optimcontrol <- list(method = "genoud", pop.size = 20, max.generations = 10)
##' critcontrol <- list(distrib = "SUR", n.points = 100)
##' omEGO3 <- GParetoptim(model = model, fn = fname, crit = "SUR", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO3$par)
##' print(omEGO3$values)
##' 
##' # Optimization 4: EMI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(nbsamp = 200)
##' omEGO4 <- GParetoptim(model = model, fn = fname, crit = "EMI", nsteps = nsteps,
##'                      lower = lower, upper = upper, optimcontrol = optimcontrol)
##' print(omEGO4$par)
##' print(omEGO4$values)
##' 
##' # graphics
##' sol.grid <- apply(expand.grid(seq(0, 1, length.out = 100),
##'                               seq(0, 1, length.out = 100)), 1, fname)
##' plot(t(sol.grid), pch = 20, col = rgb(0, 0, 0, 0.05), xlim = c(0, 1),
##'      ylim = c(-2, 10), xlab = expression(f[1]), ylab = expression(f[2]))
##' plotGPareto(res = omEGO1, add = TRUE,
##'             control = list(pch = 20, col = "blue", PF.pch = 17,
##'                            PF.points.col = "blue", PF.line.col = "blue"))
##' text(omEGO1$values[,1], omEGO1$values[,2], labels = 1:nsteps, pos = 3, col = "blue")
##' plotGPareto(res = omEGO2, add = TRUE,
##'             control = list(pch = 20, col = "green", PF.pch = 17,
##'                            PF.points.col = "green", PF.line.col = "green"))
##' text(omEGO2$values[,1], omEGO2$values[,2], labels = 1:nsteps, pos = 3, col = "green")
##' plotGPareto(res = omEGO3, add = TRUE,
##'             control = list(pch = 20, col = "red", PF.pch = 17,
##'                            PF.points.col = "red", PF.line.col = "red"))
##' text(omEGO3$values[,1], omEGO3$values[,2], labels = 1:nsteps, pos = 3, col = "red") 
##' plotGPareto(res = omEGO4, add = TRUE,
##'             control = list(pch = 20, col = "orange", PF.pch = 17,
##'                            PF.points.col = "orange", PF.line.col = "orange"))
##' text(omEGO4$values[,1], omEGO4$values[,2], labels = 1:nsteps, pos = 3, col = "orange")
##' points(response.grid[,1], response.grid[,2], col = "black", pch = 20)
##' legend("topright", c("EHI", "SMS", "SUR", "EMI"), col = c("blue", "green", "red", "orange"),
##'  pch = rep(17,4))
##'  
##'  
##' # Post-processing
##' plotGPareto(res = omEGO1, UQ_PF = TRUE, UQ_PS = TRUE, UQ_dens = TRUE)
##' 
##' ################################################
##' # NOISY PROBLEMS
##' ################################################
##' set.seed(25468)
##' library(DiceDesign)
##' d <- 2 
##' nsteps <- 3
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(refPoint = c(1, 10))
##' 
##' n.grid <- 21
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' n.init <- 30
##' design <- maximinESE_LHS(lhsDesign(n.init, d, seed = 42)$design)$design
##' 
##' fit.models <- function(u) km(~., design = design, response = response[, u],
##'                              noise.var=design.noise.var[,u])
##' 
##' # Test 1: EHI, constant noise.var
##' noise.var <- c(0.1, 0.2)
##' funnoise1 <- function(x) {ZDT3(x) + sqrt(noise.var)*rnorm(n=d)}
##' response <- t(apply(design, 1, funnoise1))
##' design.noise.var <- matrix(rep(noise.var, n.init), ncol=d, byrow=TRUE)
##' model <- lapply(1:d, fit.models)
##' 
##' omEGO1 <- GParetoptim(model = model, fn = funnoise1, crit = "EHI", nsteps = nsteps,
##'                       lower = lower, upper = upper, critcontrol = critcontrol,
##'                       reinterpolation=TRUE, noise.var=noise.var, optimcontrol = optimcontrol)
##' plotGPareto(omEGO1)
##' 
##' # Test 2: EMI, noise.var given by fn
##' funnoise2 <- function(x) {list(ZDT3(x) + sqrt(0.05 + abs(0.1*x))*rnorm(n=d), 0.05 + abs(0.1*x))}
##' temp <- funnoise2(design)
##' response <- temp[[1]]
##' design.noise.var <- temp[[2]]
##' model <- lapply(1:d, fit.models)
##' 
##' omEGO2 <- GParetoptim(model = model, fn = funnoise2, crit = "EMI", nsteps = nsteps,
##'                       lower = lower, upper = upper, critcontrol = critcontrol,
##'                       reinterpolation=TRUE, noise.var="given_by_fn", optimcontrol = optimcontrol)
##' plotGPareto(omEGO2)
##' 
##' # Test 3: SMS, functional noise.var
##' funnoise3 <- function(x) {ZDT3(x) + sqrt(0.025 + abs(0.05*x))*rnorm(n=d)}
##' noise.var <- function(x) return(0.025 + abs(0.05*x))
##' response <- t(apply(design, 1, funnoise3))
##' design.noise.var <- t(apply(design, 1, noise.var))
##' model <- lapply(1:d, fit.models)
##' 
##' omEGO3 <- GParetoptim(model = model, fn = funnoise3, crit = "SMS", nsteps = nsteps,
##'                            lower = lower, upper = upper, critcontrol = critcontrol,
##'                            reinterpolation=TRUE, noise.var=noise.var, optimcontrol = optimcontrol)
##' plotGPareto(omEGO3)
##' 
##' # Test 4: SUR, fastfun, constant noise.var
##' noise.var <- 0.1
##' funnoise4 <- function(x) {ZDT3(x)[1] + sqrt(noise.var)*rnorm(1)}
##' cheapfn <- function(x) ZDT3(x)[2]
##' response <- apply(design, 1, funnoise4)
##' design.noise.var <- rep(noise.var, n.init)
##' model <- list(km(~., design = design, response = response, noise.var=design.noise.var))
##' 
##' omEGO4 <- GParetoptim(model = model, fn = funnoise4, cheapfn = cheapfn, crit = "SUR", 
##'                       nsteps = nsteps, lower = lower, upper = upper, critcontrol = critcontrol,
##'                       reinterpolation=TRUE, noise.var=noise.var, optimcontrol = optimcontrol)
##'  plotGPareto(omEGO4)
##'                             
##'  # Test 5: EMI, fastfun, noise.var given by fn
##'  funnoise5 <- function(x) {
##'    if (is.null(dim(x))) x <- matrix(x, nrow=1)
##'    list(apply(x, 1, ZDT3)[1,] + sqrt(abs(0.05*x[,1]))*rnorm(nrow(x)), abs(0.05*x[,1]))
##'  }
##'  
##'  cheapfn <- function(x) {
##'    if (is.null(dim(x))) x <- matrix(x, nrow=1)
##'    apply(x, 1, ZDT3)[2,]
##'  }
##'  
##'  temp <- funnoise5(design)
##'  response <- temp[[1]]
##'  design.noise.var <- temp[[2]]
##'  model <- list(km(~., design = design, response = response, noise.var=design.noise.var))
##'  
##'  omEGO5 <- GParetoptim(model = model, fn = funnoise5, cheapfn = cheapfn, crit = "EMI", 
##'                        nsteps = nsteps, lower = lower, upper = upper, critcontrol = critcontrol,
##'                        reinterpolation=TRUE, noise.var="given_by_fn", optimcontrol = optimcontrol)
##'  plotGPareto(omEGO5)               
##'  
##'  # Test 6: EHI, fastfun, functional noise.var
##'  noise.var <- 0.1
##'  funnoise6 <- function(x) {ZDT3(x)[1] + sqrt(abs(0.1*x[1]))*rnorm(1)}
##'  noise.var <- function(x) return(abs(0.1*x[1]))
##'  cheapfn <- function(x) ZDT3(x)[2]
##'  response <- apply(design, 1, funnoise6)
##'  design.noise.var <- t(apply(design, 1, noise.var))
##'  model <- list(km(~., design = design, response = response, noise.var=design.noise.var))
##'  
##'  omEGO6 <- GParetoptim(model = model, fn = funnoise6, cheapfn = cheapfn, crit = "EMI", 
##'                        nsteps = nsteps, lower = lower, upper = upper, critcontrol = critcontrol,
##'                        reinterpolation=TRUE, noise.var=noise.var, optimcontrol = optimcontrol)
##'  plotGPareto(omEGO6)    
##' }

GParetoptim <- function (model, fn, ..., cheapfn = NULL, crit = "SMS", nsteps, lower, upper, type = "UK", cov.reestim=TRUE,
                         critcontrol = NULL, noise.var = NULL, reinterpolation = NULL,
                         optimcontrol = list(method = "genoud", trace = 1), ncores = 1){
  
  n <- nrow(model[[1]]@X)
  d <- model[[1]]@d
  
  if(is.null(optimcontrol$method)) optimcontrol$method <- "genoud"
  if(is.null(critcontrol$threshold)) critcontrol$threshold <- 1e-4
  if(is.null(critcontrol$distance)) critcontrol$distance <- "euclidean"
  if(is.null(optimcontrol$trace)) optimcontrol$trace <- 1
  # if(is.null(critcontrol)) critcontrol <- list()
  
  ## Check that no model is noisy
  if(is.null(noise.var)) {
    noise.var <- colSums(matrix(Reduce(cbind, lapply(model, slot, "noise.var")), ncol=length(model)))
    if(sum(noise.var)==0) noise.var <- NULL
  }
  
  if(is.null(noise.var)) reinterpolation <- FALSE
  if(is.null(reinterpolation)) reinterpolation <- TRUE
  
  if (is.null(noise.var) && reinterpolation) {
    cat("The reinterpolation option only works on noisy problems \n")
    reinterpolation <- FALSE
  }
  
  fn <- match.fun(fn)
  
  if(is.null(optimcontrol$method))
    optimcontrol$method <- "genoud"
  
  if (optimcontrol$method=="discrete") {
    if (is.null(optimcontrol$candidate.points)) {
      cat("Error in optimcontrol setting: candidate.points must be given for method=discrete \n")
      return(NULL)
    }
  }
  
  # Build fastfun if necessary
  if (!is.null(cheapfn)) {
    cheapfn <- match.fun(cheapfn)
    fastobs <- apply(model[[1]]@X, 1, cheapfn)
    fastmod <- fastfun(fn = cheapfn, design = model[[1]]@X, response = fastobs)
    model[[length(model)+1]] <- fastmod
  } else {
    Y.new.cheap <- NULL
    noise.new.cheap <- NULL
  }
  n.obj <- length(model)
  
  if (typeof(noise.var) == "double" && length(noise.var)==1) noise.var <- rep(noise.var, n.obj)
  
  if (length(model) < 2 ){
    cat("Error in model definition: 'model' must be a list of 2 or more km models \n")
    return(NULL)
  }
  
  if (length(model) >= 3 && crit=="EHI"){
    cat("Analytical Hypervolume EI only works with 2 objectives; SAA approximation used. \n")
  }
  if(length(model) > 3 && crit == "SUR"){
    cat("SUR is available for 2 or 3 objectives \n")
    return(NULL)
  }
  
  # Regroup all observations
  observations <- Reduce(cbind, lapply(model, slot, "y"))
  paretoFront <- unique(t(nondominated_points(t(observations))))
  
  if (reinterpolation) {
    # Reinterpolation models
    build.optim.model <- function(model) {
      # if (!.hasSlot(cheapfn, "covariance")){
      if(class(model) == "fastfun"){
        return(model)
      } else {
        observations.denoised <- predict(model, newdata=model@X, checkNames=FALSE, type=type)$mean
        
        uniqueX <- !duplicated(model@X)
        optim.model <- try(km(formula=model@trend.formula, design=model@X[uniqueX,], response=observations.denoised[uniqueX],
                              covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
                              coef.var=model@covariance@sd2,control=model@control), silent = TRUE)
        # If km crashes: add nugget to the interpolating model
        if(typeof(optim.model) == "character") {
          optim.model <- try(km(formula=model@trend.formula, design=model@X, response=observations.denoised[uniqueX],
                                covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
                                coef.var=model@covariance@sd2, nugget=1e-8,control=model@control), silent = TRUE)
        }
        if (typeof(optim.model) == "character")  return(NULL)
        else                                    return(optim.model)
      }
    }
    reinterp.model <- lapply(model, FUN=build.optim.model)
    if (any(unlist(lapply(reinterp.model, typeof))=="NULL")) {
      cat("Unable to build an initial reinterpolating model \n")
      return(NULL)
    }
    
    # Denoised observations
    observations.denoised <- Reduce(cbind, lapply(reinterp.model, slot, "y"))
    paretoFront <- unique(t(nondominated_points(t(observations.denoised))))
    
  } else {
    if(!is.null(noise.var)) { # noisy case
      # Denoised observations
      pred <- predict_kms(model, newdata=model[[1]]@X, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
      observations.denoised <- t(pred$mean)
      paretoFront <- unique(t(nondominated_points(t(observations.denoised))))
    }
  }
  
  if(optimcontrol$trace > 0){
    cat("----------------------------\n")
    cat("Starting optimization with : \n The criterion", crit, "\n The solver",  optimcontrol$method, "\n")
    cat("----------------------------\n")
    cat("Ite / Crit / New x / New y \n")
  }
  
  noRef <- FALSE #TRUE if a refPoint is provided
  if((crit == "SMS" | crit =="EHI") & is.null(critcontrol$refPoint)){
    if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
    noRef <- TRUE
    # estimatedRef <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed
    PF_range <- apply(paretoFront, 2, range)
    estimatedRef <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
    critcontrol$refPoint <- estimatedRef 
    if((crit == "SMS" | crit =="EHI") & optimcontrol$trace > 0)
      cat("No refPoint provided, ", signif(critcontrol$refPoint, 3), "used \n")
  }
  
  ## Change of seed for optimization unless set for genoud (for now by setting default values)
  ## NOTE: might be interesting to modify the defaults in crit_optimizer
  if(is.null(optimcontrol$unif.seed))
    changeSeed1 <- TRUE
  if(is.null(optimcontrol$int.seed))
    changeSeed2 <- TRUE
  
  #### Main loop starts here ################
  for (i in 1:nsteps) {
    # Denoised observations + optim model
    if (i>1) {
      
      ## Compute current Pareto front
      paretoFront <- unique(t(nondominated_points(t(observations))))
      
      if (reinterpolation) {
        reinterp.model <- lapply(model, FUN=build.optim.model)
        if (any(unlist(lapply(reinterp.model, typeof))=="NULL")) {
          if(optimcontrol$trace > 0){
            cat("Unable tobuild a reinterpolating model at iteration ", i, "- optimization stopped \n")
            cat("Last model returned \n")
          }
          return(list(par=model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE], values=NULL, nsteps = i-1, lastmodel = model, 
                      observations.denoised=observations.denoised))
        }
        observations.denoised <- Reduce(cbind, lapply(reinterp.model, slot, "y"))
        ## Compute current denoised Pareto front
        paretoFront <- unique(t(nondominated_points(t(observations.denoised)))) 
      } else {
        if (!is.null(noise.var)) {
          pred <- predict_kms(model, newdata=model[[1]]@X, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
          observations.denoised <- t(pred$mean)
          ## Compute current denoised Pareto front
          paretoFront <- unique(t(nondominated_points(t(observations.denoised))))
        }
      }
    }
    
    ## Change the refPoint if none is provided, if necessary
    if(noRef & (crit == "SMS" | crit =="EHI")){
      PF_range <- apply(paretoFront, 2, range)
      if(any(estimatedRef != matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj))){
        estimatedRef <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
        critcontrol$refPoint <- estimatedRef
        if(optimcontrol$trace > 0)
          cat("refPoint changed, ", signif(critcontrol$refPoint, 3), "used \n")
      }
    }
    
    ## Change the seeds for genoud to avoid selecting always the same initial values
    if(optimcontrol$method == "genoud" & changeSeed1){
      optimcontrol$unif.seed <- sample.int(1e9, 1)
    }
    
    if(optimcontrol$method == "genoud" & changeSeed2){
      optimcontrol$int.seed <- sample.int(1e9, 1)
    }
    
    # observations removed (could be reintegrated)
    if(reinterpolation){
      sol <- try(crit_optimizer(crit = crit, model = reinterp.model, lower = lower, upper = upper, 
                                optimcontrol = optimcontrol, type = type, paretoFront = paretoFront, 
                                critcontrol = c(critcontrol, "nsteps.remaining" = nsteps-i), ncores = ncores))
    }else{
      sol <- try(crit_optimizer(crit = crit, model = model, lower = lower, upper = upper, 
                                optimcontrol = optimcontrol, type = type, paretoFront = paretoFront, 
                                critcontrol = c(critcontrol, "nsteps.remaining" = nsteps-i), ncores = ncores))
    }
    
    if (typeof(sol) == "character") {
      if(optimcontrol$trace > 0){
        cat("Unable to maximize criterion at iteration ", i, "- optimization stopped \n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model[[1]]@X[-n,, drop = FALSE]
        values <- observations[-n, , drop = FALSE]
      }
      
      if (is.null(noise.var)) {
        return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
      } else {
        return(list(par=par, values=values, nsteps = i-1, lastmodel = model, observations.denoised=observations.denoised))
      }
    }
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fn(as.numeric(sol$par), ...))
    
    if (is.null(noise.var)) {
      newnoise.var <- NULL
    } else {
      if (typeof(noise.var) == "closure") {
        newnoise.var <- noise.var(X.new)
      } else if (typeof(noise.var) == "double") {
        newnoise.var <- noise.var
      } else if (noise.var =="given_by_fn") {
        newnoise.var <- Y.new[[2]]
        Y.new <- Y.new[[1]]
      }
      
      ## To avoid problems with update.km with 0 variance
      newnoise.var <- pmax(sqrt(.Machine$double.eps), newnoise.var)
    }
    
    if (!is.null(cheapfn)) {
      Y.new.cheap <- try(cheapfn(as.numeric(sol$par)))
      noise.new.cheap <- rep(0, length(Y.new.cheap))
    }
    
    if (typeof(Y.new) == "character" || (!is.null(cheapfn) && typeof(Y.new.cheap) == "character")) {
      if(optimcontrol$trace > 0){
        cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", X.new, "\n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model[[1]]@X[-n,, drop=FALSE]
        values <- observations[-n, , drop = FALSE]
      }
      if (is.null(noise.var)) {
        return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
      } else {
        return(list(par=par, values=values, nsteps = i-1, lastmodel = model, observations.denoised=observations.denoised))
      }
      
    }
    
    Y.new <- c(Y.new, Y.new.cheap)
    
    if(optimcontrol$trace > 0) cat( i, "/", signif(sol$val,3), "/", signif(X.new,3), "/", signif(Y.new,3), "\n", sep = "\t")
    
    # Remove new observation from integration points if discrete case is used
    if (optimcontrol$method=="discrete") {
      optimcontrol$candidate.points <- optimcontrol$candidate.points[-sol$index,,drop=FALSE]
    }
    
    # Update models
    observations <- rbind(observations, Y.new)
    
    my.update <- function(u, model) {
      try(update(object = model[[u]], newX = X.new, newy=Y.new[u], newX.alreadyExist=FALSE,
                 newnoise.var = newnoise.var[u], 
                 cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
    }
    newmodel <- lapply(1:length(model), my.update, model = model)
    
    for (j in 1:n.obj) {
      if (typeof(newmodel[[j]]) == "character" && cov.reestim) {
        if(optimcontrol$trace > 0)
          cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", j, "\n")
        newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newnoise.var = newnoise.var[j], 
                                    newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      }
      
      if (typeof(newmodel[[j]]) == "character" && is.null(noise.var)) {
        if(optimcontrol$trace > 0)
          cat("Update with old hyperparameter values failed for model ", j, ", try with nugget estimation \n")
        if(length(model[[j]]@control) == 0) tmpcontrol <- list(control = list(trace = FALSE)) else tmpcontrol <- model[[j]]@control
        newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newnoise.var = NULL, 
                                    newX.alreadyExist=FALSE, cov.reestim = TRUE, nugget.reestim = TRUE,
                                    kmcontrol = tmpcontrol), silent = TRUE)
        # after running with old hyperparameters the slot control is NULL, hence trace is set to TRUE
      }
      
      if (typeof(newmodel[[j]]) == "character") {
        if(optimcontrol$trace > 0){
          cat("Unable to udpate kriging model ", j, " at iteration", i-1, "- optimization stopped \n")
          cat("lastmodel ", j, " is the model at iteration", i-1, "\n")
          cat("par and values contain the ",i , "th observation \n \n")
        }
        # allX.new <- rep(NA, d)
        if (i > 1) allX.new <- rbind(model[[1]]@X[(n+1):(n+i-1),, drop=FALSE], X.new) else allX.new <- X.new
        
        if (is.null(noise.var)) {
          return(list(
            par    = allX.new,
            values = observations[(n+1):(n+i),, drop=FALSE],
            nsteps = i, 
            lastmodel = model))
        } else {
          return(list(
            par    = allX.new,
            values = observations[(n+1):(n+i),, drop=FALSE],
            nsteps = i, 
            lastmodel = model,
            observations.denoised = observations.denoised))
        }
        
      } 
    }
    
    ## If update successful for all models
    model <- newmodel
    
  }
  if (!is.null(noise.var)) {
    pred <- predict_kms(model, newdata=model[[1]]@X, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    observations.denoised <- t(pred$mean)
  }
  
  if(optimcontrol$trace > 0) cat("\n")
  #### End of main loop ################
  
  res <- list(
    par=model[[1]]@X[(n+1):(n+nsteps),, drop=FALSE], 
    values=observations[(n+1):(n+nsteps),, drop=FALSE], 
    nsteps=nsteps, 
    lastmodel=model)
  
  if (!is.null(noise.var)) res <- c(res, list(observations.denoised=observations.denoised))
  
  return(res)
}
