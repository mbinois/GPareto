## ' Sample Average Approximation for several multi-objective Expected Improvement with respect to the
## ' current Pareto front. To avoid numerical instabilities, the new point is evaluated only if it is not too close to an existing observation.
## ' 
## ' @title SAA MO Expected Improvement with m objectives
## ' 
## ' @param x a vector representing the input for which one wishes to calculate EHI,
## ' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
## ' @param critcontrol List with six arguments: 
## '        \code{nb_samp} number of random samples from the posterior distribution;
## '        \code{seed} seed used for the random samples;
## '        \code{type} choice of multi-objective improvement function: "maximin" or "hypervolume";
## '        \code{refPoint} (optional) reference point for Hypervolume Expected Improvement,
## '        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} and \code{distance} are available to avoid numerical issues occuring when adding points too close to the existing ones.
## ' @param type "SK" or "UK" (by default), depending whether uncertainty related to trend estimation 
## '        has to be taken into account. 
## ' @param paretoFront matrix corresponding to the Pareto Front (one output per column).
## ' @return The SAA hypervolume improvement at \bold{x}.
## ' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim from which \code{EHI} is an extension.
## ' @export
## ' @importFrom MASS mvrnorm
## ' @references Svenson, J. D., & Santner, T. J. (2010). Multiobjective Optimization of Expensive Black-Box
## ' Functions via Expected Maximin Improvement. Technical Report, 
## ' 
## ' @examples
## ' #---------------------------------------------------------------------------
## ' # SAAmEI surface associated with the "P1" problem at a 15 points design
## ' #---------------------------------------------------------------------------
## ' \donttest{
## ' set.seed(25468)
## ' library(DiceDesign)
## ' 
## ' n_var <- 2 
## ' f_name <- "P1" 
## ' n.grid <- 101
## ' test.grid <- expand.grid(seq(0, 1,, n.grid), seq(0, 1,, n.grid))
## ' n_appr <- 15 
## ' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
## ' response.grid <- t(apply(design.grid, 1, f_name))
## ' Front_Pareto <- t(nondominated_points(t(response.grid)))
## ' mf1 <- km(~., design = design.grid, response = response.grid[,1])
## ' mf2 <- km(~., design = design.grid, response = response.grid[,2])
## ' 
## ' SAAmEI_grid <- apply(test.grid, 1, SAA_mEI, model = list(mf1, mf2),
## '                      critcontrol = list(nb_samp = 20, type = "maximin"))
## ' 
## ' filled.contour(seq(0, 1,, n.grid), seq(0, 1,, n.grid), matrix(SAAmEI_grid, n.grid),
## '                main = "Expected Maximin Improvement", xlab = expression(x[1]),
## '                ylab = expression(x[2]), color = terrain.colors, nlevels = 50,
## '                plot.axes = {axis(1); axis(2);
## '                             points(design.grid[,1], design.grid[,2],pch = 21,bg = "white")
## '                             }
## '               )
## ' }

#' Beta test: support for batch EHI
#' @noRd
SAA_mEI <- function(x, model, batch = FALSE,
                    critcontrol = list(nb.samp=100, seed = 42, type = "maximin", refPoint = NULL),
                    type = "UK", paretoFront = NULL){
  
  if (is.null(critcontrol$seed)) critcontrol$seed <- 42
  if (is.null(critcontrol$type)) critcontrol$type <-"maximin"
  if (is.null(critcontrol$nb.samp)) critcontrol$nb.samp <- 100
  
  n.obj <- length(model)
  d <- model[[1]]@d
  #x.new <- matrix(x, 1, d)
  if (!is.matrix(x)) x <- matrix(x, 1, d)
  n.candidates <- nrow(x)
  
  if(critcontrol$type == "hypervolume" & batch & n.candidates < 2) warning("Batch with one element considered.")
  
  if(is.null(paretoFront)){
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    paretoFront <- t(nondominated_points(t(observations)))
  }
  
  refPoint <- critcontrol$refPoint
  
  nb.samp <- critcontrol$nb.samp
  if(is.null(nb.samp)){
    nb.samp <- 50
  }
  
  seed <- critcontrol$seed
  if(is.null(seed)){
    seed <- 42
  }
  
  if(critcontrol$type == "hypervolume"){
    if (is.null(refPoint)){
      if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.1
      # refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
      PF_range <- apply(paretoFront, 2, range)
      refPoint <- matrix(PF_range[2,] + (PF_range[2,] - PF_range[1,]) * critcontrol$extendper , 1, n.obj)
      cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
    }
  }
  
  # mu    <- rep(NaN, n.obj)
  # sigma <- rep(NaN, n.obj)
  # for (i in 1:n.obj){    
  #   pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  #   mu[i]    <- pred$mean
  #   sigma[i] <- pred$sd
  # }
  pred <- predict_kms(model, newdata=x, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = batch)
  mu <- t(pred$mean)
  sigma <- t(pred$sd)
  
  ## A new x too close to the known observations could result in numerical problems
  if(!batch){
    check <- checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)
    idxOk <- which(!check)
    Res <- rep(-1,n.candidates)
    if(length(idxOk) == 0) return(Res)
  }
  
  # Set seed to have a deterministic function to optimize
  # see http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else
    oldseed <- NULL
  
  set.seed(seed)
  
  if(batch){
    Samples <- array(0, dim = c(nb.samp, n.candidates, n.obj))
    for (j in 1:n.obj) Samples[,,j] <- mvrnorm(n = nb.samp, mu[,j], 
                                               pred$cov[[j]] + diag(sqrt(.Machine$double.eps),n.candidates)) #, tol=sqrt(.Machine$double.eps))
  }else{
    Samples <- matrix(0,nb.samp*length(idxOk),n.obj)
    cpt <- 1
    for (i in idxOk){
      Samples[(1+(cpt-1)*nb.samp):(cpt*nb.samp),] <- mvrnorm(n = nb.samp, mu[i,], diag(sigma[i,]))
      cpt <- cpt+1
    }
  }
  
  
  if (!is.null(oldseed)) 
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)
  
  if(batch){
    if(critcontrol$type == "hypervolume"){
      H0 <- dominated_hypervolume(points = t(paretoFront), ref = refPoint)
      ImprovementSamples <- apply(Samples, 1, Hypervolume_improvement, front = paretoFront, refPoint = refPoint, H0 = H0)
    }else{stop("Batch EMI not supported.")}
    return(mean(ImprovementSamples))
    
  }else{
    if(critcontrol$type == "hypervolume"){
      H0 <- dominated_hypervolume(points = t(paretoFront), ref = refPoint)
      ImprovementSamples <- Hypervolume_improvement_vec(points = Samples, front = paretoFront, refPoint = refPoint, HypInitial = H0)
    } 
    else ImprovementSamples <- apply(Samples, 1, Maximin_Improvement, front = paretoFront, refPoint = refPoint)
    
    for (i in 1:length(idxOk)) Res[idxOk[i]] <- mean(ImprovementSamples[(1+(i-1)*nb.samp):(i*nb.samp)])
    
    return(Res)
  }
}

Hypervolume_improvement <- function(pts, front, refPoint, H0 = NULL){
  if(is.null(H0)) H0 <- dominated_hypervolume(points = t(front), ref = refPoint)
  # if(!is_dominated(t(rbind(points,front)))[1:nrow(points)]){
  Hi <- dominated_hypervolume(t(rbind(pts,front)), refPoint)
  # }
  return(Hi - H0)
}

contribution_to_front <- function(point, front, refPoint){
  return(dominated_hypervolume(t(rbind(point,front)),refPoint))
}

Hypervolume_improvement_vec <- function(points, front, refPoint, HypInitial = NULL){
  Hi <- rep(0,nrow(points))
  if(is.null(HypInitial)) HypInitial <- dominated_hypervolume(t(front),refPoint)
  idxND <- which(nonDomSet(points,front))
  Hi[idxND] <- apply(points[idxND,,drop=FALSE], 1, contribution_to_front, front = front, refPoint = refPoint) - HypInitial
  #for (i in idxND) Hi[i] <- dominated_hypervolume(t(rbind(points[i,],front)),refPoint)
  #Hi <- Hi-HypInitial
  #for (i in idxND) Hi[i] <- hypervolume_contribution(t(rbind(points[i,],front)),refPoint)[1]
  return(Hi)
}

# ref_point added to have the same arguments than Hypervolume_improvement
Maximin_Improvement <- function(point, front, refPoint){
  Em <- 0
  if(!is_dominated(t(rbind(point,front)))[1]){
    Em <- -Inf
    for(i in 1:nrow(front)){
      # for(j in 1:ncol(front)){
      tmp <- min(Inf, point - front[i,])
      # }
      Em <- max(Em, tmp)
    }
  }
  return(-Em)
}





