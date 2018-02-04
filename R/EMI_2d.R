###### Semi-analytical formulas in the bi-objective case

## transpose paretoFront for now

EMI_2d <- function(x, model, type = "UK", critcontrol = NULL, paretoFront = NULL, r = 0){
  
  if(is.null(paretoFront)){
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    paretoFront <- t(nondominated_points(t(observations)))
  }
  
  if (is.unsorted(paretoFront[,1])){
    paretoFront <- paretoFront[order(paretoFront[,1]),]
  }
  
  ## A new x too close to the known observations could result in numerical problems
  if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    return(-1)
  }
  
  
  tmp <- data.frame(t(x))
  colnames(tmp) <- colnames(model[[1]]@X)
  # a <- predict(object = model[[1]], newdata = tmp, type = type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  # b <- predict(object = model[[2]], newdata = tmp, type = type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  # y <- c(a$mean, b$mean)
  # s <- c(a$sd, b$sd)
  
  pred <- predict_kms(model, newdata=tmp, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  y <- as.numeric(pred$mean)
  s <- as.numeric(pred$sd)
  
  res <- sum(apply(expand.grid(1:2, 1:nrow(paretoFront)), 1, EMI_ij, paretoFront = paretoFront, y = y, s = s, r = r))

  if(is.nan(res)) return(-1)
  
  # res <- 0
  # for(i in 1:2){
  #   for(j in 1:ncol(paretoFront)){
  #     res <- res + Int1(i, j, paretoFront, y, s, r) + Int2(i, j, paretoFront, y, s, r) + Int3(i, j, paretoFront, y, s, r)
  #   }
  # }
  return(res)
}

EMI_ij <- function(ij, paretoFront, y, s, r){
  return(Int1(ij[1], ij[2], paretoFront, y, s, r) + Int2(ij[1], ij[2], paretoFront, y, s, r) + Int3(ij[1], ij[2], paretoFront, y, s, r))
}


k_S <- function(a){
  if(a == 1){
    return(2)
  }
  if(a == 2){
    return(1)
  }
}


h_S <- function(i, j){
  if(i == 1){
    return(j - 1);
  }
  if(i == 2){
    return(j + 1);
  }
}

d_S <- function(i, j, P, y){
  if(i == 1 & j == (nrow(P)+1)){
    return(Inf)
  }
  if(i == 2 & j == 0){
    return(Inf)
  }
  return(P[j,i] - y[i])
}


u_S <- function(i, j, P, y, s){
  pnorm(d_S(i, j, P, y)/s[i])
}


v_S <- function(i, j, P, y, s, r){
  if(j == 0){
    return(Inf)
  }
  if(j ==(nrow(P)+1)){
    return(Inf)
  }
  return(y[i]/s[i]^2 + (P[j,2] - d_S(k_S(i), j, P, y) + r * s[k_S(i)] * y[i] / s[i]) / sqrt((1-r^2)*s[k_S(i)]^2)) 
}


q_S <- function(i, j, P, y, s, r){
  (1 - r^2) * s[k_S(i)]^2 * s[i]^2 / (s[i]^2 + s[k_S(i)]^2 - 2 * r * s[i] * s[k_S(i)])
}


Int1 <- function(i, j, P, y, s, r){
  return(
    s[i] * dnorm(d_S(i,j, P, y) / s[i]) * (
      pnorm((-d_S(k_S(i), j, P, y) + r * s[k_S(i)] * d_S(i,j, P, y) / s[i]) / sqrt((1 - r^2) * s[k_S(i)]^2))
      - pnorm((-d_S(k_S(i), h_S(i, j), P, y) + r * s[k_S(i)] * d_S(i,j, P, y) / s[i]) / sqrt((1 - r^2) * s[k_S(i)]^2))
    )
  )
}


Int2 <- function(i, j, P, y, s, r){
  aa <- q_S(i, j, P, y, s, r)
  
  a1 <- sqrt(aa) * (s[i] - s[k_S(i)] * r) / sqrt(2 * pi * (1 - r^2) * s[k_S(i)]^2)
  a2 <- exp(-1/2 * (y[i]^2 / s[i]^2 + (P[j,i] - d_S(k_S(i), j, P, y) + r * s[k_S(i)] * y[i] / s[i])^2 / sqrt(1 - r^2) * s[k_S(i)]^2)) * exp(1/2 * aa * v_S(i, j, P, y, s, r)^2) * pnorm((P[j,i] - aa * v_S(i, j, P, y, s, r)) / sqrt(aa))
  a3 <- exp(-1/2 * (y[i]^2 / s[i]^2 + (P[j,i] - d_S(k_S(i), h_S(i,j), P, y) + r * s[k_S(i)] * y[i] / s[i])^2 / sqrt(1-r^2)*s[k_S(i)]^2)) * exp(1/2 * aa * v_S(i, h_S(i, j), P, y, s, r)^2) * pnorm((P[j,i] - aa * v_S(i, h_S(i,j), P, y, s, r)) / sqrt(aa))
  
  if(is.nan(a2))
    a2 <- 0
  
  if(is.nan(a3)){
    a3 <- 0
  }
  
  return(a1 * (a2 - a3))
}


ftmp1 <- function(w, i, j, P, y, s, r){
  return(
    pnorm((d_S(i, j, P, y) - d_S(k_S(i), j, P, y) + (s[k_S(i)]*r - s[i])*qnorm(w))/sqrt((1-r^2)*s[k_S(i)]^2))
  )
}

ftmp2 <- function(w, i, j, P, y, s, r){
  return(
    pnorm((d_S(i, j, P, y) - d_S(k_S(i), h_S(i, j), P, y) + (s[k_S(i)]*r - s[i])*qnorm(w))/sqrt((1-r^2)*s[k_S(i)]^2))
  )
}

Int3 <- function(i, j, P, y, s, r){
  a1 <- d_S(i, j, P, y)
  b1 <- u_S(i, j, P, y, s)
  
  if(b1 < 1e-5){
    return(0)
  }
  
  a2 <- integrate(ftmp1, 0, b1, i=i, j=j, P=P, y=y, s=s, r=r, stop.on.error = FALSE)$value
  a3 <- integrate(ftmp2, 0, b1, i=i, j=j, P=P, y=y,s=s, r=r, stop.on.error = FALSE)$value
  
  return(a1 * (a2 - a3))
}


