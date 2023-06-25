# Function 1: Calculate gradient and Hessian matrix
gradient <- function(x, y, delta, beta, sigma, pi) {
  er <- (y - x %*% beta) / sigma
  exper <- exp(er)
  d.beta <- x * c(exper - delta) / pi / sigma #First Derivative of beta
  d.sig <- (er * exper - delta * er - delta) / pi / sigma
  return(cbind(d.sig, d.beta))
}

hessian <- function(x, y, delta, beta, sigma, pi) {
  er <- (y - x %*% beta) / sigma
  exper <- exp(er)
  d2.beta <- -t(x * c(exper)  / pi) %*% x / (sigma ^ 2)
  d2.sig <- sum((delta + 2 * er * (delta - exper) - (er ^ 2) 
                 * exper) / pi) / (sigma ^ 2)
  d2.sigbeta <- colSums((x * c(delta - exper * (1 + er)))) / (sigma ^ 2)
  return(cbind(c(d2.sig, d2.sigbeta), 
               rbind(t(d2.sigbeta), d2.beta)))
}

# Function 2: Calculate MLE for weibull regression model
#### initial values
#### set all betas except intercept as zero
osmacWei.est <- function(x, y, delta, pi, 
                      control = list(maxit = 1000, tol = 1e-5, 
                                     init = c(1, rep(0, dim(x)[2])),
                                     step = 0.8)) {
  r <- nrow(x)
  p <- ncol(x)
  beta <- control$init[-1]
  sigma <- control$init[1]
  msg <- "NA"
  for (i in 1L:control$maxit) {
    beta.old <- beta
    sigma.old <- sigma
    # update beta
    er <- (y - x %*% beta) / sigma
    exper <- exp(er)
    #First Derivative of beta
    d.beta <- colSums(x * c(exper - delta) / pi / sigma) * control$step
    #Second Derivative of beta
    d2.beta <- -t(x * c(exper)  / pi) %*% x / (sigma ^ 2)
    tryCatch({
      upd.beta <- NA
      upd.beta <- solve(d2.beta, -d.beta)
    },
    error = function(e) {
      cat("Error in the ", i, "'s loop: ", 
          conditionMessage(e),".\n",  sep = "")})
    if (is.na(upd.beta[1])) {
      msg <- "Algorithm not converges."
      beta <- sigma <- i <- NA
      return(list(scale = 0, para = 0, message = msg))
    }
    beta <- beta + upd.beta
    # update sigma
    er <- (y - x %*% beta) / sigma
    exper <- exp(er)
    # First Derivative of sigma
    d.sig <- sum((er * exper - delta * er - delta) / pi / sigma) * control$step
    # Second Derivative of sigma
    d2.sig <- sum((delta + 2 * er * (delta - exper) - (er ^ 2)
                   * exper) / pi) / (sigma ^ 2)
    upd.sig <- -d.sig / d2.sig
    sigma <- sigma + upd.sig
    tol.theta <- sqrt(sum(upd.sig ^ 2 + upd.beta ^ 2))
    if (tol.theta <= control$tol) {
      msg <- "Algorithm converges."
      break
    }
    if (i == control$maxit) {
      warning("Maxium iteration reached")
      msg <- "Algorithm not converges."
    }
  }
  return(list(scale = sigma, para = beta, message = msg))
}

# Function 3: Calculate SSP for different method.
osmacWei.ssp <- function (x, y, delta, n.pilot,
                          method, control = list(maxit = 1000)) {
  n <- nrow(x)
  if (method == "uniform") {
    return(list(ssps = rep(1 / n, n), 
                index.pilot = sample(n, n.pilot, TRUE)))
  }else{
    # Step 1: pilot estimate(use uniform ssp).
    for (i in 1L:control$maxit) {
      pr.pilot <- rep(1 / n, n.pilot)
      ind.pilot <- sample(n, n.pilot, T)
      x.pilot <- x[ind.pilot, ]
      y.pilot <- y[ind.pilot]
      delta.pilot <- delta[ind.pilot]
      mle.pilot <- osmacWei.est(x.pilot, y.pilot, delta.pilot, pr.pilot)
      if (mle.pilot$message == "Algorithm converges.") break
      if (i == control$maxit) warning("Maxium iteration reached 
                                    and you failed to get a pilot estimate.")
    }
    # Step 2: Calculating SSPs
    beta.pilot <- mle.pilot$para
    sig.pilot <- mle.pilot$scale
    g <- gradient(x, y, delta, beta.pilot, sig.pilot, 1)
    if (method == "optL") {
      g.norm <- sqrt(rowSums(g ^ 2))
      pr.mvc <- g.norm / sum(g.norm) * (1 - 0.2) + 0.2 / n
      return(list(ssps = pr.mvc, index.pilot = ind.pilot))
    }else if (method == "optA"){
      # Use pilot sample to calculate Hessian
      M <- hessian(x.pilot, y.pilot, delta.pilot, beta.pilot, sig.pilot, 1)
      M.inv <- solve(M)
      M.mse <- sqrt(colSums((M.inv %*% t(g)) ^ 2))
      pr.mmse <- M.mse / sum(M.mse) * (1 - 0.2) + 0.2 / n
      return(list(ssps = pr.mmse, index.pilot = ind.pilot))
    }
  }
}


# Function 4: Obtain estimated coefficients and its standard error for 
# selected subsampling method
osmacWei.fit <- function(x, y, delta, n.pilot, n.sample, method, se = TRUE,
                         pilot = TRUE) {
  if (method == "full"){
    ind <- c(1:length(y))
    sec.ssp <- 1
    est <- osmacWei.est(x[ind, ], y[ind], delta[ind], sec.ssp)
    if (est$message == "Algorithm not converges."){
      return("Fail to get converging results.")
    }
  }else{
    # Get subsampling probabilities
    SSP <- osmacWei.ssp(x, y, delta, n.pilot, n.sample, method)
    
    # Select subsample and calculate the estimator. Stop when converges.
    for (i in 1:1000) {
      ind.sample <- sample(length(y), n.sample, prob = SSP$ssps, replace = TRUE)
      if(pilot){
        ind <- c(ind.sample, SSP$index.pilot)
        sec.ssp <- c(SSP$ssps[ind.sample], rep((1 / length(y)), n.pilot))
      }else{
        ind <- ind.sample
        sec.ssp <- SSP$ssps[ind]
      }
      est <- osmacWei.est(x[ind, ], y[ind], delta[ind], sec.ssp)
      if (est$message == "Algorithm converges."){
        break
      }
    }
    if (i == 1000) {
      return("Fail to get converging results.")
    }
  }
  
  # Get standard error
  if(se) {
    g <- gradient(x[ind, ], y[ind], delta[ind], est$para, est$scale, sec.ssp) 
    M <- hessian(x[ind, ], y[ind], delta[ind], est$para, est$scale, sec.ssp)
    Vc <- t(g) %*% g 
    N <- solve(M)
    Vx <- N %*% Vc %*% N
    se <- sqrt(diag(Vx))
  }else se <- NA
  out <- cbind(c(est$scale, est$para), se)
  colnames(out) <- c("Est", "SE")
  rownames(out) <- c("Scale", colnames(x))
  return(out)
}

