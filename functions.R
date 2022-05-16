# Function 1: Calculating gradient and Hessian matrix
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
fitWeiReg <- function(x, y, delta, pi, 
                      control = list(maxit = 1000, tol = 1e-5, 
                                     init = c(1, rep(0, dim(x)[2]))), ...) {
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
    d.beta <- colSums(x * c(exper - delta) / pi / sigma)
    #Second Derivative of beta
    d2.beta <- -t(x * c(exper)  / pi) %*% x / (sigma ^ 2)
    tryCatch({
      upd.beta <- NA
      upd.beta <- solve(d2.beta, -d.beta) * 0.8
    },
    error = function(e) {
      cat("Error in the ", i, "'s loop: ", conditionMessage(e),".\n",  sep = "")})
    if (is.na(upd.beta[1])) {
      msg <- "Algorithm not converges. "
      beta <- sigma <- i <- NA
      return(list(scale = 0, para = 0, message = msg, iteration = i))
    }
    beta <- beta + upd.beta
    # update sigma
    er <- (y - x %*% beta) / sigma
    exper <- exp(er)
    #First Derivative of sigma
    d.sig <- sum((er * exper - delta * er - delta) / pi / sigma)
    #Second Derivative of sigma
    d2.sig <- sum((delta + 2 * er * (delta - exper) - (er ^ 2)
                   * exper) / pi) / (sigma ^ 2)
    upd.sig <- -d.sig / d2.sig * 0.8
    sigma <- sigma + upd.sig
    tol.theta <- sqrt(sum(upd.sig ^ 2 + upd.beta ^ 2))
    if (tol.theta <= control$tol) {
      msg <- "Algorithm converges."
      break
    }
    if (i == control$maxit) warning("Maxium iteration reached")
  }
  return(list(scale = sigma, para = beta, 
              message = msg, iteration = i))
}

# Function 3: Calculate SSP for different method.
weibull.ssp <- function (x, y, delta, n.pilot, n.sample, 
                         method, control = list(maxit = 1000)) {
  n <- nrow(x)
  if (method == "uniform") {
    pr.uni <- rep(1 / n, n)
    return(list(pr = pr.uni, index.pilot = sample(n, n.pilot, TRUE), 
                it.pilot = rep(NA, 2)))
  }else{
    # Step 1: pilot estimate(use uniform ssp)
    for (i in 1L:control$maxit) {
      pr.pilot <- rep(1 / n, n.pilot)
      ind.pilot <- sample(n, n.pilot, T)
      x.pilot <- x[ind.pilot, ]
      y.pilot <- y[ind.pilot]
      delta.pilot <- delta[ind.pilot]
      mle.pilot <- fitWeiReg(x.pilot, y.pilot, delta.pilot, pr.pilot)
      if (mle.pilot[[3]] == "Algorithm converges.") break
      if (i == control$maxit) warning("Maxium iteration reached 
                                    and you failed to get a pilot estimate.")
    }
    # Step 2: Calculating SSPs
    beta.pilot <- mle.pilot$para
    sig.pilot <- mle.pilot$scale
    g <- gradient(x, y, delta, beta.pilot, sig.pilot, 1)
    if (method == "mVc") {
      g.norm <- sqrt(rowSums(g ^ 2))
      pr.mvc <- g.norm / sum(g.norm) * (1 - 0.2) + 0.2 / n
      return(list(pr = pr.mvc, index.pilot = ind.pilot, 
                  it.pilot = c(i, mle.pilot$iteration)))
    }else{
      # Use pilot sample to calculate Hessian
      M <- hessian(x.pilot, y.pilot, delta.pilot, beta.pilot, sig.pilot, 1)
      M.inv <- solve(M)
      M.mse <- sqrt(colSums((M.inv %*% t(g)) ^ 2))
      pr.mmse <- M.mse / sum(M.mse) * (1 - 0.2) + 0.2 / n
      return(list(pr = pr.mmse, index.pilot = ind.pilot, 
                  it.pilot = c(i, mle.pilot$iteration)))
    }
  }
}

# Function 4: Generate the subsample (including pilot sample)
gen.subsample <- function(x, y, delta, n.sample, pr){
  n <- length(y)
  pi <- pr$pr; ind.pilot <- pr$index.pilot
  ind.sample <- sample(length(y), n.sample, prob = pi, replace = TRUE)
  ind <- c(ind.sample, ind.pilot)
  pi.ssp <- c(pi[ind.sample], rep(1 / n, length(ind.pilot)))
  return(cbind(pi.ssp, delta[ind], y[ind], x[ind, ]))
}

# Function 5: Calculate the standard error
se.ssp <- function(x, y, delta, beta, sigma, pr) {
  g <- gradient(x, y, delta, beta, sigma, pr) 
  M <- hessian(x, y, delta, beta, sigma, pr)
  Vc <- t(g) %*% g 
  N <- solve(M)
  Vx <- N %*% Vc %*% N
  return(sqrt(diag(Vx)))
}

# Function 6: Obtain estimated coefficients and its standard error for 
# selected subsampling method
esti.coe <- function(x, y, delta, n.pilot, n.sample, method) {
  n <- n.sample + n.pilot
  t1 <- system.time(pr <- weibull.ssp(x, y, delta, n.pilot, n.sample, method))[1]
  t2 <- system.time(data.ssp <- gen.subsample(x, y, delta, n.sample, pr))[1]
  t3 <- system.time(ssp.res <- fitWeiReg(data.ssp[, -c(1:3)], data.ssp[, 3], 
                                         data.ssp[, 2], data.ssp[, 1]))[1]
  se <- se.ssp(data.ssp[, -c(1:3)], data.ssp[, 3], data.ssp[, 2], 
               ssp.res$para, ssp.res$scale, data.ssp[, 1])
  coe <- c(ssp.res$scale, ssp.res$para)
  time1 <- c(t1 + t2 + t3, rep(0, length(coe) - 1))
  it.main <- c(ssp.res$iteration, rep(0, length(coe) - 1))
  it.pilot <- c(pr$it.pilot, rep(0, length(coe) - 2))
  return(cbind(coe, se, it.main, it.pilot, time1))
}

# Function 7: calculate estimated parameters 
fitWeiReg.ssp <- function(x, y, delta, n.pilot, n.sample) {
  n <- length(y)
  p <- dim(x)[2] + 1
  # result of mVc method
  res.mVc <- esti.coe(x, y, delta, n.pilot, n.sample, "mVc")
  # result of mMse method
  res.mMse <- esti.coe(x, y, delta, n.pilot, n.sample, "mMse")
  # result of uniform method
  res.uni <- esti.coe(x, y, delta, n.pilot, n.sample, "uniform")
  res.all <- cbind(res.mVc[, -5], res.mMse[, -5], res.uni[, -5])
  colnames(res.all) <- c("coe.mVc", "se.mVc", "it.mVc", "it.pilot.mVc", 
                         "coe.mMse", "se.mMse", "it.mMse", "it.pilot.mMse", 
                         "coe.uni", "se.uni", "it.uni", "it.pilot.uni")
  r.name <- c("sigma", "intercept")
  for (i in 1 : (p - 2)) {
    r.name <- c(r.name, paste("beta", i, sep = ""))
  }
  rownames(res.all) <- r.name
  time.all <- c(res.mVc[1, 5], res.mMse[1, 5], res.uni[1, 5], 
                rep(0, nrow(res.all) - 3))
  out <- cbind(res.all, time.all)
  return(out)
}
