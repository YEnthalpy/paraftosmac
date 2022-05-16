# Function 1: Function for generating mvnorm and mvt data
genData <- function(n, xdist, beta, sigma, crate, Sigma = 0, mu = 0, df = 0) {
  # Generate data when covariates follow multivariate normal distribution
  if (xdist == "mvnorm") {
    x <- cbind(rep(1, n), MASS::mvrnorm(n, mu, Sigma))
  }
  # Generate data when covariates follow multivariate t distribution
  if (xdist == "mvt") {
    x <- cbind(rep(1, n), mvtnorm::rmvt(n, sigma = Sigma * (df - 2) / df, 
                                        df = df))
  }
  stime <- log(rweibull(n, 1 / sigma, exp(x %*% beta)))
  shape <- seq(0, 5, 0.1)
  ctime <- log(sapply(shape, function(t){rweibull(n, shape = 1 / sigma, scale = t)}))
  etime <- sapply(c(1:ncol(ctime)), function(t){pmin(stime, ctime[, t])})  
  delta <- sapply(c(1:ncol(ctime)), function(t){as.numeric(stime <= etime[, t])})
  rates <- 1 - colSums(delta)/n
  cindex <- min(which(abs(crate - rates) == min(abs(crate - rates))))
  return(list(covariates = x, response = etime[, cindex], delta = delta[, cindex]))
}

# Funciton 2: Function for summarizing MSE
mySummary.mse <- function(sim, ent.fit) {
  coe.ent.fit <- c(ent.fit$scale, ent.fit$para)
  emp.mse <- rowMeans(colSums((sim[, c("coe.mVc", "coe.mMse", "coe.uni"), ] 
                               - coe.ent.fit) ^ 2))
  est.mse <- rowMeans(colSums(sim[, c("se.mVc", "se.mMse", "se.uni"), ] ^ 2))
  emp.bias <- colSums((rowMeans(sim[, c("coe.mVc", "coe.mMse", "coe.uni"), ], 
                                dim = 2) - coe.ent.fit) ^ 2)
  emp.var <- emp.mse - emp.bias
  out <- cbind(emp.mse, est.mse, emp.var)
  rownames(out) <- c("mVc", "mMse", "uni")
  return(out)
}

# Function 3: Function for summarizing coverage rate
mySummary.converage <- function(sim, ent.fit) {
  coe.ent.fit <- c(ent.fit$scale, ent.fit$para)
  beta1.true <- coe.ent.fit[3]
  ci.beta.1 <- t(sim["beta1", seq(14, 19, 1), ])
  converage <- rep(0, 3)
  names(converage) <- c("mVc", "mMse", "uni")
  for (i in 1 : nrow(ci.beta.1)) {
    for (j in 1 : 3) {
      if (ci.beta.1[i, 2 * j - 1] < beta1.true && 
          beta1.true <  ci.beta.1[i, 2 * j]) {
        converage[j] <- converage[j] + 1
      }
    }
  }
  return(converage/nrow(ci.beta.1))
}

# Function 4: Function for summarizing computation time
mySummary.time <- function(sim) {
  time.res <- sim[, "time.all", ][c(1, 2, 3), ]
  rownames(time.res) <- c("mVc", "mMse", "uni")
  rowSums(time.res)
}

# Function 5: Function for summarizing iteration
mySummart.ite <- function(sim) {
  it.res <- sim[, c("it.mVc", "it.mMse", "it.uni"), ][c(1:3), ,]
  out1 <- rowMeans(it.res, dim = 2)
  rownames(out1) <- c("it.all", "it.beta", "it.sig")
  it.pilot <- sim[, c("it.pilot.mVc", "it.pilot.mMse"),  ][c(1:4), ,]
  out2 <- rowMeans(it.pilot, dim = 2)
  rownames(out2) <- c("it.pilot", "it.all", "it.beta", "it.sig")
  return(list(out1, out2))
}

# Function 6: Function for summarizing all results together
mySummary.all <- function(sim, ent.fit) {
  out2 <- sapply(1:length(sim), function(n){
    mySummary.mse(sim[[n]], ent.fit)})
  out3 <- sapply(1:length(sim), function(n)
  {mySummary.converage(sim[[n]], ent.fit)})
  out4 <- sapply(1:length(sim), function(n)
  {mySummary.time(sim[[n]])})
  out <- rbind(out2, out3, out4)
  rownames(out) <- c("emp.mVc", "emp.mMse", "emp.uni", "est.mVc", "est.mMse", 
                     "est.uni", "var.mVc", "var.mMse", "var.uni",
                     "con.mVc", "con.mMse", "con.uni", "time.mVc",
                     "time,mMse", "time.uni")
  colnames(out) <- c(1000, 2000, 3000, 4000)
  return(out)
}

# Function 7: Function for one replicate
do1rep <- function(dat, n.pilot, n.sample) {
  x <- dat$covariates; y <- dat$response; delta <- dat$delta
  fit <- fitWeiReg.ssp(x, y, delta, n.pilot, n.sample)
  con.int <- matrix(NA, nrow = nrow(fit), ncol = 6)
  colnames(con.int) <- c("lower.mVc", "upper.mVc", "lower.mMse",
                         "upper.mMse", "lower.uni", "upper.uni")
  for (i in 1 : 3) {
    con.int[, 2 * i - 1] <- fit[, 4 * i - 3] - 1.96 * fit[, 4 * i - 2]
    con.int[, 2 * i] <- fit[, 4 * i - 3] + 1.96 * fit[, 4 * i - 2]
  }
  return(cbind(fit, con.int))
}

# Function 8: Function for several replicates
my.replicate <- function(n, dat, n.pilot, n.sample) {
  k0 <- 1; k1 <- 1
  out <- array(NA, dim = c(ncol(dat$covariates) + 1, 19, n))
  while(k0 <= n) {
    skip <- FALSE
    tryCatch(
      temp <- do1rep(dat, n.pilot, n.sample),
      error = function(e) skip <<- TRUE
    )
    k1 <- k1 + 1
    if (skip) next
    out[, , k0] <- temp
    k0 <- k0 + 1
  }
  rownames(out) <- rownames(temp)
  colnames(out) <- colnames(temp)
  return(out)
}

# Function 9: Function of getting simulation results for a given data set
simu <- function(n, xdist, beta, sigma, crate, n.sample, n.pilot, n.rep, 
                 Sigma, df = 0, mu = 0) {
  dat <- genData(n, xdist, beta, sigma, crate, Sigma, df = df, mu = mu)
  x <- dat$covariates; y <- dat$response; delta <- dat$delta
  sim <- lapply(n.sample, function(x){my.replicate(n.rep, dat, n.pilot, x)})
  names(sim) <- n.sample
  sim[["ent.fit"]] <- fitWeiReg(x, y, delta, rep(1/n, n))
  sim[["Data"]] <- dat
  return(sim)
}
