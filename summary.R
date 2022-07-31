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
  colnames(x) <- c("Intercept", paste0("beta", 1:(ncol(x) - 1)))
  stime <- rweibull(n, 1 / sigma, exp(x %*% beta))
  ctime <- sapply(seq(0.01, 10, 0.05), function(t){runif(n, max = t)})
  etime <- sapply(c(1:ncol(ctime)), function(t){pmin(stime, ctime[, t])})  
  delta <- sapply(c(1:ncol(ctime)), function(t){as.numeric(stime <= etime[, t])})
  rates <- 1 - colSums(delta)/n
  cindex <- min(which(abs(crate - rates) == min(abs(crate - rates))))
  return(list(covariates = x, response = etime[, cindex], 
              delta = delta[, cindex], xdist = xdist, crate = crate,
              sigma = sigma))
}


# Funciton 2: Function for summarizing empirical and estiamted MSE

summ <- function(method, n.s, xdist = NA, crate = NA, sigma = NA,
                 real = FALSE) {
  if(real){
    full.coe <- full.fit[, 1]
    res <- eval(parse(text = paste0("real", ".", method, ".", n.s)))
  }else{
    full <- coe.full[[which(dat.xdist == xdist & dat.crate == crate 
                            & dat.sig == sigma)]]
    full.coe <- full[, 1]
    res <- eval(parse(text = paste0(xdist, crate, ".", sigma,
                                    ".", method, ".", n.s)))
  }
  tmp <- array(unlist(res), dim = c(nrow(res[[1]]), ncol(res[[1]]), length(res)))
  emp.mse <- sum(rowMeans((tmp[, 1, ] - full.coe) ^ 2))
  est.mse <- sum(rowMeans(tmp[, 2, ] ^ 2))
  atime <- sum(tmp[1, 3, ])
  if(real){
    return(data.frame(value = c(est.mse, emp.mse, atime), 
                      type = c("Estimated", "Empirical", "Time"),
                      method = rep(method, 3), n.s = rep(n.s, 3)))
  }else{
    return(data.frame(value = c(est.mse, emp.mse, atime), 
                      type = c("est.mse", "emp.mse", "Time"),
                      method = rep(method, 3), xdist = rep(xdist, 3),
                      crate = rep(crate, 3), sigma = rep(sigma, 3),
                      n.s = rep(n.s, 3)))
  }
}


