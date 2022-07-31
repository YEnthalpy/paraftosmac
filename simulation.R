source("functions.R")
source("summary.R")

library(ggplot2)
library(tikzDevice)

# 1. Generate Simulation Data

Sigma <- matrix(data = 0.5, nrow = 7, ncol = 7); diag(Sigma) <- 1

set.seed(101)
dat.1 <- mapply(
  function(xdist, crate, sigma){
    data = genData(10000, xdist, rep(0.01, 8), sigma, crate, 
                   Sigma = Sigma, mu = rep(0, 7), df = 5)}, 
  xdist = c(rep("mvnorm", 6), rep("mvt", 6)),
  crate = rep(c(0.25, 0.5, 0.75), 4), 
  sigma = rep(c(rep(1, 3), rep(2, 3)), 2),
  SIMPLIFY = FALSE, USE.NAMES = FALSE)

# 2. Simulation

## 2.1. Generate seeds for 1000 replicates

RNGkind("L'Ecuyer-CMRG")
set.seed(928)
seeds <- list(.Random.seed)
for (i in 2:1000) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}

## 2.2. Do the calculation

for (i in 1:length(dat.1)) {
  dat <- dat.1[[i]]
  for (j in c(1000, 2000, 3000, 4000)) {
    for (k in c("mMSE", "mVc", "uniform")) {
      eval(parse(text = paste0(dat$xdist, dat$crate, ".", 
                               dat$sigma,".", k, ".", j,'<- 
      lapply(seeds, function(seed, dat, n.s, meds){
          assign(".Random.seed", seed, envir = .GlobalEnv)
          x <- dat$covariates
          y <- log(dat$response)
          delta <- dat$delta
          t <- system.time(res <- osmacWei.fit(x, y, delta, 1000, n.s, meds))
          return(cbind(res, "Time" = c(t[1], rep(NA, (nrow(res) - 1)))))},
        dat = dat, n.s = j, meds = k)')))
    }
  }
}


## 2.3. Get results based on full Data.

coe.full <- lapply(dat.1, function(dat){
  x <- dat$covariates
  y <- log(dat$response)
  delta <- dat$delta
  osmacWei.fit(x, y, delta, 0, 0, "full")})

# 3. Draw plots

## 3.1. Summary the results

dat.xdist <- sapply(dat.1, function(dat){dat$xdist})
dat.crate <- sapply(dat.1, function(dat){dat$crate})
dat.sig <- sapply(dat.1, function(dat){dat$sigma})


summ.dat <- data.frame()
for(xdist in c("mvnorm", "mvt")){
  for(crate in c(0.25, 0.5, 0.75)){
    for(sigma in c(1, 2)){
      for(method in c("mMSE", "mVc", "uniform")){
        for(n.s in seq(1000, 4000, 1000)){
          summ.dat <- rbind(summ.dat, 
                            summ(method, n.s, xdist = xdist, 
                                 crate = crate, sigma = sigma))
        }
      }
    }
  }
}


cr.labs <- c("0.25" = "$c_r=0.25$", "0.50" = "$c_r=0.50$", 
             "0.75" = "$c_r=0.75$")
dist.labs <- c("mvnorm" = "Normal", "mvt" = "T5")

## 3.2. Plots for the situation: \sigma = 1

### 3.2.1. Compare MSE for different SSPs

tikz('mse_sig1.tex', standAlone = TRUE ,width = 12, height = 9)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}"),
        "\\usepackage{amsfonts}")

ggplot(summ.dat[summ.dat$sigma == 1 & summ.dat$type == "emp.mse", ], 
       mapping = aes(x = n.s, y = value, color = method, group = method)) +
  facet_grid(crate~xdist, scales = "free_y", 
             labeller = labeller(xdist = dist.labs, crate = cr.labs)) +
  geom_point(aes(shape = method), size = 4) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank()) +
  ylab("MSE") +
  xlab("$r$")
dev.off()
tools::texi2dvi('mse_sig1.tex', pdf = T)

### 3.2.2. Compare empirical and estimated MSE

tikz('mse_est_emp_sig1.tex', standAlone = TRUE ,width = 12, height = 9)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}"),
        "\\usepackage{amsfonts}")
ggplot(summ.dat[summ.dat$sigma == 1 & summ.dat$type %in% c("emp.mse", "est.mse")
                & summ.dat$method == "mVc", ], 
       mapping = aes(x = n.s, y = value, color = type, group = type)) +
  facet_grid(crate~xdist, scales = "free_y", 
             labeller = labeller(xdist = dist.labs, crate = cr.labs)) +
  geom_point(aes(shape = type), size = 4.5) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank()) +
  ylab("MSE") +
  xlab("$r$")
dev.off()
tools::texi2dvi('mse_est_emp_sig1.tex', pdf = T)

## 3.3. Plots for the situation: \sigma = 2

### 3.3.1. Compare MSE for different SSPs

tikz('mse_sig2_SuppInfo.tex', standAlone = TRUE ,width = 12, height = 9)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}"),
        "\\usepackage{amsfonts}")

ggplot(summ.dat[summ.dat$sigma == 2 & summ.dat$type == "emp.mse", ], 
       mapping = aes(x = n.s, y = value, color = method, group = method)) +
  facet_grid(crate~xdist, scales = "free_y", 
             labeller = labeller(xdist = dist.labs, crate = cr.labs)) +
  geom_point(aes(shape = method), size = 4) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank()) +
  ylab("MSE") +
  xlab("$r$")
dev.off()
tools::texi2dvi('mse_sig2_SuppInfo.tex', pdf = T)

### 3.3.2. Compare empirical and estimated MSE

tikz('mse_est_emp_sig2_SuppInfo.tex', standAlone = TRUE ,width = 12, height = 9)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}"),
        "\\usepackage{amsfonts}")
ggplot(summ.dat[summ.dat$sigma == 2 & summ.dat$type %in% c("emp.mse", "est.mse")
                & summ.dat$method == "mVc", ], 
       mapping = aes(x = n.s, y = value, color = type, group = type)) +
  facet_grid(crate~xdist, scales = "free_y", 
             labeller = labeller(xdist = dist.labs, crate = cr.labs)) +
  geom_point(aes(shape = type), size = 4.5) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank()) +
  ylab("MSE") +
  xlab("$r$")
dev.off()
tools::texi2dvi('mse_est_emp_sig2_SuppInfo.tex', pdf = T)



