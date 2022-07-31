load("lympDat.RData")
source("functions.R")
source("summary.R")

library(ggplot2)
library(ggpubr)
library(tikzDevice)



# 1. Clean the data
data2 <- lymp[, c("AGE_DX", "SEX", "Nonwhite", "DATE_yr", "survtime", 
                  "alldeath", "DATE_yr")]
## Delete missing values
data2 <- data2[-which(data2$AGE_DX==999), ]
## Because the unit of survival time is month
## Since the survival time cannot be 0 in numerical computation.
## We consider 0 survival time is rounded by 0.25. 
data2$survtime[which(data2$survtime==0)] <- 0.25
## Scale and center the covariates.
data2$AGE_DX <- scale(data2$AGE_DX, scale = TRUE, center = TRUE)
data2$DY <- scale(data2$DATE_yr, scale = TRUE, center = TRUE)
## Redefine categorical variables.
data2$Nonwhite <- ifelse(data2$Nonwhite=="TRUE", 1, 0)
data2$SEX <- ifelse(data2$SEX=="2", 1, 0)
## Construct interactions between Age and categorical variables.
data2$AN <- data2$AGE_DX * data2$Nonwhite
data2$AS <- data2$AGE_DX * data2$SEX
## Define covariate, response and censoring indicator.

real.dat <- list(covariates = as.matrix(cbind(rep(1, nrow(data2)), 
                               data2[, c("AGE_DX", "SEX", "Nonwhite",
                                         "AN", "AS", "DY")])), 
             response = data2$survtime, delta = data2$alldeath)
colnames(real.dat$covariates) <- c("Intercept", "Age", "Male", "Nonwhite",
                                   "Age*Nonwhite", "Age*Male", "Year")

# 2. Calaculate the true MLE
time.full <- system.time(
  full.fit <- osmacWei.fit(real.dat$covariates, 
                           log(real.dat$response), 
                           real.dat$delta, 0, 0, "full"))[1]

# 3. Calculate estimations based on different SSPs.
## Pilot subsample size r0 = 500.
## Step twp subsample size %in% c(1000, 2000, 3000, 4000).

RNGkind("L'Ecuyer-CMRG")
set.seed(928)
seeds <- list(.Random.seed)
for (i in 2:1000) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}

for (i in seq(1000, 4000, 1000)) {
  for (k in c("mMSE", "mVc", "uniform")) {
    eval(parse(text = paste0("real", ".", k, ".", i,'<- 
    lapply(seeds, function(seed, dat, n.s, meds){
        assign(".Random.seed", seed, envir = .GlobalEnv)
        t <- system.time(res <- osmacWei.fit(dat$covariates, log(dat$response), 
                                           dat$delta, 500, n.s, meds))
        return(cbind(res, "Time" = c(t[1], rep(NA, (nrow(res) - 1)))))},
      dat = real.dat, meds = k, n.s = i)')))
  }
}

# 4. Derive empirical standard error of full 
#    data MLE by Bootstrap method.

boots <- lapply(seeds, function(seed, dat){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  osmacWei.fit(dat$covariates, log(dat$response), 
               dat$delta, 0, length(dat$response), "uniform")},
  dat = real.dat)

# 4. Draw plots

## 4.1. Summary the results

full.coe <- full.fit[, 1]
tmp.boots <- array(unlist(boots), 
                   dim = c(nrow(boots[[1]]), ncol(boots[[1]]), 
                           length(boots)))
emp.mse.boots <- sum(rowMeans((tmp.boots[, 1, ] - full.coe) ^ 2))
est.mse.boots <- sum(rowMeans(tmp.boots[, 2, ] ^ 2))
summ.dat <- data.frame(value = rep(c(emp.mse.boots, est.mse.boots), 4), 
                       type = rep(c("Estimated", "Empirical"), 4),
                       method = rep("Full", 8), 
                       n.s = c(rep(1000, 2), rep(2000, 2), 
                               rep(3000, 2), rep(4000, 2)))

for(method in c("mMSE", "mVc", "uniform")){
  for(n.s in seq(1000, 4000, 1000)){
    summ.dat <- rbind(summ.dat, summ(method, n.s, real = TRUE))
  }
}

summ.dat$method <- factor(summ.dat$method, 
                          levels = c("mMSE", "mVc", "uniform", "Full"))

## 4.2. MSE of differents SSPs

mse.methods <- ggplot(summ.dat[summ.dat$type == "Empirical", ], 
       mapping = aes(x = n.s, y = value, color = method, group = method)) +
  geom_point(aes(shape = method), size = 4.5) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(c(0, 0.23)) +
  ylab("MSE") +
  xlab("$r$")

## 4.3. Compare empirical and estimated MSE

se.types <- ggplot(summ.dat[summ.dat$method == "mVc", ], 
       mapping = aes(x = n.s, y = value, color = type, group = type)) +
  geom_point(aes(shape = type), size = 4.5) +
  geom_line(size = 1.1) +
  theme(text = element_text(size = 25), legend.position = "bottom",
        legend.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylim(c(0, 0.23)) +
  ylab("") +
  xlab("$r$")

tikz('mse_realdata.tex', standAlone = TRUE ,width = 13,height = 5)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}"),
        "\\usepackage{amsfonts}")
ggarrange(mse.methods, NULL, se.types, common.legend = FALSE,
          legend = "bottom", widths = c(1, 0.05, 1), nrow = 1)

dev.off()
tools::texi2dvi('mse_realdata.tex', pdf = T)

