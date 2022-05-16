load("lympDat.RData")
source("functions.R")
source("summary.R")
library(ggplot2)
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
## Define covariate `x`, response `y` and censoring indicator `delta`.
x <- as.matrix(cbind(rep(1, nrow(data2)), data2[, c("AGE_DX", "SEX", "Nonwhite",
                                                    "AN", "AS", "DY")]))
y <- log(data2$survtime); delta <- data2$alldeath
dat <- list(covariates = x, response = y, delta = delta)

# 2. Calaculate the true MLE
ent.fit <- fitWeiReg(x, y, delta, 1)

# Calculate estiamted MSE by different SSPs over 1000 replicates.
## Pilot subsample size r0 = 500.
## Step twp subsample size %in% c(1000, 2000, 3000, 4000).
set.seed(928)
res.real.500 <- lapply(c(1000, 2000, 3000, 4000), function(n){
  my.replicate(1000, dat, 500, n)})

# 3. Derive empirical standard error of full data MLE by Bootstrap method.
boots <- function(x, y, delta) {
  n <- nrow(x)
  index <- sample(n, n, replace = TRUE)
  temp <- fitWeiReg(x[index, ], y[index], delta[index], 1)
  out <- c(temp[[1]], temp[[2]])
  return(out)
}

boots.rep <- function (n, x, y, delta) 
{
  k0 <- 1
  out <- matrix(data = NA, nrow = n, ncol = ncol(x) + 1)
  while (k0 <= n) {
    skip <- FALSE
    tryCatch(temp <- boots(x, y, delta), error = function(e) skip <<- TRUE)
    if (skip) 
      next
    out[k0, ] <- temp
    k0 <- k0 + 1
  }
  colnames(out) <- c("Scale", "Intercept", "Age", "Sex", "Nonwhite", 
                     "Age*Nonwhite", "Age*Sex", "DY")
  return(out)
}

real.boots <- boots.rep(1000, x, y, delta)

# 4. Draw plot for showing empirical MSE derived by different methods
summ.real <- mySummary.all(res.real.500, ent.fit)
plot.mse <- my.expand(summ.real)
fit <- c(ent.fit[[1]], ent.fit[[2]])
mse.boots <- sum(rowMeans((t(real.boots) - fit) ^ 2))
var.boots <- sum(rowMeans((t(real.boots) - colMeans(real.boots)) ^ 2))

## Compare MSE for different SSPs
options(tikzLatexPackages 
        =c(getOption("tikzLatexPackages"),"\\usepackage{amsmath}"))
tikz("plot5.tex", width = 8, height = 4)
plot.dat1 <- plot.mse[which(plot.mse$Methods %in% 
                              c("emp.mVc", "emp.mMse", "emp.uni", 
                                "var.mVc", "var.mMse", "var.uni")), ]
plot.dat1$type[which(plot.dat1$Methods %in% 
                       c("var.mVc", "var.mMse", "var.uni"))] <- "Var"
plot.dat1$type[which(plot.dat1$Methods %in% 
                       c("emp.mVc", "emp.mMse", "emp.uni"))] <- "MSE"
plot.dat1$Methods[which(plot.dat1$Methods %in% 
                          c("emp.mVc", "var.mVc"))] <- "mVc"
plot.dat1$Methods[which(plot.dat1$Methods %in% 
                          c("emp.mMse", "var.mMse"))] <- "mMSE"
plot.dat1$Methods[which(plot.dat1$Methods %in% 
                          c("emp.uni", "var.uni"))] <- "uniform"
plot.dat1$r <- factor(plot.dat1$r, levels = unique(plot.dat1$r))
temp1 <- data.frame(MSE = rep(mse.boots, 4), r = seq(1000, 4000, 1000), 
                    Methods = rep("zBootstrap", 4), type = rep("MSE", 4))
temp2 <- data.frame(MSE = rep(var.boots, 4), r = seq(1000, 4000, 1000), 
                    Methods = rep("zBootstrap", 4), type = rep("Var", 4))
plot.dat1 <- rbind(plot.dat1, temp1, temp2)

plot1 <- ggplot(plot.dat1, aes(x = r, y = MSE)) + 
  geom_line(aes(group = Methods, colour = Methods), size = 1.1) + 
  geom_point(aes(group = Methods, colour = Methods, shape = Methods), size = 4.5) + 
  facet_wrap(~type) + 
  theme(text = element_text(size = 18), legend.position = "bottom", 
        legend.title = element_blank()) +
  xlab("$r$")
plot(plot1)
dev.off()

# 5. Generating the table


