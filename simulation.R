source("functions.R")
source("summary.R")
library(ggplot2)
library(tikzDevice)

# 1. Initiate settings
n <- 100000; beta <- rep(0.01, 8); crate <- c(0.25, 0.50, 0.75)
n.sample <- c(1000, 2000, 3000, 4000); n.pilot <- 1000

Sigma <- matrix(data = 0.5, nrow = 7, ncol = 7); diag(Sigma) <- 1
sigma.1 <- 1; sigma.2 <- 2

set.seed(101)

# 2. Simulate for different data sets
## sigma = 1
mvt.1 <- lapply(crate, function(x){
  simu(n, xdist = "mvt", beta, sigma = sigma.1, x, n.sample, n.pilot, n.rep = 1000, 
       Sigma = Sigma, df = 3)})
mvn.1 <- lapply(crate, function(x){
  simu(n, xdist = "mvnorm", beta, sigma = sigma.1, x, n.sample, n.pilot, n.rep = 1000, 
       Sigma = Sigma, mu = rep(0, 7))})
## sigma = 2
mvt.2 <- lapply(crate, function(x){
  simu(n, xdist = "mvt", beta, sigma = sigma.2, x, n.sample, n.pilot, n.rep = 1000, 
       Sigma = Sigma, df = 3)})
mvn.2 <- lapply(crate, function(x){
  simu(n, xdist = "mvnorm", beta, sigma = sigma.2, x, n.sample, n.pilot, n.rep = 1000, 
       Sigma = Sigma, mu = rep(0, 7))})

# 3. Draw plots for the simulation study

names(mvt.2) <- c("0.25", "0.50", "0.75")
names(mvt.1) <- c("0.25", "0.50", "0.75")
names(mvn.2) <- c("0.25", "0.50", "0.75")
names(mvn.1) <- c("0.25", "0.50", "0.75")

expand.1 <- function(sim) {
  out <- data.frame()
  for (i in 1:length(sim)) {
    temp.1 <- my.expand(mySummary.all(sim[[i]][-c(5, 6)], 
                                      sim[[i]][[5]]))
    censoring.rate <- rep(names(sim)[i], nrow(temp.1))
    name <- deparse(substitute(sim))
    sig <- rep(unlist(strsplit(name, ""))[5], nrow(temp.1))
    Dist <- rep(name, nrow(temp.1))
    temp.2 <- cbind(temp.1, censoring.rate, Dist, sig)
    out <- rbind(out, temp.2)
  }
  return(out)
}

expand.2 <- function(sim){
  out <- data.frame()
  for (i in 1:length(sim)) {
    out <- rbind(out, eval(parse(text = paste(
      "expand.1(", sim[i], ")", sep = ""))))
  }
  return(out)
}

plot.simu <- function(sim.full) {
  # Compare MSE for different SSPs
  options(tikzLatexPackages 
          =c(getOption("tikzLatexPackages"),"\\usepackage{amsmath}"))
  tikz("plot2.tex", width = 12, height = 9)
  sim.mse1 <- sim.full[which(sim.full$Methods%in%c("emp.mVc", "emp.mMse", "emp.uni")), ]
  sim.mse1$Methods[which(sim.mse1$Methods=="emp.mVc")] <- "mVc"
  sim.mse1$Methods[which(sim.mse1$Methods=="emp.mMse")] <- "mMSE"
  sim.mse1$Methods[which(sim.mse1$Methods=="emp.uni")] <- "uniform"
  sim.mse1$r <- factor(sim.mse1$r, levels=unique(sim.mse1$r))
  dist.labs <- c(mvn.1 = "Normal", mvn.2 = "Normal", mvt.1 = "T5"
                 , mvt.2 = "T5")
  cr.labs <- c("$c_r=0.25$", "$c_r=0.50$", "$c_r=0.75$")
  names(cr.labs) <- c("0.25", "0.50", "0.75")
  plot1 <- ggplot(sim.mse1, aes(x = r, y = MSE)) + 
    geom_line(aes(group = Methods, colour = Methods), size = 1.1) +
    geom_point(aes(group = Methods, colour = Methods, shape = Methods), size = 4) +
    facet_grid(rows = censoring.rate ~ Dist, scales="free_y", labeller = 
                 labeller(Dist = dist.labs, censoring.rate = cr.labs)) +
    theme(text = element_text(size = 20), legend.position = "bottom",
          legend.title = element_blank()) +
    xlab("$r$")
  plot(plot1)
  dev.off()
  
  # Compare empirical and estimated MSE
  options(tikzLatexPackages 
          =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
  tikz("plot3.tex", width = 12, height = 9)
  sim.mse2 <- sim.full[which(sim.full$Methods%in%c("emp.mVc", "est.mVc")), ]
  sim.mse2$Methods[which(sim.mse2$Methods=="emp.mVc")] <- "Empirical"
  sim.mse2$Methods[which(sim.mse2$Methods=="est.mVc")] <- "Estimated"
  sim.mse2$r <- factor(sim.mse2$r, levels=unique(sim.mse2$r))
  plot2 <- ggplot(sim.mse2, aes(x = r, y = MSE)) + 
    geom_line(aes(group = Methods, colour = Methods), size = 1.1) +
    geom_point(aes(group = Methods, shape = Methods, colour = Methods), size = 4.5) +
    facet_grid(rows = censoring.rate ~ Dist, scales="free_y", labeller = 
                 labeller(Dist = dist.labs, censoring.rate = cr.labs)) +
    theme(text = element_text(size = 20), legend.position = "bottom",
          legend.title = element_blank()) +
    xlab("$r$")
  plot(plot2)
  dev.off()
  
  # plot for coverage rate
  options(tikzLatexPackages 
          =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
  tikz("plot4.tex", width = 12, height = 9)
  sim.cov <- sim.full[which(sim.full$Methods%in%c("con.mVc", "con.mMse", "con.uni")), ]
  sim.cov$Methods[which(sim.cov$Methods=="con.mVc")] <- "mVc"
  sim.cov$Methods[which(sim.cov$Methods=="con.mMse")] <- "mMSE"
  sim.cov$Methods[which(sim.cov$Methods=="con.uni")] <- "uniform"
  sim.cov$r <- factor(sim.cov$r, levels=unique(sim.cov$r))
  plot3 <- ggplot(sim.cov, aes(x = r, y = MSE)) + 
    geom_line(aes(group = Methods, colour = Methods), size = 1.1) +
    geom_point(aes(group = Methods, shape = Methods, colour = Methods), size = 4) +
    geom_hline(aes(yintercept = 0.95), 
               colour="#990000", linetype="dashed") +
    facet_grid(rows = censoring.rate ~ Dist, labeller = 
                 labeller(Dist = dist.labs, censoring.rate = cr.labs)) +
    ylim(0.8, 1) + 
    theme(text = element_text(size = 20), legend.position = "bottom",
          legend.title = element_blank()) +
    ylab("Coverage Rate") +
    xlab("$r$")
  plot(plot3)
  dev.off()
}

dat <- c("mvn.1", "mvn.2", "mvt.1", "mvt.2")
sim.full <- expand.2(dat)
sim.sig1 <- sim.full[which(sim.full$sig==1), ]
sim.sig2<- sim.full[which(sim.full$sig==2), ]
## Draw plots for sigma = 1 (shown in the main file)
plot.simu(sim.sig1)

## Draw plots for sigma = 2 (shown in the supplement file)
plot.simu(sim.sig2)
