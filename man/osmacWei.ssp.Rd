\name{osmacWei.ssp}
\alias{osmacWei.ssp}
\title{Calculate the Optimal Subsample Probabilities}
\usage{
osmacWei.ssp(x, y, delta, n.pilot,
             method, control = list(maxit = 1000))
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{n.pilot}{the input pilot subsample size.}

\item{method}{the type of optimal subsampling probabilities (SSPs). `optA` for the
                `optA` SSPs; `optL` for the `optL` SSPs; `uniform` for the
                uniform SSPs.}
}

\value{
\item{ssp}{the subsampling probabilities.}

\item{index.pilot}{the index of pilot subsample in the full sample.}
}

\description{
This function calculate the optimal subsampling probabilities for the full sample.
The inputs are full sample covariate matrix x, response vector y and censoring
indicator delta. Moreover, the pilot subsample size `r_0` needs to be specified. For
more arguments, please read the details section.
It returns a list with two elements: ssp, the subsampling probabilities;
index.pilot, the index of pilot subsample in the full sample.
}
\examples{
library(paraftosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Get the optimal SSPs based on A-optimality
set.seed(101)
ssps <- osmacWei.ssp(x, y, delta, n.pilot = 500, method = "optA")
}
\keyword{osmacWei.ssp}
