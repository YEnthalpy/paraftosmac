\name{osmacWei.est}
\alias{osmacWei.est}
\title{Calculate the Subsample Estimator}
\usage{
osmacWei.est(x, y, delta, pi, 
             control = list(maxit = 1000, tol = 1e-5, 
                            init = c(1, rep(0, dim(x)[2])),
                            step = 0.8))
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{pi}{the subsampling probabilities.}

\item{control}{the list of controling variables, including the initial value
               ,the tolerance, the maximum number of iterations and the step size.}
}

\value{
\item{scale}{the estimates of the scale parameter}

\item{par}{the estimates of the regression coefficients}

\item{message}{the message of converging}

}

\description{
This function calculate the subsampling estiamtor for the subsample
covariate matrix x, response vector y, censoring indicator delta and the 
subsampling probabilities. The model is the Weibull AFT model.
It returns a list with three elements: 
scale, the estimates of the scale parameter; 
par, the estimates of regression coefficients; message, the converging message.
}
\examples{
library(paraftosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Sampling with replacement based on uniform subsampling method
# Subsample size 500
set.seed(101)
ind <- sample(n, 500, replace = TRUE)

# Get the subsample estimator
osmacWei.est(x[ind, ], y[ind], delta[ind], rep(1/n, 500))
}
\keyword{osmacWei.est}
