\name{osmacWei.fit}
\alias{osmacWei.fit}
\title{The Two-step Procedure}
\usage{
osmacWei.fit(x, y, delta, n.pilot, n.sample, method, se = TRUE,
             pilot = TRUE)
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{n.pilot}{the input pilot subsample size.}

\item{n.sample}{the input second-step subsample size.}

\item{method}{the type of optimal subsampling probabilities (SSPs). `optA` for the
                `optA` SSPs; `optL` for the `optL` SSPs; `uniform` for the
                uniform SSPs.}
              
\item{se}{the indicator of calculating standard errors. The default is `TRUE`.}
           
\item{pilot}{the indicator of whether combining the pilot and the second-step subsample.}
}

\value{
\item{Est}{the first column of the output which is the resultant estimator.}

\item{SE}{the second column of the output which is the estimated standard errors.}
}

\description{
This function apply the two-step method for the full sample
covariate matrix x, response vector y and censoring indicator delta with user
specified pilot subsample size and second-step subsample size.
For more arguments, please read the details section.
It returns a matrix with two columns: Est, the first column which is the
subample estimator; SE, the second column which is the estimated standard errors.
}
\examples{
library(semirkosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Get the subsample estimator using A-optimal SSPs when \code{r0 = 500} and \code{r = 2000}.
set.seed(101)
osmacWei.fit(x, y, delta, n.pilot = 500, n.sample = 2000, 
            method = "optA", se = "TRUE", pilot = "TRUE")
}
\keyword{osmacWei.fit}
