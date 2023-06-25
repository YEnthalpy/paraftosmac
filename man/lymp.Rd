\docType{data}
\name{lymp}
\alias{lymp}
\title{The Lymphoma Data: Data Used in the Example Code}
\format{
  The dataset \code{lymp} contains 159,149 lymphoma patients diagnosed 
  from 1973 to 2012 and the censoring rate is 58.3\%. All risk factors
  are centered and scaled so that the optimal SSPs will not dominate by
  some risk factors with large means and standard errors. The interaction
  is defined after centering and scaling.
}
\description{
  Estimate the survival time of Lymphoma patients
}
\details{
  The variables are as follows:
    
    \itemize{
      \item{Age: continous variable, the age in the unit of year}
      \item{Nonwhite: the nonwhite race indicator (1 =nonwhite)}
      \item{Male: male indicator (1 = male)}
      \item{Diagnostic Year: continuous variable, the year when patients are diagnosed}
      \item{Age*Nonwhite: interaction between Age and Nonwhite}
      \item{Age*Male: interaction between Age and Male}
    }
}
\keyword{datasets}