\name{housing}
\alias{housing}
\docType{data}
\title{
Homeless Data
}
\description{
Housing status for 362 severely mentally ill homeless subjects measured at baseline and at three follow-up times.
}
\usage{data(housing)}
\format{
  A data frame with 1448 observations on the following 4 variables:
  \describe{
    \item{\code{id}}{Subject identifier variable.}
    \item{\code{y}}{Housing status response, coded as (0) for street living, (1) for community living and (2) for independent housing.}
    \item{\code{time}}{Measurement occasions variable, coded as (0) at baseline, (1) at 6 months, (2) at twelve months and (3) at twenty four months.}
    \item{\code{sec}}{Section 8 rent certificate indicator.}
  }
}
\source{
Hedeker, D. (2003). A mixed-effects multinomial logistic regression model. 
\emph{Statistics in Medicine}, 22, 1433--1446.
}
\examples{
data(housing)
str(housing)
}
\keyword{datasets}
