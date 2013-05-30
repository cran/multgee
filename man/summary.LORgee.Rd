\name{summary.LORgee}
\alias{summary.LORgee}
\title{
Summary Method Function for Objects of Class `LORgee'
}
\description{
Returns a summary method for GEE marginal models for correlated multinomial responses using \code{ordLORgee} and \code{nomLORgee}.
}
\usage{
\method{summary}{LORgee}(object, \dots)
}
\arguments{
  \item{object}{
a fitted LORgee object.
}
  \item{\dots}{
arguments to be passed to or from other methods.
}
}
\value{
Returns an object of the class "\code{LORgee}". This has components:
\item{title}{title for the GEE model.}
\item{version}{the current version of the GEE solver.}
\item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values.}
\item{link}{the marginal link function.}
\item{local.odds.ratios}{the `working' local odds ratios structure name along with the estimates.}
\item{call}{the matched call.}
\item{residuals}{the residuals of the marginal regression model. The \code{j}-th column corresponds to the \code{j}-th response category.}
\item{niter}{The number of iterations the algorithm needed for convergences.}
\item{pvalue}{the p-value based on a Wald test that no covariates are statistically significant.}
}
\author{
Anestis Touloumis
}
\keyword{internal}
