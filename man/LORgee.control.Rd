\name{LORgee.control}
\alias{LORgee.control}
\title{
Control For The GEE Solver
}
\description{
Control function for the GEE solver in the \link{nomLORgee} and \link{ordLORgee} functions.
}
\usage{
LORgee.control(tolerance = 0.001, maxiter = 15, verbose = FALSE, TRACE = FALSE)
}
\arguments{
  \item{tolerance}{
positive convergence tolerance. The algorithm converges when the pairwise maximum of the absolute relative difference in parameter estimates is less than or equal to \code{tolerance}.
}
  \item{maxiter}{
positive integer that indicates the maximum number of iterations in the Fisher-scoring iterative procedure. 
}
  \item{verbose}{
logical that indicates if output should be printed for each iteration.
}
  \item{TRACE}{
logical that indicates if the parameter estimates and the convergence criterion for each iteration should be saved.
}
}
\author{
Anestis Touloumis
}
\seealso{
\link{nomLORgee} and \link{ordLORgee}.
}
\examples{
data(arthritis)
fitmod <- ordLORgee(y~factor(trt)+factor(baseline)+factor(time), id = "id", 
                    repeated = "time", data = arthritis)
\dontrun{A one-step GEE estimator}
fitmod1 <- update(fitmod, control = LORgee.control(maxiter=1))
coef(fitmod)
coef(fitmod1)
}
