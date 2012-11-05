\name{nomLORgee}
\alias{nomLORgee}
\title{
Marginal Models For Correlated Nominal Multinomial Responses
}
\description{
Solving the generalized estimating equations for correlated nominal multinomial responses assuming a baseline category logit model for the marginal probabilities.
}
\usage{
nomLORgee(formula = formula, data = data, id = id, repeated = repeated, 
          bstart = NULL, LORstr = "time.exch", LORem = "3way", LORterm = NULL, 
          add = 0, homogeneous = TRUE, control = LORgee.control(), 
          ipfp.ctrl = ipfp.control(), IM = "solve")
}
\arguments{
  \item{formula}{
a formula expression as for other regression models for multinomial responses. An intercept term must be included.
}
 \item{data}{
a mandatory data frame that should include the variables provided in the \code{formula}, \code{id} and \code{repeated} arguments. 
}
  \item{id}{
a vector that identifies the subjects.
}
  \item{repeated}{
a vector that identifies the order of the observations within each subject.
}
  \item{bstart}{
a vector that includes an initial estimate for the marginal regression parameter vector.
}
 \item{LORstr}{
a character string that indicates the local odds ratios structure. Options include "\code{independence}", "\code{time.exch}", "\code{RC}" or "\code{fixed}".
}
  \item{LORem}{
a character string that indicates if the marginalized local odds ratios structure is estimated simultaneously ("\code{3way}") or seperately at each level pair of \code{repeated} ("\code{2way}").
}
  \item{LORterm}{
a matrix that contains the desired local odds ratios structure. It should be used when \code{LORstr} is "\code{fixed}".
}
  \item{add}{
a positive constant to be added at each cell of the full marginalized contingency table in the presence of zero observed counts.
}
  \item{homogeneous}{
a logical that indicates homogeneous score parameters when \code{LORstr} is "\code{time.exch}" or "\code{RC}".
}
  \item{control}{
a vector that specifies the control variables for the GEE solver.
}
  \item{ipfp.ctrl}{
a vector that specifies the control variables for the \code{ipfp} function.
}
  \item{IM}{
a character string that indicates the method used for inverting a matrix. Options include "\code{solve}", "\code{qr.solve}" or "\code{cholesky}".
}
}
\details{
The \code{data} must be provided in a subject level or equivalently in `long' format. See details about the `long' format in the \link{reshape} function.

A term of the form \code{offset(expression)} is allowed in the \code{formula}.

The \code{id} and the \code{repeated} do not need to be pre-sorted. Instead the function reshapes \code{data} in an ascending order of \code{id} and \code{repeated}.

The default set for the response categories is \eqn{1,\ldots,I}, where \eqn{I>2} is the maximum observed response category. If otherwise, the function recodes the observed response categories onto this set. 

The default set for the levels of \code{repeated} is \eqn{1,\ldots,T}, where \eqn{T} is the number of observed levels. If otherwise, the function recodes the observed levels onto this set. 

The \eqn{I}-th response category is treated as baseline. 

The linear predictor is of the form 
\deqn{\beta_{0j} +\beta^{'}_j x_{it}}
where \eqn{\beta_{0j}} is the \eqn{j}-th intercept, \eqn{\beta_j} is the \eqn{j}-th response category specific parameter vector and \eqn{x_{it}} is the covariate vector for the \eqn{i}-th subject at the \eqn{t}-th level of \code{repeated}.

The \code{LORterm} argument must be an \eqn{L} x \eqn{I^2} matrix, where \eqn{L} is the number of level pairs of \code{repeated}. These are ordered as \eqn{(1,2), (1,3), ...,(1,T), (2,3),...,(T-1,T)} and the rows of \code{LORterm} are supposed to preserve this order. Each row is assumed to contain the vectorized form of a probability table that satisfies the desired local odds ratios structure.

}
\value{
Returns an object of the class "\code{LORgee}". Generic \link{coef}, \link{summary}, \link{print}, \link{fitted} and \link{residuals} methods are available. The \code{pvalue of the Null model} corresponds to the hypothesis \eqn{H_0: \beta_1=...=\beta_{I-1}=0} based on the Wald test statistic.  
}
\references{
Touloumis, A., Agresti, A. and Kateri, M. (2012). GEE for multinomial responses using a local odds ratios parameterization. \emph{Submitted}.
}
\author{
Anestis Touloumis
}
\seealso{
For an ordinal response scale use \link{ordLORgee}. 
}
\examples{
data(housing)
data <- housing
fitmod <- nomLORgee(y~factor(time)*sec, id="id",repeated="time",data=data)
summary(fitmod) 
coef(fitmod)
}