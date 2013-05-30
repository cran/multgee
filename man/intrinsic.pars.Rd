\name{intrinsic.pars}
\alias{intrinsic.pars}
\title{
Intrinsic Parameters Estimation 
}
\description{
Utility function to assess the underlying association pattern.
}
\usage{
intrinsic.pars(response, id, repeated, ncategories, rscale = "ordinal")
}
\arguments{
  \item{response}{
a vector that indicates the response variables.
}
  \item{id}{
a vector that identifies the clusters.
}
  \item{repeated}{
a vector that identifies the order of the observations within each cluster.
}
  \item{ncategories}{
a positive constant that indicates the number of observed response categories.
}
  \item{rscale}{
a character string that indicates the nature of the response scale. Options include "\code{ordinal}" or "\code{nominal}". 
}
}
\details{
Simulation studies in Touloumis et al. (2013) suggested that if the range of the intrinsic parameter estimates is small then simple local odds ratios structures should adequately approximate the association pattern. Otherwise more complicated structures should be employed.

The intrinsic parameters are estimated under the heterogeneous linear-by-linear association model (Agresti, 2013) for ordinal response categories and under the RC-G(1) model (Becker and Clogg, 1989) with homogeneous score parameters for nominal response categories. 
}
\value{
Returns a numerical vector with the estimated intrinsic parameters.
}
\author{
Anestis Touloumis
}
\references{
Agresti, A. (2013). \emph{Categorical Data Analysis}. New York: John Wiley and Sons, Inc., 3rd Edition.

Becker, M. and Clogg, C. (1989). Analysis of sets of two-way contingency tables using association models. \emph{Journal of the American Statistical Association}, \bold{84}, 142-151.

Touloumis, A., Agresti, A. and Kateri, M. (2013). GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}.
}
\seealso{
\link{nomLORgee} and \link{ordLORgee}.
}
\examples{
data(arthritis)
intrinsic.pars(arthritis$y,arthritis$id,arthritis$time,5)
## The intrinsic parameters do not differ much. The 'uniform' local odds ratios
## structure might be a good approximation for the association pattern.

set.seed(1)
data(housing)
intrinsic.pars(housing$y,housing$id,housing$time,3,rscale="nominal")
## The intrinsic parameters vary. The 'RC' local odds ratios structure
## might be a good approximation for the association pattern.
}
