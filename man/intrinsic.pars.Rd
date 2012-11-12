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
a vector that identifies the subjects.
}
  \item{repeated}{
a vector that identifies the order of the observations within each subject.
}
  \item{ncategories}{
a positive constant that indicates the maximum number of the observed response categories.
}
  \item{rscale}{
a character string that indicates the nature of the response scale. Options include "\code{ordinal}" or "\code{nominal}". 
}
}
\details{
Simulation studies in Touloumis, Agresti and Kateri (2012) suggested that if the range of the intrinsic parameter estimates is small then simple local odds ratios structures should adequately approximate the association pattern.

The intrinsic parameters are estimated under the heterogeneous uniform association model for ordinal response categories and under the RC-G model with homogeneous score parameters for nominal response categories. 
}
\value{
Returns a vector with the estimated intrinsic parameters.
}
\author{
Anestis Touloumis
}
\references{
Touloumis, A., Agresti, A. and Kateri, M. (2012). GEE for multinomial responses using a local odds ratios parameterization. \emph{Submitted}.
}
\seealso{
\link{nomLORgee} and \link{ordLORgee}
}
\examples{
data(arthritis)
intrinsic.pars(arthritis$y,arthritis$id,arthritis$time,5)
## The intrinsic parameters do not differ much. The 'uniform' local odds ratios
## structure might be a good approximation for the association pattern.

set.seed(1)
data(housing)
intrinsic.pars(housing$y,housing$id,housing$time,3,rscale="nominal")
## The intrinsic parameters seem to vary. The 'RC' local odds ratios structure
## might be a good approximation for the association pattern
}
