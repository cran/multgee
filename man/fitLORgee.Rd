\name{fitLORgee}
\alias{fitLORgee}
\title{
Internal Fit Function
}
\description{
Internal Fit Function.
}
\usage{
fitLORgee(Y, X_mat, coeffs, ncategories, id, repeated, offset, link,
          LORterm, marpars, ipfp.ctrl, control, IM, LORem = LORem,
           LORstr = LORstr, add)
}
\arguments{
  \item{Y}{
response vector.
}
  \item{X_mat}{
design matrix.
}
  \item{coeffs}{
initial coefficients.
}
  \item{ncategories}{
number of categories.
}
  \item{id}{
id variable.
}
  \item{repeated}{
repeated variable.
}
  \item{offset}{
offset variable.
}
  \item{link}{
link function.
}
  \item{LORterm}{
interaction term.
}
  \item{marpars}{
marginal parameter for counts.
}
  \item{ipfp.ctrl}{
control function for the ipfp.
}
  \item{control}{
control function for the LOR functions.
}
  \item{IM}{
inversion method.
}
  \item{LORem}{
local odds ratio estimation model.
}
  \item{LORstr}{
local odds ratio estimation structure.
}
  \item{add}{
adding constant.
}
}
\value{
A list that produces the appropriate fitting values.
}
\references{
Later on.
}
\author{
Anestis Touloumis
}
\keyword{ internal }