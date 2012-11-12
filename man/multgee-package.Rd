\name{multgee-package}
\alias{multgee}
\docType{package}
\title{
A GEE Solver For Correlated Nominal Or Ordinal Multinomial Responses
}
\description{
A generalized estimating equations (GEE) solver for fitting marginal regression models with correlated nominal or ordinal multinomial responses using a local odds ratios parameterization for the association structure.
}
\details{
The package contains two core functions; \link{ordLORgee} for an ordinal response scale, and \link{nomLORgee} for a nominal response scale.

The following arguments are required in both functions: (i) a data frame (\code{data}), (ii) a model formula (\code{formula}), (iii) a subject identifier variable (\code{id}) and (iv) a vector that identifies the order of observations within each subject (\code{repeated}). Note that \code{data} must contain \code{id}, \code{repeated} and any variable present in \code{formula}.

Options for the marginal model in \link{ordLORgee} include cumulative link models or an adjacent category logit model, determined by the argument \code{link}. A marginal baseline category logit model is offered in \link{nomLORgee}. See the help files of \link{nomLORgee} and \link{ordLORgee} for the form of the linear predictor.

The association structure among the correlated multinomial responses is expressed via local odds ratios (Touloumis, Agresti and Kateri, 2012). The estimating procedure for the local odds ratios can be summarized as follows: For each level pair of the \code{repeated} variable, the responses are aggregated across subjects to form a square marginalized contingency table. 
Treating these tables as independent, the argument \code{LORem} indicates whether a model for counts is fitted either seperately to each marginalized contingency table ("\code{2way}") or simultaneously to the full marginalized contingency table ("\code{3way}"). The model applied to the counts is the RC-G model proposed by Becker and Clogg (1989) and it belongs to the family of association models developed by Goodman (1981). 
The argument \code{LORstr} determines the form of the local odds ratios structure. 

If the underlying association pattern does not change dramatically across the level pairs of \code{repeated} then a simple form for the local odds ratios should be a sufficient approximation for the association structure. To assess this, one might compare the intrinsic parameters of an RC-G model. See the utility function \link{intrinsic.pars} for more details.

Alternatively a user-defined local odds ratios structure can be provided if \code{LORstr} is "\code{fixed}". In this case, the utility function \link{matrixLOR} is useful in constructing the required \code{LORterm}.
}

\author{
Anestis Touloumis
Maintainer: Anestis Touloumis <anestis@ebi.ac.uk>
}

\references{
Becker, M. and Clogg, C. (1989). Analysis of sets of two-way contingency tables using association models. \emph{Journal of the American Statistical Association} 84, 142-151.

Goodman, L. (1985). The analysis of cross-classified data having ordered and or/unordered categories: Association models, correlation models, and asymmetry models for contingency tables with or without missing entries. \emph{The Annals of Statistics} 13, 10-69.

Touloumis, A., Agresti, A. and Kateri, M. (2012). GEE for multinomial responses using a local odds ratios parameterization. \emph{Submitted}.
}

\seealso{
\link{nomLORgee} and \link{ordLORgee}.
}

\keyword{ package }

\examples{
data(arthritis)
data <- arthritis
fitord <- ordLORgee(y~factor(time)+factor(trt)+factor(baseline),
          id="id",repeated="time",data=data)
summary(fitord) 

data(housing)
data <- housing
fitnom <- nomLORgee(y~factor(time)*sec, id="id",repeated="time",data=data)
summary(fitnom) 
}
