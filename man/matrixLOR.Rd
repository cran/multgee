\name{matrixLOR}
\alias{matrixLOR}
\title{
Creating A Probability Matrix With Specified Local Odds Ratios
}
\description{
Utility function to create a square probability matrix with a specified local odds ratios structure.
}
\usage{
matrixLOR(x)
}
\arguments{
  \item{x}{
a square matrix with positive entries that describes the desired local odds ratios matrix.
}
}
\value{
Returns a square probability matrix that satisfies the local odds ratios structure described by \code{x}.
}
\section{Warning}{
Caution is needed for local odds ratios close to zero.
}
\details{
This function is designed to ease the construction of a fixed local odds ratios structure in the functions \link{nomLORgee} or \link{ordLORgee}. See also the example.
}
\author{
Anestis Touloumis
}
\seealso{
\link{nomLORgee} and \link{ordLORgee}.
}
\examples{
\dontrun{Illustrating the construction of a "fixed" local odds ratio structure 
         using the arthritis dataset. Here, we assume a uniform local odds ratios 
         structure equal to 2 for each time pair.}
\dontrun{Create the uniform local odds ratios structure.}
lorterm <- matrixLOR(matrix(2,4,4))
\dontrun{Create the LORterm argument.}
lorterm <- c(lorterm)
lorterm <- matrix(c(lorterm),3,25,TRUE)
data(arthritis)
fitmod.fixed <- ordLORgee(y~factor(trt)+factor(time)+factor(baseline), id="id",
                 repeated="time", data=arthritis, LORstr="fixed", LORterm=lorterm)
fitmod.fixed
}

