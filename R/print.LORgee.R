print.LORgee <-
function (x, ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nNumber of clusters:", x$max.id, "\n")
    cat("Maximum cluster size:", max(x$clusz), "\n")
    cat("\nNumber of categories:", max(x$categories), "\n")
    cat("Number of observations:", x$nobs, "\n")
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nNumber of iterations:", x$convergence$niter, "\n")
    cat("Algorithm converged:", x$convergence$conv, "\n")
    if (6 <= nrow(x$odds.ratio$theta)) {
        cat("\nLocal Odds Ratio Structure[1:6,1:6]\n")
        print(x$odds.ratio$theta[1:6, 1:6])
    }
    else {
        cat("\nLocal Odds Ratio Structure[1:6,1:6]\n")
        print(x$odds.ratio$theta)
    }
    if (!is.null(x$pvalue)) 
        cat("\npvalue of Null model:", round(x$pvalue, 4), "\n")
}

