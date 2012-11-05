nomLORgee <-
function (formula = formula, data = data, id = id, repeated = repeated, 
    bstart = NULL, LORstr = "time.exch", LORem = "3way",
    LORterm = NULL, add = 0, homogeneous = TRUE, 
    control = LORgee.control(), ipfp.ctrl = ipfp.control(), IM = "solve") 
{
    restricted <- NULL
    LORstrs <- c("independence", "time.exch", "RC", "fixed")
    icheck <- as.integer(match(LORstr, LORstrs, -1))
    if (icheck < 1) {
        stop("unknown odds ratio structure")
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        LORem <- NULL
    } else {
        if (LORstr == "RC") 
            LORem <- "2way"
        if (LORem != "2way" & LORem != "3way") 
            stop("'LORem' must be '2way' or '3way'")
    }
    if (LORstr == "time.exch" | LORstr == "RC") {
        if (!is.logical(homogeneous)) 
            stop("'homogeneous' must be 'TRUE' or 'FALSE'")
    }
    else {
        homogeneous <- NULL
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        add <- NULL
    }
    else {
        if (!is.numeric(add) | add < 0) 
            stop("'add' must be >=0")
    }
    ipfp.ctrl <- ipfp.ctrl
    control <- control
    verbose <- control$verbose
    IMs <- c("cholesky", "solve", "qr.solve" )
    icheck <- as.integer(match(IM, IMs, -1))
    if (icheck < 1) 
        stop("unknown method for inverting a matrix")
    if (missing(data)) 
        stop("Dataframe not found")
    if (!is.data.frame(data)) 
        data <- data.frame(data)
    rownames(data) <- 1:nrow(data)
    m <- model.frame(formula, data)
    nonmissrows <- as.numeric(rownames(m))
    data <- data[nonmissrows, ]
    if (missing(id)) 
        stop("ID variable not determined")
    id <- as.character(id)
    if (match(id, names(data), 0) == 0) 
        stop("ID variable not found in the dataframe")
    if (missing(repeated)) 
        stop("'repeated' not found")
    repeated <- as.character(repeated)
    if (match(repeated, names(data), 0) == 0) 
        stop("'repeated' not found in 'data'")
    data <- data[order(data[, id], data[, repeated], na.last = NA), 
        ]
    rownames(data) <- seq.int(nrow(data))
    id <- as.numeric(factor(data[, id]))
    repeated <- as.numeric(factor(data[, repeated]))
    dummy <- split(repeated, id)
    if (all(unlist(lapply(dummy, length)) != unlist(lapply(lapply(dummy, 
        unique), length)))) 
        stop("'repeated' does not have unique values per 'id'")
    Y <- as.numeric(factor(model.response(model.frame(formula, 
        data = data))))
    ncategories <- nlevels(factor(Y))
    if (ncategories <= 2) 
        stop("The response should have more than 2 categories")
    if (LORstr != "independence" & LORstr != "fixed") {
        data.model <- datacounts(Y, id, repeated, ncategories)
        marpars <- mmpar(LORem, LORstr, max(data.model$tp), homogeneous)
        LORem <- marpars$LORem
        LORstr <- marpars$LORstr
        LORterm <- fitmm(data.model, marpars, homogeneous, 
            NULL, add)
    }
    id <- rep(id, each = ncategories - 1)
    repeated <- rep(repeated, each = ncategories - 1)
    link <- "bcl"
    if (is.null(bstart)) {
        coeffs <- vglm(formula, multinomial(refLevel = ncategories), 
            data = data)
        coeffs <- c(matrix(coef(coeffs), ncol = ncategories - 
            1, byrow = TRUE))
        if (!is.numeric(coeffs)) 
            stop("Please insert initial values")
        if (verbose) {
            cat("\nGEE FOR NOMINAL MULTINOMIAL RESPONSES\n")
            cat("\nrunning 'vglm' function to get initial regression estimates\n")
            print(matrix(coeffs, ncol = 1, dimnames = list(1:length(coeffs), 
                "Initial.Values")))
        }
    }
    data <- data.frame(lapply(data, function(x) rep(x, each = ncategories - 
        1)))
    Intercept <- rep.int(seq(ncategories - 1), length(id)/(ncategories - 
        1))
    data <- cbind(data, Intercept)
    m <- model.frame(update(formula, ~factor(Intercept) + . - 
        1), data = data)
    Terms <- attr(m, "terms")
    Y <- as.numeric(factor(model.response(m)))
    Y <- as.numeric(Y == Intercept)
    Xinit_mat <- model.matrix(Terms, m)
    xxnames <- colnames(Xinit_mat)
    X_mat <- Xinit_mat[, 1:(ncategories - 1)]
    if (ncol(Xinit_mat) != (ncategories - 1)) {
        X_inter <- Xinit_mat[, (1:(ncategories - 1))]
        for (i in ncategories:ncol(Xinit_mat)) X_mat <- cbind(X_mat, 
            X_inter * Xinit_mat[, i])
    }
    X_mat <- matrix(X_mat, ncol = ncol(X_mat), dimnames = NULL)
    X_mat <- X_mat[, c(matrix(seq(ncol(X_mat)), ncol = ncategories - 
        1, byrow = TRUE))]
    offset <- model.offset(m)
    if (length(offset) <= 1) 
        offset <- rep(0, length(Y))
    if (length(offset) != length(Y)) 
        stop("'offset' and 'y' not same length")
    offset <- as.double(offset)
    if (!is.null(bstart)) {
        coeffs <- as.numeric(bstart)
        if (length(coeffs) != ncol(X_mat)) 
            stop("'bstart' and 'beta' differ in length")
        if (verbose) {
            cat("\nGEE FOR NOMINAL MULTINOMIAL RESPONSES\n")
            cat("\nuser's initial regression estimate\n")
            print(matrix(coeffs, ncol = 1, dimnames = list(1:length(coeffs), 
                "Initial.Values")))
        }
    }
    fitmod <- fitLORgee(Y, X_mat, coeffs, ncategories, id, repeated, offset, 
    link, LORterm, marpars, ipfp.ctrl, control, IM, LORem = LORem, 
    LORstr = LORstr, add)
    fit <- list()
    fit$title <- "GEE FOR NOMINAL MULTINOMIAL RESPONSES"
    fit$version <- "version 1.0 modified 04-11-2012"
    fit$link <- c("Baseline-category Logit")
    fit$odds.ratio <- list()
    fit$odds.ratio$structure <- LORstr
    fit$odds.ratio$model <- LORem
    fit$odds.ratio$homogeneous <- homogeneous
    fit$odds.ratio$theta <- fitmod$theta
    fit$terms <- attr(m, "terms")
    fit$contrasts <- attr(model.matrix(Terms, m), "contrasts")
    fit$nobs <- nrow(data)
    fit$convergence <- list()
    fit$convergence$niter <- fitmod$iter
    fit$convergence$criterion <- fitmod$crit[fitmod$iter]
    fit$convergence$conv <- fitmod$conv
    xxnames[1:(ncategories - 1)] <- paste("beta0", 1:(ncategories - 
        1), sep = "")
    xnames <- xxnames[1:(ncategories - 1)]
    if (length(xxnames) > length(xnames)) {
        for (i in 1:((ncol(X_mat) - ncategories + 1)/(ncategories - 
            1))) xnames <- c(xnames, paste(xxnames[i + ncategories - 
            1], 1:(ncategories - 1), sep = ":"))
    }
    xnames <- xnames[c(matrix(seq(ncol(X_mat)), ncol = ncategories - 
        1, byrow = TRUE))]
    fit$coefficients <- fitmod$beta_mat[, fitmod$iter + 1]
    names(fit$coefficients) <- xnames
    fit$linear.predictors <- matrix(fitmod$linear.predictor, 
        ncol = ncategories - 1, byrow = TRUE)
    rownames(fit$linear.predictors) <- 1:nrow(fit$linear.predictors)
    colnames(fit$linear.predictors) <- 1:(ncategories - 1)
    fitted.values <- fitmod$fitted.values
    fitted.values.1 <- matrix(fitted.values, ncol = ncategories - 
        1, byrow = TRUE)
    fitted.values.2 <- 1 - rowSums(fitted.values.1)
    fitted.values <- cbind(fitted.values.1, fitted.values.2)
    rownames(fitted.values) <- 1:nrow(fitted.values.1)
    colnames(fitted.values) <- 1:ncategories
    fit$fitted.values <- fitted.values
    fit$residuals <- matrix(fitmod$residuals, ncol = ncategories - 
        1, byrow = TRUE)
    rownames(fit$residuals) <- 1:nrow(fit$residuals)
    colnames(fit$residuals) <- 1:(ncategories - 1)
    y <- Y
    y <- apply(matrix(y, ncol = ncategories - 1, byrow = TRUE), 
        1, function(x) which(x == 1))
    y <- as.numeric(y)
    y[is.na(y)] <- ncategories
    fit$y <- y
    fit$id <- id
    fit$max.id <- max(unique(id))
    fit$clusz <- unlist(lapply(split(id, id), length))/(ncategories - 
        1)
    fit$robust.variance <- fitmod$robust
    dimnames(fit$robust.variance) <- list(xnames, xnames)
    fit$naive.variance <- fitmod$naive
    dimnames(fit$naive.variance) <- list(xnames, xnames)
    fit$xnames <- xnames
    fit$categories <- ncategories
    fit$occasions <- sort(unique(repeated))
    fit$gee.control <- control
    fit$ipfp.control <- ipfp.ctrl
    fit$inverse.method <- IM
    fit$adding.constant <- add
    if (control$TRACE) {
        fit$trace <- list()
        fit$trace$coeffs <- fitmod$beta_mat
        fit$trace$crit <- fitmod$crit
    }
    if (length(xxnames) == (ncategories - 1)) 
        fit$pvalue <- NULL
    else {
        dummy <- seq(1, length(xxnames), ncategories - 1)
        waldts <- fit$coefficients[-dummy] %*% solve((fit$robust.variance)[-dummy, 
            -dummy])
        waldts <- waldts %*% fit$coefficients[-dummy]
        fit$pvalue <- 1 - pchisq(waldts, length(xxnames) - length(dummy))
    }
    fit$call <- match.call()
    class(fit) <- "LORgee"
    fit
}