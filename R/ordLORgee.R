ordLORgee <-
function (formula = formula, data = data, id = id, repeated = repeated, 
    link = "logistic", bstart = NULL, LORstr = "category.exch", 
    LORem = "3way", LORterm = NULL, add = 0, homogeneous = TRUE, 
    restricted = FALSE, control = LORgee.control(), 
    ipfp.ctrl = ipfp.control(), IM = "solve") 
{
    options(contrasts=c("contr.treatment", "contr.poly"))
    LORstrs <- c("independence", "uniform", "category.exch", 
        "time.exch", "RC","fixed")
    icheck <- as.integer(match(LORstr, LORstrs, -1))
    if (icheck < 1) {
        stop("unknown odds ratio structure")
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        LORem <- NULL
    }
    else if (LORstr == "category.exch") {
        LORem <- "3way"
    }
    else {
        if (LORstr == "RC") 
            LORem <- "2way"
        if (LORem != "2way" & LORem != "3way") 
            stop("'LORem' must be '2way' or '3way'")
    }
    if (LORstr == "time.exch" | LORstr == "RC") {
        if (!is.logical(homogeneous)) 
            stop("'homogeneous' must be 'TRUE' or 'FALSE'")
        if (!is.logical(restricted)) 
            stop("'restricted' must be 'TRUE' or 'FALSE'")
        restricted <- if (!restricted) NULL else TRUE
    }
    else {
        homogeneous <- restricted <- NULL
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
    IMs <- c("cholesky", "solve","qr.solve")
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
        LORterm <- fitmm(data.model, marpars,homogeneous, 
            restricted, add)
    }
    id <- rep(id, each = ncategories - 1)
    repeated <- rep(repeated, each = ncategories - 1)
    link <- as.character(link)
    link <- switch(link, logistic = "logit", probit = "probit", 
        cloglog = "cloglog", cauchit = "cauchit", acl = "acl")
    if (is.null(link)) 
        stop("'link' must be \"logistic\", \"probit\", \"cloglog\", \"cauchit\" or \"acl\"")
    if (is.null(bstart)) {
        if (link == "acl") family <- acat(reverse = TRUE, parallel = TRUE)
        if (link == "logit") family <- cumulative("logit", parallel = TRUE)
        if (link == "probit") family <- cumulative("probit", parallel = TRUE)
        if (link == "cloglog") family <- cumulative("cloglog", parallel = TRUE)
        if (link == "cauchit") family <- cumulative("cauchit", parallel = TRUE)
        capture.output(coeffs <- vglm(formula, data = data, family = family), 
            file = "NUL")
        coeffs <- as.numeric(as.vector(coef(coeffs)))
        if (!is.numeric(coeffs)) 
            stop("Please insert initial values")
        if (verbose) {
            cat("\nGEE FOR ORDINAL MULTINOMIAL RESPONSES\n")
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
    X_mat <- model.matrix(Terms, m)
    xnames <- colnames(X_mat)
    X_mat <- matrix(X_mat, ncol = ncol(X_mat), dimnames = NULL)
    if (link == "acl") {
        dummy <- ncategories - 1
        dummy.matrix <- diagmod(rep.int(1, dummy))
        dummy.matrix[upper.tri(dummy.matrix)] <- 1
        X_mat[, 1:dummy] <- kronecker(rep.int(1, nrow(X_mat)/dummy), 
            dummy.matrix)
        if (dummy != ncol(X_mat)) {
            X_mat[, -c(1:dummy)] <- X_mat[, -c(1:dummy)] * rep(dummy:1, 
                nrow(X_mat)/dummy)
        }
    }
    offset <- model.offset(m)
    if (length(offset) <= 1) 
        offset <- rep(0, length(Y))
    if (length(offset) != length(Y)) 
        stop("'offset' and 'y' not same length")
    offset <- as.double(offset)
    if (!is.null(bstart)) {
        coeffs <- as.numeric(bstart)
        if (length(coeffs) != ncol(X_mat)) 
            stop("Starting values and parameters vector differ in length")
        if (any(diff(coeffs[1:(ncategories - 1)]) < 0)) 
            stop("cutpoints are not increasing")
        if (verbose) {
            cat("\nGEE FOR ORDINAL MULTINOMIAL RESPONSES\n")
            cat("\nuser's initial regression estimate\n")
            print(matrix(coeffs, ncol = 1, dimnames = list(1:length(coeffs), 
                "Initial.Values")))
        }
    }
    fitmod <- fitLORgee(Y, X_mat, coeffs, ncategories, id, repeated, offset, 
    link, LORterm, marpars, ipfp.ctrl, control, IM, LORem = LORem, 
    LORstr = LORstr, add)
    fit <- list()
    fit$title <- "GEE FOR ORDINAL MULTINOMIAL RESPONSES"
    fit$version <- "version 1.1 modified 10-11-2012"
    fit$link <- if (link == "acl") 
        paste("Adjacent Category Logit")
    else paste("Cumulative", link, sep = " ")
    fit$odds.ratio <- list()
    fit$odds.ratio$structure <- LORstr
    fit$odds.ratio$model <- LORem
    fit$odds.ratio$homogeneous <- homogeneous
    fit$odds.ratio$restricted <- restricted
    fit$odds.ratio$theta <- fitmod$theta
    fit$terms <- attr(m, "terms")
    fit$contrasts <- attr(model.matrix(Terms, m), "contrasts")
    fit$nobs <- nrow(data)
    fit$convergence <- list()
    fit$convergence$niter <- fitmod$iter
    fit$convergence$criterion <- fitmod$crit[fitmod$iter]
    fit$convergence$conv <- fitmod$conv
    xnames[1:(ncategories - 1)] <- paste("beta0", 1:(ncategories - 
        1), sep = "")
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
    fit$call <- match.call()
    if (control$TRACE) {
        fit$trace <- list()
        fit$trace$coeffs <- fitmod$beta_mat
        fit$trace$crit <- fitmod$crit
    }
    if (length(xnames) == (ncategories - 1)) 
        fit$pvalue <- NULL
    else {
        dummy <- 1:(ncategories - 1)
        waldts <- fit$coefficients[-dummy] %*% solve((fit$robust.variance)[-dummy, 
            -dummy])
        waldts <- waldts %*% fit$coefficients[-dummy]
        fit$pvalue <- 1 - pchisq(waldts, length(xnames) - length(dummy))
    }
    class(fit) <- "LORgee"
    fit
}

