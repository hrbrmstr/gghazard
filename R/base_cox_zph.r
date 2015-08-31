#' Graphical Test of Proportional Hazards
#'
#' Displays a graph of the scaled Schoenfeld residuals, along with a smooth curve.
#'
#' This is an augmented version of the \code{survival::plot.cox.zph} function
#' that allows for specifying separate colors for the smooth fit line and
#' confidence interval lines
#'
#' @param x	result of the \code{cox.zph} function.
#' @param resid	a logical value, if \code{TRUE} the residuals are included on
#'        the plot, as well as the smooth fit.
#' @param se a logical value, if \code{TRUE}, confidence bands at two standard
#'        errors will be added.
#' @param df the degrees of freedom for the fitted natural spline, \code{df=2}
#'        leads to a linear fit.
#' @param nsmo number of points used to plot the fitted spline.
#' @param var	 the set of variables for which plots are desired. By default,
#'        plots are produced in turn for each variable of a model.
#'        Selection of a single variable allows other features to be added to
#'        the plot, e.g., a horizontal line at zero or a main title.\cr
#'        \cr
#'        This has been superseded by a subscripting method; see the example below.
#' @param smoothcol color for the smooth fit line
#' @param confcol color for the upper/lower confidence interval lines
#' @param ...	additional graphical arguments passed to the \code{plot} function.
#' @note Load this pacakge AFTER the \code{survival} pacakge.
#' @export
#' @examples
#' library(survival)
#' vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) +
#'               karno + age, data=veteran, x=TRUE)
#'
#' plot(cox.zph(vfit))
plot.cox.zph <- function (x, resid = TRUE, se = TRUE, df = 4,
                          nsmo = 40, var,
                          smoothcol="black", confcol=smoothcol, ...)  {
    xx <- x$x
    yy <- x$y
    d <- nrow(yy)
    df <- max(df)
    nvar <- ncol(yy)
    pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
    temp <- c(pred.x, xx)
    lmat <- splines::ns(temp, df = df, intercept = TRUE)
    pmat <- lmat[1:nsmo, ]
    xmat <- lmat[-(1:nsmo), ]
    qmat <- qr(xmat)
    if (qmat$rank < df)
        stop("Spline fit is singular, try a smaller degrees of freedom")
    if (se) {
        bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
        xtx <- bk %*% t(bk)
        seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
    }
    ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
    if (missing(var))
        var <- 1:nvar
    else {
        if (is.character(var))
            var <- match(var, dimnames(yy)[[2]])
        if (any(is.na(var)) || max(var) > nvar || min(var) <
            1)
            stop("Invalid variable requested")
    }
    if (x$transform == "log") {
        xx <- exp(xx)
        pred.x <- exp(pred.x)
    }
    else if (x$transform != "identity") {
        xtime <- as.numeric(dimnames(yy)[[1]])
        indx <- !duplicated(xx)
        apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx),
            length = 17)[2 * (1:8)])
        temp <- signif(apr1$y, 2)
        apr2 <- approx(xtime[indx], xx[indx], temp)
        xaxisval <- apr2$y
        xaxislab <- rep("", 8)
        for (i in 1:8) xaxislab[i] <- format(temp[i])
    }

    for (i in var) {
        y <- yy[, i]
        yhat <- pmat %*% qr.coef(qmat, y)
        if (resid)
            yr <- range(yhat, y)
        else yr <- range(yhat)
        if (se) {
            temp <- 2 * sqrt(x$var[i, i] * seval)
            yup <- yhat + temp
            ylow <- yhat - temp
            yr <- range(yr, yup, ylow)
        }
        if (x$transform == "identity")
            plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i],
                ...)
        else if (x$transform == "log")
            plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i],
                log = "x", ...)
        else {
            plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i],
                axes = FALSE, ...)
            axis(1, xaxisval, xaxislab)
            axis(2)
            box()
        }
        if (resid)
            points(xx, y)
        lines(pred.x, yhat, col=smoothcol)
        if (se) {
            lines(pred.x, yup, lty = 2, col=confcol)
            lines(pred.x, ylow, lty = 2, col=confcol)
        }
    }
}
