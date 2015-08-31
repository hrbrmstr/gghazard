#' Create data structures from survival hazard Cox Regression
#'
#' This function generates data to be used with \code{ggplot2} and is custom-
#' tailored for use with \code{gg_cox_zpn}. It returns a list the length of
#' \code{x} with names matching the row names of the output of \code{cox.zph}.\cr
#' \cr
#' Each \code{list} element has two \code{data.frames}, \code{resid} and
#' \code{pred} matching the residuals and fit+confidence interval for that element
#' along with \code{ylab} (which matches the list element name) and a boolean
#' \code{log} which indicates whether a log scale should be used.
#'
#' @param x result of the \code{cox.zph} function.
#' @param df the degrees of freedom for the fitted natural spline, \code{df=2}
#'         leads to a linear fit.
#' @param nsmo number of points used to plot the fitted spline.
#' @return \code{list} with residuals and fit+confidence intervals for each
#'         component of \code{x}
#' @note Designed to be used with gg_cox_zph but can be used separately.
#' @export
#' @examples
#' library(survival)
#' library(ggplot2)
#' vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) +
#'               karno + age, data=veteran, x=TRUE)
#'
#' gg_cox_zph(fortify(cox.zph(vfit)))
fortify.cox.zph <- function(x, df=4, nsmo=40) {

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

  bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
  xtx <- bk %*% t(bk)
  seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)

  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])

  var <- 1:nvar

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

  lapply(var, function(i) {

    y <- yy[, i]
    yhat <- pmat %*% qr.coef(qmat, y)

    temp <- 2 * sqrt(x$var[i, i] * seval)
    yup <- yhat + temp
    ylow <- yhat - temp

    list(resid=data.frame(x=xx, y=y),
         pred=data.frame(x=pred.x, y=yhat, ymin=ylow, ymax=yup),
         ylab=ylab[i],
         log=x$transform == "log")

  }) -> out

  names(out) <- ylab

  out

}

#' Graphical Test of Proportional Hazards
#'
#' Displays a graph of the scaled Schoenfeld residuals, along with a smooth curve.
#'
#' @param x	result of the \code{fortify.cox.zph} function.
#' @param resid_col color for the residuals points
#' @param resid_shp point shape for the residuals points
#' @param pred_col color for the prediction fit line
#' @param pred_minmax_col color for the confidence interval lines
#' @param pre_minmax_linetype line type for the confidence interval lines
#' @param ribbon_col color for the shaded confidence interval region
#' @param ribbon_alpha alpha value for the confidence interval region
#' @return \code{grid} object
#' @export
#' @examples
#' library(survival)
#' library(ggplot2)
#' vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) +
#'               karno + age, data=veteran, x=TRUE)
#'
#' gg_cox_zph(fortify(cox.zph(vfit)))
gg_cox_zph <- function(x,
                       resid_col="black", resid_shp=1,
                       pred_col=resid_col, pred_minmax_col=pred_Col,
                       pred_minmax_linetype="dashed",
                       ribbon_col=pred_Col, ribbon_alpha=0.25) {

  lapply(x, function(czplot) {

    gg <- ggplot()
    gg <- gg + geom_ribbon(data=czplot$pred,
                           aes(x, ymax=ymax, ymin=ymin),
                           alpha=0.25)
    gg <- gg + geom_line(data=czplot$pred, aes(x, ymin), linetype="dashed")
    gg <- gg + geom_line(data=czplot$pred, aes(x, ymax), linetype="dashed")
    gg <- gg + geom_line(data=czplot$pred, aes(x, y))
    gg <- gg + geom_point(data=czplot$resid, aes(x, y), shape=1)
    if (czplot$log) gg <- gg + scale_x_log10()
    gg <- gg + labs(x="Time", y=czplot$ylab)
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.grid=element_blank())
    gg

  }) -> czplots

  do.call(grid.arrange, czplots)

}
