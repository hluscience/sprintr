#' Build a regression model with hierarchically constrained pairwise interactions from \code{hierNet} package.
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param family response type. Default is ’gaussian’. The other choice is ’binomial’ for logistic regression
#' @param nlam Number of values of lambda to be tried.
#' @param nfolds Number of folds in cross-validation. Default value is 10.
#' @param num_keep Number of variables to keep in the screening phase
#' @param strong Indicator to specify strong hierarchy (TRUE) or weak hierarchy (FALSE). Default FALSE.
#' @param ... other arguments to be passed to the \code{hierNet} calls.
#'
#' @return An object of S3 class "\code{cv.hierNet}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{family}}{The response \code{family} used in \code{cv.hierNet}.}
#'   \item{\code{strong}}{The \code{strong} parameter passed into \code{cv.hierNet}.}
#'   \item{\code{a0}}{Intercept value.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @import hierNet, Matrix
#' @export
cv.hierNet <- function(x, y, family = "gaussian", nlam = 20, nfolds = 10, strong = TRUE, ...){
  
  n <- nrow(x)
  p <- ncol(x)
  
  # fit a logistic path of hierNet models over different values of the regularization parameter. 
  if(family == "gaussian"){
    fit <- hierNet::hierNet.path(x, y, nlam = nlam, strong = strong, ...)
  }else if(family == "binomial"){
    fit <- hierNet::hierNet.logistic.path(x, y, nlam = nlam, strong = strong, ...)
  }
  
  # use cross-validation to estimate the regularization parameter for hierNet
  fitcv <- hierNet::hierNet.cv(fit, x, y, nfolds = nfolds, ...)
  
  # the index of best lambda
  lambda_best <- which.min(fitcv$cv.err)
  
  # grab coefficient estimate
  # overall main effect estimated coefficients are bp-bn
  # bp: "positive part" main effect
  # pn: "negative part" main effect;
  main <- fit$bp[, lambda_best] - fit$bn[, lambda_best]
  main.ix <- seq(p)[main!=0]
  main.coef <- main[main!=0]
  
  # th: Matrix of estimated interaction coefficients, of dimension p-by-p. 
  # Note: when output from hierNet is printed, th is symmetrized (set to (th+t(th))/2) for simplicity.
  inter.mt <- (t(fit$th[, , lambda_best])+fit$th[, , lambda_best])/2
  inter.vc <- inter.mt[lower.tri(inter.mt, diag=T)] 
  inter.ix <- Matrix::which(Matrix::tril(inter.mt!=0), arr.ind = TRUE)[,c("col","row")]
  inter.coef <- inter.vc[inter.vc!=0]
  
  # we now combine estimated main effects and interactions coefficients
  coef.cb <- c(main.coef, inter.coef)
  ix.cb <- rbind(cbind(rep(0,sum(main!=0)), main.ix), inter.ix)
  compact <- cbind(ix.cb, coef.cb)
  colnames(compact) <- c("index_1", "index_2", "coefficient")
  rownames(compact) <- NULL
  a0 <- as.numeric(fit$b0[lambda_best])
  
  out <- list(n = n,
              p = p,
              family = family, 
              strong = strong,
              a0 = a0,
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}