#' Build a quadratic regression model with regularization algorithm under marginality principle from \code{RAMP} package.
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param family response type. Default is ’gaussian’. The other choice is ’binomial’ for logistic regression
#' @param hier whether to specify strong or weak heredity. Default is ’Strong’.
#' @param ... other arguments to be passed to the \code{RAMP} calls.
#'
#' @return An object of S3 class "\code{ramp.fit}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{family}}{The response \code{family} used in \code{ramp.fit}.}
#'   \item{\code{hier}}{The \code{hier} parameter passed into \code{ramp.fit}.}
#'   \item{\code{a0}}{Intercept value.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @import RAMP
#' @export
ramp.fit <- function(x, y, family = "gaussian", hier = "Strong", ...){
  
  n <- nrow(x)
  p <- ncol(x)
  
  # call RAMP algorithm for logistic regression model fitting
  fit <- RAMP::RAMP(X=x, y=y, family = family, hier = hier, ...)
  
  # grab coefficient estimate
  a0 <- fit$a0
  # main effects
  main <- fit$mainInd
  if(is.null(main)){
    compact.m <- NULL
  }else{
    main.ix <-  matrix(cbind(rep(0,length(main)), main), ncol=2)
    main.coef <- fit$beta.m
    compact.m <- matrix(cbind(main.ix, main.coef), ncol=3)
  }
  
  # interactions
  inter <- fit$interInd
  if(is.null(inter)){
    compact.i <- NULL
  }else{
    inter.ix <- matrix(unlist(regmatches(inter, gregexpr("[[:digit:]]+", inter))),ncol=2,byrow = T)
    inter.coef <- fit$beta.i
    compact.i <- matrix(cbind(inter.ix, inter.coef),ncol=3)
  }
  
  # we now combine estimated main effects and interactions coefficients
  compact <- matrix(as.numeric(rbind(compact.m, compact.i)), ncol=3)
  colnames(compact) <- c("index_1", "index_2", "coefficient")
  rownames(compact) <- NULL
  
  out <- list(n = n,
              p = p,
              family = family, 
              hier = hier,
              a0 = a0,
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}