#' Learn a quadratic regression model via hierarchical group-lasso regularization from \code{glinternet} package. 
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param family A character string describing the target variable: "gaussian" for continuous (the default), "binomial" for logistic.
#' @param ... other arguments to be passed to the \code{glinternet} calls.
#'
#' @return An object of S3 class "\code{cv.glinternet}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{family}}{The response \code{family} used in \code{cv.glinternet}.}
#'   \item{\code{a0}}{Intercept value.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @import glinternet
#' @export
cv.glinternet <- function(x, y, family = "gaussian", ...){
  
  n <- nrow(x)
  p <- ncol(x)
  
  # conduct cross validation for glinternet and returns a value of lambda.
  # numLevels: Number of levels for each variable, of length p. Set to 1 for continuous variables.
  fit <- glinternet::glinternet.cv(x, y, numLevels = rep(1, p), family = family, ...)
  
  # grab estimated coefficients
  a0 <- fit$betahat[[1]][1]
  main.ix <- coef(fit)$mainEffects$cont
  main.n <- length(main.ix)
  inter.ix <- coef(fit)$interactions$contcont
  inter.n <- nrow(inter.ix)
  ix <- matrix(rbind(cbind(rep(0,main.n),main.ix),inter.ix),ncol=2)
  
  # main effects
  main.coef <- NULL
  if(is.null(main.ix)){
    main.coef <- NULL
  }else{
    for(i in seq(main.n)){main.coef <- append(main.coef, coef(fit)$mainEffectsCoef$cont[[i]])}
  }
  
  # interactions
  inter.coef <- NULL
  if(is.null(inter.ix)){
    inter.coef <- NULL
  }else{
    for(i in seq(inter.n)){inter.coef <- append(inter.coef, coef(fit)$interactionsCoef$contcont[[i]])}
  }
  
  # we now combine estimated main effects and interactions coefficients
  compact <- cbind(ix, c(main.coef,inter.coef))
  colnames(compact) <- c("index_1", "index_2", "coefficient")
  
  out <- list(n = n,
              p = p,
              family = family,
              a0 = a0,
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}