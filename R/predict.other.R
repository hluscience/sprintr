#' Calculate prediction from a \code{other} object.
#'
#' @param object a fitted \code{other} object. (Notice: not used for ordinal logistic fit object)
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the fit \code{object}.
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.other <- function(object, newdata, ...){

  # input check
  stopifnot(ncol(newdata) == object$p)

  idx <- object$compact[, 1:2, drop = FALSE]
  # selected indices for main effects
  idxm <- idx[idx[, 1] == 0, 2]
  # selected index pairs for interactions
  idxi <- idx[idx[, 1] != 0, , drop = FALSE]

  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)

  if(nrow(idxi) == 1)
    xint <- matrix(xm[, idxi[, 1]] * xm[, idxi[, 2]], ncol = 1)
  else
    xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]

  design <- cbind(newdata[, idxm], xint)

  mu <- as.numeric(object$a0 + design %*% object$compact[, 3])

  return(mu)
}
