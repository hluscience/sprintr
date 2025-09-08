#' Oracle method
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param family Either a character string representing one of the built-in exponential families, including "gaussian", "poisson", and "binomial". Default is "gaussian".
#' @param num_keep Number of variables to keep in the screening phase
#' @param ... other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}
#'
#' @return An object of S3 class "\code{sis_lasso}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{a0}}{Intercept value.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 30
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- oracle_lasso(x = x, y = y)
#'
#' @import glmnet
#' @export

oracle_lasso <- function(x, y, family = "gaussian", num_keep = NULL, type.measure = "deviance", ...){

  # x is the unstandardized main effects
  p <- ncol(x)
  n <- nrow(x)
  # Total number of square effects and pairwise interactions
  q <- (p^2 + p) / 2

  # Set default value for num_keep if not specified
  if(is.null(num_keep))
    num_keep <- ceiling(n / log(n))

  # Ensure num_keep does not exceed the number of possible interactions
  num_keep <- min(num_keep, q)

  # Step 1: screen interactions
  # Standardize the main effects
  x <- myscale(x)
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")

  # idx contains indices for square effects and pairwise interactions
  idx <- rbind(cbind(seq(p), seq(p)), t(combn(p, 2)))

  # Initialize an empty matrix to store interaction strengths
  interaction_strengths <- matrix(NA, nrow = q, ncol = 3)

  # Iterate through each interaction index
  for (i in 1:q) {

    # Construct the design matrix
    xx <- myscale(x[, idx[i, 1]] * x[, idx[i, 2]])
    design <- cbind(x, xx)
    col_mean_cb <- c(col_mean, attr(x = xx, which = "scaled:center"))
    col_sd_cb <- c(col_sd, attr(x = xx, which = "scaled:scale"))

    # Fit cv.glmnet
    fit <- glmnet::cv.glmnet(x = design, y = y,
                             family = family,
                             type.measure = type.measure,
                             intercept = TRUE,
                             standardize = FALSE, ...)

    # Extract coefficients at the best lambda
    ibest <- which.min(fit$cvm)
    beta <- as.numeric(fit$glmnet.fit$beta[, ibest])
    beta <- beta / col_sd_cb

    # Store interaction strengths
    interaction_strengths[i, ] <- c(idx[i, 1], idx[i, 2], abs(beta[p + 1]))

  }

  # Select the top 'num_keep' interactions
  ordered_interactions <- interaction_strengths[order(interaction_strengths[, 3], decreasing = TRUE), ]
  top_interactions <- ordered_interactions[1:num_keep, ]
  top_interactions <- top_interactions[top_interactions[,3] != 0, ]

  # Step 2: Fit cv.glmnet with all main effects and selected interactions
  # Construct design matrix
  xx <- myscale(x[, top_interactions[, 1]] * x[, top_interactions[, 2]])
  col_mean_cb <- c(col_mean, attr(x = xx, which = "scaled:center"))
  col_sd_cb <- c(col_sd, attr(x = xx, which = "scaled:scale"))
  design <- cbind(x, xx)

  # Fit cv.glmnet
  fit <- glmnet::cv.glmnet(x = design, y = y, family = family,
                           intercept = FALSE,
                           standardize = FALSE,
                           type.measure = type.measure, ...)

  # Extract coefficients at the best lambda
  coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
  coef <- coef / col_sd_cb
  a0 <- as.numeric(fit$glmnet.fit$a0[which.min(fit$cvm)] - crossprod(col_mean_cb, coef))

  # Create a compact representation of the selected variables
  idx_all <- rbind(cbind(rep(0, p), seq(p)), top_interactions[, 1:2])
  compact <- cbind(idx_all[which(coef != 0), , drop = FALSE], coef[coef != 0])
  rownames(compact) <- NULL
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  out <- list(n = n,
              p = p,
              a0 = a0,
              compact = compact,
              call = match.call())

  class(out) <- "other"
  return(out)
}
