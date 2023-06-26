
# functions for monotone splines

#' Generate monotone (I-spline) spline basis
#'
#' @param y outcome vector
#' @param quantiles quantiles of y to use as spline basis knots
#' @param degree spline degree
#' @param knots optional; manually specified spline basis knots
#'
#' @return I-spline basis matrix
get_ispline_basis <- function(y,
                              quantiles = c(0.25, 0.5, 0.75),
                              degree = 3,
                              knots = NULL) {
  if (is.null(knots)) {
    splines2::iSpline(
      y,
      knots = quantile(y, probs = quantiles),
      degree = degree,
      intercept = TRUE,
      Boundary.knots = range(y, na.rm=T))
  } else {
    splines2::iSpline(
      y,
      knots = knots[-c(1, length(knots))], # no boundary knots here
      degree = degree,
      intercept = TRUE,
      Boundary.knots = knots[c(1, length(knots))])
  }
}


#' Generate M-spline spline basis
#'
#' @param y outcome vector
#' @param quantiles quantiles of y to use as spline basis knots
#' @param degree spline degree
#' @param knots optional; manually specified spline basis knots
#'
#' @return M-spline basis matrix
get_mspline_basis <- function(y,
                              quantiles = c(0.25, 0.5, 0.75),
                              degree = 3,
                              knots = NULL) {
  if (is.null(knots)) {
    splines2::mSpline(
      y,
      knots = quantile(y, probs = quantiles),
      degree = degree,
      intercept = TRUE,
      Boundary.knots = range(y, na.rm=T))
  } else {
    splines2::mSpline(
      y,
      knots = knots[-c(1, length(knots))],
      degree = degree,
      intercept = TRUE,
      Boundary.knots = knots[c(1, length(knots))])
  }
}


#' Get I-spline coefficients from y ~ basis spline fit
#'
#' @param y outcome vector
#' @param basis I-spline basis
#' @param nonneg restrict coefficients to be non-negative?
#'
#' @return monotone spline coefficient vector (intercept, coefs)
get_ispline_coefs <- function(y,
                              basis,
                              nonneg = T) {
  ll <- ifelse(nonneg, 0, -Inf)
  mod <- glmnet::glmnet(
    x = basis,
    y = y,
    lambda = 0,         # no lasso penalty
    lower.limits = ll,  # run non-negative least squares
    intercept = T)   # keep intercept term
  return(c(as.vector(mod$a0), as.vector(mod$beta)))
}

# get all knots (boundary and interior quantile) for a given df
get_knots <- function(df, outcome=y, num_knots = 3) {
  outcome <- dplyr::pull(df, {{outcome}})
  if (is.null(num_knots)) {
    num_knots <- 3
  }

  # evenly space num_knots quantiles
  quants <- seq(0, 1, length.out = num_knots + 2)
  quants <- quants[-c(1, length(quants))]

  # store spline basis knots
  int_knots <- quantile(outcome, quants, na.rm=T)

  # interior knots may end up on the boundary due to repeated outcomes,
  #  so we need to do some manual adjustment
  if (min(int_knots) == min(outcome) | max(int_knots) == max(outcome)) {
    warning("Interior knots placed on boundary: manually adjusting...")
    min_knot <- min(int_knots)
    max_knot <- max(int_knots)
    inc <- sd(outcome)

    int_knots[1] <- min_knot+0.1*inc
    int_knots[length(int_knots)] <- max_knot-0.1*inc
  }

  # return all knots (interior and boundary)
  return(c(min(outcome), int_knots, max(outcome)))
}

