#' Bayesian inference of TK model with Stan
#'
#' @export
#' @param stanTKdata List of Data require for computing
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
fitTK <- function(stanTKdata, ...) {
  out <- rstan::sampling(stanmodels$TK, data = stanTKdata, ...)
  return(out)
}