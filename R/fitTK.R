#' Bayesian inference of TK model with Stan
#'
#' @export
#' @param stanTKdata List of Data require for computing
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
fitTK <- function(stanTKdata, ...) {
  stanfit <- rstan::sampling(stanmodels$TK, data = stanTKdata, ...)
  out <- list(stanTKdata = stanTKdata, stanfit = stanfit)
  # class(out) <- append("fitTK", class(out))
  return(out)
}
