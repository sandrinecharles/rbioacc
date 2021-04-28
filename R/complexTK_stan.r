#' Bayesian inference of simple model with Stan
#'
#' @export
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#'
#'
complexTK_stan <- function(stanTKdata, ...) {
  out <- rstan::sampling(stanmodels$complexTK, data = stanTKdata, ...)
  return(out)
}





