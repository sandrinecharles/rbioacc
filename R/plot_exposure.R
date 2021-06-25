#' Plot exposure profile
#' 
#' 
#' @export
#' 
plot_exposure <- function(object){
  
  col_exposure <- .index_col_exposure(object)
  col_exp_time <- c(col_exposure, match("time", colnames(object)))
  
  col_names <- colnames(object[,col_exposure])
  ls_exp <- lapply(col_exposure, function(i) object[,i])
  
  df <- data.frame(
    time = rep(object$time, length(col_names)),
    route = rep(col_names, each = nrow(object)),
    exposure = do.call("c",ls_exp)
  )

  plt <- ggplot(data = df) + 
    theme_classic() +
    geom_line(aes(x = time, y = exposure)) +
    facet_wrap(~route, scale = "free")
  
  return(plt)
}
