

tmm_norm_median_ref <- function(x, ...) {
  
  median_column <- 
    apply(x, 
          MARGIN = 1, 
          median)
  
  x <- cbind(x, median_column)
  
  norm_data <- tmm(x, refColumn = dim(x)[2], ...)
  
  norm_data[, -(dim(x)[2])]
}


