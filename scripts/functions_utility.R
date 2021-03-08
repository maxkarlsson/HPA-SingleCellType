
savepath <- 
  function(savename) { 
    wd <- getwd()
    result_folder <- paste0(str_extract(wd, ".*Blood cell deconvolution"), "/results/", Sys.Date())
    
    dir.create(result_folder, showWarnings = FALSE)
    
    paste0(result_folder, "/", savename)
  }

omega_sq <- function(aov_in, neg2zero=T){
  aovtab <- summary(aov_in)[[1]]
  n_terms <- length(aovtab[["Sum Sq"]]) - 1
  output <- rep(-1, n_terms)
  SSr <- aovtab[["Sum Sq"]][n_terms + 1]
  MSr <- aovtab[["Mean Sq"]][n_terms + 1]
  SSt <- sum(aovtab[["Sum Sq"]])
  for(i in 1:n_terms){
    SSm <- aovtab[["Sum Sq"]][i]
    DFm <- aovtab[["Df"]][i]
    output[i] <- (SSm-DFm*MSr)/(SSt+MSr)
    if(neg2zero & output[i] < 0){output[i] <- 0}
  }
  output <- c(output, 1 - sum(output))
  names(output) <- c(rownames(aovtab)[1:n_terms], "Residuals")
  
  return(output)
}

multispread <- function(df, key, value) {
  
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}
