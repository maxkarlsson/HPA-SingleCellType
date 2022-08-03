
savepath <- 
  function(savename) { 
    wd <- getwd()
    result_folder <- paste0(str_extract(wd, ".*HPA-SingleCellType"), "/results/", Sys.Date())
    
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


rotate_coords <- 
  function(x, y, rotate_angle, rotate_center = c(0, 0)) {
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             # Angle
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)) + rotate_angle,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata
  }

shrink_rotation_coords <- 
  function(x, y, shrink_angle, rotate_center = c(0, 0)) {
    
    
    shrink_factor <- 
      1 - shrink_angle / (2 * pi)
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
             
             # Shrink angle
             # angle = case_when(x == 0 & y == 0 ~ 0, 
             #                   quadrant %in% c(1,4) ~ angle + shrink_angle,
             #                   quadrant %in% c(2,3) ~ angle - shrink_angle),
             
             angle = angle * shrink_factor,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata 
  }


calculate_retina_cut_angle <- 
  function(clust) {
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)
    
    g_edgepoints <- 
      g$data %>% 
      as_tibble() %>% 
      filter(height == 0) %>% 
      left_join(cutree(clust, k = 2) %>% 
                  enframe("label", 
                          "cluster"),
                by = "label") %>% 
      mutate(angle = calculate_coord_angle(x, y))
    
    expand_grid(node1 = g_edgepoints$.ggraph.index,
                node2 = g_edgepoints$.ggraph.index) %>% 
      left_join(g_edgepoints %>% 
                  select(node1 = .ggraph.index,
                         angle1 = angle, 
                         cluster1 = cluster),
                by = "node1") %>% 
      left_join(g_edgepoints %>% 
                  select(node2 = .ggraph.index,
                         angle2 = angle, 
                         cluster2 = cluster),
                by = "node2") %>% 
      filter(cluster1 == 1,
             cluster2 == 2) %>% 
      group_by_all() %>% 
      mutate(dist = c(angle1 - angle2,
                      (angle1 - 2 * pi) - angle2,
                      angle1 - (angle2 - 2 * pi),
                      (angle1 - 2 * pi) - (angle2 - 2 * pi)) %>% 
               abs() %>% 
               min()) %>% 
      ungroup() %>% 
      arrange(dist) %>% 
      slice(1:2) %>% 
      mutate(cut_angle = (angle1 + angle2) / 2) %>% 
      pull(cut_angle) %>% 
      {2 * pi - .}
  }


calculate_coord_angle <- 
  function(x, y, rotate_center = c(0, 0)) {
    tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp))) %>% 
      pull(angle)
  }

