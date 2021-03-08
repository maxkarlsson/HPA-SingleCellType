# ----- Classification plots -----

all_class_plots <- 
  function(class_table, savename) {
    # n_genes barplot & pie
    spec_dist_barplots <- class_spec_dist_barplot(class_table)
    spec_dist_barplots
    
    ggsave(savepath(paste(savename, "class n_genes per dist barplot.pdf")), spec_dist_barplots[[1]], width = 5, height = 4)
    ggsave(savepath(paste(savename, "class n_genes per spec barplot.pdf")), spec_dist_barplots[[2]], width = 5, height = 4)
    
    spec_dist_piecharts <- class_spec_dist_piechart(class_table)
    spec_dist_piecharts
    
    ggsave(savepath(paste(savename, "class n_genes spec piechart.pdf")), spec_dist_piecharts[[1]], width = 6, height = 6)
    ggsave(savepath(paste(savename, "class n_genes dist piechart.pdf")), spec_dist_piecharts[[2]], width = 6, height = 6)
    
    
    # n_genes chord
    pdf(savepath(paste(savename, "class comb n_genes chord.pdf")), width = 6, height = 6, useDingbats = F)
    class_chord_plot(class_table)
    dev.off()
    class_chord_plot(class_table)
    
    # n_genes tissue
    n_genes_tissue_spec_plot <- 
      class_tissue_n_enriched_barplot(class_table) 
    ggsave(savepath(paste(savename, "class n_genes enriched per tissue barplot.pdf")), n_genes_tissue_spec_plot, width = 5, height = 6)
    
    n_genes_tissue_dist_plot <-
      class_tissue_n_expressed_barplot(class_table)
    ggsave(savepath(paste(savename, "class n_genes detected per tissue barplot.pdf")), n_genes_tissue_dist_plot, width = 5, height = 6)
  }



class_spec_dist_piechart <- function(class_table) {
  spec_plot <-
    class_table %>%
    group_by(spec_category) %>%
    summarise(n_genes = n()) %>%
    mutate(spec_category = factor(spec_category, levels = spec_category_levels),
           perc = paste0("(", round(100 * n_genes / sum(n_genes), digits = 1), "%)")) %>%

    ggplot(aes("", n_genes, fill = spec_category)) +
    geom_bar(stat = "identity",
             show.legend = F,
             color = "white",
             size = 1)+
    geom_text(aes(y = n_genes, label = paste(n_genes, perc, sep = "\n")),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+
    geom_text(aes(x = 1.5, y = n_genes, label = spec_category),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+

    coord_polar("y", start = 0)+
    scale_fill_manual(values = gene_category_pal) +
    theme_void()


  ### Distribution
  dist_plot <-
    class_table %>%
    group_by(dist_category) %>%
    summarise(n_genes = n()) %>%
    mutate(dist_category = factor(dist_category, levels = dist_category_levels),
           perc = paste0("(", round(100 * n_genes / sum(n_genes), digits = 1), "%)")) %>%

    ggplot(aes("", n_genes, fill = dist_category)) +
    geom_bar(stat = "identity",
             show.legend = F,
             color = "white",
             size = 1)+
    geom_text(aes(y = n_genes, label = paste(n_genes, perc, sep = "\n")),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+
    geom_text(aes(x = 1.5, y = n_genes, label = dist_category),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+

    coord_polar("y", start = 0)+
    scale_fill_manual(values = gene_category_pal) +
    theme_void()

  list(spec_plot, dist_plot)
}

class_elevated_bar_plot <- function(class_table){

  class_table %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(tissue = names)} %>%
    gather(key = "Classification", value = "Number of genes", -tissue) %>%
    mutate(tissue = factor(tissue, levels = rev(unique(tissue[order(mapply(tissue, FUN = function(x) sum(`Number of genes`[tissue == x & Classification %in% c("Tissue enriched","Celltype enriched",
                                                                                                                                                               "Group enriched","Tissue enhanced","Celltype enhanced")])))]))),
           Classification = factor(Classification, levels = c("Not detected in any tissues","Not detected in any celltypes","Not detected in this tissue","Not detected in this celltype","Mixed in this tissue", "Mixed in this celltype","Expressed in all tissues","Expressed in all celltypes","Tissue enhanced", "Celltype enhanced","Group enriched","Tissue enriched", "Celltype enriched")),
           Classification = gsub(pattern = translate_categories, replacement = names(translate_categories), Classification)) %>%
    filter(Classification %in% c("Tissue enriched","Celltype enriched",
                                 "Group enriched","Tissue enhanced","Celltype enhanced")) %>%
    ggplot(aes(tissue, `Number of genes`, fill = Classification))+
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "",values = cat2.cols)+
    simple_theme+
    xlab("")+
    ylab("Number of genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = c(0.8, 0.8))

}


class_spec_dist_barplot <- function(class_table) {
  cat <-
    class_table %>%
    group_by(spec_category, dist_category) %>%
    summarise(n_genes=n()) %>%
    ungroup() %>%
    mutate(spec_category = factor(spec_category, levels = spec_category_levels),
           dist_category = factor(dist_category, levels = dist_category_levels))


  list(ggplot(cat, aes(x=dist_category, y=n_genes, fill=spec_category)) +
         geom_bar(stat = "identity")+
         ggtitle("Distribution category")+
         scale_fill_manual(values = gene_category_pal, name = "Specificity category")+
         ylab("Number of genes")+
         xlab("")+
         simple_theme+
         theme(axis.text.x = element_text(angle=60, hjust=1)),

       ggplot(cat, aes(x=spec_category, y=n_genes, fill=dist_category)) +
         geom_bar(stat = "identity")+
         ggtitle("Specificity category")+
         scale_fill_manual(values = gene_category_pal, name = "Distribution category")+
         ylab("Number of genes")+
         xlab("")+
         simple_theme+
         theme(axis.text.x = element_text(angle=60, hjust=1)))
}

chord_classification <- function(from, to, sizes, grid.col, groups, plot.order, size_labels = F){
  require(circlize)

  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)


  tb <-
    tibble(from, to, sizes)

  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }

  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }

  circos.par(gap.after = gap.after.par)

  chord <-
    tb %>%
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05,
                 preAllocateTracks = 1,
                 order = plot.order)

  if(size_labels) {
    for(i in 1:nrow(chord)) {
      value <- chord$value[i]
      if(is.null(value)) value <- chord$value1[i]
      x1 <- chord$x1[i] - value / 2
      x2 <- chord$x2[i] - value / 2
      to_ <- chord$cn[i]
      from_ <- chord$rn[i]
      circos.text(x = x1, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = from_, niceFacing = T)
      circos.text(x = x2, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = to_, niceFacing = T)
    }
  }



  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")

    adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)

    circos.segments(x0 = mean(xlim), x1 = mean(xlim),
                    y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment,
                    sector.name)


    circos.text(mean(xlim), mean(ylim) + adjustment, sector.name, niceFacing = TRUE, facing = "bending")
  }, bg.border = NA)

  circos.clear()
}


chord_with_title <- function(from, to, sizes, grid.col, groups, plot.order, titles, 
                             from_labels, to_labels, size_labels = F, label_style = "horisontal", gap_adj = 1){
  require(circlize)
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  
  tb <-
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2 * gap_adj)
    } else {
      gap.after.par <- c(gap.after.par, 15 * gap_adj)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2 * gap_adj)
  } else {
    gap.after.par <- c(gap.after.par, 15 * gap_adj)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>%
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05,
                 preAllocateTracks = 1,
                 order = plot.order)
  
  if(size_labels) {
    for(i in 1:nrow(chord)) {
      value <- chord$value[i]
      if(is.null(value)) value <- chord$value1[i]
      x1 <- chord$x1[i] - value / 2
      x2 <- chord$x2[i] - value / 2
      to_ <- chord$cn[i]
      from_ <- chord$rn[i]
      circos.text(x = x1, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = from_, niceFacing = T)
      circos.text(x = x2, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = to_, niceFacing = T)
    }
  }
  
  
  
  group_coordinates <- 
    sapply(plot.order, FUN = function(si) get.cell.meta.data("xlim",sector.index = si)) %>% 
    t() %>% 
    as_tibble(rownames = "sector_name") %>% 
    mutate(group = groups) %>% 
    group_by(group) %>% 
    mutate(x = cumsum(max.data)) %>%
    
    mutate(group_n = row_number(), 
           group_first = group_n == 1) 
  
  
  if(label_style == "horisontal") {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.name <- get.cell.meta.data("sector.index")
      sector_label <- c(from_labels, to_labels)[match(sector.name, c(from, to))]
      sector.index <- get.cell.meta.data("sector.numeric.index")
      
      adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)
      
      circos.segments(x0 = mean(xlim), x1 = mean(xlim),
                      y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment,
                      sector.name)
      
      
      circos.text(mean(xlim), mean(ylim) + adjustment, sector_label, niceFacing = TRUE, facing = "bending")
      
      sector_i <- match(sector.name, group_coordinates$sector_name)
      
      if(group_coordinates$group_first[sector_i]) {
        
        group_coordinates_lim <- 
          group_coordinates %>% 
          filter(group == group_coordinates$group[sector_i]) %>% 
          filter(group_n == max(group_n) | group_n == min(group_n) )
        
        sector_indices <- which(unique(group_coordinates_lim$group) == groups)
        
        group_n <- unique(group_coordinates_lim$group)
        
        max_gr <- max(group_coordinates_lim$x)
        
        titl_x <-
          max_gr / 2
        
        titl_lab <- titles[match(group_n, unique(group_coordinates$group))]
        
        circos.text(titl_x, 1.2, titl_lab, niceFacing = TRUE, facing = "bending", cex = 2)
      }
    }, bg.border = NA)
  } else if(label_style == "vertical") {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.name <- get.cell.meta.data("sector.index")
      sector_label <- c(from_labels, to_labels)[match(sector.name, c(from, to))]
      sector.index <- get.cell.meta.data("sector.numeric.index")
      
      
      
      circos.text(mean(xlim), min(ylim), sector_label, niceFacing = TRUE, facing = "clockwise", adj = c(0, 0))
      
      sector_i <- match(sector.name, group_coordinates$sector_name)
      
      if(group_coordinates$group_first[sector_i]) {
        
        group_coordinates_lim <- 
          group_coordinates %>% 
          filter(group == group_coordinates$group[sector_i]) %>% 
          filter(group_n == max(group_n) | group_n == min(group_n) )
        
        sector_indices <- which(unique(group_coordinates_lim$group) == groups)
        
        group_n <- unique(group_coordinates_lim$group)
        
        max_gr <- max(group_coordinates_lim$x)
        
        titl_x <-
          max_gr / 2
        
        titl_lab <- titles[match(group_n, unique(group_coordinates$group))]
        
        circos.text(titl_x, 1.2, titl_lab, niceFacing = TRUE, facing = "bending", cex = 2)
      }
    }, bg.border = NA)
  }
  
  
  
  circos.clear()
}



class_chord_plot <- function(class_table) {

  class_table %>%
    group_by(spec_category, dist_category) %>%
    summarise(n_genes=n()) %>%
    ungroup() %>%
    mutate(spec_category = case_when(spec_category == "not detected" ~ "not detected ",
                                     T ~ spec_category)) %$%
    chord_classification(from = spec_category,
                         to = dist_category,
                         sizes = n_genes,
                         grid.col = gene_category_pal,
                         groups = c(rep(1, 5), rep(2, 5)),
                         plot.order = c(c(spec_category_levels[-5], "not detected "),
                                        dist_category_levels),
                         size_labels = T)

}

class_tissue_n_enriched_barplot <- function(class_table) {

  class_table_temp <-
    class_table %>%
    select(gene, spec_category, enriched_tissues) %>%
    separate_rows(enriched_tissues, sep = ";")


  tissues_not_in_plot <- with(tissue_mapping, consensus_tissue_name[which(!consensus_tissue_name %in% class_table_temp$enriched_tissues)])

  if(length(tissues_not_in_plot) != 0) warning(paste0("These tissues are not in the plot: ",
                                                      paste(tissues_not_in_plot, collapse = ", ")))

  class_table_temp %>%
    filter(!is.na(enriched_tissues)) %>%
    group_by(enriched_tissues, spec_category) %>%
    summarise(n_genes = n()) %>%
    ungroup() %>%
    mutate(enriched_tissues = factor(enriched_tissues, levels = {group_by(., enriched_tissues) %>%
        summarise(n_genes = sum(n_genes)) %$%
        enriched_tissues[order(n_genes)]})) %>%
    ggplot(aes(enriched_tissues, n_genes, fill = spec_category)) +
    geom_col(width = 0.9, color = "white", size = 0.1) +
    simple_theme +
    scale_fill_manual(values = gene_category_pal, name = "Specificity") +
    coord_flip() +
    ggtitle("Number of genes enriched per tissue") +
    xlab("Tissue") +
    ylab("Number of genes") +
    scale_y_continuous(position = "bottom")
}

class_tissue_n_enriched_barplot_dendro <- function(class_table, dendro, pal, width = 1) {
  
  class_table_temp <-
    class_table %>%
    select(gene, spec_category, enriched_tissues) %>%
    separate_rows(enriched_tissues, sep = ";")
  
  
  tissues_not_in_plot <- with(tissue_mapping, consensus_tissue_name[which(!consensus_tissue_name %in% class_table_temp$enriched_tissues)])
  
  if(length(tissues_not_in_plot) != 0) warning(paste0("These tissues are not in the plot: ",
                                                      paste(tissues_not_in_plot, collapse = ", ")))
  
  
  dendr <- dendro_data(dendro)
  
  
  dendro_plot_data <- 
    left_join(dendr$segments, 
              dendr$labels, 
              by = c("x" = "x", "yend" = "y")) 
  
  dendro_plot <- 
    dendro_plot_data %>%
    ggplot() +
    geom_segment(aes(x=y, y=x, xend=yend, yend=xend, group = label))+
    geom_rect(aes(xmin=0, ymin=x + 0.5, 
                  xmax=-width, ymax=xend - 0.5, 
                  fill = label), 
              show.legend = F) +
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal)+
    scale_x_reverse(expand = expand_scale(mult = 0.25))+
    theme(axis.text.y = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), units = "mm"), 
          panel.background = element_blank()) 

  
  
  grid.arrange(dendro_plot, 
               class_table_temp %>%
                 filter(!is.na(enriched_tissues)) %>%
                 group_by(enriched_tissues, spec_category) %>%
                 summarise(n_genes = n()) %>%
                 ungroup() %>%
                 mutate(enriched_tissues = factor(enriched_tissues, levels = dendr$labels$label)) %>%
                 ggplot(aes(enriched_tissues, n_genes, fill = spec_category)) +
                 geom_col(width = 0.95, color = "white", size = 0.1) +
                 simple_theme +
                 scale_fill_manual(values = gene_category_pal, name = "Specificity") +
                 coord_flip() +
                 xlab("Tissue") +
                 ylab("Number of genes") +
                 scale_y_continuous(position = "bottom", expand = c(0,0)) + 
                 
                 theme(axis.text.y = element_text(hjust = 0.5), legend.position = c(0.7, 0.5)),
               nrow = 1, 
               widths = c(0.5, 1))
  
  
}

class_tissue_n_expressed_barplot <- function(class_table) {

  class_table_temp <-
    class_table %>%
    select(gene, dist_category, tissues_detected) %>%
    separate_rows(tissues_detected, sep = ";")


  tissues_not_in_plot <- with(tissue_mapping, consensus_tissue_name[which(!consensus_tissue_name %in% class_table_temp$tissues_detected)])

  if(length(tissues_not_in_plot) != 0) warning(paste0("These tissues are not in the plot: ",
                                                      paste(tissues_not_in_plot, collapse = ", ")))

  class_table_temp %>%
    filter(!is.na(tissues_detected)) %>%
    group_by(tissues_detected, dist_category) %>%
    summarise(n_genes = n()) %>%
    ungroup() %>%
    mutate(tissues_detected = factor(tissues_detected, levels = {group_by(., tissues_detected) %>%
        summarise(n_genes = sum(n_genes)) %$%
        tissues_detected[order(n_genes)]})) %>%
    ggplot(aes(tissues_detected, n_genes, fill = dist_category)) +
    geom_col(width = 0.9, color = "white", size = 0.1) +
    simple_theme +
    scale_fill_manual(values = gene_category_pal, name = "Distribution") +
    coord_flip() +
    ggtitle("Number of genes detected per tissue") +
    xlab("Tissue") +
    ylab("Number of genes") +
    scale_y_continuous(position = "bottom")
}

## Retinagram

map_colors_to_edge <- function(edge_data, color_mapping, label_col, color_col, mean_color = F, default_color = "gray80") {
  mapped_colors <-
    edge_data %>%
    as_tibble() %>%
    left_join(color_mapping %>%
                select(label = label_col,
                       color = color_col),
              by = c("node2.label" = "label")) %>%
    mutate(radius = xend^2 + yend^2,
           edge.id = as.character(edge.id)) %>%
    arrange(-radius)

  mapped_colors %>%
    filter(!near(radius, 1)) %$%
    sapply(edge.id,
           FUN = function(edge_id_) {
             xend_ = mapped_colors$xend[which(mapped_colors$edge.id == edge_id_)]
             yend_ = mapped_colors$yend[which(mapped_colors$edge.id == edge_id_)]

             new_color <-
               mapped_colors %>%
               filter(near(x, xend_) & near(y, yend_)) %$%
               ifelse(mean_color,
                      colorRampPalette(color)(3)[2],
                      ifelse(length(unique(color)) == 1,
                             unique(color), default_color))


             new_color_column <- mapped_colors$color
             new_color_column[which(mapped_colors$edge.id == edge_id_)] <- new_color

             mapped_colors$color <<- new_color_column
             NULL


           })


  mapped_colors
}

circular_dendrogram_retinastyle <-
  function(clust, color_mapping, label_col, color_col, mean_color = F) {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)

    dendrogram <-
      clust %>%
      as.dendrogram()



    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)

    edge_data <- get_edges()(g$data)

    edge_data_colors <- map_colors_to_edge(edge_data, color_mapping, label_col, color_col, mean_color = mean_color)


    g +
      scale_edge_width(range = c(1.5, 6))+
      geom_edge_diagonal(aes(edge_color = as.character(edge_data$edge.id),
                             edge_width = 1 - sqrt(xend^2 + yend^2)),
                         strength = 0.8,
                         show.legend = F) +
      scale_edge_color_manual(values = edge_data_colors %$%
                                set_names(c(color, "gray80"), c(edge.id, "")))  +
      g$data %>%
      filter(label != "") %>%
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi,
                                x < 0 ~ 360 - asin(y) * 180 / pi)) %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = "label") %>%
                {geom_node_text(data = .,
                                aes(label = label),
                                angle = .$degree,
                                hjust = ifelse(.$x < 0, 1, 0),
                                vjust = 0.5,
                                size = 3)}  +
      scale_x_continuous(expand = expand_scale(c(0.25, 0.25))) +
      scale_y_continuous(expand = expand_scale(c(0.25, 0.25))) +

      coord_fixed() +
      theme_void()
  }

circular_dendrogram_retinastyle_2 <-
  function(clust, color_mapping, label_col, color_col, 
           scale_expansion = c(0.25, 0.25), text_size = 3, width_range = c(1.5, 6), 
           arc_strength = 0.8, default_color = "gray80") {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)
    
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    
    
    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)
    # 
    # g +
    #   geom_edge_fan(data = edge_data %>%
    #                    mutate(hghl = edge.id == 99),
    #                  aes(label = edge_id, color = as.factor(rank_radius)),
    #                  width =4) +
    #   geom_node_text(aes(label = label))
    
    edge_data <- 
      get_edges()(g$data) %>%
      as_tibble() %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = c("node2.label" = "label")) %>%
      mutate(radius = xend^2 + yend^2,
             edge.id = as.character(edge.id)) %>%
      arrange(-radius) %>%
      mutate(edge_id = as.character(row_number()),
             rank_radius = unclass(factor(-radius)),
             x_m = round(x, 10),
             y_m = round(y, 10),
             xend_m = round(xend, 10),
             yend_m = round(yend, 10)) 
    
    
    edge_id_colors <- 
      edge_data %>% 
      filter(!is.na(color)) %$%
      set_names(color, edge_id)
    
    
    for(rank_rad in 2:max(edge_data$rank_radius)) {
      edge_id_colors_new <- 
        left_join(edge_data %>%
                    select(edge_id, radius, xend_m, yend_m, rank_radius) %>%
                    filter(rank_radius == rank_rad),
                  edge_data %>%
                    select(edge_id, radius, x_m, y_m, rank_radius) %>%
                    filter(rank_radius < rank_rad),
                  by = c("xend_m" = "x_m", "yend_m" = "y_m")) %>%
        left_join(enframe(edge_id_colors),
                  by = c("edge_id.y" = "name")) %>%
        group_by(edge_id.x) %>% 
        summarise(color = ifelse(n_distinct(value) == 1 & any(value != default_color), 
                                 as.character(unique(value)),
                                 default_color)) %$%
        set_names(color, edge_id.x)
      edge_id_colors <- 
        c(edge_id_colors, edge_id_colors_new)
    }
    
    
    
    g +
      scale_edge_width(range = width_range)+
      geom_edge_diagonal(data = edge_data,
                         aes(edge_color = as.character(edge_id),
                             edge_width = 1 - sqrt(xend^2 + yend^2)),
                         strength = arc_strength,
                         show.legend = F) +
      
      
      scale_edge_color_manual(values = edge_id_colors)  +
      g$data %>%
      filter(label != "") %>%
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi,
                                x < 0 ~ 360 - asin(y) * 180 / pi)) %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = "label") %>%
      {geom_node_text(data = .,
                      aes(label = label),
                      angle = .$degree,
                      hjust = ifelse(.$x < 0, 
                                     1, 
                                     0),
                      vjust = 0.5,
                      size = text_size)}  +
      scale_x_continuous(expand = expand_scale(scale_expansion)) +
      scale_y_continuous(expand = expand_scale(scale_expansion)) +
      
      coord_fixed() +
      theme_void()
  }

circular_dendrogram_retinastyle_3 <-
  function(clust, color_mapping, label_col, color_col, 
           scale_expansion = c(0.25, 0.25), text_size = 3, width_range = c(1.5, 6), 
           arc_strength = 0.8, default_color = "gray80", rotate_circle = 0, shrink_circle = 1,
           flip_text = F, text_vjust = 0.5) {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)
    
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    # graph_layout <- 
    #   create_layout(dendrogram, layout = 'dendrogram', circular = T)
    # 
    # 
    # new_layout <-
    #   graph_layout %>% 
    #   as_tibble() %>%
    #   mutate(hyp = sqrt(x^2 + y^2),
    #          angle = ifelse(x == 0 & y == 0, 0, atan(y/x)),
    #          angle_just = angle / 2,
    #          x = cos(angle_just) * hyp,
    #          y = sin(angle_just) * hyp) 
    #   ggplot(aes(x,y)) +
    #   geom_point()
    # 
    #   
    #   g <-
    #     ggraph(graph_from_data_frame(new_layout), layout = new_layout)
      g <-
        ggraph(dendrogram, layout = 'dendrogram', circular = T)
    
      # g +
    #   geom_edge_fan(data = edge_data %>%
    #                    mutate(hghl = edge.id == 99),
    #                  aes(label = edge_id, color = as.factor(rank_radius)),
    #                  width =4) +
    #   geom_node_text(aes(label = label))
    
      
      k1 = tan(rotate_circle)
      k2 = -1/tan(rotate_circle)
      
      if(!(is.finite(k1) & is.finite(k2))) {
        k1 <- 0
        k2 <- 0
      }
      
    ggplot() + 
      geom_point(data = tibble(x = runif(300, -1, 1),
                               y = runif(300, -1, 1)) %>%
                   mutate(quadrant = case_when(x >= y / k2 & y >= x * k1 ~ 1,
                                               x < y / k2 & y >= x * k1 ~ 2,
                                               x < y / k2 & y < x * k1 ~ 3,
                                               x >= y / k2 & y < x * k1 ~ 4)),
                 aes(x, y, color = as_factor(quadrant))) +
      
      geom_abline(slope = c(k1, 
                            k2)) +
      geom_hline(yintercept = 0) + 
      geom_vline(xintercept = 0) +
      # geom_abline(slope = c(-1, 1)) +
      geom_path(data = {
        tt <- seq(0,2*pi,length.out = 100)
        xx <- 1 * cos(tt)
        yy <- 1 * sin(tt)
        data.frame(x = xx, y = yy)
      },
      aes(x, y),
      inherit.aes = F) +
      coord_fixed()
    
    edge_data <- 
      get_edges()(g$data) %>%
      as_tibble() %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = c("node2.label" = "label")) %>%
      mutate(radius = xend^2 + yend^2,
             edge.id = as.character(edge.id)) %>%
      arrange(-radius) %>%
      mutate(edge_id = as.character(row_number()),
             rank_radius = unclass(factor(-radius)),
             x_m = round(x, 10),
             y_m = round(y, 10),
             xend_m = round(xend, 10),
             yend_m = round(yend, 10),
             
             # First pass - rotate circle
             # quadrant = case_when(x >= y / k2 & y >= x * k1 ~ 1,
             #                      x < y / k2 & y >= x * k1 ~ 2,
             #                      x < y / k2 & y < x * k1 ~ 3,
             #                      x >= y / k2 & y < x * k1 ~ 4),
             # quadrant_end = case_when(xend >= yend / k2 & yend >= xend * k1 ~ 1,
             #                          xend < yend / k2 & yend >= xend * k1 ~ 2,
             #                          xend < yend / k2 & yend < xend * k1 ~ 3,
             #                          xend >= yend / k2 & yend < xend * k1 ~ 4),
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             quadrant_end = case_when(xend >= 0 & yend >= 0 ~ 1,
                                      xend < 0 & yend >= 0 ~ 2,
                                      xend < 0 & yend < 0 ~ 3,
                                      xend >= 0 & yend < 0 ~ 4),
             
             hyp = sqrt(x^2 + y^2),
             hyp_end = sqrt(xend^2 + yend^2),
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)) + rotate_circle,
             angle_end = case_when(xend == 0 & yend == 0 ~ 0, 
                                   quadrant_end %in% 1:2 ~ acos(xend/hyp_end),
                                   quadrant_end %in% 3:4 ~ 2 * pi - acos(xend/hyp_end)) + rotate_circle,
             
             x = cos(angle) * hyp,
             y = sin(angle) * hyp, 
             xend = cos(angle_end) * hyp_end,
             yend = sin(angle_end) * hyp_end,
             
             #Second pass - shrink circle
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             quadrant_end = case_when(xend >= 0 & yend >= 0 ~ 1,
                                      xend < 0 & yend >= 0 ~ 2,
                                      xend < 0 & yend < 0 ~ 3,
                                      xend >= 0 & yend < 0 ~ 4),
             
             hyp = sqrt(x^2 + y^2),
             hyp_end = sqrt(xend^2 + yend^2),
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
             angle_end = case_when(xend == 0 & yend == 0 ~ 0, 
                                   quadrant_end %in% 1:2 ~ acos(xend/hyp_end),
                                   quadrant_end %in% 3:4 ~ 2 * pi - acos(xend/hyp_end)),
             
             
             angle_just = angle / shrink_circle,
             angle_end_just = angle_end / shrink_circle,
             
             x = cos(angle_just) * hyp,
             y = sin(angle_just) * hyp, 
             xend = cos(angle_end_just) * hyp_end,
             yend = sin(angle_end_just) * hyp_end)
    
    ####
    # b <- 
    #   a %>%
    #   mutate(edge_id = as.character(row_number()),
    #          rank_radius = unclass(factor(-radius)),
    #          x_m = round(x, 10),
    #          y_m = round(y, 10),
    #          xend_m = round(xend, 10),
    #          yend_m = round(yend, 10),
    #          xend_ori = xend,
    #          yend_ori = yend,
    #          x_ori = x,
    #          y_ori = y) %>% 
    #   select(x, y, xend, yend, x_ori, y_ori, xend_ori, yend_ori)
    # 
    # c <- b %>%
    #   mutate(quadrant = case_when(x >= 0 & y >= 0 ~ 1,
    #                               x < 0 & y >= 0 ~ 2,
    #                               x < 0 & y < 0 ~ 3,
    #                               x >= 0 & y < 0 ~ 4),
    #          quadrant_end = case_when(xend >= 0 & yend >= 0 ~ 1,
    #                                   xend < 0 & yend >= 0 ~ 2,
    #                                   xend < 0 & yend < 0 ~ 3,
    #                                   xend >= 0 & yend < 0 ~ 4),
    #          
    #          hyp = sqrt(x^2 + y^2),
    #          hyp_end = sqrt(xend^2 + yend^2),
    #          
    #          angle = case_when(x == 0 & y == 0 ~ 0, 
    #                            quadrant %in% 1:2 ~ acos(x/hyp),
    #                            quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
    #          angle_end = case_when(xend == 0 & yend == 0 ~ 0, 
    #                                quadrant_end %in% 1:2 ~ acos(xend/hyp_end),
    #                                quadrant_end %in% 3:4 ~ 2 * pi - acos(xend/hyp_end)),
    #          # angle = ifelse(x == 0 & y == 0, 0, atan(y/x)),
    #          # angle_end = ifelse(xend == 0 & yend == 0, 0, atan(yend/xend)),
    #          
    #          angle_just = angle / 2,
    #          angle_end_just = angle_end / 2,
    #          
    #          x = cos(angle_just) * hyp,
    #          y = sin(angle_just) * hyp, 
    #          xend = cos(angle_end_just) * hyp_end,
    #          yend = sin(angle_end_just) * hyp_end)
    # c
    # 
    # c %>% 
    #   ggplot(aes(x_ori, y_ori, xend = xend_ori, yend = yend_ori)) +
    #   geom_segment()
    # 
    # c %>% 
    #   ggplot(aes(x, y, xend = xend, yend = yend)) +
    #   geom_segment()
    # 
    # c %>% 
    #   ggplot(aes(x_ori, y_ori, xend = x, yend = y)) +
    #   geom_segment(arrow = arrow())
    # 
    # c %>% 
    #   ggplot(aes(xend_ori, yend_ori, xend = xend, yend = yend)) +
    #   geom_segment(arrow = arrow())
    ####
    
    edge_id_colors <- 
      edge_data %>% 
      filter(!is.na(color)) %$%
      set_names(color, edge_id)
    
    
    for(rank_rad in 2:max(edge_data$rank_radius)) {
      edge_id_colors_new <- 
        left_join(edge_data %>%
                    select(edge_id, radius, xend_m, yend_m, rank_radius) %>%
                    filter(rank_radius == rank_rad),
                  edge_data %>%
                    select(edge_id, radius, x_m, y_m, rank_radius) %>%
                    filter(rank_radius < rank_rad),
                  by = c("xend_m" = "x_m", "yend_m" = "y_m")) %>%
        left_join(enframe(edge_id_colors),
                  by = c("edge_id.y" = "name")) %>%
        group_by(edge_id.x) %>% 
        summarise(color = ifelse(n_distinct(value) == 1 & any(value != default_color), 
                                 as.character(unique(value)),
                                 default_color)) %$%
        set_names(color, edge_id.x)
      edge_id_colors <- 
        c(edge_id_colors, edge_id_colors_new)
    }
    
    
    node_data <- 
      g$data %>% 
      # as_tibble() %>%
      mutate(hyp = sqrt(x^2 + y^2),
             # First pass - rotate circle
             
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)) + rotate_circle,
             
             x = cos(angle) * hyp,
             y = sin(angle) * hyp, 
             
             #Second pass - shrink circle
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             hyp = sqrt(x^2 + y^2),
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
             
             angle_just = angle / shrink_circle,
             
             x = cos(angle_just) * hyp,
             y = sin(angle_just) * hyp)  
    
    # node_data %>%
    #   as_tibble() %>%
    #   ggplot() +
    #   geom_point(aes(x,y)) +
    #   geom_segment(data = edge_data,
    #                aes(x, y, xend = xend, yend = yend))
    # 
    # node_data %>%
    #   as_tibble() %>%
    #   ggplot(aes(x,y, xend = x_ori, yend = y_ori)) +
    #   geom_segment()
    # 
    # test_data <-
    #   node_data %>% 
    #   as_tibble %>% 
    #   filter(height == 0) %>%
    #   select(x, y, x_ori, y_ori, quadrant, angle, angle_just) %>% 
    #   mutate(degrees = angle * 180 / pi, degrees_just = angle_just * 180 / pi) %>% 
    #   slice(1, 2, 8, 21) %>%
    #   arrange(quadrant)
    # 
    # test_data %>% 
    #   ggplot(aes(x_ori,y_ori, xend = x, yend = y)) + 
    #   geom_segment(color = "red", arrow = arrow()) +
    #   geom_text(aes(label = round(degrees))) +
    #   geom_text(aes(x, y, label = round(degrees_just))) +
    #   geom_hline(yintercept = 0) + 
    #   geom_vline(xintercept = 0) + 
    #   geom_abline(slope = c(-1, 1)) + 
    #   geom_path(data = {
    #     tt <- seq(0,2*pi,length.out = 100)
    #     xx <- 1 * cos(tt)
    #     yy <- 1 * sin(tt)
    #     data.frame(x = xx, y = yy)
    #   },
    #             aes(x, y), 
    #             inherit.aes = F) +
    #   coord_fixed()
    
      
    text_degree <-
      ifelse(flip_text,
             180,
             0)
      
    g +
      scale_edge_width(range = width_range)+
      geom_edge_diagonal(data = edge_data,
                         aes(edge_color = as.character(edge_id),
                             edge_width = 1 - sqrt(xend^2 + yend^2)),
                         strength = arc_strength,
                         show.legend = F) +
      
      
      scale_edge_color_manual(values = edge_id_colors)  +
      node_data %>%
      filter(label != "") %>%
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi + text_degree,
                                x < 0 ~ 360 - asin(y) * 180 / pi + text_degree)) %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = "label") %>%
      {geom_node_text(data = .,
                      aes(label = label),
                      angle = .$degree,
                      hjust = case_when(.$x < 0 & !flip_text ~ 1,
                                        .$x >= 0 & !flip_text ~ 0,
                                        .$x < 0 & flip_text ~ 0,
                                        .$x >= 0 & flip_text ~ 1),
                      vjust = text_vjust,
                      size = text_size)}  +
      scale_x_continuous(expand = expansion(scale_expansion)) +
      scale_y_continuous(expand = expansion(scale_expansion)) +
      
      coord_fixed() +
      theme_void()
  }

pairwise_alluvial_plot <- 
  function(data, var1, var2, cat1, cat2, cat_levels, cat_names, pal) {
    
    alluv_1 <- 
      data %>% 
      select(var1 = var1, 
             var2 = var2, 
             cat1 = cat1, 
             cat2 = cat2) %>% 
      ungroup() %>%
      mutate(row_n = row_number()) %>%
      gather(cat_type, cat, -var1, -var2, -row_n) %>%
      mutate(cat = factor(cat, levels = cat_levels),
             cat_type = case_when(cat_type == "cat1" ~ cat_names[1],
                                  cat_type == "cat2" ~ cat_names[2]),
             cat_type = factor(cat_type, 
                               levels = cat_names)) %>%
      
      ggplot(aes(x = cat_type, stratum = cat, alluvium = row_n,
                 y = 1,
                 fill = cat)) +
      scale_x_discrete(expand = c(.1, .1), position = "top") +
      geom_flow(show.legend = F) +
      geom_stratum(show.legend = F, color = NA) + 
      scale_fill_manual(values = pal) + 
      # theme_void() +
      theme(axis.text.x = element_text(size = 18, face = "bold"),
            axis.text.y = element_blank(), 
            axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.title = element_blank())
    
    alluv_1 <- 
      alluv_1 +
      
      # Flow label
      geom_text(stat = "flow", 
                aes(label = 1),
                nudge_x = ifelse(ggplot_build(alluv_1)$data[[1]]$x == 1, 
                                 1/6,
                                 -1/6), 
                hjust = ifelse(ggplot_build(alluv_1)$data[[1]]$x == 1, 
                               0,
                               1), 
                size = 2) +
      
      # Stratum label
      geom_text(stat = "stratum", 
                aes(label = 1), 
                size = 2,
                nudge_x = ifelse(ggplot_build(alluv_1)$data[[2]]$x == 1, 
                                 -1/6 - 1/100,
                                 1/6 + 1/100),
                angle = ifelse(ggplot_build(alluv_1)$data[[2]]$x == 1, 
                               90,
                               -90), 
                vjust = 0) + 
      
      geom_text(stat = "stratum", 
                aes(label = cat), 
                size = 3)
    
    
    alluv_1
  }

pairwise_subcat_alluvial_plot <- 
  function(data, var1, var2, cat1, cat2, cat_levels, subcat, cat_names, weight = 1, pal, cat_pal) {
    
    
      
      
    
    ####
    
    
    
    if(is.character(weight)) {
      alluv_data <-
        data %>% 
        select(var1 = var1, 
               var2 = var2, 
               cat1 = cat1, 
               cat2 = cat2, 
               subcat = subcat, 
               weight = weight) 
    } else {
      alluv_data <-
        data %>% 
        
        select(var1 = var1, 
               var2 = var2, 
               cat1 = cat1, 
               cat2 = cat2, 
               subcat = subcat) %>%
        mutate(weight = weight)
    }
    
    alluv_1 <- 
      alluv_data %>%
      mutate(cat1 = paste(cat1, subcat),
             cat2 = paste(cat2, subcat)) %>% 
      ungroup() %>%
      mutate(row_n = row_number()) %>%
      gather(cat_type, cat, -var1, -var2, -row_n, -subcat, -weight) %>%
      mutate(cat = factor(cat, levels = c(cat_levels)),
             cat_type = case_when(cat_type == "cat1" ~ cat_names[1],
                                  cat_type == "cat2" ~ cat_names[2]),
             cat_type = factor(cat_type, 
                               levels = cat_names)) %>%
      
      ggplot(aes(x = cat_type, stratum = cat, alluvium = row_n,
                 y = weight,
                 fill = cat)) +
      scale_x_discrete(expand = c(.1, .1), position = "top") +
      geom_flow(show.legend = F) +
      geom_stratum(show.legend = F, color = NA) + 
      scale_fill_manual(values = pal) + 
      # theme_void() +
      theme(axis.text.x = element_text(size = 18, face = "bold"),
            axis.text.y = element_blank(), 
            axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.title = element_blank())
    
    
    flow_data <- 
      ggplot_build(alluv_1)$data[[1]]
    
    stratum_data <- 
      ggplot_build(alluv_1)$data[[2]]
    
    condensed_stratum_data <- 
      stratum_data %>% 
      mutate(label = gsub(data[[subcat]] %>% 
                            unique() %>% 
                            paste(collapse = "$| ") %>% 
                            paste0(" ", ., "$"), 
                          "", 
                          stratum)) %>% 
      group_by(label, x) %>%
      summarise(ymin = min(ymin), 
                ymax = max(ymax), 
                xmin = min(xmin), 
                xmax = max(xmax),
                y = (ymin + ymax) / 2) %>% 
      left_join(cat_pal %>% 
                  enframe(), 
                by = c("label" = "name"))
    
    alluv_1 <- 
      alluv_1 +
      
      # Stratum label
      geom_text(stat = "stratum",
                aes(label = 1), 
                size = 2,
                nudge_x = ifelse(stratum_data$x == 1, 
                                 -1/6 - 1/100,
                                 1/6 + 1/100),
                angle = ifelse(stratum_data$x == 1, 
                               90,
                               -90), 
                vjust = 0) + 
      
      geom_rect(data = condensed_stratum_data,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = F, 
                fill = NA, 
                size = 0.1,
                color = condensed_stratum_data$value) +
      
      geom_text(data = condensed_stratum_data,
                aes(x = x, y = y, label = label), 
                size = 3, 
                inherit.aes = F) + 
      
      # Flow label
      geom_text(stat = "flow", 
                aes(label = 1),
                nudge_x = ifelse(flow_data$x == 1, 
                                 1/6,
                                 -1/6), 
                hjust = ifelse(flow_data$x == 1, 
                               0,
                               1), 
                size = 2) 
    
    
    
    alluv_1
  }

multi_alluvial_plot <- 
  function(data, vars, chunk_levels, pal, color_by = c(1, 3, 3)) {
    
    selvars = vars
    
    if(!is.null(names(vars))) {
      vars = names(vars)
    }
    
    alluv_1 <-
      data %>%
      ungroup() %>%
      select(selvars) %>% 
      ungroup() %>%
      mutate(row_n = row_number()) %>%
      gather(bar, chunk, -row_n) %>%
      left_join(tibble(bar = vars, 
                       color_vars = color_by), 
                by = "bar") %>% 
      group_by(row_n) %>%
      mutate(chunk_color = chunk[match(vars[color_vars], bar)]) %>% 
      ungroup() %>%
      
      mutate(chunk = factor(chunk, levels = chunk_levels),
             bar = factor(bar, levels = vars)) %>%
      
      
      ggplot(aes(x = bar, stratum = chunk, alluvium = row_n,
                 y = 1)) +
      
      geom_flow(aes(fill = chunk_color), 
                show.legend = F) +
      geom_stratum(aes(fill = chunk), 
                   show.legend = F, color = NA) +
      
      scale_x_discrete(expand = c(.1, .1), position = "top") +
      scale_fill_manual(values = pal) + 
      
      
      theme(axis.text.x = element_text(size = 18, face = "bold"),
            axis.text.y = element_blank(), 
            axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.title = element_blank())
    
    
    
    
    flow_data <-
      ggplot_build(alluv_1)$data[[1]] %>%
      as_tibble() %>%
      {
        if("side" %in% names(.)) {
          .
        } else{
          mutate(.,
                 side = case_when(contact == "back" ~ "start",
                                  contact == "front" ~ "end"))
        }}
    
    
    stratum_data <- 
      ggplot_build(alluv_1)$data[[2]]
    
    flow_data_labels <-
      flow_data %>% 
      as_tibble() %>% 
      
      select(x, stratum, group, side, ymin, ymax) %>% 
      pivot_wider(names_from = side, values_from = c(x, stratum, ymin, ymax)) %>%
      
      mutate_at(c("x_end", "ymax_end", "ymin_end", "x_start", "ymax_start", "ymin_start"), as.numeric) %>% 
      group_by(stratum_start, stratum_end, x_start, x_end) %>%
      summarise(y_end = (min(ymin_end) + max(ymax_end)) / 2, 
                y_start = (min(ymin_start) + max(ymax_start)) / 2, 
                size = max(ymax_start) - min(ymin_start))
    
    alluv_1 <- 
      alluv_1 +
      geom_text(data = flow_data_labels,
                aes(x = x_start + 1/6,
                    y = y_start, 
                    label = size), 
                inherit.aes = F, 
                size = 3, 
                hjust = 0) +
      geom_text(data = flow_data_labels,
                aes(x = x_end - 1/6,
                    y = y_end, 
                    label = size), 
                inherit.aes = F, 
                size = 3, 
                hjust = 1) +
      
      # Stratum label
      
      geom_text(data = stratum_data,
                aes(x = x, 
                    y = y,
                    label = paste(stratum, 
                                  ymax - ymin, sep = "\n")), 
                size = 4, 
                inherit.aes = F)
    
    
    alluv_1
  }


## Specialized plots - do not reuse:

classification_network_plot <- 
  function(data_, title = "") {
    enrichment_whole_group <- 
      data_ %>% 
      select(gene, spec_category, enriched_tissues) %>% 
      filter(spec_category %in% c("tissue enriched", "group enriched")) %>% 
      group_by(enriched_tissues, spec_category) %>% 
      summarise(n_genes = n()) %>%
      ungroup()
    
    net_data <-
      enrichment_whole_group %>% 
      mutate(all_enriched_tissues = enriched_tissues) %>%
      separate_rows(enriched_tissues, sep = ";") %>%
      group_by(enriched_tissues, spec_category) %>% 
      mutate(rank = rank(-n_genes, ties.method = "min")) %>%
      group_by(all_enriched_tissues) %>%
      mutate(any_low_rank = any(rank <= 2)) %>%
      ungroup() %>%
      
      filter((spec_category == "tissue enriched" | any_low_rank | n_genes >= 8) &
               (spec_category == "tissue enriched" | n_genes >= 2)) %>% 
      mutate(edge_id = paste("enriched:", all_enriched_tissues)) %>% 
      arrange(n_genes)
    
    net_edges <- 
      net_data %$% 
      tibble(node1 = enriched_tissues, node2 = edge_id, n = n_genes) %>% 
      unique()
    
    g <-
      net_edges %>%
      graph_from_data_frame(directed = FALSE) %>%
      ggraph(layout = "kk") 
    
    link_map <- 
      net_edges %>% 
      gather(node, id, -(3:4)) %>% 
      mutate(tissue_node = node == "node1", 
             color_id = case_when(tissue_node ~ id, 
                                  grepl(";", id) ~ "Group enriched",
                                  !grepl(";", id) ~ "Tissue enriched"), 
             label = ifelse(tissue_node, color_id, n)) %>%
      select(n, node, id, tissue_node, color_id, label) %>%
      unique()
    
    
    edge_data <- get_edges()(g$data)
    node_data <- 
      get_nodes()(g$data) %>% 
      as_tibble() %>%
      left_join(link_map, 
                by = c("name" = "id")) 
    
    g + 
      geom_edge_arc(aes(width = n), 
                    color = "gray", 
                    strength = 0, 
                    show.legend = F) + 
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_width_continuous(range = c(1, 3)) +
      
      geom_node_point(data = node_data  %>%
                        filter(!tissue_node),
                      aes(size = log(n), 
                          fill = color_id), 
                      stroke = 1,
                      # size = 10,
                      shape = 21,
                      show.legend = F)+
      geom_node_point(data = node_data %>%
                        filter(tissue_node),
                      aes(fill = gsub(paste0(" (", paste(LETTERS, collapse = "|"), ")$"), "", color_id)), 
                      stroke = 1,
                      size = 20,
                      shape = 21,
                      show.legend = F)+
      geom_node_text(data = node_data,
                     aes(label = label),
                     size = 4) +
      ggtitle(title) +
      scale_size_continuous(range = c(5, 10)) +
      scale_fill_manual(values = c(cell_type_palette, gene_category_pal)) +
      
      theme_void()
    
  }

classification_network_plot_2 <- 
  function(class_table, gene_col, spec_col, enriched_col, spec_filter, pal, savename, enriched_sep = ";", 
           node_filter_rank = 2, node_filter_show_cat = c("tissue enriched"), 
           node_filter_min = 2, node_filter_n_show = 8, scale_factor = 1) {
    
    enrichment_table <- 
      class_table %>% 
      select(gene = gene_col, 
             spec = spec_col, 
             enriched = enriched_col) %>% 
      filter(spec %in% spec_filter) %>% 
      group_by(enriched, spec) %>% 
      summarise(n_genes = n()) %>%
      ungroup()
    
    net_data <-
      enrichment_table %>% 
      mutate(all_enriched = enriched) %>%
      separate_rows(enriched, sep = enriched_sep) %>%
      group_by(enriched, spec) %>% 
      mutate(rank = rank(-n_genes, ties.method = "min")) %>%
      group_by(all_enriched) %>%
      mutate(any_low_rank = any(rank <= node_filter_rank)) %>%
      ungroup() %>%
      
      filter((spec == node_filter_show_cat | any_low_rank | n_genes >= node_filter_n_show) &
               (spec == node_filter_show_cat | n_genes >= node_filter_min)) %>% 
      mutate(edge_id = paste("enriched:", all_enriched)) %>% 
      arrange(n_genes)
    
    net_edges <- 
      net_data %$% 
      tibble(node1 = enriched, node2 = edge_id, n = n_genes) %>% 
      unique()
    
    g <-
      net_edges %>%
      graph_from_data_frame(directed = FALSE) %>%
      ggraph(layout = "kk") 
    
    link_map <- 
      net_edges %>% 
      gather(node, id, -(3)) %>% 
      mutate(tissue_node = node == "node1", 
             color_id = case_when(tissue_node ~ id, 
                                  grepl(";", id) ~ "Group enriched",
                                  !grepl(";", id) ~ "Tissue enriched"), 
             label = ifelse(tissue_node, color_id, n)) %>%
      select(n, node, id, tissue_node, color_id, label) %>%
      unique()
    
    
    edge_data <- get_edges()(g$data)
    node_data <- 
      get_nodes()(g$data) %>% 
      as_tibble() %>%
      left_join(link_map, 
                by = c("name" = "id")) 
    
    
    fig <- 
      g + 
      geom_edge_arc(aes(width = n), 
                    color = "gray", 
                    strength = 0, 
                    alpha = 0.5,
                    show.legend = F) + 
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_width_continuous(range = c(1, 3)) +
      
      geom_node_point(data = node_data  %>%
                        filter(!tissue_node),
                      aes(size = log(n), 
                          fill = color_id), 
                      stroke = 1,
                      # size = 10,
                      shape = 21,
                      show.legend = F)+
      geom_node_point(data = node_data %>%
                        filter(tissue_node),
                      aes(fill = color_id), 
                      stroke = 1,
                      size = 20 * scale_factor,
                      shape = 21,
                      show.legend = F)+
      geom_node_text(data = node_data,
                     aes(label = label),
                     size = 4 * scale_factor) +
      scale_size_continuous(range = c(5, 10) * scale_factor) +
      scale_fill_manual(values = pal) +
      
      theme_void()
    
    ## ----- Save
    
    
    cyto_summary <- 
      net_edges %>% 
      mutate(category = ifelse(!grepl(enriched_sep, node2), "Tissue enriched", "Group enriched"), 
             node_id = unclass(factor(node2)),
             node1 = node1, 
             n_sqrt = sqrt(n), 
             str_len = str_length(node1)) %>%
      select(category, node1, node2, node_id, n, n_sqrt, str_len) 
    
    cyto_summary %>% 
      write_delim(savepath(paste(savename, "cytoscape nodes summary.txt")), delim = "\t")
    
    bind_rows(cyto_summary %>% 
                left_join(pal %>% 
                            enframe("node1", "color")) %>% 
                select(node_id = node1, 
                       color) %>% 
                unique() %>%
                mutate(node_type = "Tissue"),
              cyto_summary %>% 
                mutate(color = case_when(category == "Tissue enriched" ~ "#e41a1c",
                                         category == "Group enriched" ~ "#FF9D00"), 
                       node_id = as.character(node_id)) %>%
                select(node_id, color) %>% 
                unique() %>%
                
                mutate(node_type = "Enrichment")) %>%
      mutate(color2 = case_when(node_type == "Enrichment" ~ color,
                                node_type == "Tissue" ~ "#D3D3D3FF"),
             color3 = case_when(node_type == "Enrichment" ~ color,
                                node_type == "Tissue" ~ "#BEBEBEFF")) %>% 
      
      write_delim(savepath(paste(savename, "cytoscape nodes color.txt")), delim = "\t") 
    
    bind_rows(cyto_summary %>% 
                select(node_id = node1) %>% 
                mutate(label = node_id) %>%
                unique(),
              cyto_summary %>% 
                mutate(node_id = as.character(node_id), 
                       label = as.character(n)) %>%
                select(node_id, label) %>% 
                unique()) %>% 
      write_delim(savepath(paste(savename, "cytoscape nodes label whole group.txt")), 
                  delim = "\t")
    
    ## ----
    
    fig
  }

single_cell_tissue_overlap_network <- 
  function(data_, tissue_class, con_tis_, title = "") {
    
    
    tissue_enriched_overlap <- 
      tissue_class %>% 
      filter(specificity_category %in% c("Tissue enriched", "Group enriched")) %>%
      mutate(enhanced_tissues = gsub(", ", "_", enhanced_tissues)) %>%
      separate_rows(enhanced_tissues, sep = ",") %>% 
      mutate(enhanced_tissues = gsub("_", ", ", enhanced_tissues)) %>%
      left_join(data_ %>%
                  filter(spec_category %in% c("tissue enriched", "group enriched", "tissue enhanced")) %>%
                  select(gene, spec_category, dist_category, enriched_tissues) %>%
                  mutate(all_enriched_cells = enriched_tissues) %>%
                  separate_rows(enriched_tissues, sep = ";"),
                by = c("ensg_id" = "gene")) %>%
      mutate(enriched_tissues = ifelse(is.na(enriched_tissues), "Not cell enriched", enriched_tissues)) %>%
      rename(enriched_cells = enriched_tissues) %>% 
      filter(enhanced_tissues == con_tis_)
    
    
    tissue_enriched_overlap %>%
      select(ensg_id, tissue_category = specificity_category, cell_category = spec_category, all_enriched_cells) %>%
      left_join(gene_info92 %>%
                  select(ensg_id, gene_name, gene_description)) %>%
      unique %>%
      mutate(cell_category = case_when(is.na(cell_category) ~ "not cell enriched",
                                       T ~ gsub("tissue", "cell", cell_category))) %>%
      write_csv(savepath(paste(title, "overlap enriched genes.csv")))
    
    net_data <-
      tissue_enriched_overlap %>% 
      group_by(enhanced_tissues, enriched_cells, all_enriched_cells) %>% 
      summarise(n = n()) %>% 
      arrange(-n) %>% 
      
      ungroup() %>%
      # filter(n > 3) %>%
      mutate(edge_id = paste("N", enhanced_tissues, all_enriched_cells),
             community = enhanced_tissues,
             cell_node = paste("CELL", enhanced_tissues, enriched_cells)) %>%
      gather(node, id, enhanced_tissues, cell_node) %>% 
      arrange(community)
    
    net_edges <- 
      net_data %$% 
      tibble(node1 = id, node2 = edge_id, n = n, community = community, cell = enriched_cells) %>% 
      unique() 
    
    unique(net_edges$node1)
    
    
    
    g <-
      net_edges %>%
      graph_from_data_frame(directed = FALSE) %>%
      ggraph(layout = "kk") 
    
    
    link_map <- 
      net_edges %>% 
      gather(node, id, 1:2) %>% 
      mutate(node_type = case_when(grepl("^CELL ", id) ~ "cell",
                                   grepl("^N ", id) ~ "n",
                                   T ~ "tissue"), 
             color_id = case_when(node_type %in% c("cell") ~ cell,
                                  node_type %in% c("tissue") ~ community, 
                                  grepl(";", id) ~ "Group enriched",
                                  !grepl(";", id) ~ "Tissue enriched"), 
             label = ifelse(node_type %in% c("cell", "tissue"), color_id, n)) %>%
      unique()
    
    
    edge_data <- 
      get_edges()(g$data) %>%
      as_tibble() %>%
      mutate(name = node1.name) 
    node_data <- 
      get_nodes()(g$data) %>% 
      as_tibble() %>%
      left_join(link_map, 
                by = c("name" = "id")) 
    
    
    g1 <- 
      g + 
      geom_edge_arc(data = edge_data,
                    aes(width = n), 
                    color = "gray", 
                    strength = 0, 
                    show.legend = F, 
                    alpha = 0.5) + 
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_width_continuous(range = c(1, 3)) +
      
      geom_node_point(data = node_data  %>%
                        filter(node_type == "n"),
                      aes(size = log(n), 
                          fill = color_id), 
                      stroke = 1,
                      # size = 10,
                      shape = 21,
                      show.legend = F)+
      geom_node_point(data = node_data %>%
                        filter(node_type != "n"),
                      aes(fill = gsub(paste0(" (", paste(LETTERS, collapse = "|"), ")$"), "", color_id)), 
                      stroke = 1,
                      size = 10,
                      shape = 21,
                      show.legend = F)+
      geom_node_text(data = node_data,
                     aes(label = label),
                     size = 2) +
      ggtitle(title) +
      scale_size_continuous(range = c(5, 10)) +
      scale_fill_manual(values = c(tissue_colors, cell_type_palette, gene_category_pal)) +
      scale_x_continuous(expand = expansion(0.5)) +
      scale_y_continuous(expand = expansion(0.5)) +
      theme_void() + 
      theme(strip.text = element_text(size = 5, face = "bold"))
    
    g1
    
    
  }

# ----- Misc plots -----

pie_plot <- 
  function(plot_data, class_col, y_col, pal, dodge = 0.1){
    plot_data_ <- 
      plot_data %>%
      rename(class = class_col, 
             y = y_col)
    
    if(!is.factor(plot_data$term)) {
      plot_data_ <- 
        plot_data_ %>%
        mutate(class = factor(class))
    }
    plot_data_ %>% 
      arrange(match(class, rev(levels(class)))) %>% 
      mutate(y_stack = cumsum(y) - y/2) %>% 
      {ggplot(., aes(1, y, fill = class, group = class, label = paste0(class, "\n", round(y * 100, 1), " %"))) +
          geom_col(show.legend = F, 
                   color = "white", 
                   width = 1) + 
          geom_segment(aes(x = 1.5 + dodge, xend = 1.5, 
                           y = y_stack, yend = y_stack), size = 0.5) +
          geom_text_repel(aes(x = 1.5 + dodge, y = y_stack), 
                          color = "black", nudge_x = dodge, 
                          segment.size = 0.5) +
          scale_fill_manual(values = pal) + 
          coord_polar("y") +
          theme_void() + 
          scale_x_continuous(expand = expand_scale(c(0,0.8)))}
  }

pie_plot_facet <- 
  function(plot_data, class_col, y_col, facet_col, pal, dodge = 0.1){
    plot_data_ <- 
      plot_data %>%
      rename(class = class_col, 
             y = y_col, 
             facet = facet_col)
    
    if(!is.factor(plot_data$term)) {
      plot_data_ <- 
        plot_data_ %>%
        mutate(class = factor(class))
    }
    
    plot_data_ %>% 
      arrange(match(class, rev(levels(class)))) %>% 
      group_by(facet) %>%
      mutate(y_stack = cumsum(y) - y/2) %>% 
      {ggplot(., aes(1, y, fill = class, group = class, label = paste0(class, "\n", round(y * 100, 1), " %"))) +
          geom_col(show.legend = F, 
                   color = "white", 
                   width = 1) + 
          geom_segment(aes(x = 1.5 + dodge, xend = 1.5, 
                           y = y_stack, yend = y_stack), size = 0.5) +
          geom_text_repel(aes(x = 1.5 + dodge, y = y_stack), 
                          color = "black", nudge_x = dodge, 
                          segment.size = 0.5) +
          scale_fill_manual(values = pal) + 
          
          facet_wrap(~facet) +
          
          coord_polar("y") +
          theme_void() + 
          scale_x_continuous(expand = expand_scale(c(0,0.8)))}
  }


scatter_cor_plot <- 
  function(plot_data, 
           val1 = "nx_human", val2 = "nx_pig", 
           lab1 = "ensg_id", lab2 = "enssscg_id", 
           color_col = "region_tissue_name", 
           pal) {
    
    plot_data_ <- 
      plot_data %>% 
      rename(val1 = val1, 
             val2 = val2, 
             lab1 = lab1, 
             lab2 = lab2, 
             color = color_col) 
    
    cor_data_ <-
      plot_data_ %>% 
      group_by(lab1, lab2) %>% 
      summarise(cor_s = cor(val1, val2, method = "spearman", use = "pairwise.complete.obs"),
                cor_p = cor(log10(val1 + 1), log10(val2 + 1), method = "pearson", use = "pairwise.complete.obs"), 
                max_val = log10(max(c(val1, val2)) + 1),
                min_val = log10(min(c(val1, val2)) + 1))
    
    plot_data_ %>%
      
      {ggplot(., aes(log10(val1 + 1), log10(val2 + 1))) + 
          geom_point(aes(color = color), 
                     show.legend = F) + 
          geom_smooth(method = "lm", color = "black", fill = NA) +
          geom_abline(slope = 1) +
          geom_text(data = cor_data_,
                    aes(label = paste0("spearman ", round(cor_s, 2), "\n",
                                       "pearson ", round(cor_p, 2)),
                        x = min_val, y = max_val), 
                    hjust = 0, vjust = 1) +
          # facet_grid("Human" + gene_name_human ~ "Pig" + gene_name_pig) + 
          facet_grid(~ lab1 + lab2) + 
          xlab(paste0("log10(", lab1, " + 1)")) + 
          ylab(paste0("log10(", lab2, " + 1)")) + 
          stripped_theme_facet + 
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
          scale_color_manual(values = pal) + 
          scale_x_continuous(limits = log10(c(min(.$val1, .$val2), max(.$val1, .$val2)) + 1)) +
          scale_y_continuous(limits = log10(c(min(.$val1, .$val2), max(.$val1, .$val2)) + 1)) +
          coord_fixed()}
  }

# ----- tissue clustering plots -----

evolutionary_tree_plot <- function(data,
                                   dist_fun = function(x) cor(x, method = 'spearman', use="pairwise.complete.obs") %>%
                                     {as.dist(1-.)},
                                   tree_fun = function(x) nj(),
                                   color_mapping,
                                   mapping_col,
                                   color_col) {

  # calculate distance
  if(is.function(dist_fun) & is.function(tree_fun)) {
    data_dist <- dist_fun(data)
  } else {
    data_dist <- dist_fun
  }

  # neighbor-joining tree estimation
  if(is.function(tree_fun)) {
    nj_tree <- tree_fun(data_dist)
  } else {
    nj_tree <- tree_fun
  }



  plot.phylo(as.phylo(nj_tree),
             type="u",
             lab4ut = "axial",
             font = 1,
             cex = 0.8,
             tip.color = color_mapping %>%
               rename(mapping_col = mapping_col,
                      color_col = color_col) %>%
               {.[match(nj_tree$tip.label, .$mapping_col),]} %$%
               color_col,
             edge.col="black",
             edge.width=2,
             show.node.label=TRUE, no.margin=TRUE,
             use.edge.length = T)

}

pca_calc <- function(data, npcs) {

  pca_res <-
    data %>%
    pca(nPcs = npcs)

  pca_stats <-
    tibble(PC = 1:npcs,
           R2cum = R2cum(pca_res))

  informative_pcs <- pca_stats$PC[which(pca_stats$R2cum > 0.95)[1]]

  pca_stats_plot <-
    pca_stats %>%
    select(PC, R2cum) %>%
    ggplot(aes(PC, R2cum)) +
    geom_point() +
    geom_line() +
    simple_theme +
    geom_vline(xintercept = informative_pcs, linetype = "dashed") +
    annotate("text",
             x = informative_pcs,
             y = 0.55,
             label = paste0("PC ", informative_pcs,
                            "\nR2 = ", round(pca_stats[informative_pcs, ]$R2cum, 3)),
             hjust = 1,
             vjust = 0)

  list(pca = pca_res,
       scores = pcaMethods::scores(pca_res),
       loadings = pcaMethods::loadings(pca_res),
       stats = pca_stats,
       pc_95 = informative_pcs,
       plot = pca_stats_plot)
}

pca_transform <- function(data, loadings) {
  scale(as.matrix(data), scale = F) %*% (as.matrix(loadings))
}

pca_score_plot <- function(pca_scores, group_mapping, mapping_col, group_col, pal, xpc = 1, ypc = 2) {

  pca_scores_mapped <-
    pca_scores %>%
    as_tibble(rownames = "mapping_col") %>%
    left_join(group_mapping,
              by = c("mapping_col" = mapping_col)) %>%
    rename(xpc = xpc + 1,
           ypc = ypc + 1)

  if(group_col == mapping_col) {
    pca_scores_mapped <-
      pca_scores_mapped %>%
      mutate(group_col = mapping_col)
  } else {
    pca_scores_mapped <-
      pca_scores_mapped %>%
      rename(group_col = group_col)
  }


  pca_scores_mapped %>%
    group_by(group_col) %>%
    mutate(mean_x = mean(xpc),
           mean_y = mean(ypc)) %>%
    ungroup() %>%

    {ggplot(., aes(xpc, ypc, color = group_col)) +
        geom_point(show.legend = F) +
        scale_color_manual(values = pal) +
        simple_theme +
        geom_segment(aes(xend = mean_x, yend = mean_y),
                     show.legend = F) +
        geom_text(data = select(., group_col, mean_x, mean_y) %>%
                    unique(),
                  aes(x = mean_x, y = mean_y, label = group_col),
                  show.legend = F) +
        xlab(paste0("PC", xpc)) +
        ylab(paste0("PC", ypc))}


}

pca_score_plot_simple <- function(pca_scores, xpc = 1, ypc = 2) {
  
  pca_scores_mapped <-
    pca_scores %>%
    as_tibble(rownames = "names") %>%
    rename(xpc = xpc + 1,
           ypc = ypc + 1)
  
  
  
  pca_scores_mapped %>%
    
    {ggplot(., aes(xpc, ypc, label = names)) +
        geom_point(show.legend = F) +
        simple_theme +
        geom_text(show.legend = F) +
        xlab(paste0("PC", xpc)) +
        ylab(paste0("PC", ypc))}
  
  
}

pca_R2_plot <- function(pca_stats, mark_pc) {
  pca_stats %>% 
    ggplot(aes(PC, R2cum)) +
    geom_point() + 
    geom_line() + 
    simple_theme +
    geom_vline(xintercept = mark_pc, linetype = "dashed") + 
    annotate("text", 
             x = mark_pc,
             y = pca_stats$R2cum[mark_pc], 
             label = paste0("PC ", mark_pc, ": R2 = ", round(pca_stats$R2cum[mark_pc], 3)),
             hjust = -0.1, 
             vjust = -0.5, 
             angle = -90) 
}

plot_dendrogram <- 
  function(clust, 
           color_mapping, label_col, color_col, pal, do_label = F) {
    
    dendr <- dendro_data(clust)
    
    
    dendro_plot_data <- 
      left_join(dendr$segments, 
                dendr$labels, 
                by = c("x" = "x", "yend" = "y")) 
    
    dendro_plot <- 
      dendro_plot_data %>%
      ggplot() +
      geom_segment(aes(x=y, y=x, xend=yend, yend=xend, group = label))+
      
      scale_color_manual(values = pal)+
      scale_fill_manual(values = pal)+
      scale_x_reverse(expand = expand_scale(mult = 0.25))+
      theme(axis.text.y = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(1,1,1,1), units = "mm"), 
            panel.background = element_blank()) 
    
    if(do_label) {
      dendro_plot + 
        geom_label(data = ggdendro::label(dendr) %>%
                    as_tibble(),
                  aes(label=label,
                      x=0,
                      y=x,
                      fill = label),
                  size = 3,
                  hjust = 0,
                  label.r = unit(0, units = "cm"),
                  nudge_x = 0.005,
                  show.legend = F)
    } else {
      dendro_plot + 
        geom_text(data = ggdendro::label(dendr) %>%
                    as_tibble(),
                  aes(label=label,
                      x=0,
                      y=x,
                      color = label),
                  size = 3,
                  hjust = 0,
                  nudge_x = 0.005,
                  show.legend = F)
    }
    
  }

umap_calc <- function(data, npcs = 2, n_epochs = 400, n_neighbors = round(sqrt(dim(data)[1]))) {
  
  data %>%
    umap(n_components = npcs,
         n_epochs = n_epochs,
         n_neighbors = n_neighbors)
  
}

umap_score_plot <- function(umap_scores, group_mapping, mapping_col, group_col, pal, xpc = 1, ypc = 2, repel = F) {
  
  umap_scores_mapped <-
    umap_scores %>%
    as_tibble(rownames = "mapping_col") %>%
    left_join(group_mapping,
              by = c("mapping_col" = mapping_col)) %>%
    rename(xpc = xpc + 1,
           ypc = ypc + 1)
  
  if(group_col == mapping_col) {
    umap_scores_mapped <-
      umap_scores_mapped %>%
      mutate(group_col = mapping_col)
  } else {
    umap_scores_mapped <-
      umap_scores_mapped %>%
      rename(group_col = group_col)
  }
  
  
  g <- 
    umap_scores_mapped %>%
    group_by(group_col) %>%
    mutate(mean_x = mean(xpc),
           mean_y = mean(ypc)) %>%
    ungroup() %>%
    
    {ggplot(., aes(xpc, ypc, color = group_col)) +
        geom_point(show.legend = F) +
        scale_color_manual(values = pal) +
        simple_theme +
        geom_segment(aes(xend = mean_x, yend = mean_y),
                     show.legend = F) +
        
        xlab(paste0("V", xpc)) +
        ylab(paste0("V", ypc))}
  
  if(repel) {
    g + 
      umap_scores_mapped %>%
      group_by(group_col) %>%
      mutate(mean_x = mean(xpc),
             mean_y = mean(ypc)) %>%
      ungroup() %>%
      {geom_text_repel(data = select(., group_col, mean_x, mean_y) %>%
                         unique(),
                       aes(x = mean_x, y = mean_y, label = group_col),
                       show.legend = F)}
  } else {
    g + 
      umap_scores_mapped %>%
      group_by(group_col) %>%
      mutate(mean_x = mean(xpc),
             mean_y = mean(ypc)) %>%
      ungroup() %>%
      {geom_text(data = select(., group_col, mean_x, mean_y) %>%
                   unique(),
                 aes(x = mean_x, y = mean_y, label = group_col),
                 show.legend = F)}
  }
  
  
}

# ----- correlation plots -----

basic_network_plot <- function(data, 
                         cor_col, 
                         var1_col, var2_col, 
                         var1_label_col, var2_label_col, 
                         var1_fill_col, var2_fill_col, 
                         fill_pal = NA, 
                         var1_color = "black",
                         var2_color = "black",
                         scale_edge = F, 
                         edge_color = "black") {
  
  if(is.na(fill_pal)) fill_pal <- rep("white", 1000)
  
  cor_links <- 
    data %>%
    do(tibble(var1 = .[,var1_col][[1]],
              var2 = .[,var2_col][[1]], 
              cor = .[,cor_col][[1]], 
              var1_label = .[,var1_label_col][[1]], 
              var2_label = .[,var2_label_col][[1]], 
              var1_fill = .[,var1_fill_col][[1]], 
              var2_fill = .[,var2_fill_col][[1]]))
  
  link_mapping <- 
    cor_links %>% 
    do(bind_rows(select(., 
                        var = var1, 
                        label = var1_label, 
                        fill = var1_fill) %>% 
                   mutate(var_n = "var1"), 
                 select(., 
                        var = var2, 
                        label = var2_label, 
                        fill = var2_fill) %>%
                   mutate(var_n = "var2"))) %>% 
    unique()

  g <-
    cor_links %>%
    graph_from_data_frame(directed = FALSE) %>%
    ggraph(layout = "kk") 
  
  edge_data <- get_edges()(g$data)
  node_data <- 
    get_nodes()(g$data) %>% 
    left_join(link_mapping, 
              by = c("name" = "var"))
  
  if (scale_edge) {
    g <-
      g + 
      geom_edge_fan(aes(width = cor, 
                        alpha = cor), 
                    color = edge_color) + 
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_width_continuous(range = c(1, 3))
    
  } else {
    g <- 
      g + 
      geom_edge_fan(color = edge_color)
  }
  
  
  g + 
    geom_node_point(data = node_data,
                    aes(fill = fill, 
                        color = var_n), 
                    stroke = 2,
                    size = 10,
                    shape = 21,
                    show.legend = F)+
    geom_node_text(data = node_data,
                   aes(label = label),
                   size = 4) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = c("var1" = var1_color, "var2" = var2_color)) +
      # scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("dodgerblue2", "white", "firebrick2")) +
      
      
    theme_void()
  
  
  
    
}

# ----- tissue grouping plots -----

tissue_connection_plot <- 
  function(mapping_table, col1, col2, col1_name, col2_name, title, pal = rep("black", 1000), na.rm = F) {
    tissue_overlap <- 
      mapping_table %>% 
      rename(col1 = col1,
             col2 = col2)
    
    if(na.rm) {
      tissue_overlap <- 
        tissue_overlap %>%
        filter(!(is.na(col1) | is.na(col2)))
    }
    tissue_overlap <- 
      tissue_overlap %>% 
      select(col1, col2) %>% 
      distinct() %>% 
      mutate(col1_i = unclass(factor(col1, levels = unique(col1))),
             col2_i = unclass(factor(col2, levels = unique(col2)))) %>% 
      group_by(col2) %>%
      mutate(col2_y = case_when(!is.na(col2) ~ mean(col1_i))) %>% 
      ungroup()
    
    
    
    tissue_overlap  %>%
      ggplot(aes(x = 1, xend = 2, 
                 y = col1_i, yend = col2_y)) + 
      annotate("text",
               x = 1:2,
               y = max(tissue_overlap$col2_y) + 1,
               label = c(col1_name, col2_name), 
               vjust = 0, 
               hjust = 1:0,
               fontface = "bold", 
               size = 4) +
      geom_segment(aes(color = col1), 
                   size = 2, 
                   alpha = 1, 
                   show.legend = F) + 
      geom_text(aes(x = 1, y = col1_i, label = col1, color = col1), 
                inherit.aes = F,
                hjust = 1, 
                show.legend = F) + 
      geom_text(aes(x = 2, y = col2_y, label = col2, color = col1), 
                inherit.aes = F, 
                hjust = 0, 
                show.legend = F) + 
      theme_void() + 
      scale_x_continuous(expand = expand_scale(3)) + 
      scale_color_manual(values = pal) + 
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5))
  }

tissue_connection_multi_plot <- 
  function(tissue_overlap, col_names, title, pal = rep("black", 1000), na.rm = F, expand_x = 0.5) {
    
    cols_ <- grep("^col\\d$", names(tissue_overlap), value = T)
    if(na.rm) {
      tissue_overlap <- 
        tissue_overlap[complete.cases(tissue_overlap),]
    }
    
    tissue_overlap <- 
      tissue_overlap %>% 
      select(cols_) %>% 
      distinct() %>% 
      mutate_at(.vars = set_names(cols_, paste0(cols_, "_i")), 
                .funs = function(x) unclass(factor(x, levels = unique(x)))) 
    
    n_dist <- 
      paste0("max(c(", paste0("tissue_overlap$", cols_, "_i", collapse = ", "), "), na.rm = T)") %>% 
      parse(text = .) %>%
      eval()
    
    tissue_overlap <- 
      tissue_overlap %>% 
      mutate_at(.vars = paste0(cols_, "_i"), 
                .funs = function(x) scales::rescale(x, c(1,n_dist))) 

             
    paste0("tissue_overlap  %>%
      mutate(", paste0(cols_, "_i = as.numeric(", cols_, "_i)", collapse = ", "), ") %>%
      ggplot() + 
      annotate('text',
               x = 1:", length(cols_), ",
               y = max(c(", paste0("tissue_overlap$", cols_, "_i", collapse = ", "), "), na.rm = T) + 1,
               label = c(col_names), 
               vjust = 0, 
               hjust = 0.5,
               fontface = 'bold', 
               size = 4) +", 
           paste0("geom_segment(data = tissue_overlap,
                   aes(x = ", 1:(length(cols_) - 1), ", xend = ", 2:(length(cols_)), ",
                       y = col", 1:(length(cols_) - 1), "_i, yend = col", 2:(length(cols_)), "_i,
                       color = col", 1:(length(cols_) - 1), "), 
                   size = 2, 
                   alpha = 1, 
                   show.legend = F) + ", collapse = "\n"), 
           paste0("geom_label(data = tissue_overlap, 
                 aes(x = ", 1:length(cols_), ", y = col", 1:length(cols_), "_i, 
                  label = col", 1:length(cols_), ", fill = col", 1:length(cols_), "),
                 inherit.aes = F,
                 hjust = 0.5,
                 show.legend = F) +", collapse = "\n"), 
           "theme_void() + 
      scale_x_continuous(expand = expand_scale(expand_x)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5))") %>%  
      parse(text = .) %>%
      eval()
    
  }

# ----- ggplot2 utility -----

label_unique <-  function (labels) {
  apply(labels, MARGIN = 1, 
        FUN = function(row_) unique(row_) %>% paste(collapse = " - ")) %>%
    list()
}

# ----- Atlas Demo Plots -----

ortholog_connection_plot <- 
  function(orth_connection_data, example_gene) {
    example_community <- 
      orth_connection_data %>% 
      filter(enssscg_id == example_gene | ensg_id == example_gene) %>% 
      pull(community) %>% 
      unique()
    
    connection_data <- 
      orth_connection_data %>% 
      filter(community == example_community) %>% 
      mutate(human_label = ifelse(is.na(human_gene_name), 
                                  ensg_id, 
                                  human_gene_name),
             pig_label = ifelse(is.na(pig_gene_name), 
                                enssscg_id, 
                                pig_gene_name),
             
             left_id = case_when(grepl("^ENSG", example_gene) ~ ensg_id,
                                 grepl("^ENSSSCG", example_gene) ~ enssscg_id),
             right_id = case_when(grepl("^ENSG", example_gene) ~ enssscg_id,
                                  grepl("^ENSSSCG", example_gene) ~ ensg_id),
             
             selected_gene = example_gene == enssscg_id | example_gene == ensg_id, 
             selected_primary = selected_gene & ((grepl("^ENSG", example_gene) & 
                                                    ortholog_class_human == "primary") |
                                                   (grepl("^ENSSSCG", example_gene) & 
                                                      ortholog_class_pig == "primary")),
             
             left_factor = unclass(factor(left_id, levels = unique(left_id[order(expression_identity)]))),
             right_factor = unclass(factor(right_id, levels = unique(right_id[order(expression_identity)]))),
             
             primary_ortholog_direction = case_when(mutual_primary ~ 0,
                                                    
                                                    grepl("^ENSG", example_gene) & 
                                                      ortholog_class_human == "primary" &
                                                      ortholog_class_pig == "secondary" ~ 1,
                                                    grepl("^ENSG", example_gene) & 
                                                      ortholog_class_human == "secondary" &
                                                      ortholog_class_pig == "primary" ~ -1,
                                                    
                                                    grepl("^ENSSSCG", example_gene) & 
                                                      ortholog_class_human == "primary" &
                                                      ortholog_class_pig == "secondary" ~ -1,
                                                    grepl("^ENSSSCG", example_gene) & 
                                                      ortholog_class_human == "secondary" &
                                                      ortholog_class_pig == "primary" ~ 1),
             
             left_y = scales::rescale(left_factor, 
                                      c(1, unique(community_size))),
             right_y = scales::rescale(right_factor, 
                                       c(1, unique(community_size))), 
             
             
             left_label = case_when(grepl("^ENSG", example_gene) ~ human_label,
                                    grepl("^ENSSSCG", example_gene) ~ pig_label),
             right_label = case_when(grepl("^ENSG", example_gene) ~ pig_label,
                                     grepl("^ENSSSCG", example_gene) ~ human_label))
    
    connection_data_primary <- 
      connection_data %>% 
      filter(selected_primary)
    
    connection_data %>%
      ggplot() + 
      
      geom_segment(aes(x = 1, xend = 2, 
                       y = left_y, 
                       yend= right_y, 
                       linetype = ortholog_confidence,
                       # size = expression_identity, 
                       color = selected_primary),
                   # size = 1,
                   inherit.aes = F, 
                   show.legend = F, 
                   
                   arrow = arrow(length = case_when(is.na(connection_data$primary_ortholog_direction) ~ 0, 
                                                    T ~ 0.1) %>% 
                                   unit("inches"),
                                 
                                 ends = with(connection_data, 
                                             case_when(primary_ortholog_direction == 0 ~ "both", 
                                                       primary_ortholog_direction == -1 ~ "first", 
                                                       primary_ortholog_direction == 1 ~ "last", 
                                                       T ~ "first")), 
                                 type = "open")) + 
      
      geom_segment(data = connection_data_primary,
                   aes(x = 1, xend = 2, 
                       y = left_y, 
                       yend= right_y, 
                       linetype = ortholog_confidence,
                       # size = expression_identity, 
                       color = selected_primary),
                   # size = 1,
                   inherit.aes = F, 
                   show.legend = F, 
                   
                   arrow = arrow(length = case_when(is.na(connection_data_primary$primary_ortholog_direction) ~ 0, 
                                                    T ~ 0.1) %>% 
                                   unit("inches"),
                                 
                                 ends = with(connection_data_primary, 
                                             case_when(primary_ortholog_direction == 0 ~ "both", 
                                                       primary_ortholog_direction == -1 ~ "first", 
                                                       primary_ortholog_direction == 1 ~ "last", 
                                                       T ~ "first")), 
                                 type = "open")) + 
      
      geom_text(data = connection_data %>% 
                  select(left_label, left_id, left_y, selected_gene) %>% 
                  group_by(left_label, left_id) %>% 
                  summarise(left_y = unique(left_y), 
                            selected_gene = any(selected_gene)),
                
                aes(x = 1, y = left_y, 
                    label = left_label, 
                    color = selected_gene), 
                show.legend = F, 
                hjust = 1.1) + 
      
      geom_text(data = connection_data %>% 
                  select(right_label, right_id, right_y, selected_primary) %>% 
                  group_by(right_label, right_id) %>% 
                  summarise(right_y = unique(right_y), 
                            selected_primary = any(selected_primary)),
                aes(x = 2, y = right_y, 
                    label = right_label, 
                    color = selected_primary), 
                show.legend = F, 
                hjust = -0.1) + 
      
      
      geom_point(data = connection_data,
                 aes(x = 1.5, y = (left_y + right_y) / 2, 
                     color = selected_gene), 
                 show.legend = F, 
                 fill = "white",
                 shape = 21,
                 size = 8) + 
      
      geom_text(data = connection_data,
                aes(x = 1.5, y = (left_y + right_y) / 2, 
                    label = round(expression_identity * 100, 1), 
                    color = selected_gene), 
                show.legend = F, 
                size = 3,
                hjust = 0.5,
                vjust = 0.45) + 
      
      
      theme_void() +
      stripped_theme_HPA + 
      
      scale_linetype_manual(values = c("1" = "solid", "0" = "dashed")) +
      scale_x_continuous(expand = expand_scale(mult = 0.9)) + 
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray")) +
      scale_size_continuous(range = c(1, 2), limits = 0:1)
  }

plot_ortholog_pair <- 
  function(pig_id, human_id) {
    
    plot_data <- 
      ortholog_connection_data %>% 
      rename(pig = 1, 
             human = 2) %>%
      filter(pig == pig_id & 
               human == human_id) %>%
      gather(species, gene, 1:2) %>% 
      select(species, gene) %>% 
      unique() %>% 
      left_join(bind_rows(pig_atlas_region %>% 
                            rename(gene = 1) %>% 
                            select(gene, region_tissue_name, nx),
                          human_atlas_region %>% 
                            rename(gene = 2) %>% 
                            select(gene, region_tissue_name, nx)),
                by = "gene") %>% 
      left_join(bind_rows(gene_mapping %>% 
                            select(gene = 1, 
                                   name = 2),
                          human_gene_mapping %>% 
                            select(gene = 1, 
                                   name = 2)),
                by = "gene") %>% 
      mutate(gene_name = case_when(name == "NULL" ~ gene, 
                                   T ~ paste(name,  "-",  gene)), 
             region_tissue_name = factor(region_tissue_name, levels = region_levels)) %>%
      left_join(comparison_table %>% 
                  select(1, 2, spec_category_human, spec_category_pig) %>% 
                  gather("species", "enrichment", -1, -2) %>% 
                  mutate(species = gsub("spec_category_", "", species),
                         gene = ifelse(species == "human", ensg_id, enssscg_id)), 
                by = c("gene", "species"))
    
    g1 <- 
      plot_data %>%
      ggplot(aes(region_tissue_name, nx, fill = region_tissue_name)) + 
      geom_col(show.legend = F) + 
      facet_wrap(paste(species, "-", enrichment) ~ gene_name, scales = "free", ncol = 1) + 
      stripped_theme_HPA + 
      ylab("NX") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), 
            strip.text.x = element_text(hjust = 0), 
            axis.title.x = element_blank()) + 
      scale_fill_manual(values = region_palette) + 
      scale_y_continuous(expand = expand_scale(c(0, 0.05)))
    
    plot_data <- 
      ortholog_connection_data %>% 
      rename(pig = 1, 
             human = 2) %>%
      
      filter(pig == pig_id & 
               human == human_id) %>%
      left_join(orth_region, 
                by = c("pig" = "enssscg_id", 
                       "human" = "ensg_id")) %>% 
      left_join(gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("pig" = "gene")) %>% 
      left_join(human_gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("human" = "gene"), 
                suffix = c("_pig", "_human")) %>% 
      mutate(gene_name_pig = case_when(name_pig == "NULL" ~ pig, 
                                       T ~ paste(name_pig,  "-",  pig)),
             gene_name_human = case_when(name_human == "NULL" ~ human, 
                                         T ~ paste(name_human,  "-",  human)))
    
    g2 <- 
      plot_data %>%
      {ggplot(., aes(log2(nx_pig + 1), log2(nx_human + 1))) + 
          geom_point(aes(color = region_tissue_name), 
                     show.legend = F) + 
          geom_smooth(method = "lm", color = "black", fill = NA) +
          geom_abline(slope = 1) +
          # facet_grid("Human" + gene_name_human ~ "Pig" + gene_name_pig) + 
          facet_grid("Human" ~ "Pig") + 
          xlab("log2(NX + 1)") + 
          ylab("log2(NX + 1)") +
          stripped_theme_HPA + 
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
          scale_color_manual(values = tissue_colors_palette_full_humanpig) + 
          scale_x_continuous(limits = log2(c(min(.$nx_pig, .$nx_human), max(.$nx_pig, .$nx_human)) + 1)) +
          scale_y_continuous(limits = log2(c(min(.$nx_pig, .$nx_human), max(.$nx_pig, .$nx_human)) + 1)) +
          coord_fixed() + 
          with(comparison_table %>% 
                 filter(enssscg_id == pig_id & 
                          ensg_id == human_id),
               ggtitle(paste0("Spearman: ", round(cor_spearman, 2), ", Pearson: ", round(cor_pearson, 2))))}
    
    
    grid.arrange(g1, g2, widths = c(2, 1))
  }

ortholog_network_plot <- 
  function(ortholog_connection_data, gene, species) {
    connection_data <- 
      ortholog_connection_data %>% 
      rename(pig = 1, 
             human = 2) %>%
      mutate(selected_species = eval(parse(text = species)))  %>%
      
      # Get edges in community of selected gene
      filter(community == .$community[match(gene, selected_species)]) %>% 
      
      mutate(n_genes = max(c(n_distinct(.$human), n_distinct(.$pig))), 
             pig_factor = unclass(factor(pig)), 
             pig_y = scales::rescale(pig_factor, 
                                     c(1, unique(n_genes))), 
             selected_gene = gene == selected_species) %>%
      group_by(pig) %>%
      mutate(human_order = 1:length(human)) %>% 
      ungroup() %>%
      arrange(pig_y, human_order) %>%
      mutate(human_factor = unclass(factor(human, levels = unique(human))),
             human_y = scales::rescale(human_factor, 
                                       c(1, unique(n_genes)))) 
    
    
    plot_connection_data <- 
      connection_data %>%
      mutate(x = 1, 
             xend = 2, 
             y = pig_y, 
             yend = human_y) %>% 
      select(x, xend, 
             y, yend, 
             selected_gene, 
             cor)
    
    
    plot_label_data <- 
      connection_data %>% 
      gather(species, ID, 1:2) %>%
      mutate(x = case_when(species == "pig" ~ 1,
                           species == "human" ~ 2), 
             y = case_when(species == "pig" ~ pig_y,
                           species == "human" ~ human_y)) %>% 
      select(x, y, ID, species, selected_gene) %>% 
      group_by(ID) %>%
      mutate(selected_gene = any(selected_gene)) %>%
      ungroup() %>%
      unique() %>%
      left_join(gene_mapping %>% 
                  select(1, 
                         pig_name = display_name) %>% 
                  filter(pig_name != "NULL"),
                by = c("ID" = "enssscg_id")) %>% 
      left_join(human_gene_mapping %>% 
                  select(1, 
                         human_name = gene_name),
                by = c("ID" = "ensg_id")) %>% 
      mutate(use_ID = case_when(species == "human" & is.na(human_name) ~ ID, 
                                species == "human" & !is.na(human_name) ~ human_name,
                                species == "pig" & is.na(pig_name) ~ ID, 
                                species == "pig" & !is.na(pig_name) ~ pig_name))
    
    plot_label_data %>%
      ggplot(aes(x, y, label = use_ID, alpha = selected_gene)) + 
      geom_text(aes(hjust = 2 - x), 
                show.legend = F) + 
      geom_segment(data = plot_connection_data, 
                   aes(x = x, xend = xend, 
                       y = y, yend= yend, 
                       linetype = selected_gene, 
                       color = cor),
                   size = 1,
                   inherit.aes = F, 
                   show.legend = F) + 
      scale_x_continuous(expand = expand_scale(mult = 0.9)) + 
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
      theme_void() +
      stripped_theme_HPA + 
      scale_color_gradient2(low = "blue", 
                            mid = "gray80", 
                            high = "red",
                            midpoint = 0, limits = c(-1, 1)) + 
      # scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray")) + 
      scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = F) 
    
  }

ortholog_heatmap <- 
  function(ortholog_connection_data, gene, species) {
    
    
    
    plot_data <- 
      ortholog_connection_data %>% 
      rename(pig = 1, 
             human = 2) %>%
      mutate(selected_species = eval(parse(text = species)))  %>%
      
      # Get edges in community of selected gene
      filter(community == .$community[match(gene, selected_species)]) %>%
      left_join(orth_region, 
                by = c("pig" = "enssscg_id", 
                       "human" = "ensg_id")) %>% 
      left_join(gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("pig" = "gene")) %>% 
      left_join(human_gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("human" = "gene"), 
                suffix = c("_pig", "_human")) %>% 
      mutate(gene_name_pig = case_when(name_pig == "NULL" ~ pig, 
                                       T ~ paste(name_pig,  "-",  pig)),
             gene_name_human = case_when(name_human == "NULL" ~ human, 
                                         T ~ paste(name_human,  "-",  human)))
    
    plot_data2 <- 
      plot_data %>% 
      select(nx_pig, nx_human, gene_name_human, gene_name_pig, region_tissue_name) %>% 
      gather(species, nx, 1:2) %>% 
      mutate(species = gsub("nx_", "", species),
             gene_name = ifelse(species == "pig", gene_name_pig, gene_name_human))  
      
    plot_data2 %>% 
      select(region_tissue_name, nx, gene_name) %>% 
      mutate(nx = log10(nx + 1)) %>%
      unique() %>%
      spread(gene_name, nx) %>% 
      column_to_rownames("region_tissue_name") %>% 
      pheatmap(color = heatmap_palette, 
               clustering_method = "ward.D2", 
               annotation_col = plot_data2 %>% 
                 select(gene_name, species) %>% 
                 unique() %>% 
                 column_to_rownames("gene_name"))
    
    
    
    
  }

ortholog_barplot <- 
  function(ortholog_connection_data, gene, species) {
    
    
     
      plot_data <- 
        ortholog_connection_data %>% 
        rename(pig = 1, 
               human = 2) %>%
        mutate(selected_species = eval(parse(text = species)))  %>%
        
        # Get edges in community of selected gene
        filter(community == .$community[match(gene, selected_species)]) %>%
        gather(species, gene, 1:2) %>% 
        select(species, gene) %>% 
        unique() %>% 
        left_join(bind_rows(pig_atlas_region %>% 
                              rename(gene = 1) %>% 
                              select(gene, region_tissue_name, nx),
                            human_atlas_region %>% 
                              rename(gene = 2) %>% 
                              select(gene, region_tissue_name, nx)),
                  by = "gene") %>% 
        left_join(bind_rows(gene_mapping %>% 
                              select(gene = 1, 
                                     name = 2),
                            human_gene_mapping %>% 
                              select(gene = 1, 
                                     name = 2)),
                  by = "gene") %>% 
        mutate(gene_name = case_when(name == "NULL" ~ gene, 
                                     T ~ paste(name,  "-",  gene)), 
               region_tissue_name = factor(region_tissue_name, levels = region_levels))
      
      lapply(unique(plot_data$gene),
             function(gn) {
               plot_data %>%
                 filter(gene == gn) %>%
                 ggplot(aes(region_tissue_name, nx, fill = region_tissue_name)) + 
                 geom_col(show.legend = F) + 
                 facet_grid(species ~ gene_name, scales = "free") + 
                 stripped_theme_HPA + 
                 ylab("NX") +
                 theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), 
                       strip.text.x = element_text(hjust = 0), 
                       axis.title.x = element_blank()) + 
                 scale_fill_manual(values = region_palette) + 
                 scale_y_continuous(expand = expand_scale(c(0, 0.05)))
             }) %>% 
        set_names(unique(plot_data$gene))
      
      
      
      
  }

ortholog_scatterplot <- 
  function(ortholog_connection_data, gene, species) {
    
    
    
    plot_data <- 
      ortholog_connection_data %>% 
      rename(pig = 1, 
             human = 2) %>%
      mutate(selected_species = eval(parse(text = species)))  %>%
      
      # Get edges in community of selected gene
      filter(community == .$community[match(gene, selected_species)]) %>%
      left_join(orth_region, 
                by = c("pig" = "enssscg_id", 
                       "human" = "ensg_id")) %>% 
      left_join(gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("pig" = "gene")) %>% 
      left_join(human_gene_mapping %>% 
                  select(gene = 1, 
                         name = 2), 
                by = c("human" = "gene"), 
                suffix = c("_pig", "_human")) %>% 
      mutate(gene_name_pig = case_when(name_pig == "NULL" ~ pig, 
                                       T ~ paste(name_pig,  "-",  pig)),
             gene_name_human = case_when(name_human == "NULL" ~ human, 
                                         T ~ paste(name_human,  "-",  human)))
      
      
    
    lapply(unique(plot_data$human),
           function(gn) {
             plot_data %>%
               filter(human == gn) %>%
               {ggplot(., aes(log2(nx_pig + 1), log2(nx_human + 1))) + 
                   geom_point(aes(color = region_tissue_name), 
                              show.legend = F) + 
                   geom_smooth(method = "lm", color = "black", fill = NA) +
                   geom_abline(slope = 1) +
                   # facet_grid("Human" + gene_name_human ~ "Pig" + gene_name_pig) + 
                   facet_grid("Human" ~ "Pig") + 
                   xlab("log2(NX + 1)") + 
                   ylab("log2(NX + 1)") +
                   stripped_theme_HPA + 
                   theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
                   scale_color_manual(values = tissue_colors_palette_full_humanpig) + 
                   scale_x_continuous(limits = log2(c(min(.$nx_pig, .$nx_human), max(.$nx_pig, .$nx_human)) + 1)) +
                   scale_y_continuous(limits = log2(c(min(.$nx_pig, .$nx_human), max(.$nx_pig, .$nx_human)) + 1)) +
                   coord_fixed()}
               
               
           }) %>% 
      set_names(unique(plot_data$human))
    
  }

ortholog_metrics <- 
  function(ortholog_connection_data, gene, species) {
    
    plot_data <- 
      ortholog_connection_data %>% 
      select(-cor) %>%
      left_join(comparison_table, 
                by = c("enssscg_id", "ensg_id")) %>% 
      rename(pig = 1, 
             human = 2) %>%
      mutate(selected_species = eval(parse(text = species)))  %>%
      
      # Get edges in community of selected gene
      filter(community == .$community[match(gene, selected_species)]) %>% 
      gather(cor_type, cor, cor_pearson, cor_spearman) %>% 
      mutate(y = unclass(factor(cor_type)), 
             cor_type = gsub("cor_", "", cor_type))
    
    
    
    lapply(unique(plot_data$human),
           function(gn) {
             plot_data %>%
               filter(human == gn) %>%
               {ggplot(.) + 
                   geom_text(aes(x = 1, 
                                 y = 1, 
                                 vjust = y,
                                 label = paste0(cor_type, ":", round(cor, 2))),
                             hjust = 0) +
                   scale_x_continuous(limits = c(1,2))+
                   stripped_theme_HPA + 
                   theme_void() + 
                   theme(plot.background = element_rect(fill ="#D9D9D9"))}
             }) %>% 
      set_names(unique(plot_data$human))
    
  }

quick_atlas_barplot <- 
  function(data, 
           tis_col = "comparison_tissue", exp_col = "nx", 
           gene_name_col = "gene_name", gene_id_col = "species_gene_id") {
    data %>%
      rename(tissue = tis_col, 
             exp = exp_col, 
             gene_name = gene_name_col, 
             species_gene_id = gene_id_col) %>%
      ggplot(aes(tissue, 
                 exp,
                 fill = tissue)) + 
      geom_col(show.legend = F) + 
      facet_wrap(~ paste(gene_name, "-", species_gene_id), ncol = 1) +
      stripped_theme_facet + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
      scale_fill_manual(values = tissue_colors_palette_full_humanpig)
  }

quick_atlas_fishboneplot <- 
  function(data, orth_connection_data, example_gene) {
    
    connection_data <- 
      orth_connection_data %>% 
      filter(enssscg_id == example_gene | ensg_id == example_gene) 
    
    selected_species <- 
      case_when(grepl("^ENSG", example_gene) ~ "human",
                grepl("^ENSSSCG", example_gene) ~ "pig")
    
    compare_species <- 
      case_when(selected_species == "human" ~ "pig",
                selected_species == "pig" ~ "human")
    
    x_order <- 
      data %>%
      filter(gene_id %in% c(connection_data$enssscg_id, connection_data$ensg_id)) %>%
      mutate(selected_gene = gene_id == example_gene) %>% 
      select(nx, comparison_tissue, selected_gene) %>%
      arrange(!selected_gene, -nx) %>% 
      filter(nx != 0) %>%
      pull(comparison_tissue) %>%
      as.character() %>%
      c(unique(as.character(data$comparison_tissue))) %>%
      unique() 
      
    plot_data <- 
      data %>% 
      filter(gene_id %in% c(connection_data$enssscg_id, 
                            connection_data$ensg_id)) %>% 
      select(1, 2, human_gene_name, comparison_tissue, nx, species) %>%
      
      filter(!is.na(comparison_tissue)) %>% 
      select(-gene_id) %>%
      
      spread(species, nx) %>% 
      filter(complete.cases(.)) %>%
      
      mutate(comparison_tissue = factor(comparison_tissue, levels = x_order), 
             x = unclass(comparison_tissue)) %>% 
      rename(upper_nx = selected_species,
             lower_nx = compare_species) %>%
      separate(mutual_id, into = c("enssscg_id", "ensg_id"), remove = F) %>%
      mutate(label = human_gene_name)
    
    plot_range <- 
      c(-1, 1) *
      max(c(plot_data$upper_nx, 
            plot_data$lower_nx))
    
    plot_data %>%
      ggplot() + 
      geom_segment(aes(x = comparison_tissue,
                       xend = comparison_tissue,
                       y = 0, 
                       yend = upper_nx,
                       color = comparison_tissue), 
                   show.legend = F, 
                   size = 4) + 
      geom_segment(aes(x = comparison_tissue,
                       xend = comparison_tissue,
                       y = 0, 
                       yend = -lower_nx, 
                       color = comparison_tissue),
                   show.legend = F, 
                   size = 4) + 
      
      geom_hline(yintercept = 0, color = "darkgrey") +
      # geom_text(aes(x = x,
      #               y = upper_nx,
      #               label = comparison_tissue), 
      #           angle = 60, 
      #           hjust = -0.1, 
      #           vjust = 0) +
      
      ylab("NX") +
      
      facet_wrap(~ label + ensg_id, ncol = 1, strip.position = "right") +
      
      scale_color_manual(values = tissue_colors_palette_full_humanpig)  +
      scale_y_continuous(limits = plot_range) +
      
      stripped_theme_facet + 
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.line.x = element_blank(),
            panel.border = element_blank(), 
            axis.text.x = element_text(angle = 60, hjust = 1))
    
    
      
  }
# ----- Not currently in use -----





tissue_distributions_plot <- function(data, expression_col, content_col, det_lim, pal, do.tissues = "all") {

  data %>%
    filter(consensus_tissue_name %in% sample(unique(consensus_tissue_name), 5)) %>%
    rename(content = content_col,
           expression = expression_col) %>%
        {ggplot(., aes(content, expression, fill = content, color = content))+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "median",
                         fun.args = c("na.rm" = T),
                         geom = "line",
                         size = 2,
                         alpha = 0.5)+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "min",
                         geom = "line",
                         inherit.aes = F,
                         size = 1,
                         alpha = 0.5)+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "max",
                         geom = "line",
                         inherit.aes = F,
                         size = 1,
                         alpha = 0.5)+
            # stat_summary(fun.y = "min", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            # stat_summary(fun.y = "max", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+


            geom_violin(draw_quantiles = 0.5, alpha = 0.5, position = "identity")+
            geom_hline(yintercept = det_lim, linetype = "dashed")+
            simple_theme+
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
            theme(axis.title = element_blank())+
            scale_fill_manual(values = pal)+
            scale_color_manual(values = pal)}

}


