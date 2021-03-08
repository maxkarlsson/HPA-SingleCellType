
# expression_col = "nx"
# tissue_col = "consensus_tissue_name" 
# gene_col = "enssscg_id"  
# enr_fold = 4
# max_group_n = 5
# det_lim = 1
# data <- pig_atlas_consensus

hpa_gene_classification <- 
  function(data, expression_col, tissue_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      select(gene = gene_col,
             expression = expression_col,
             tissue = tissue_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    
    n_groups <- length(unique(data_$tissue))
  
    gene_class_info <- 
      data_ %>%
      group_by(gene) %>%
      summarise(
        
        # Gene expression distribution metrics
        mean_exp = mean(expression, na.rm = T),
        min_exp = min(expression, na.rm = T),
        max_exp = max(expression, na.rm = T), 
        max_2nd = sort(expression)[length(expression)-1],
        
        # Expression frequency metrics
        n_exp = length(which(expression >= det_lim)),
        frac_exp = n_exp/length(expression[!is.na(expression)])*100,
        
        # Limit of enhancement metrics
        lim = max_exp/enr_fold, 
        
        exps_over_lim = list(expression[which(expression >= lim & expression >= det_lim)]),
        n_over = length(exps_over_lim[[1]]), 
        mean_over = mean(exps_over_lim[[1]]),
        min_over = ifelse(n_over == 0, NA,
                          min(exps_over_lim[[1]])),
        
        max_under_lim = max(expression[which(expression < min_over)], det_lim*0.1),
        
        
        exps_enhanced = list(which(expression/mean_exp >= enr_fold & expression >= det_lim)),
        
        
        
        
        # Expression patterns
        enrichment_group = paste(sort(tissue[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(tissue[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(tissue[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        tissues_not_detected = paste(sort(tissue[which(expression < det_lim)]), collapse=";"),
        tissues_detected = paste(sort(tissue[which(expression >= det_lim)]), collapse=";")) 
      
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are tissue enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "tissue enriched", 
                                  
                                  # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in tissues fold times more than the mean are tissue enhance
                                  n_enhanced > 0 ~ "tissue enhanced", 
                                  
                                  # Genes expressed with low tissue specificity
                                  T ~ "low tissue specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 31 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "tissue enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "tissue enhanced" ~ max_exp/mean_exp)) 
    
      
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "tissue enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_enriched,
                                    spec_category == "tissue enhanced" ~ n_enhanced)) %>%
      select(gene, 
             spec_category, 
             dist_category, 
             spec_score,
             n_expressed = n_exp, 
             fraction_expressed = frac_exp,
             max_exp = max_exp,
             enriched_tissues,
             n_enriched,
             n_na = n_na,
             tissues_not_detected,
             tissues_detected) 
      
    
    
  }	

hpa_gene_classification_multi_sample <- 
  function(data, expression_col, tissue_col, gene_col, sample_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      select(gene = gene_col,
             expression = expression_col,
             tissue = tissue_col, 
             sample = sample_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    if(any(is.na(data_$sample))) stop("NAs in sample column")
    
    n_tissues <- length(unique(data_$tissue))
    n_samples <- length(unique(data_$sample))
    
    gene_wise_info <-
      data_ %>%
      group_by(gene) %>%
      mutate(max_exp = sample == sample[order(expression, decreasing = T)][1],
             max_2nd = sample == sample[which(tissue != tissue[which(max_exp)][1])][order(expression[which(tissue != tissue[which(max_exp)][1])], decreasing = T)][1],
             
             
             lim = max(max(expression) / enr_fold, det_lim),
             
             expressed = expression >= det_lim,
             enriched = expression >= lim & expressed,
             enhanced = expression/mean(expression) >= enr_fold & expressed,
             
             min_over_lim = sample == sample[which(enriched)][order(expression[which(enriched)], decreasing = F)][1],
             
             max_under_lim = sample == sample[which(!enriched)][order(expression[which(!enriched)], decreasing = T)][1]) %>% 
      arrange(gene, -expression) %>%
      mutate(min_over_lim = ifelse(is.na(min_over_lim), F, T),
             max_under_lim = ifelse(is.na(max_under_lim), F, T))
    
    gene_categories <- 
      gene_wise_info %>%
      summarise(n_sample_exp = length(which(expressed)),
                n_tissue_exp = n_distinct(tissue[which(expressed)]),
                
                frac_sample_exp = n_sample_exp / n_samples,
                frac_tissue_exp = n_tissue_exp / n_tissues,
                
                n_sample_enr = length(which(enriched)),
                n_tissue_enr = n_distinct(tissue[which(enriched)]),
                
                n_sample_enh = length(which(enhanced)),
                n_tissue_enh = n_distinct(tissue[which(enhanced)]),
                
                mean_enriched = mean(expression[which(enriched)]),
                mean_expression = mean(expression),
                
                enriched_samples = paste(sort(sample[which(enriched)]), collapse = ";"),
                enriched_tissues = paste(unique(sort(tissue[which(enriched)])), collapse = ";"),
                
                enhanced_samples = paste(sort(sample[which(enhanced)]), collapse = ";"),
                enhanced_tissues = paste(unique(sort(tissue[which(enhanced)])), collapse = ";"),
                
                tissue_enriched_score = expression[max_exp] / max(expression[max_2nd], det_lim*0.1),
                group_enriched_score = ifelse(any(max_under_lim),
                                              mean_enriched / expression[max_under_lim],
                                              0), 
                tissue_enhanced_score = expression[max_exp] / mean_expression, 
                
                spec_category = case_when(n_sample_exp == 0 ~ "not detected", 
                                          
                                          # Genes with expression fold times more than anything else are tissue enriched
                                          n_tissue_enr == 1 ~ "tissue enriched", 
                                          
                                          # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                          n_tissue_enr > 1 & 
                                            n_tissue_enr <= max_group_n ~ "group enriched", 
                                          
                                          # Genes with expression in tissues fold times more than the mean are tissue enhance
                                          n_tissue_enh > 0 ~ "tissue enhanced", 
                                          
                                          # Genes expressed with low tissue specificity
                                          T ~ "low tissue specificity"), 
                
                
                dist_category = case_when(frac_tissue_exp == 100 ~ "detected in all",
                                          frac_tissue_exp >= 31 ~ "detected in many",
                                          n_tissue_exp > 1 ~ "detected in some",
                                          n_tissue_exp == 1 ~ "detected in single",
                                          n_tissue_exp == 0 ~ "not detected")) %>%
      mutate(spec_score = case_when(spec_category == "tissue enriched" ~ 
                                      tissue_enriched_score,
                                    spec_category == "group enriched" ~ 
                                      group_enriched_score, 
                                    spec_category == "tissue enhanced" ~ 
                                      tissue_enhanced_score, 
                                    T ~ 0)) 
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enriched_tissues,
                                          spec_category == "tissue enhanced" ~ enhanced_tissues),
             enriched_samples = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enriched_samples,
                                          spec_category == "tissue enhanced" ~ enhanced_samples),
             n_tissues_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_tissue_enr,
                                            spec_category == "tissue enhanced" ~ n_tissue_enr),
             n_samples_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_sample_enr,
                                            spec_category == "tissue enhanced" ~ n_sample_enr)) %>%
      select(gene, 
             spec_category, 
             dist_category, 
             spec_score,
             n_samples_expressed = n_sample_exp, 
             n_tissues_expressed = n_tissue_exp, 
             fraction_samples_expressed = frac_sample_exp,
             fraction_tissues_expressed = frac_tissue_exp,
             
             enriched_samples,
             enriched_tissues,
             n_samples_enriched = n_sample_enr,
             n_tissues_enriched = n_tissue_enr) 
    
    
    
  }	

calc_gene_correlations <- 
  function(data, var1, var2, val1, val2, cor_method = "spearman", p_adjust_method = "BH", alternative = "two.sided") {
    
    data_ <- 
      data %>%
      rename(var1 = var1, 
             var2 = var2, 
             val1 = val1, 
             val2 = val2)
      
    
    data_ %>%
      group_by(var1, var2) %>% 
      do(if(cor_method == "pearson") {
        cor.test(.$val1, .$val2, method = cor_method, alternative = alternative) %$%
          tibble(pval = p.value, 
                 cor = estimate, 
                 lo_confint = conf.int[1], 
                 hi_confint = conf.int[2])
      } else if(cor_method == "spearman") {
        cor.test(.$val1, .$val2, method = cor_method, alternative = alternative) %$%
          tibble(pval = p.value, 
                 cor = estimate)
        }) %>% 
      ungroup() %>%
      mutate(padj = p.adjust(pval, method = p_adjust_method), 
             significant = padj <= 0.05, 
             log10P = -log10(padj)) %>% 
      arrange(padj) %>% 
      set_colnames(c(var1, var2, colnames(.)[-c(1, 2)]))
  }

calc_gene_distance <- 
  function(data, var1, var2, val1, val2) {
    
    data_ <- 
      data %>%
      rename(var1 = var1, 
             var2 = var2, 
             val1 = val1, 
             val2 = val2)
    
    
    data_ %>%
      group_by(var1, var2) %>% 
      summarise(dist = rbind(val1, val2) %>% 
                  dist() %>% 
                  as.numeric(), 
                mean_var1 = mean(val1), 
                mean_var2 = mean(val2), 
                common_mean = mean(c(val1, val2))) %>% 
      ungroup() %>%
      arrange(dist) %>% 
      set_colnames(c(var1, var2, "dist", paste0("mean_", c(var1, var2)), "common_mean"))
  }



make_ortholog_net <- 
  function(orthologs) {
    edges <- 
      orthologs %>%
      select(from = 1, 
             to = 2) %>% 
      mutate(weight = 1)
    
    nodes <- 
      orthologs %$%
      c(enssscg_id,
        ensg_id) %>% 
      unique()
    
    
    orth_net <- 
      graph_from_data_frame(d = edges, 
                            vertices = nodes, 
                            directed = F)
    
    ortholog_communities <- 
      components(orth_net) %$% 
      left_join(enframe(membership, 
                        "node", "community"), 
                enframe(csize, 
                        "community", "community_size"), 
                by = "community")
    
    list(edges = edges, 
         nodes = nodes, 
         network = orth_net,
         communities = ortholog_communities)
  }








