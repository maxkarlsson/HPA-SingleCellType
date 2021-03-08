#!/usr/bin/env Rscript

# required_packages <-
#   c("tidyverse", "pheatmap", "viridis", "magrittr")
# 
# install.packages(required_packages[!required_packages %in% rownames(installed.packages())], repos = "https://cloud.r-project.org")

require(tidyverse)
require(pheatmap)
require(magrittr)

# ------------ INPUTS ------------
# Gene info

# Format of geneinfo
# # A tibble: 19,670 x 2
# ensg_id         gene_name
# <chr>           <chr>    
# 1 ENSG00000000003 TSPAN6   
# 2 ENSG00000000005 TNMD     
# 3 ENSG00000000419 DPM1     
# 4 ENSG00000000457 SCYL3    
# 5 ENSG00000000460 C1orf112 
# 6 ENSG00000000938 FGR      
# 7 ENSG00000000971 CFH      
# 8 ENSG00000001036 FUCA2    
# 9 ENSG00000001084 GCLC     
# 10 ENSG00000001167 NFYA     
# # ... with 19,660 more rows

gene_info92 <- 
  read_delim("geninfo_92.tsv", delim = "\t") %>% 
  select(1:2)

# Cell markers

# Format of cell markers
# # A tibble: 237 x 4
# tissue   cell_group       cell_type            marker 
# <chr>    <chr>            <chr>                <chr>  
# 1 _General Adipose tissue   Adipocytes           LIPE   
# 2 _General Adipose tissue   Adipocytes           PLIN1  
# 3 _General Adipose tissue   Adipocytes           FABP4  
# 4 _General Adipose tissue   Adipocytes           SLC7A10
# 5 _General Adipose tissue   Adipocytes           PLIN4  
# 6 _General Epithelial cells neuroendocrine cells CHGA   
# 7 _General Epithelial cells neuroendocrine cells SCG2   
# 8 _General Epithelial cells neuroendocrine cells SCG5   
# 9 _General Epithelial cells neuroendocrine cells PCSK1N 
# 10 _General Epithelial cells neuroendocrine cells SCG3   
# # ... with 227 more rows

# uppsala_cell_markers <-
#   read_delim("data/meta/20200320 Cell type markers.tsv", delim = "\t") %>%
#   select(tissue = 1,
#          cell_group = 2,
#          cell_type = 3,
#          marker = 4)

cell_markers <-
  read_delim("Cell type markers.tsv", delim = "\t") %>%
  select(tissue = 1,
         cell_group = 2,
         cell_type = 3,
         marker = 4) %>% 
  left_join(gene_info92,
            by = c("marker" = "gene_name"))




# Cell type palette

cell_type_palette <- 
  read_delim("cell_type_palette.tsv", delim = "\t") %$% 
  set_names(value, name)

# Data

# Format of data
# # A tibble: 7,474,600 x 4
# dataset              ensg_id         cluster_id norm_count
# <chr>                <chr>           <chr>           <dbl>
# 1 Axillary lymph nodes ENSG00000121410 Cluster-0      5.51  
# 2 Axillary lymph nodes ENSG00000148584 Cluster-0      0     
# 3 Axillary lymph nodes ENSG00000175899 Cluster-0      8.77  
# 4 Axillary lymph nodes ENSG00000166535 Cluster-0      0     
# 5 Axillary lymph nodes ENSG00000184389 Cluster-0      0     
# 6 Axillary lymph nodes ENSG00000128274 Cluster-0     62.7   
# 7 Axillary lymph nodes ENSG00000118017 Cluster-0      0     
# 8 Axillary lymph nodes ENSG00000094914 Cluster-0     11.0   
# 9 Axillary lymph nodes ENSG00000081760 Cluster-0      7.95  
# 10 Axillary lymph nodes ENSG00000114771 Cluster-0      0.0459
# # ... with 7,474,590 more rows

cluster_norm_count <- 
  read_delim("cluster_exp.tab", delim = "\t") 
  
  
# cluster_annotation %>% 
#   select(dataset) %>% 
#   unique %>% 
#   mutate(B = c(NA, sort(unique(cluster_norm_count$dataset)), NA), 
#          C = paste0("'", dataset, "' = '", B, "'")) %>% 
#   pull(C) %>% 
#   paste0(collaspe = ",\n") %>% 
#   cat



# Cluster annotation

# Format of cluster annotation
# # A tibble: 380 x 5
# dataset_id dataset              cluster_id n_cells cell_type        
# <dbl> <chr>                <chr>        <dbl> <chr>            
# 1          1 Axillary lymph nodes Cluster-0     5335 Endothelial cells
# 2          1 Axillary lymph nodes Cluster-1     3843 Endothelial cells
# 3          1 Axillary lymph nodes Cluster-2     2971 Endothelial cells
# 4          1 Axillary lymph nodes Cluster-3     2787 Endothelial cells
# 5          1 Axillary lymph nodes Cluster-4     2357 Endothelial cells
# 6          1 Axillary lymph nodes Cluster-5     1879 Endothelial cells
# 7          1 Axillary lymph nodes Cluster-6     1639 Endothelial cells
# 8          1 Axillary lymph nodes Cluster-7     1309 Endothelial cells
# 9          1 Axillary lymph nodes Cluster-8      672 Endothelial cells
# 10          1 Axillary lymph nodes Cluster-9      489 Endothelial cells
# # ... with 370 more rows

cluster_annotation <- 
  read_delim("Cluster annotation.tsv", delim = "\t") %>% 
  select(dataset_id = 1,
         dataset = 2, 
         cluster_id = 3,
         cell_type = 4) %>%
  filter(!is.na(dataset_id)) %>%
  mutate(row_id = paste(dataset_id, str_extract(cluster_id, "\\d*$"), cell_type)) 



# ------------  Format data ------------

  
heatmap_palette <-
  viridis::inferno(20, direction = -1)

dataset_map <-
  c('Axillary lymph nodes' = 'Lymphatic system',
    'Breast' = 'Breast',
    'Colon' = 'Colon',
    'Colon 2' = 'Colon',
    'Eyes' = 'Eye',
    'Eyes macula' = 'Eye',
    'Eyes peripheral' = 'Eye',
    'Head and neck lymph nodes' = 'Lymphatic system',
    'Heart' = 'Muscle',
    'Ileum' = 'Ileum',
    'Kidney' = 'Kidney',
    'Liver' = 'Liver',
    'Liver hep- CD45-' = 'Liver',
    'Liver hep- CD45+' = 'Liver',
    'Lung' = 'Lung',
    'Muscle' = 'Muscle',
    'NK cells blood' = 'Lymphatic system',
    'NK cells bone marrow' = 'Lymphatic system',
    'PBMCs' = 'Lymphatic system',
    'Placenta' = 'Placenta',
    'Placenta blood' = 'Placenta',
    'Prostate' = 'Prostate',
    'Prostate 2' = 'Prostate',
    'Prostate 3' = 'Prostate',
    'Rectum' = 'Rectum',
    'Testis' = 'Testis',
    'Testis 2' = 'Testis',
    'Testis 3' = 'Testis') %>%
  enframe("dataset", "tissue")


heatmap_data <- 
  cluster_norm_count %>% 
  
  filter(ensg_id %in% cell_markers$ensg_id) %>% 
  left_join(cluster_annotation,
            by = c("dataset", "cluster_id")) %>% 
  left_join(gene_info92) 
# ------------  Format data ------------

# All markers and data 

row_annotation <- 
  cluster_annotation %>% 
  select(row_id, cell_type, dataset) %>%
  as.data.frame() %>%
  set_rownames(.$row_id) %>%
  {.[,-1]}

col_annotation <-
  cell_markers %>%
  select(cell_type, marker) %>%
  mutate(yes = "yes") %>%
  unique() %>%
  spread(cell_type, yes, fill = "no") %>%
  as.data.frame() %>%
  set_rownames(.$marker) %>%
  {.[,-1]}

colors_annotation <-
  col_annotation %>% 
  names() %>%
  {set_names(rep(list(c('yes' = 'red', 'no' = 'white')), length(.)), .)}

colors_annotation[["cell_type"]] <- cell_type_palette



plot_data <- 
  heatmap_data %>%
  select(row_id, gene_name, norm_count) %>%
  group_by(gene_name) %>% 
  mutate(max_count = max(norm_count),
         norm_count = norm_count / max_count) %>%
  ungroup() %>% 
  filter(max_count > 0) %>%
  select(-max_count) %>%
  spread(row_id, norm_count) %>% 
  column_to_rownames("gene_name") %>% 
  t() 


# write_delim(x = row_annotation, "row.tsv", delim = "\t")
# write_delim(x = col_annotation, "col.tsv", delim = "\t")
# write_delim(x = colors_annotation, "col.tsv", delim = "\t")

svg("Cluser marker heatmap.svg", width = 30, height = 40)
plot_data %>%
  
  pheatmap(clustering_method = "ward.D2", 
           color = heatmap_palette,
           border_color = NA,
           cutree_rows = 1,
           cutree_cols = 1,
           annotation_row = row_annotation,
           annotation_col = col_annotation,
           annotation_colors = colors_annotation,
           # filename = "Cluster marker heatmap.pdf",
           # width = 30, height = 40,
            
                                
           annotation_legend = F)

dev.off()

unique_tissues <-
  dataset_map %>%
  pull(tissue) %>%
  unique()


lapply(unique_tissues,
       function(tis_) {
         cat(tis_)
         tis_cell_markers <-
           cell_markers %>%
           filter(tissue %in% c("_General", tis_)) %>%
           left_join(gene_info92 %>%
                       select(gene_name, ensg_id))

         tis_col_annotation <-
           col_annotation[sort(unique(tis_cell_markers$marker)),] %>%
           {.[, sapply(., n_distinct) > 1]}

         heatmap_data %>%

           filter(grepl(tis_, tissue)) %>%
           filter(ensg_id %in% tis_cell_markers$ensg_id) %>%

           group_by(gene_name) %>%
           mutate(norm_count = norm_count / max(norm_count)) %>%
           filter(!is.na(norm_count)) %>%
           select(gene_name, row_id, norm_count) %>%
           spread(row_id, norm_count) %>%

           column_to_rownames("gene_name") %>%

           t() %>%
           pheatmap(clustering_method = "ward.D2",
                    color = heatmap_palette,
                    border_color = NA,
                    annotation_row = row_annotation,
                    annotation_col = tis_col_annotation,
                    annotation_colors = colors_annotation,
                    annotation_legend = F,
                    main = tis_,
                    filename = paste(tis_, "marker heatmap.pdf"),
                    width = 16, height = 12)

       })
