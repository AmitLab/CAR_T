  
  
  library(dplyr)
  #library(patchwork)
  library(ggplot2)
  library(Seurat)
  
  library(devtools)
  library(infercnv)
  library(Matrix)
  
  setwd("/home/labs/amit/weitan/projects/yonatan/")
  project_name <- "Apheresis"
  seqtype = "10x"
  species_id = "Hs"
  

sample_order <- read.csv("knn_classifier_sample_order.csv")[[1]]  ## a column with names of samples in the order they should be plotted

reference_library = "SH_HC"

seurat_b <- readRDS("./data/seurat.B.ds300.rds")  

min_cells_per_patient = 30

annot_fn <- "CNV_annotations_file.txt"

  ncells_counts <- 
    seurat_b@meta.data %>% 
    group_by(sample_name) %>%
    dplyr::count(name = "ncells") %>% 
    arrange(ncells) %>% 
    as.data.frame() 
  
  samples.keep <- ncells_counts %>% dplyr::filter(ncells >= min_cells_per_patient) %>% pull(sample_name)
  seurat_b <- subset(seurat_b, sample_name %in% samples.keep)

  low_n_samples <- setdiff(sample_order, unique(seurat_b@meta.data$sample_name))
  sample_order <- setdiff(sample_order,low_n_samples)
  sample_order <- rev(sample_order)
  
  sample_rename_cnv <- function(sname, index){
    prefix <- stringr::str_split(sname, pattern = "-")[[1]][1]
    suffix <- stringr::str_split(sname, pattern = "-")[[1]][2]
    cnv_name <- paste(letters[index],prefix,suffix, sep = "-")
  }
  
  sample_names_cnv = c()
  for (i in seq(1,length(sample_order))){
    sample_names_cnv <- c(sample_names_cnv, sample_rename_cnv(sample_order[length(sample_order)-i+1],i))
  }
  sample_name_df <- data.frame("sample_name" = sample_order,"sample_name_cnv" = rev(sample_names_cnv))
  seurat_b@meta.data <- seurat_b@meta.data %>% 
          tibble::rownames_to_column(var = "cell_id") %>% 
          left_join(sample_name_df) %>% 
          tibble::column_to_rownames(var = "cell_id")
          
  annot_data <- seurat_b@meta.data %>% 
                tibble::rownames_to_column(var = "cell_id") %>% 
                select(cell_id, sample_name_cnv, response_binary) %>% distinct()  %>% as.data.frame() %>% 
                rowwise() %>% 
                mutate(group = case_when(
                        response_binary == "HC" ~ "HC_B",
                        response_binary %in% c("NR", "R") ~ paste("malignant", sample_name_cnv, sep = "_"),
                        TRUE ~ "unknown")) %>% 
                filter(group != "unknown") %>% 
                select(cell_id, group)
  
  rownames(annot_data) <- NULL


  write.table(annot_data, 
            file = annot_fn, 
            row.names =  FALSE, sep = "\t")    
                  
  options(scipen = 100)
  
  counts <- GetAssayData(seurat_b, layer = "counts")
  print(dim(counts))
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                      annotations_file=annot_fn,
                                      delim="\t",
                                      gene_order_file="~/data/CNV/hg38_gencode_v27.txt",
                                      ref_group_names=c("HC_B")
                                      ) 
  
  
   analysis_mode_var = "subclusters" 
   cluster_by_groups_var = TRUE 
   
   infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=0.1 works well for 10x Genomics , cutoff=1 works well for Smart-seq2
                               analysis_mode = analysis_mode_var,
                               tumor_subcluster_partition_method = "qnorm", # default is leiden
                               cluster_by_groups=cluster_by_groups_var, 
                               denoise=TRUE,
                               HMM=TRUE,
                               k_nn=25) 
  
   
  
