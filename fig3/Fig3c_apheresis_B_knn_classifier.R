
library(dplyr)
library(patchwork)
library(ggplot2)
library(metacell)


## metacell based transcriptional classifier
##############################################################################

compute_cell_proximity <- function(id, method = c('spearman'), geneset_id = NULL, numi = 1500){
  
  marker_gset <- scdb_gset(geneset_id) 
  geneset <- names(marker_gset@gene_set)

  ##normalization by downsampling nUMI
  mat <- scdb_mat(id)
  good_genes <- mat@genes
  clean_mat <- scm_sub_mat(scmat = mat, genes = good_genes)
  numi = scm_which_downsamp_n(mat)
  print(paste0("cell cell correlations were computed after downsampling to ", numi, " UMIs" ))
  dsmat <- .downsamp(as.matrix(clean_mat@mat), n = numi)
  expMat <- log2(t(as.matrix(dsmat))+1)
  ds <- expMat[,geneset]
  # dim(ds)
  # 18611  3573
  
  if(method == 'spearman'){
    cell_cell_cor <- tgstat::tgs_cor(x = t(ds), spearman = TRUE)
    return(cell_cell_cor)
  } else {
    print("No proximinity measure declared")
    return()
  }
  
}


## identifying closest k neighbors per cell
compute_nbh <- function(id, r, samples, ref_samples = NULL, k, min_cells_sample = 10){
  
  mat <- scdb_mat(id)
  sample_metadata <- mat@cell_metadata %>% 
    tibble::rownames_to_column(var = "cell_id") %>% 
    dplyr::filter(sample_name %in% samples) %>% 
    dplyr::select(cell_id, sample_name, disease_status) 
  slist <- split(sample_metadata$cell_id, f = sample_metadata$sample_name)
  
  if(!is.null(ref_samples)){
    ref_sample_metadata <- mat@cell_metadata %>% 
      tibble::rownames_to_column(var = "cell_id") %>% 
      dplyr::filter(sample_name %in% ref_samples) %>% 
      dplyr::select(cell_id, sample_name, disease_status) 
    
    ref_slist <- split(ref_sample_metadata$cell_id, f = ref_sample_metadata$sample_name)
  }
  
  sample_cell_distances <- list()
  for (s in samples) {
    
    s_cells <- slist[[s]]
    s_cells_used <- s_cells[s_cells %in% rownames(r)]
    
    if(is.null(ref_samples)){
      other_cells <- sample_metadata %>% dplyr::filter(!(cell_id %in% s_cells)) %>% pull(cell_id)  
    } else{
      other_cells <- ref_sample_metadata %>% pull(cell_id)  
    }
    
    other_cells_used <- other_cells[other_cells %in% rownames(r)]
    r_sample <- r[s_cells_used, other_cells_used] 
    
    ## get top k neighbours per cell
    ncells <- length(s_cells_used)
    if (ncells > min_cells_sample){
      k_mat <- matrix(nrow = ncells, ncol = k)
      for (i in seq(1:ncells)){
        cell <- s_cells_used[i]
        top_cells <- sort(r_sample[cell,], decreasing =  TRUE) %>% head(k) %>% as.numeric()
        k_mat[i,] <- top_cells
      }
      rownames(k_mat) <- s_cells_used
      sample_cell_distances[[s]] <- k_mat
    }
  }
  
  return(sample_cell_distances)
  
}


test_nbh_significance <- function(ref_nbh, group_nbh) {
  
  ref_mean <- mean(unlist(ref_nbh))
  ref_sd <- sd(unlist(ref_nbh))
  #norm_hist_p <- 
  
  res <- list()
  for(sample in names(group_nbh)){
    sample_mat <- group_nbh[[sample]]
    sample_pvals <- c()
    for (cell in rownames(sample_mat)){
      sample_mean <- mean(sample_mat[cell,])
      sample_sd <-   sd(sample_mat[cell,])
      q = (sample_mean - ref_mean) / sample_sd # z-score 
      pval = pnorm(q, mean = 0, sd = 1, lower.tail = TRUE) # p value from z-score  
      sample_pvals <- c(sample_pvals, pval)
    }
    res[[sample]] <- setNames(object = sample_pvals, nm = rownames(sample_mat))
  }
  
  return(res)
}

## 
summarize_nbh_significance <- function(test_res){
  
  res_df <- lapply(test_res, function(x) {as.data.frame(x) %>% dplyr::rename(pvalue = x) %>% mutate(sample = names(x))})   
  
  data <- data.frame()
  for (s in names(res_df)) {
    s_df <- as.data.frame(x = res_df[[s]]) %>% mutate(sample = s) %>% tibble::rownames_to_column(var = "cell_id"); colnames(s_df) <- c("cell_id", "pvalue", "sample")
    s_df[["FDR"]] <- p.adjust(s_df$pvalue, method = "BH", n = length(s_df$pvalue)) 
    data <- rbind(data, s_df)
  }
  
  data <- data %>% mutate(logFDR = log10(FDR))
  return(data)
}



classify_cell_state <- function(id, method = "spearman", k = 100, min_cells_sample = 10, healthy_mc){
  
  mat <- scdb_mat(id)
  mc <- scdb_mc(id)
  mc_assignments <- mc@mc[rownames(mat@cell_metadata)]
  mat@cell_metadata$mc <- mc_assignments
  mat@cell_metadata$cell <- rownames(mat@cell_metadata)
  
  ## compute cell-cell correlations
  geneset_nm = paste0(id, "_lateral_filtered") # metacell variable gene set
  cell_cell_cor <- compute_cell_proximity(id, method = 'spearman', geneset_id = geneset_nm)
  
  ## distribution of cell correlations for healthy population
  df <- mat@cell_metadata[rownames(cell_cell_cor), c("disease_status","sample_name", "cell", "mc")] %>% 
        dplyr::filter(disease_status == "Healthy", mc %in% healthy_mc)
  healthy_samples <- df$sample_name %>% unique()
  healthy_cells <- df$cell
  healthy_cor = cell_cell_cor[healthy_cells, healthy_cells]
  healthy_samples_nbh <- compute_nbh(id, r =  healthy_cor, samples = healthy_samples, k = k)
  
  ## distribution of cell correlations for patient vs healthy populations
  patient_samples <- mat@cell_metadata %>% dplyr::filter(!(sample_name %in% healthy_samples)) %>% group_by(sample_name) %>% dplyr::count() %>% filter(n >= min_cells_sample) %>% pull(sample_name) %>% unique()
  patient_samples_nbh <- compute_nbh(id, r = cell_cell_cor, samples = patient_samples, ref_samples = healthy_samples, k = k)
  
  ## testing cell neighborhood distribution using z-score statistics  
  patient_nbh_sig_res <- test_nbh_significance(ref_nbh = healthy_samples_nbh, group_nbh = patient_samples_nbh)
  patient_nbh_sig_summary <- summarize_nbh_significance(patient_nbh_sig_res)
  healthy_nbh_sig_res <- test_nbh_significance(ref_nbh = healthy_samples_nbh, group_nbh = healthy_samples_nbh)
  healthy_nbh_sig_summary <- summarize_nbh_significance(healthy_nbh_sig_res)
  
  ## summarizing naive comparisons 
  md <- mat@cell_metadata %>% dplyr::select(sample_name,disease_status) %>% distinct() %>% tibble::remove_rownames()
  
  nbh_sig_summary <- rbind(healthy_nbh_sig_summary, patient_nbh_sig_summary) %>% 
    dplyr::rename(sample_name=sample) %>% 
    left_join(md) %>% 
    as.data.frame()
  
  nbh_sig_order_df <- nbh_sig_summary %>% 
    group_by(sample_name, disease_status) %>% 
    dplyr::summarise(median_logFDR = median(-logFDR)) %>% 
    arrange(disease_status,median_logFDR) %>% 
    as.data.frame() 
  
  nbh_sig_order <- nbh_sig_order_df %>% pull(sample_name)
  nbh_sig_summary$sample_name <- factor(nbh_sig_summary$sample_name, levels = nbh_sig_order)
  nbh_sig_summary <- nbh_sig_summary %>% arrange(sample_name)
  return(nbh_sig_summary)
}


## Main 
##################################################################################

outdir_name <- id <-  "apheresis_B_base_clean"  
work_dir <- paste0("./saved_work/", outdir_name)
figdir <- paste0("./figures/", outdir_name)
anndir <- paste0("./annotations/", outdir_name)
dir.create(work_dir, recursive = TRUE, showWarnings = TRUE)
scdb_init(work_dir, force_reinit=T)
dir.create(figdir, showWarnings = TRUE, recursive = TRUE)
scfigs_init(figdir)
dir.create(anndir, showWarnings = TRUE, recursive = TRUE)
disease_status_name2color <- c("Healthy" = "darkblue", "Malignant" = "darkred")

reference_mc = c(seq(1,4), seq(8,11)) # healthy reference cells are sampled from naive and memory metacells
nbh_sig_summary = classify_cell_state (id, method = "spearman", k = 100, min_cells_sample = 10,  healthy_mc = reference_mc) 

outfile <- paste0(anndir, "/knn_classifier.csv")
write.csv(nbh_sig_summary, file = outfile, row.names = FALSE)  

#outfile <- paste0(anndir, "/knn_classifier_sample_order.csv")
#write.csv(nbh_sig_order, file = outfile, row.names = FALSE)
  