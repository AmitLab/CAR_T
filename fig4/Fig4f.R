

library(dplyr)
library(patchwork)
library(ggplot2)


create_prediction_heatmap <- function(sample_feature_file, metadata_file, prediction_file, features = NULL, outfile){
  
  df <- read.csv(sample_feature_file)
  if (!("sample_id" %in% colnames(df))){
    df <- df %>% dplyr::rename(sample_id = X)
  }
  
  if (!is.null(features)){
    df1 <- df %>% dplyr::select(all_of(c("sample_id", features))) %>% tibble::column_to_rownames("sample_id") %>% as.data.frame()
  } else {
    df1 <- df %>% select(-c("B_category","response_3m","LDH_prior_tx")) %>% as.data.frame()  # Product
  }
  
  obj_md <- read.csv(metadata_file)
  
  rename_samples <- function (sid) {
    dict = list("NOV"="Tis","GIL"="Axi")
    return(unlist(lapply(sid, function(x) {
      p = strsplit(x, split = "_")[[1]][1]
      p_id = strsplit(x,split = "_")[[1]][2]
      paste(dict[p],p_id,sep="-")
    })))
  }
  
  obj_md$sample_id <- rename_samples(obj_md$sample_id)
  
  prediction_probs <- read.csv(prediction_file) %>% rename(sample_id = X)
  
  annotation_row_df = obj_md %>% select(sample_id, response_3m) %>% 
    dplyr::rename(response = response_3m) %>% 
    left_join(prediction_probs) %>% 
    dplyr::select(sample_id, response, prediction_scores) %>% 
    tibble::column_to_rownames(var = "sample_id") 
  annotation_row_df <- annotation_row_df[!is.na(annotation_row_df$prediction_scores),]
  
  df1 <- scale(df1, scale = TRUE) 
  df1[df1 > 3] <- 3
  df1[df1 < -3] <- -3
  response_status_name2color <- setNames(c("#b30000", "#71cceb"), c("NR","R"))
  breaksList = seq(-3, 3, by = 0.1)
  cols = colorRampPalette(c("blue", "white", "red"))(length(breaksList))
  prob_breaksList = seq(0, 1, by = 0.05)
  prob_cols = colorRampPalette(c("#b30000", "white", "#71cceb"))(length(prob_breaksList))
  
  pheatmap::pheatmap(df1[,features], scale = "none", cluster_rows = TRUE, annotation_row = annotation_row_df, 
                     annotation_colors = list(response = response_status_name2color, 
                                              predicted = response_status_name2color,
                                              prediction_scores = prob_cols
                     ), breaks = breaksList, color = cols, width = 6, height = 10, treeheight_row = 25, treeheight_col = 25) 
  
}




sample_abundance_file <- "./Fig4_data/prediction/partial_features.csv"
metadata_file <- "./Fig4_data/sample_response_duration.csv"
prediction_file = "./Fig4_data/prediction/partial_features_prediction.csv"
features_use <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB.CD16.Mono","B_category","B","CD4.T","CD8.T","Myeloid","NKT","Tregs","CD14.Mono","CD16.Mono","mon_ratio")
create_prediction_heatmap(sample_abundance_file, metadata_file, prediction_file, features = features_use)

