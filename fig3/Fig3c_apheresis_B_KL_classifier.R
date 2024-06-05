
library(Seurat)
library(dplyr)
library(ggplot2)

calculate_malignant_b_fraction <- function(obj = NULL, 
                                           min_chain_exp = 0, 
                                           n_informative_cells = 10, 
                                           percent_informative = 0.85,
                                           kl_chain_ratio = 4) {
  
  data <- FetchData(obj, vars = c("IGKC","IGLC2"))
  
  ## set minimum cell chain expression cutoff
  if(min_chain_exp > 0){
    data <- apply(data, 2, function(x) ifelse(x > min_chain_exp, x, 0))
  }
  
  ## checking whether a cell expresses no chain, a single chain or both chains
  bl <- table(as.integer(apply(data, 1, function(x) length(which(x>0))))) %>% as.list()
  if (is.null(bl$`1`)) {
    return(list())  
  }
  if (is.null(bl$`2`)) {
    bl$`2` = 0
  } 
  if (is.null(bl$`0`)) {
    bl$`0` = 0
  }
  
  total_b <- dim(data)[1]
  informative_b <- bl$`1`  
  
  if(informative_b > n_informative_cells) {
    
    IGKC_restricted <- which(apply(data, 1, function(x) {length(which(x["IGKC"]>0))==1 & length(which(x["IGLC2"]>0))==0}))
    informative_fraction_IGKC <- length(IGKC_restricted)/informative_b
    is_IGKC_restricted <- ifelse(informative_fraction_IGKC > percent_informative, TRUE, FALSE)
    IGLC2_restricted <- which(apply(data, 1, function(x) {length(which(x["IGKC"]>0))==0 & length(which(x["IGLC2"]>0))==1}))
    informative_fraction_IGLC2 <- length(IGLC2_restricted)/informative_b
    is_IGLC2_restricted <- ifelse(informative_fraction_IGLC2  > percent_informative, TRUE, FALSE)
    is_malignant <- ifelse(is_IGKC_restricted | is_IGLC2_restricted, TRUE, FALSE)
    
    if(informative_b/total_b < 0.5){
      is_malignant = "UNDEF"
    }
    
    return(list(is_malignant = is_malignant, 
                informative_cells = informative_b,
                total_cells = total_b,
                iIGKC = informative_fraction_IGKC,
                iIGLC2 = informative_fraction_IGLC2))
  } else {
    return(list(is_malignant = "UNDEF",  informative_cells= NA, total_cells=NA, iIGKC = NA, iGLC2 = NA))
  }
}

###############################################################################################


seurat_b_ds <- readRDS("./seurat.B.ds300.clean.rds"))

ncells_counts <- 
  seurat_b_ds@meta.data %>% 
  group_by(sample_id) %>%
  dplyr::count(name = "ncells") %>% 
  arrange(ncells) %>% 
  as.data.frame() 
samples.keep <- ncells_counts %>% dplyr::filter(ncells >= 10) %>% pull(sample_id)
seurat <- subset(seurat_b_ds, sample_id %in% samples.keep)


#FeatureScatter(seurat, feature1 = "IGKC", feature2 = "IGLC2", pt.size = 0.5)
min_chain_exp = 2.5 ## cutoff for chain expression 

b_cell_annot <- list() 
samples <- unique(seurat$sample_name)
for (s in samples) {
  print(s)
  subset <- subset(seurat, sample_name == s)
  b_cell_annot[[s]] <- tryCatch({
    calculate_malignant_b_fraction(subset, min_chain_exp = min_chain_exp) %>% as.data.frame()
  },  error = function(cond) {
    message(conditionMessage(cond))
    TRUE
  })
}
b_cell_annot <- plyr::ldply(b_cell_annot) %>% dplyr::rename(sample_name = `.id`)

## Figure 3c - KL
md <- seurat@meta.data %>% 
      tibble::rownames_to_column(var = "cell_id") %>% 
      select(sample_name,  disease_status) %>% 
      mutate(patient_alias = sample_name) %>% 
      distinct()

ggdata <- b_cell_annot %>% 
  select(sample_name, iIGKC, iIGLC2, is_malignant) %>% 
  dplyr::rename(IGKC=iIGKC, IGLC2=iIGLC2) %>% 
  tidyr::pivot_longer(-c(sample_name, is_malignant), names_to = "Chain", values_to = "fraction" ) %>% 
  left_join(md)
ggdata$is_malignant <- factor(ggdata$is_malignant, levels = c("UNDEF", FALSE, TRUE))
ggdata$disease_status <- factor(ggdata$disease_status, levels = c("Healthy", "Malignant"))

ggplot(ggdata, aes(x = sample_name, y = Chain, color = is_malignant, size = fraction)) +
  geom_point(aes(colour = is_malignant)) +
  scale_colour_manual(values = list("TRUE" = "red" , "FALSE" = "black", "UNDEF" = "gray")) +
  xlab("Sample id") +
  ylab("Chain") +
  theme_bw() +
  ggpubr::rotate_x_text(45)



