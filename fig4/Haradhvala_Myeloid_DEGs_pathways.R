########################################################################
#  DEG and pathway analysis Myeloid cells - Haradhvala et al. Nature Medicine
########################################################################
# load packages
library(Seurat)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(stringr)

# AggregateExpression for pseudobulk and Deseq2
########################################################################
# load Seurat object
Seurat_Myeloid <- readRDS("/Haradhvala_pbmc/seurat_myeloid.annot.rds")

# Use aggregate Expression to pseudobulk counts based on response-patient-celltype
Seurat_Myeloid_pseudo <- AggregateExpression(Seurat_Myeloid, assays = "RNA", return.seurat = T, group.by = c("response_binary", "sample_id", "cell_type"))

#add new metadata column with cell_type and orig.ident
metadata <- Seurat_Myeloid_pseudo@meta.data
metadata  <- metadata  %>%
  mutate(cell_type = gsub(".*_", "", rownames(.)),  # Extract cell type from row names
         cell_type_response = paste(cell_type, orig.ident, sep = "_"))  # Combine cell type and orig.ident
Seurat_Myeloid_pseudo@meta.data <- metadata

# set new idents (response_patient_cell_type)
Idents(Seurat_Myeloid_pseudo) <- Seurat_Myeloid_pseudo@meta.data$cell_type_response

# Pseudobulk with DESeq2 and save excel sheets
output <- "/Output"
celltypes <- unique(Seurat_Myeloid$cell_type)  

# Initialize an empty list to store results
results_list <- list()

# Perform analysis for each cell type
for (i in celltypes) { 
  bulk.Mcells <- FindMarkers(object = Seurat_Myeloid_pseudo, 
                             ident.1 = paste0(i, "_R"), 
                             ident.2 = paste0(i, "_NR"),
                             test.use = "DESeq2")
  
  # Add a new column indicating the comparison
  bulk.Mcells$comparison <- paste0(i, "_R_vs_NR")
  # Add a new column with gene names
  bulk.Mcells$gene_names <- rownames(bulk.Mcells)
  # Save the dataframe to CSV
  write.csv(bulk.Mcells, file = paste0(output, i, "DESeq2_pseudobulk_R_vs_NR.csv"))
  
  # Store the dataframe in the list
  results_list[[i]] <- bulk.Mcells
}

# Combine all dataframes into one big dataframe
all_results_DEseq2 <- do.call(rbind, results_list)

# Save the combined dataframe to CSV
write.csv(all_results_DEseq2, file = paste0(output, "Myeloid_cells_Haradhvala_pseudobulk_results_DEseq2.csv"))


# HEATMAPS OF TOP 30 genes, p val < 0.05
########################################################################
#get expression matrix
expr_matrix <- GetAssayData(object = Seurat_Myeloid_pseudo, layer = "scale.data")
expr_df <- as.data.frame(expr_matrix)

# Define the cell types
cell_types <- unique(Seurat_Myeloid$cell_type)

# Iterate over each cell type
for (cell_type in cell_types) {
  # Extract the matching columns for the current cell type
  matching_columns <- colnames(expr_df)[endsWith(colnames(expr_df), cell_type)]
  
  # Subset expr_df to include only the specified cell type
  subset_expr_df <- expr_df[, matching_columns]
  
  # Select a subset of genes
  Mono_sig_genes <- subset(all_results_DEseq2, comparison == paste(cell_type, "_R_vs_NR", sep = "") & p_val < 0.05 &
                             pct.1 >= 0.2 & 
                             pct.2 >= 0.2)$gene_names
  Mono_sig_genes <- head(Mono_sig_genes, 30)
  
  # Subset the expression dataframe for selected genes
  subset_expr_df <- subset_expr_df[Mono_sig_genes, ]
  
  # Define breaks for color scale
  breaks <- seq(-2, 2, length.out = 101)
  
  # adapt colnames
  colnames_old <- colnames(subset_expr_df)
  colnames_new <- sub("_[^_]*$", "", colnames_old)
  colnames(subset_expr_df) <- colnames_new
  
  # Create a heatmap using pheatmap
  pm <- pheatmap(subset_expr_df, 
                 scale = "row",  # Scale rows
                 cluster_rows = TRUE,  # Cluster rows
                 cluster_cols = FALSE,  # Do not cluster columns
                 color = colorRampPalette(c("royalblue4", "white", "red4"))(100),# Define color palette
                 breaks = breaks,  # Specify the breaks for the color scale
                 main = paste("Heatmap of Gene Expression for", cell_type))
  
  # Save the heatmap
  ggsave(paste0(cell_type, "_heatmap_top30.pdf"), pm, width = 8, height = 7)
}


# Geneset enrichment with FGSEA Multilevel
########################################################################
# read in file containing lists of genes for each pathway
# downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

setwd("..data/")
hallmark_pathway <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

setwd("/home/labs/amit/pascale/Projects/CART/Haradhvala_Nat_Med/Myeloid_pathways/")

# Get unique CellTypes R vs NR
cellTypes <- unique(all_results_DEseq2$comparison)

# Create an empty dataframe to store all results
all_results_df <- data.frame()

# Loop over each CellType for the current Comparison
for (i in cellTypes) {
  
  # Subset the data for the current CellType
  subset_df <- subset(all_results_DEseq2, comparison == i)
  subset_df2 <- subset(subset_df, p_val  < 0.05)
  
  # Order the data by avg_log2FC and  pvalue for fgseaMultilevel
  geneList <- subset_df2 %>%
    filter(!is.na(gene_names), !is.na(avg_log2FC)) %>%  #remove genes and log2FC which is zero
    mutate(ranking_metric = -log10(p_val)*sign(avg_log2FC)) %>%   #multiply p value with -1 if avg_log2FC is negative, 0 if 0 and +1 if positive 
    group_by(gene_names) %>% 
    summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
    arrange(-ranking_metric) %>% # sort descending (important!)
    tibble::deframe() # convert to named vector

  # Perform fgsea analysis
  fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                   stats = geneList,
                                   minSize = 5,
                                   maxSize = Inf,
                                   nPermSimple = 1000)
  
  # Add column for cell type to the results
  fgsea_results$CellType <- rep(i, nrow(fgsea_results))
  
  # Append the results to the overall dataframe
  all_results_df <- rbind(all_results_df, fgsea_results)
  
}


# Save the combined results dataframe to a CSV file
all_results_df$leadingEdge <- sapply(all_results_df$leadingEdge, paste, collapse = ",")
write.csv(all_results_df, "pathways_all_results_genes_Haradhvala_5genes_min_ranking.csv", row.names = FALSE)

# Create a dotplot with all pathways
# add the count of the leading edge as new column to plot later
all_results_df[, leadingEdge_count := sapply(strsplit(leadingEdge, ","), length)]

# Define significance threshold and add a new column indicating significance
significance_threshold <- 0.05
all_results_df[, significant := padj < significance_threshold]

#subset on pathways with leading edge of at least 5
all_results_df_leadingEdge_5 <- subset(all_results_df, leadingEdge_count>=5)


# Module score patient level analysis for selected pathways
########################################################################
#  calculate Z scores for selected pathways (significantly changed btw R and NR in our own dataset)
selected_pathways <- unique(all_results_df_leadingEdge_5$pathway)[c(8,7,5,21,19,20,17,2,1)]
all_results_df_leadingEdge_5_subset <- all_results_df_leadingEdge_5[pathway %in% selected_pathways, ]

# add celltype as column
all_results_df_leadingEdge_5_subset$cell_type <- gsub("_.*", "", all_results_df_leadingEdge_5_subset$CellType)

# select all cell types which have significant pathways
celltypes <- unique(all_results_df_leadingEdge_5_subset$cell_type)

all_data_long <- data.frame()

for (cellsubset in celltypes) {
  subset_leadingEdge <- subset(all_results_df_leadingEdge_5_subset, cell_type == cellsubset)
  Pseudo_subset <- subset(Seurat_Myeloid_pseudo, cell_type == cellsubset)
  
  pathways_with_leading_edge <- list()
  
  for (i in 1:nrow(subset_leadingEdge)) {
    pathway <- subset_leadingEdge$pathway[i]
    leading_edge_genes <- strsplit(subset_leadingEdge$leadingEdge[i], ",")[[1]]
    
    pathways_with_leading_edge[[pathway]] <- leading_edge_genes
  }
  
  for (pathway_name in names(pathways_with_leading_edge)) {
    genes <- pathways_with_leading_edge[[pathway_name]]
    
    Pseudo_subset <- AddModuleScore(object = Pseudo_subset, 
                                    features = list(genes = genes), 
                                    name = paste0("pathway_", pathway_name), 
                                    assay = "RNA",
                                    scale = TRUE)
  }
  
  data <- Pseudo_subset@meta.data
  # response <- rownames(data)
  data$response_3m <- data$orig.ident
  data$patient <- rownames(Pseudo_subset@meta.data)
  data$patient <- str_extract(data$patient, "(?<=_)[A-Za-z0-9-]+(?=_)")
  
  data_long <- tidyr::pivot_longer(data, cols = starts_with("pathway_"), 
                                   names_to = "pathway", values_to = "zscore")
  
  plot1 <- ggplot(data_long, aes(x = response_3m, y = zscore, fill = response_3m)) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha = 1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(color = response_3m, fill = response_3m), size = 0.8, shape = 21) +  
    scale_fill_manual(values = c("#9d1c20", "#81c4db")) +
    scale_color_manual(values = c("black", "black")) +  
    labs(title = paste0(cellsubset, ", Pathway Module Scores by Response"), x = "Response", y = "Zscore") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.position = "none",
          aspect.ratio = 1,
          panel.spacing = unit(0.5, "cm")) +
    facet_wrap(~ pathway, scales = "free_y", ncol = 3, as.table = TRUE) +
    ylim(-0.3, 0.5)
  
  print(plot1)
  ggsave(paste0("Zscore_plots_", cellsubset, ".pdf"), plot = plot1, width = 8, height = 4)
  
  pathway_tests <- data_long %>%
    group_by(pathway) %>%
    summarize(p_value = t.test(zscore ~ response_3m)$p.value)
  
  write.csv2(pathway_tests, file = paste0("T_test_pathway_zscore_", cellsubset, ".csv"))
  
  all_data_long <- rbind(all_data_long, data_long)
}

write.csv(all_data_long, "Haradhvala_selected_pathways_modulescore_long.csv", row.names = FALSE)

# rearange, clean and add new column
Zscore_myeloid_df <- transform(all_data_long,
                               cleaned_pathways = gsub("^pathway_|\\d+$", "", pathway),
                               celltype_pathway = paste(cell_type, gsub("^pathway_|\\d+$", "", pathway), sep = "_"),
                               Patient = sub("Patient([0-9])$", "Patient0\\1", sub("-.*$", "", patient))
)


# Add Product from Seurat_Myeloid@meta.data
matching_indices <- match(Zscore_myeloid_df$Patient, Seurat_Myeloid@meta.data$Patient)
Zscore_myeloid_df$Product <- Seurat_Myeloid@meta.data$Product[matching_indices]

# order by pathway and save
Zscore_myeloid_df_order <- Zscore_myeloid_df[order(Zscore_myeloid_df$pathway), ]
Zscore_myeloid_df_save <- Zscore_myeloid_df_order[,c(2,3,4,6,7,8,10,11)]

write.csv(Zscore_myeloid_df_save, "Myeloid_Modulescore_all_pathways5genes.csv", row.names = FALSE)

