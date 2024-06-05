########################################################################
#  DEG and pathway analysis NK/T cells
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
# load NKT cell Seurat object 
Seurat_NKT <- readRDS("..data/seurat.nkt.annot.rds")

# remove healthy control group for DEG analysis
Seurat_NKT_noHC <- subset(Seurat_NKT, subset = response_binary != "HC")

# Use aggregate Expression to pseudobulk counts based on response-patient-celltype
Seurat_NKT_pseudo <- AggregateExpression(Seurat_NKT_noHC, assays = "RNA", return.seurat = T, group.by = c("response_3m", "sample_id", "cell_type"))

#add new metadata column with cell_type and orig.ident
metadata <- Seurat_NKT_pseudo@meta.data
metadata  <- metadata  %>%
  mutate(cell_type = gsub(".*_", "", rownames(.)),  # Extract cell type from row names
         cell_type_response = paste(cell_type, orig.ident, sep = "_"))  # Combine cell type and orig.ident
Seurat_NKT_pseudo@meta.data <- metadata

# set new idents (response_patient_cell_type)
Idents(Seurat_NKT_pseudo) <- Seurat_NKT_pseudo@meta.data$cell_type_response

# Pseudobulk with DESeq2 and save excel sheets
output <- "../Output/"
 
# Initialize an empty list to store results
results_list <- list()

# Perform analysis for each cell type
celltypes <- unique(Seurat_NKT_noHC$cell_type) 
for (i in celltypes) { 
  bulk.Tcells <- FindMarkers(object = Seurat_NKT_pseudo, 
                             ident.1 = paste0(i, "_R"), 
                             ident.2 = paste0(i, "_NR"),
                             test.use = "DESeq2")
  
  # Add a new column indicating the comparison
  bulk.Tcells$comparison <- paste0(i, "_R_vs_NR")
  # Add a new column with gene names
  bulk.Tcells$gene_names <- rownames(bulk.Tcells)
  # Save the dataframe to CSV
  write.csv2(bulk.Tcells, file = paste0(output, i, "DESeq2_pseudobulk_R_vs_NR.csv"))
  
  # Store the dataframe in the list
  results_list[[i]] <- bulk.Tcells
}

# Combine all dataframes into one big dataframe
all_results_DEseq2 <- do.call(rbind, results_list)

# Save the combined dataframe to CSV
write.csv2(all_results_DEseq2, file = paste0(output, "Tcells_all_pseudobulk_results_DEseq2.csv"))


# HEATMAPS OF TOP 30 genes, p val < 0.05
########################################################################
#get expression matrix
expr_matrix <- GetAssayData(object = Seurat_NKT_pseudo, layer = "scale.data")
expr_df <- as.data.frame(expr_matrix)

# Define the cell types
cell_types <- unique(Seurat_NKT_noHC$cell_type)

# Iterate over each cell type
for (cell_type in cell_types) {
  # Extract the matching columns for the current cell type
  matching_columns <- colnames(expr_df)[endsWith(colnames(expr_df), cell_type)]
  
  # Subset expr_df to include only the specified cell type
  subset_expr_df <- expr_df[, matching_columns]
  
  # Select a subset of genes
  sig_genes <- subset(all_results_DEseq2, comparison == paste(cell_type, "_R_vs_NR", sep = "") & p_val < 0.05 &
                             pct.1 >= 0.2 & 
                             pct.2 >= 0.2)$gene_names
  sig_genes <- head(sig_genes, 30)
  
  # Subset the expression dataframe for selected genes
  subset_expr_df <- subset_expr_df[sig_genes, ]
  
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

# Get unique CellTypes R vs NR
cellTypes <- unique(all_results_DEseq2$comparison)

# Create an empty dataframe to store all results
all_results_df <- data.frame()

# Loop over each CellType for the current Comparison
for (i in cellTypes) {
  
  # Subset the data for the current CellType and pval < 0.05
  subset_df <- subset(all_results_DEseq2, comparison == i)
  subset_df2 <- subset(subset_df, p_val  < 0.05)
  
  # Order the data by avg_log2FC and pvalue for fgseaMultilevel
  geneList <- subset_df2 %>%
    filter(!is.na(gene_names), !is.na(avg_log2FC)) %>%  #remove genes and log2FC which is zero
    mutate(ranking_metric = -log10(p_val)*sign(avg_log2FC)) %>%   #multiply p value with -1 if avg_log2FC is negative, 0 if 0 and +1 if positive 
    group_by(gene_names) %>% 
    summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
    arrange(-ranking_metric) %>% # sort descending 
    tibble::deframe() # convert to named vector
  head(geneList)
  tail(geneList)

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
write.csv(all_results_df, "pathways_all_results_genes_5genes_min_ranking.csv", row.names = FALSE)

# Create a dotplot with all pathways
# add the count of the leading edge as new column to plot later
all_results_df[, leadingEdge_count := sapply(strsplit(leadingEdge, ","), length)]

# Define significance threshold and add a new column indicating significance
significance_threshold <- 0.05
all_results_df[, significant := padj < significance_threshold]

#subset on pathways with leading edge of at least 5 and at least one occurence of padj < 0.05
all_results_df_leadingEdge_5 <- subset(all_results_df, leadingEdge_count>=5)
all_results_df_leadingEdge_5_sig <- subset(all_results_df_leadingEdge_5, padj<0.05)
significant_pathways <- unique(all_results_df_leadingEdge_5_sig $pathway)
all_results_df_leadingEdge_5_sig <- all_results_df_leadingEdge_5[pathway %in% significant_pathways, ]

write.csv(all_results_df_leadingEdge_5_sig, "all_results_df_leadingEdge_5_sig.csv", row.names = FALSE)

# exclude non-relevant pathways for Dotplot
exclude <- c("HALLMARK_MYC_TARGETS_V2", "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_ESTROGEN_RESPONSE_EARLY")
all_results_df_leadingEdge_5_sig_sub <- all_results_df_leadingEdge_5_sig[!(all_results_df_leadingEdge_5_sig$pathway %in% exclude), ]


# plot and save dotplot
pt2 <- ggplot(all_results_df_leadingEdge_5_sig_sub, aes(x = CellType, y = pathway, color = NES, size = leadingEdge_count)) + 
  geom_point(data = subset(all_results_df_leadingEdge_5_sig_sub, padj < significance_threshold),  # Subset data for significant points
             aes(fill = NES, size = leadingEdge_count), colour = "red", pch = 21) +  # Specify aesthetics for significant points
  geom_point(data = subset(all_results_df_leadingEdge_5_sig_sub, padj >= significance_threshold),  # Subset data for non-significant points
             aes(fill = NES, size = leadingEdge_count), colour = "black", pch = 21) +  # Specify aesthetics for non-significant points
  scale_fill_gradient2(low ="#9d1c20" , mid = "grey95", high = "#81c4db") + 
  scale_size(range = c(2, 10)) +  # Adjust the range for the dot size
  xlab(" ") + ylab("pathway") + ggtitle("Tcell cells - pathways") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(size = guide_legend(title = "Leading Edge Count"))  # Include leadingEdge_count in legend with a title

# Show and save the plot
print(pt2)
ggsave("Pathway_dotplot_Tcells.pdf", pt2, width = 8, height = 8)

# Module score patient level analysis
########################################################################
#select all cell types which have significant pathways
all_results_df_leadingEdge_5_subset_sign <- subset(all_results_df_leadingEdge_5_sig, significant == "TRUE")
all_results_df_leadingEdge_5_subset_sign$cell_type <- gsub("_.*", "", all_results_df_leadingEdge_5_subset_sign$CellType)
cell_types <- unique(all_results_df_leadingEdge_5_subset_sign$cell_type)

#create empty dataframe
all_data_long <- data.frame()

# iterate over each celltype
for (cellsubset in celltypes) {
  # Subset the leading edge results and Seurat object for the current cell type
  subset_leadingEdge <- subset(all_results_df_leadingEdge_5_subset_sign, cell_type == cellsubset)
  Pseudo_subset <- subset(Seurat_NKT_pseudo, cell_type == cellsubset)
  
  # Initialize a list to store pathways and their leading edge genes
  pathways_with_leading_edge <- list()
  
  # Extract the leading edge genes for each significant pathway
  for (i in 1:nrow(subset_leadingEdge)) {
    pathway <- subset_leadingEdge$pathway[i]
    leading_edge_genes <- strsplit(subset_leadingEdge$leadingEdge[i], ",")[[1]]
    
    pathways_with_leading_edge[[pathway]] <- leading_edge_genes
  }
  
  # Add module scores for each pathway based on leading edge genes
  for (pathway_name in names(pathways_with_leading_edge)) {
    genes <- pathways_with_leading_edge[[pathway_name]]
    
    Pseudo_subset <- AddModuleScore(object = Pseudo_subset, 
                                    features = list(genes = genes), 
                                    name = paste0("pathway_", pathway_name), 
                                    assay = "RNA",
                                    scale = TRUE)
  }
  
  # Extract metadata and transform for plotting
  data <- Pseudo_subset@meta.data
  response <- rownames(data)
  data$response_3m <- data$orig.ident
  data$patient <- rownames(Pseudo_subset@meta.data)
  data$patient <- str_extract(data$patient, "(?<=_)[A-Za-z0-9-]+(?=_)")
  
  # Transform data to long format for ggplot
  data_long <- tidyr::pivot_longer(data, cols = starts_with("pathway_"), 
                                   names_to = "pathway", values_to = "zscore")
 
# Create a boxplot for pathway Z-scores by response
 plot1 <- ggplot(data_long, aes(x = pathway, y = zscore, fill = response_3m)) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha = 1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(color = response_3m, fill = response_3m), size = 0.8, shape = 21) +  
    scale_fill_manual(values = c("#9d1c20", "#81c4db")) +
    scale_color_manual(values = c("black", "black")) +  
    labs(title = paste0(cellsubset, ", Pathway Z-Scores by Response"), x = "Response", y = "Module score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.position = "none",
          aspect.ratio = 1,
          panel.spacing = unit(0.5, "cm")) +
    ylim(-0.9, 2)
  
 # Print and save plots
  print(plot1)
  ggsave(paste0("Module_score_plots_", cellsubset, ".pdf"), plot = plot1, width = 4, height = 8)
  
  # Perform t-tests for each pathway and save the results
  pathway_tests <- data_long %>%
    group_by(pathway) %>%
    summarize(p_value = t.test(zscore ~ response_3m)$p.value)
  
  write.csv2(pathway_tests, file = paste0("T_test_pathway_zscore_", cellsubset, ".csv"))
  
  # Combine all data into a single data frame
  all_data_long <- rbind(all_data_long, data_long)

}


# Select the top 10 significant modules
selected_pathways <- all_results_df_leadingEdge_5_subset_sign %>%
  arrange(padj) %>%
  slice(1:10)
selected_pathways <- selected_pathways %>%
  mutate(celltype_pathway = paste(cell_type, pathway, sep = "_"))
Top10_celltype_pathways <- selected_pathways$celltype_pathway

#select top 10 significant modules in the Zscore dataset
Zscore_Tcell_df <- as.data.frame(all_data_long)

# Clean pathway names
Zscore_Tcell_df <- Zscore_Tcell_df %>%
  mutate(cleaned_pathways = str_remove_all(pathway, "^pathway_|\\d$"),
         celltype_pathway = paste(cell_type, cleaned_pathways, sep = "_"))

# Filter top 10 significant pathways
Zscore_Tcell_df_top10 <- Zscore_Tcell_df %>%
  filter(celltype_pathway %in% Top10_celltype_pathways) %>%
  arrange(pathway) %>%
  select(cell_type, orig.ident, response_3m, patient, pathway, zscore, celltype_pathway)

# Save the filtered dataset as CSV file
write.csv(Zscore_Tcell_df_top10, "Tcells_zscore_top5_padj_pathways.csv", row.names = FALSE)






