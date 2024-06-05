########################################################################
#  DEG and pathway analysis Myeloid cells
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
# load myeloid Seurat object
Seurat_Myeloid <- readRDS("..data/seurat_myeloid.annot.rds")

# remove healthy control group for DEG analysis
Seurat_Myeloid_noHC <- subset(Seurat_Myeloid, subset = response_binary != "HC")

# Use aggregate Expression to pseudobulk counts based on response-patient-celltype
Seurat_Myeloid_pseudo <- AggregateExpression(Seurat_Myeloid_noHC, assays = "RNA", return.seurat = T, group.by = c("response_3m", "sample_id", "cell_type"))

#add new metadata column with cell_type and orig.ident
metadata <- Seurat_Myeloid_pseudo@meta.data
metadata  <- metadata  %>%
  mutate(cell_type = gsub(".*_", "", rownames(.)),  # Extract cell type from row names
         cell_type_response = paste(cell_type, orig.ident, sep = "_"))  # Combine cell type and orig.ident
Seurat_Myeloid_pseudo@meta.data <- metadata

# set new idents (response_patient_cell_type)
Idents(Seurat_Myeloid_pseudo) <- Seurat_Myeloid_pseudo@meta.data$cell_type_response

# Pseudobulk with DESeq2 and save excel sheets
output <- "../Output"

# Initialize an empty list to store results
results_list <- list()

# Perform analysis for each cell type
celltypes <- unique(Seurat_Myeloid_noHC$cell_type)  
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
  write.csv2(bulk.Mcells, file = paste0(output, i, "DESeq2_pseudobulk_R_vs_NR.csv"))
  
  # Store the dataframe in the list
  results_list[[i]] <- bulk.Mcells
}

# Combine all dataframes into one big dataframe
all_results_DEseq2 <- do.call(rbind, results_list)

# Save the combined dataframe to CSV
write.csv2(all_results_DEseq2, file = paste0(output, "Myeloid_cells_all_pseudobulk_results_DEseq2.csv"))

# HEATMAPS OF TOP 30 genes, p val < 0.05
########################################################################
#get expression matrix
expr_matrix <- GetAssayData(object = Seurat_Myeloid_pseudo, layer = "scale.data")
expr_df <- as.data.frame(expr_matrix)

# Define the cell types
cell_types <- unique(Seurat_Myeloid_noHC$cell_type)

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
  ggsave(paste0(cell_type, "_heatmap_top50.pdf"), pm, width = 8, height = 7)
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

# plot and save dotplot
significance_threshold <- 0.05
pt2 <- ggplot(all_results_df_leadingEdge_5_sig , aes(x = CellType, y = pathway, color = NES, size = leadingEdge_count)) + 
  geom_point(data = subset(all_results_df_leadingEdge_5_sig, padj < significance_threshold),  # Subset data for significant points
             aes(fill = NES, size = leadingEdge_count), colour = "red", pch = 21) +  # Specify aesthetics for significant points
  geom_point(data = subset(all_results_df_leadingEdge_5_sig, padj >= significance_threshold),  # Subset data for non-significant points
             aes(fill = NES, size = leadingEdge_count), colour = "black", pch = 21) +  # Specify aesthetics for non-significant points
  scale_fill_gradient2(low = "#9d1c20", mid = "grey95", high = "#81c4db") + 
  scale_size(range = c(2, 10)) +  # Adjust the range for the dot size
  xlab(" ") + ylab("pathway") + ggtitle("Myeloid cells - pathways") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(size = guide_legend(title = "Leading Edge Count"))  # Include leadingEdge_count in legend with a title

# Show the plot
print(pt2)
ggsave("Pathway_dotplot_Myeloidcells.pdf", pt2, width = 8, height = 8)

# remove Myogenesis and plot again
all_results_df_leadingEdge_5_sig_no_Myogenesis <- subset(all_results_df_leadingEdge_5_sig, pathway != "HALLMARK_MYOGENESIS")
pt3 <- ggplot(all_results_df_leadingEdge_5_sig_no_Myogenesis , aes(x = CellType, y = pathway, color = NES, size = leadingEdge_count)) + 
  geom_point(data = subset(all_results_df_leadingEdge_5_sig_no_Myogenesis, padj < significance_threshold),  # Subset data for significant points
             aes(fill = NES, size = leadingEdge_count), colour = "red", pch = 21) +  # Specify aesthetics for significant points
  geom_point(data = subset(all_results_df_leadingEdge_5_sig_no_Myogenesis, padj >= significance_threshold),  # Subset data for non-significant points
             aes(fill = NES, size = leadingEdge_count), colour = "black", pch = 21) +  # Specify aesthetics for non-significant points
  # scale_fill_viridis_c(option = "mako") +  # Use magma color palette
  scale_fill_gradient2(low = "#9d1c20", mid = "grey95", high = "#81c4db") + 
  scale_size(range = c(2, 10)) +  # Adjust the range for the dot size
  xlab(" ") + ylab("pathway") + ggtitle("Myeloid cells - pathways") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(size = guide_legend(title = "Leading Edge Count"))  # Include leadingEdge_count in legend with a title

# Show and save the plot
print(pt3)
ggsave("Pathway_dotplot_Myeloidcells_no_myogenesis.pdf", pt3, width = 8, height = 8)


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
  subset_leadingEdge <- subset(all_results_df_leadingEdge_5_subset_sign, cell_type == cellsubset)
  Pseudo_subset <- subset(Seurat_Myeloid_pseudo, cell_type == cellsubset)
  
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
  plot1 <- ggplot(data_long, aes(x = response_3m, y = zscore, fill = response_3m)) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha = 1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(color = response_3m, fill = response_3m), size = 0.8, shape = 21) +  
    scale_fill_manual(values = c("#9d1c20", "#81c4db")) +
    scale_color_manual(values = c("black", "black")) +  
    labs(title = paste0(cellsubset, ", Pathway Z-Scores by Response"), x = "Response", y = "Zscore") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.position = "none",
          aspect.ratio = 1,
          panel.spacing = unit(0.5, "cm")) +
    facet_wrap(~ pathway, scales = "free_y", ncol = 3, as.table = TRUE) +
    ylim(-0.9, 2)
  
  # Print and save plots
  print(plot1)
  ggsave(paste0("Zscore_plots_", cellsubset, ".pdf"), plot = plot1, width = 8, height = 4)
  
  # Perform t-tests for each pathway and save the results
  pathway_tests <- data_long %>%
    group_by(pathway) %>%
    summarize(p_value = t.test(zscore ~ response_3m)$p.value)
  
  write.csv2(pathway_tests, file = paste0("T_test_pathway_zscore_", cellsubset, ".csv"))
  
  # Combine all data into a single data frame
  all_data_long <- rbind(all_data_long, data_long)
}

write.csv(all_data_long, "all_data_long.csv", row.names = FALSE)

# plot selected CD16 mono pathways for figure
all_data_long_df <- as.data.frame(all_data_long)
selected_figure <- c("pathway_HALLMARK_INFLAMMATORY_RESPONSE1", "pathway_HALLMARK_TNFA_SIGNALING_VIA_NFKB1", "pathway_HALLMARK_INTERFERON_GAMMA_RESPONSE1")
all_data_long_CD16 <- subset(all_data_long_df, cell_type == "CD16 Mono" & pathway %in%  selected_figure)

point_fill_colors <- c("#9d1c20", "#81c4db")
bp2 <-  ggplot(all_data_long_CD16 , aes(x = pathway, y = zscore, fill = response_3m)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
             aes(color = response_3m, fill = response_3m), size = 3, shape = 21) +  # Set shape to 21 for filled circles
  scale_fill_manual(values = c("#9d1c20", "#81c4db")) +
  scale_color_manual(values = c("black", "black")) +  # Set border color to black
  ylim(-0.15, 0.4) +
  labs(title = "Zscore CD16 Monos", x = "Pathways", y = "Zscore") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp2
ggsave(("selected_CD16_Zscore_Figure.pdf"), plot = bp2, width = 6, height = 10)

# Select the top 10 significant modules
selected_pathways <- all_results_df_leadingEdge_5_subset_sign %>%
  arrange(padj) %>%
  slice(1:10)
selected_pathways <- selected_pathways %>%
  mutate(celltype_pathway = paste(cell_type, pathway, sep = "_"))
Top10_celltype_pathways <- selected_pathways$celltype_pathway

#select top 10 significant modules in the Zscore dataset
Zscore_Mcell_df <- as.data.frame(all_data_long)

# Clean pathway names
Zscore_Mcell_df <- Zscore_Mcell_df %>%
  mutate(cleaned_pathways = str_remove_all(pathway, "^pathway_|\\d$"),
         celltype_pathway = paste(cell_type, cleaned_pathways, sep = "_"))

# Filter top 10 significant pathways
Zscore_Mcell_df_top10 <- Zscore_Mcell_df %>%
  filter(celltype_pathway %in% Top10_celltype_pathways) %>%
  arrange(pathway) %>%
  select(cell_type, orig.ident, response_3m, patient, pathway, zscore, celltype_pathway)

# Save the filtered dataset as CSV file
write.csv(Zscore_Mcell_df_top10, "Myeloid_zscore_top5_padj_pathways.csv", row.names = FALSE)

