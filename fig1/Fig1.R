

library(dplyr)
library(ggplot2)
library(Seurat)

clinical_response_status_name2color <- setNames(c("pink","#90ee90", "#71cceb" ), c("NR","PR", "CR"))
relapse_name2color <- setNames(c("#b30000","#336699"), c("NR","R"))
disease_status_name2color <- c("Healthy" = "darkblue", "Malignant" = "darkred")

## Fig 1b - swimmer plot 
################################################
dat <- read.csv("sample_response_duration.csv")
dat$sample_name = factor(dat$sample_name, levels=dat$sample_name[order(dat$Product, dat$n_months_to_event)])

dat.m = dat %>% mutate(n_months_to_event = ifelse(n_months_to_event >= 1, 1, 0),
                       n_days_to_event = ifelse(n_days_to_event >= 30, 30, 0))
n_days_response_cutoff = 90
intercept_n_days_to_event = n_days_response_cutoff
intercept_n_months_to_event = n_days_response_cutoff/30

p <- ggplot(dat, aes(sample_name, n_months_to_event)) +
      geom_bar(stat="identity", aes(fill = response_3m), width=0.5) +
      geom_point(data=dat.m, aes(sample_name, n_months_to_event, colour=clinical_response_1m), size=3) +
      geom_segment(data=dat %>% filter((Relapse_progression == "R") & (last_event != "death")), 
                   aes(x=sample_name, xend=sample_name, y=n_months_to_event + 0.2, yend=n_months_to_event + 1.5, linetype="arrow"), 
                   pch=25, linewidth=0.5, arrow=arrow(type="closed", length=unit(0.1,"in"))) +
      geom_segment(data=dat %>% filter((Relapse_progression == "NR") & (last_event == "relapse")), 
                   aes(x=sample_name, xend=sample_name, y=n_months_to_event + 0.2, yend=n_months_to_event + 0.8, color = Relapse_progression), 
                   pch=100, size=1) + 
      geom_hline(yintercept=intercept_n_months_to_event, linetype="dashed", color = "black")+
      coord_flip() +
      scale_fill_manual(values = relapse_name2color)+
      scale_color_manual(values = clinical_response_status_name2color)+
      scale_y_continuous(limits=c(-1,40), breaks=seq(0,40,2)) +
      labs(fill="Patient response\n(3 months)", colour="1 Month PET/CT",  
           x="", y = "Months to event") +
      theme_bw() +
      theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.text.x=element_text(angle=90,hjust=1,vjust=0.5) 
      ) 



## Fig 1c - UMAP
################################################
 
seurat <- readRDS(file= paste0("seurat.annot.rds") )

cell_annot_colormap = list(
  "B" = "#ff2500",
  "CD14 Mono" = "#bd8480",
  "CD16 Mono" = "#7C0A02",
  "Inter Mono" = "#a3534d", 
  "DCs" = "#e06666",
  "pDC" = "gray", 
  "Tregs" = "#008c5f", 
  "MAIT" =	"#07d0c7",
  "CD4.Naive" = "#B0C399", 
  "CD4.CM" = "#72A400", 
  "CD4.EM" = "#557b00",
  "CD8.naive" =  "#99c7d1",
  "CD8.CM" = "#328fa3", 
  "CD8.EM" = "#00677e",
  "Cycling" = "#004858", 
  "NK.bright" = "#251188",
  "NK.dim" = "#bdb7db"	
)
cell_types_levels <- c(
  "CD4.Naive", "CD4.CM",  "CD4.EM", "Tregs" , 
  "CD8.naive", "CD8.CM", "MAIT", "CD8.EM",  "Cycling", "NK.dim", "NK.bright",
  "B", "CD14 Mono", "Inter Mono", "CD16 Mono","DCs", "pDC" )
cell_annot_colormap = cell_annot_colormap[cell_types_levels]
label_colors <- as.character(unlist(cell_annot_colormap[names(cell_annot_colormap)]))
Idents(seurat) <- seurat$cell_type_annot
DimPlot(seurat, reduction = "umap", 
             label = F, label.box = F,repel = T,
             label.size = 3, cols = label_colors) 


## Fig 1e 
################################################################

library(sccomp)
seurat <- readRDS(file= paste0("seurat.annot.rds") )

cells.use <- seurat@meta.data %>% 
  tibble::rownames_to_column(var = "cell_id") %>% 
  dplyr::filter(product %in% c("Axi","Tis"), response_3m %in% c("NR","R")) %>% 
  pull(cell_id)
seurat_response <- subset(seurat, cells = cells.use)


for (p in c("Tis","Axi")){
  
  print(p)
  obj <- subset(seurat_response, product == p)
  set.seed(167534)
  res <- obj %>% droplevels() |>
    sccomp_glm( 
      formula_composition = ~ response_3m, 
      .sample =  sample_name, 
      .cell_group = cell_group_l2,
      bimodal_mean_variability_association = TRUE,
      test_composition_above_logit_fold_change = 0.2, #reducing the threshold imrpoves outcomes
      check_outliers = FALSE,
      cores = 4
    )
  res <- res[, colnames(res)[1:length(res)-1]] %>% as.data.frame() %>% mutate(Product = p) %>% rename(cell_type = cell_group_l2)
  #write.csv(res, file = ...)
  
}

