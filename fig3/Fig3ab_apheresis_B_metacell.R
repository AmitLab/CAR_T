

library(dplyr)
library(patchwork)
library(ggplot2)
library(metacell)

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


## Metacell construction
#########################################################################################################

newumis <- Matrix::readMM("B.ds300.raw.counts.mtx")
metadata = read.csv("B.ds300.raw.metadata.csv") %>% 
          tibble::column_to_rownames(var = "X")
colnames(newumis) <- read.csv("B.ds300.raw.Cells.csv")[[1]]
rownames(newumis)  <- read.csv("B.ds300.raw.Genes.csv")[[1]]

amat = tgScMat(newumis, stat_type = "umi", cell_metadata = metadata)
scdb_add_mat(id, amat)

mat <- scdb_mat(id) 

mcell_plot_umis_per_cell(id)
min_umis_per_cell <- 500
species = "Hs"
mt.fraction.cutoff <- 0.15
T_vm = 0.1
T_tot = 150
T_top3 = 3

n_edges = 100 # K edges per cell (100)
min_mc_size = 20 # 20
n_resamp = 500 #  500
dsamp_graph = TRUE # FALSE if number of UMIs limits analysis (1000)
alpha = 2 # when alpha is larger fewer clusters are generated. (2)

# remove bad cells
mat = scdb_mat(id)
nms = mat@genes
umis <- mat@mat

umicount = colSums(as.matrix(umis))
cell_stats = mat@cell_metadata

mt_genes = grep("^MT-", rownames(umis), v=T)
if(length(mt_genes)>0){
  mt_count = colSums(as.matrix(umis[ mt_genes, ]))
  bad_cells = unique(names(which(((mt_count / umicount) > mt.fraction.cutoff) | (umicount < min_umis_per_cell))))  
} else {
  bad_cells = unique(names(which(umicount < min_umis_per_cell)))  
}

if(length(bad_cells) > 1){
  mcell_mat_ignore_cells(id, id, ig_cells = bad_cells)
}

# remove bad_genes
if(species  == "Hs"){
  ig_genes = c(grep("^IGHV",nms,v=T), grep("^IGK", nms, v=T), grep("^IGL", nms, v=T))
  bad_genes = unique(c(grep("MIR|-PS|RPL|RPS|JCHAIN|FTL1|MT-|MTMR|MTND|MALAT1|XIST|LINC|RP11-|orf|^AF\\d|^AC\\d|^AL\\d", mat@genes, v=T)))
  bad_genes <- c(bad_genes, ig_genes)
} else if (species == "Mm"){
  ig_genes = c(grep("^Igj", nms, v=T), grep("^Igh",nms,v=T), grep("^Igk", nms, v=T), grep("^Igl", nms, v=T))
  bad_genes = unique(c(grep("Gm[0-9].|Mir|-ps|Rpl|Rps|Jchain|Ftl1|mt-|Hsp", mat@genes, v=T), ig_genes))
}

v <- setNames(rep(1,length(bad_genes)), nm = bad_genes)
g <- gset_new_gset(v, "bad_genes")
scdb_add_gset(gset = g, id=paste0(id, "_bad_genes"))
mcell_mat_ignore_genes(new_mat_id=id, mat_id=id, bad_genes, reverse=F) 

## feature gene selection
mcell_add_gene_stat(gstat_id=id, mat_id=id, force=T)
mcell_gset_filter_varmean(gset_id = id, gstat_id=id, T_vm=T_vm, force_new=T) 
mcell_gset_filter_cov(gset_id = id, gstat_id=id, T_tot=T_tot, T_top3=T_top3) 
genes_gn <- scdb_gset(id)
message(paste0("number of features:", length(genes_gn@gene_set), " genes"))
mcell_plot_gstats(gstat_id=id, gset_id=id)

## lateral geneset to filter
ext_lateral_genes_nm = "EA" #early response / stress genes
if(!is.null(ext_lateral_genes_nm)){
  if (ext_lateral_genes_nm == "EA") {
    ext_lateral_genes = c("FOS","JUN","DUSP1","JUNB","CD69")
  }
  lateral_genes <- setNames(rep(1,length(ext_lateral_genes)), nm = ext_lateral_genes)
  lateral_gset = gset_new_gset(sets = lateral_genes, desc = ext_lateral_genes_nm)
  scdb_add_gset(ext_lateral_genes_nm, lateral_gset)
}

if(!is.null(ext_lateral_genes_nm)) {
  
  lateral_gset <- scdb_gset(id = ext_lateral_genes_nm)
  marker_gset = scdb_gset(id) 
  marker_gset = gset_new_restrict_gset(gset = marker_gset, filt_gset = lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
  scdb_add_gset(paste0(id, "_lateral_filtered"), marker_gset)
  mcell_add_cgraph_from_mat_bknn(mat_id=id, 
                                 gset_id = paste0(id, "_lateral_filtered"), 
                                 graph_id=id,
                                 K=n_edges,
                                 dsamp=dsamp_graph) 
} 


mcell_coclust_from_graph_resamp(
  coc_id=id, 
  graph_id=id,
  min_mc_size=min_mc_size, 
  p_resamp=0.75, 
  n_resamp=n_resamp)

mcell_mc_from_coclust_balanced(
  coc_id=id, 
  mat_id= id,
  mc_id=id, 
  K=n_edges, 
  min_mc_size=min_mc_size,
  alpha=alpha)


## Fig 3b - metacell heatmap
################################################################################

b_naive <- c("TCL1A", "FCER2", "PLPP5", "IGHD",  "IL4R",  "CXCR4", "BACH2", "TSPAN13")
b_memory <- c("CRIP1","AIM2","ANXA2", "TNFRSF13B","ITGB1","CD27","CD82","IGHG1","GPR183" )
b_other <- c("CD19")
b_tumor <- c("BCL2","TCF4","FKBP5","CARD11", "JAK2")
b_markers <- c(b_naive,b_memory,b_other,b_tumor)
fig_name = paste0(figdir, "/", id, ".Fig3b.png")
mcell_mc_plot_marks(mc_id=id, mat_id=id, add_metadata = "sample_name", gene_list = b_markers, fig_fn = fig_name  )


## Metacell projection 
################################################################################
mat = scdb_mat(id)
mcell_mc2d_force_knn(mc2d_id=id, mc_id=id, graph_id=id)
tgconfig::set_param("mcell_mc2d_height", 800, "metacell")
tgconfig::set_param("mcell_mc2d_width", 800, "metacell")
fig_name = paste0(figdir, "/", id, ".Fig3a.png")
mcell_mc2d_plot(mc2d_id = id, fig_fn = fig_name, show_mcid = TRUE)



##### Metacell projection with sample ids
###############################################################################################

mat <- scdb_mat(id)
mc2d = scdb_mc2d(id)

mat@cell_metadata$sc_x = mc2d@sc_x
mat@cell_metadata$sc_y = mc2d@sc_y

data <- mat@cell_metadata %>%   
  tibble::rownames_to_column(var = "cell_id") %>% 
  mutate(sample_alias = case_when(
    startsWith(sample_name, "HC") ~ "Healthy",
    TRUE ~ sample_name
  ))  %>%
  select(cell_id, sc_x, sc_y, disease_status, sample_name, sample_alias) %>% 
  tibble::column_to_rownames(var = "cell_id")


## sample label locations
xynames <- setNames(object = c('sc_x','sc_y'), nm = c("x","y"))
groups <- unique(x = data[,"sample_alias"])
labels.loc <- lapply(X = groups, FUN = function(group) {
  data.use <- data[data[, "sample_alias"] == group, , drop = FALSE]
  data.medians <- as.data.frame(x = t(x = apply(X = data.use[, xynames, drop = FALSE], MARGIN = 2, FUN = median, na.rm = TRUE)))
  data.medians[, "sample_alias"] <- group
  return(data.medians)
})
labels.loc <- lapply(X = labels.loc, FUN = function(x) {
    group.data <- data[as.character(x = data[, "sample_alias"]) == as.character(x["sample_alias"]), ]
    nearest.point <- RANN::nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1, 2)]), k = 1)$nn.idx
    furthest.point <- RANN::nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1, 2)]), k = dim(group.data)[1])$nn.idx[dim(group.data)[1]]
    x[1:2] <- group.data[nearest.point, 1:2]
    return(x)
  })
labels.loc <- do.call(what = "rbind", args = labels.loc)


# Meta cell cover of a given cell graph (or more generally of a scRNA data matrix)
mc = scdb_mc(id) 
mc_compute_fraction <- function (mc, scmat, f) 
{
  tb = table(scmat@cell_metadata[names(mc@mc), f], mc@mc)
  n_bc = matrix(tb, ncol = dim(tb)[2])
  rownames(n_bc) = dimnames(tb)[[1]]
  return(n_bc)
}
mc_f <- mc_compute_fraction(mc, mat, f = "disease_status")
mc_f <- apply(mc_f, 1, function(x) x/sum(x)) * sum(mc_f)  # class normalization
mc_f_ratio <- apply(mc_f, 1, function(x) {x/sum(x)})
mc_df <- data.frame(mc_x = mc2d@mc_x, mc_y = mc2d@mc_y) %>% tibble::rownames_to_column(var = "mc_id")
mc_df <- cbind(mc_df,t(mc_f_ratio))

sc_coord <- data
rownames(sc_coord) <- NULL
sc_coord$source = "sc"
mc_coord <- mc_df %>% dplyr::rename(sc_x = mc_x, sc_y = mc_y) %>% dplyr::select(sc_x, sc_y)
mc_coord$disease_status <- "UNDEF"
mc_coord$sample_name <- NA
mc_coord$sample_alias <- NA
mc_coord$source = "mc"
df <- rbind(sc_coord, mc_coord)

plot_edges = TRUE
if (plot_edges) {
  min_edge_l = 0
  edge_w = 1
  short_edge_w = 0
  fr = mc2d@graph$mc1
  to = mc2d@graph$mc2
  dx = mc2d@mc_x[fr] - mc2d@mc_x[to]
  dy = mc2d@mc_y[fr] - mc2d@mc_y[to]
  f = sqrt(dx * dx + dy * dy) > min_edge_l
  grps = seq(1,length(mc2d@graph$mc1))
  df_from <- cbind(x = mc2d@mc_x[fr], y = mc2d@mc_y[fr], grp = grps)
  df_to <- cbind(x = mc2d@mc_x[to], y = mc2d@mc_y[to], grp = grps)
  segments_df <- rbind(df_from,df_to) %>% as.data.frame()
  rownames(segments_df) <- NULL
}


## Fig 3a
library(ggrepel)
library(scales)
pal <- colorRampPalette(c("darkblue", "white","darkred"))(11)
p <- ggplot(df %>% dplyr::filter(source != "mc"), aes(x=sc_x, y=sc_y)) + #, color = disease_status color="#636363",
  geom_point(alpha = 1, color = "gray") +
  geom_line(data = segments_df, aes(x, y, group = grp), colour = "#636363", alpha = 0.75, linetype="dashed") + #
  geom_point(data = mc_df, mapping = aes(x = mc_x, y = mc_y, color = Malignant ), alpha = 1, na.rm = FALSE, show.legend = TRUE,size = 12.5) +
  scale_colour_gradientn(colours = pal ) +
  geom_text_repel(data = labels.loc, mapping = aes_string(x = xynames["x"], y = xynames["y"], label = "sample_alias"),
                  na.rm = FALSE, show.legend = TRUE, color="black",  nudge_x = 0, nudge_y = 0, size = 7, max.overlaps = 15 ) 
  #theme(legend.position="none") # 

plot(p)


## Fig 3a disease/healthy subsets
p1 <- ggplot(data, aes(x=sc_x, y=sc_y)) +
  geom_point(color = "gray") +  
  geom_point(data = data %>% dplyr::filter(disease_status == "Healthy"), aes(x=sc_x, y=sc_y) ,color = as.character(disease_status_name2color["Healthy"]))  +
  theme(legend.position="none")
p2 <- ggplot(data, aes(x=sc_x, y=sc_y)) +
  geom_point(color = "gray") +  
  geom_point(data = data %>% dplyr::filter(disease_status == "Malignant"), aes(x=sc_x, y=sc_y) ,color = as.character(disease_status_name2color["Malignant"]))  +
  theme(legend.position="none")

plot(p1+p2)


