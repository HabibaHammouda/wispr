library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(tidyverse)
library(reshape2)
library(S4Vectors)
library(png)
library(RColorBrewer)
library(igraph)
library(reticulate)
library(Giotto)
require(data.table)
library(usethis) 
library(stringr)
usethis::edit_r_environ()

#load spatial data

data_dir1 <- "filtered_feature_bc_matrix/"
list.files(data_dir1) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix1 <- Read10X(data.dir = data_dir1)

coor_ST_mouse <- read.csv("tissue_positions_list.csv", header = FALSE) ##x- y- coordinates of the spots
rownames(coor_ST_mouse) <- coor_ST_mouse$V1
coor_ST_mouse <- coor_ST_mouse[(coor_ST_mouse$V2=="1"),]

#align coordinates
match_n <- match(colnames(expression_matrix1),rownames(coor_ST_mouse))
coor_ST_mouse_reor <- coor_ST_mouse[match_n,]


write_csv(cbind((as.data.frame(coor_ST_mouse_reor))), "tissue_positions_list_aligned.csv")
coor_ST_mouse2 <- read.table(paste0("tissue_positions_list_aligned.csv"))


counts.meta_ST.mouse <-CreateSeuratObject(counts = expression_matrix1, assay = "RNA")
counts.meta_ST.mouse <- AddMetaData(object = counts.meta_ST.mouse, metadata = coor_ST_mouse$V5, col.name = "X_coor")
counts.meta_ST.mouse <- AddMetaData(object = counts.meta_ST.mouse, metadata = coor_ST_mouse$V6, col.name = "Y_coor")

###read scRNA SEQ

counts <- fread("20220926133842848471.cnt_data.dat", check.names=FALSE)
header <- read.table(paste0("20220926133842848471.cnt_data_col1.tsv"))
counts <- t(counts)
colnames(counts) <- header$V1

counts <-counts[,-1]

counts_meta <- read.table(paste0("20220926133842848471.mta_data.tsv")) #read filtered meta data
rownames(counts_meta) = counts_meta$V1
counts_meta <- dplyr::select(counts_meta, -1)

counts.metasc <-CreateSeuratObject(counts = counts)

#save rds
saveRDS(counts.metasc, "seurat_counts.rds")
counts.metasc <- readRDS("seurat_counts.rds")

colnames(counts.metasc) == rownames(counts_meta)
#set idents
counts.metasc <- AddMetaData(object = counts.metasc, metadata = counts_meta$V2, col.name = "cellname")
counts.metasc <- SetIdent(counts.metasc, value = "cellname")
Idents(counts.metasc) <- factor(x = Idents(counts.metasc), levels = sort(levels(counts.metasc)))

counts.metasc2 <- counts.metasc[-1,]
#HVG selection and dimensionality reduction
counts.metasc2 <- FindVariableFeatures(counts.metasc2, nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 10), dispersion.cutoff = c(0.5, Inf))

counts.metasc2 <- ScaleData(counts.metasc2, vars.to.regress = c("nFeature_RNA"), features = VariableFeatures(counts.metasc2))
counts.metasc2 <- RunPCA(counts.metasc2, features = VariableFeatures(counts.metasc2))
counts.metasc2 <- RunUMAP(counts.metasc2, reduction = "pca", n.neighbors = 5, features = VariableFeatures(object = counts.metasc2))
counts.metasc2 <- RunTSNE(counts.metasc2, perplexity = 75, use.pca = TRUE,  features = VariableFeatures(object = counts.metasc2))

#plotting the clusters
color2 <- DiscretePalette(30, palette = "glasbey")
p1 <- DimPlot(counts.metasc2, reduction = "umap", pt.size = 1.5,  group.by = "cellname")
p1

#find intersecting genes between sc and st datasets
inter <- intersect(rownames(counts.metasc), rownames(counts.meta_ST.mouse))
counts.metasc0 <- counts.metasc[inter,]
spat.st2 <- counts.meta_ST.mouse[inter,]

cell_types <- as.factor(counts_meta$V2)

counts.metasc0@meta.data$seurat_clusters <- cell_types

counts.metasc0 <- SetIdent(counts.metasc0, value = "seurat_clusters")
Idents(counts.metasc0) <- factor(x = Idents(counts.metasc0), levels = sort(levels(counts.metasc0)))

#rename the cluster ids
df <- as.data.frame(counts.metasc0$seurat_clusters)
head(df)
colnames(df) <- col
majors <- df %>% extract(`counts.metasc0$seurat_clusters`, into = c("Main", "Num"), "(.*)_([^_]+)$")
counts.metasc0@meta.data$cells <- majors$Main

counts.metasc2 <- ScaleData(counts.metasc0, vars.to.regress = c("nFeature_RNA"), features = topgenes_gini$genes)
counts.metasc2 <- RunPCA(counts.metasc2, features = topgenes_gini$genes)
counts.metasc2 <- RunUMAP(counts.metasc2, reduction = "pca", n.neighbors = 20, features = topgenes_gini$genes, slot = "counts")
counts.metasc2 <- RunTSNE(counts.metasc2, perplexity = 75, use.pca = TRUE,  features = topgenes_gini$genes)
color2 <- DiscretePalette(8, palette = "glasbey")
p1 <- DimPlot(counts.metasc2, reduction = "umap", pt.size = 1.5,  group.by = "cells", cols = color2)
p1

#path to python
Sys.setenv(RETICULATE_PYTHON = "/Users/alperenerdogan/Library/r-miniconda/envs/giotto_env/bin/pythonw")

reticulate::virtualenv_create("/Users/alperenerdogan/Library/r-miniconda/envs/giotto_env/bin/", python = "/Users/alperenerdogan/Library/r-miniconda/envs/giotto_env/bin/pythonw")

reticulate::use_virtualenv("/Users/alperenerdogan/Library/r-miniconda/envs/giotto_env/bin/")

py_available()

my_python_path = "/Users/alperenerdogan/Library/r-miniconda/envs/giotto_env/bin/pythonw"

instrs = createGiottoInstructions(python_path = my_python_path)

sc_matrix<-counts.metasc0@assays$RNA@counts
sc_lable<-counts.metasc0$cellname %>% as.factor()

sc_matrix2 <- createGiottoObject(raw_exprs = sc_matrix,instructions = instrs)
sc_matrix2norm <- normalizeGiotto(gobject = sc_matrix2)
sc_matrix2norm@cell_metadata$leiden_clus<-as.character(sc_lable)
sc_matrix2norm@cell_metadata$major_clus <- majors$Main

gini_markers_subclusters = findMarkers_one_vs_all(gobject = sc_matrix2norm,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus')


topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']

sc_norm_exp <- 2^(sc_matrix2norm@norm_expr)-1
ExprSubset<-sc_norm_exp[as.character(topgenes_gini$genes),]

Sig<-NULL
for (i in as.character(unique(sc_lable))){
  Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(sc_lable==i)]))))
}
colnames(Sig)<-unique(sc_lable)

sc_matrix2norm <- runPCA(gobject = sc_matrix2norm, feats_to_use = topgenes_gini$genes, scale_unit = FALSE, center = FALSE, ncp = 20)
# Run UMAP
sc_matrix2norm <- runUMAP(gobject = sc_matrix2norm, expression_values = c("normalized"), n_neighbors = 20, dimensions_to_use = 1:15, n_components = 3, n_threads = 6)

#3D plot
plotUMAP_2D(gobject = sc_matrix2norm, cell_color = 'major_clus', point_size = 3,
            show_center_label = F, save_param = list(save_name = 'brain_majorclust_UMAP_3D'))

df.m <- reshape2::melt(Sig, id.vars = NULL)
df.m_vascular70 <- df.m[(df.m$Var1=="Acsbg1"),]
deduped.data <- unique( df.m_vascular70[, 1:3 ] )

max(deduped.data$value)

####circular plotting of genes per clusters

label_data <- deduped.data
rownames(label_data)
label_data$id <- c(1:nrow(label_data))
label_data$id
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
angle
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #
label_data

# Start the plot
p <- ggplot(label_data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("#D55E00", 1)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20,70) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=value+10, label=Var2, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data$angle, inherit.aes = FALSE ) 

p


###end circular plot

st_matrix<-spat.st2@assays$RNA@counts %>% as.matrix()
st_matrix <- as.data.frame(st_matrix)
st_matrix2 <- st_matrix
sort(colSums(st_matrix2), decreasing = FALSE)

st_matrix2 <- createGiottoObject(raw_exprs = st_matrix2,instructions = instrs)
st_matrix2norm <- normalizeGiotto(gobject = st_matrix2)
all(is.finite(st_matrix2norm@norm_expr))

spot_giotto <- calculateHVG(gobject = st_matrix2norm)
gene_metadata = fDataDT(spot_giotto)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
spot_giotto <- runPCA(gobject = spot_giotto, genes_to_use = featgenes, scale_unit = F)
signPCA(spot_giotto, genes_to_use = featgenes, scale_unit = F)
spot_giotto <- createNearestNetwork(gobject = spot_giotto, dimensions_to_use = 1:10, k = 10)

spot_giotto <- Giotto::doLeidenCluster(gobject = spot_giotto, resolution = 0.2, n_iterations = 1000)
gini_markers_subclusters_st = findMarkers_one_vs_all(gobject = spot_giotto,
                                                     method = 'gini',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus')

topgenes_gini_st = gini_markers_subclusters_st[, head(.SD, 100), by = 'cluster']
st_norm_exp_st <- 2^(st_matrix2norm@norm_expr)-1
ExprSubset_st<-st_norm_exp_st[as.character(topgenes_gini_st$genes),]
st_lable <- spot_giotto@cell_metadata$leiden_clus
allmarkers_unique <- unique(cbind(c(topgenes_gini$genes, topgenes_gini_st$genes)))

ExprSubset_st2<-st_norm_exp_st[as.character(allmarkers_unique),]

Sig_st2<-NULL
for (i in as.character(unique(unique(spot_giotto@cell_metadata$leiden_clus)))){
  Sig_st2<-cbind(Sig_st2,(apply(ExprSubset_st2,1,function(y) mean(y[which(st_lable==i)]))))
}

write_csv(cbind(rownames(Sig_st2), as.data.frame(Sig_st2)), "mouse_STSeq_WISpR_gene_vs_spatclusts_visium.dat")
write_csv(cbind((as.data.frame(spot_giotto@cell_metadata$leiden_clus))), "mouse_STSeq_WISpR_spot_vs_spatclusts2_visium.dat")


spat2 <- spat.st2[allmarkers_unique, ]
sc_matrix_2 <- counts.metasc0[allmarkers_unique, ]


AE <- AverageExpression(object =sc_matrix_2, features = allmarkers_unique, 
                        return.seurat = FALSE, group.by = "cellname", slot = "counts", verbose = TRUE)

match <- match(rownames(AE$RNA),rownames(spat2))
spat.st2_AE25 <- spat2[match,]

rownames(AE$RNA)==rownames(spat.st2_AE25)

write_csv(cbind(rownames(AE$RNA), as.data.frame(AE$RNA)), "mouse_scRNASeq_WISpR_visium.dat")
write_csv(cbind(rownames(spat.st2_AE25@assays$RNA@counts), (as.data.frame(spat.st2_AE25@assays$RNA@counts))), "mouse_STSeq_WISpR_visium.dat")
