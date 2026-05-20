library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(tidyverse)
library(reshape2)
library(S4Vectors)
library(png)
library(RColorBrewer)

counts_meta <- read.csv("counts.csv", quote = "", sep = "\t") #read filtered meta data
rownames(counts_meta) = counts_meta$cell
counts_meta <- dplyr::select(counts_meta, -1)
head(counts_meta)

singlets <- readRDS(file = ("Heart/singlets.rds"))
singlets <- RunPCA(singlets, features = VariableFeatures(object = singlets))
singlets@meta.data$seurat_clusters <- counts_meta$origident
Idents(singlets) <- factor(x = Idents(singlets), levels = sort(levels(singlets)))
singlets <- FindVariableFeatures(singlets, nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 10), dispersion.cutoff = c(0.5, Inf))

singlets <- RunUMAP(singlets, reduction = "pca", n.neighbors = 5, features = VariableFeatures(object = singlets))
singlets <- RunTSNE(singlets, perplexity = 100, use.pca = FALSE,  features = VariableFeatures(object = singlets))

color2 <- DiscretePalette(17, palette = "glasbey")
p1 <- DimPlot(singlets, reduction = "umap", pt.size = 2, label = FALSE, group.by = "seurat_clusters", cols = color2)
p2 <- DimPlot(singlets, reduction = "tsne", pt.size = 1.5, group.by = "seurat_clusters", cols = color2)


p1 & theme(text = element_text(face = "bold"),
         axis.text.x=element_text(angle=45, hjust=1, size=25),
         axis.title = element_text(size=20,face="bold"),
         axis.title.y.right = element_text(size = 25),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.line = element_line(size=2))

write_csv(cbind(rownames(p1$data), p1$data), "Heart_sc_clustering_UMAP_Fig1a.csv")



Sys.setenv(RETICULATE_PYTHON = "/Applications/anaconda3/bin/python3")
library(reticulate)

reticulate::virtualenv_create("/Applications/anaconda3/bin/", python = "/Applications/anaconda3/bin/python3")

reticulate::use_virtualenv("/Applications/anaconda3/bin/")
py_available()


my_python_path= "/Applications/anaconda3/bin/python3"
library(Giotto)
instrs = createGiottoInstructions(python_path = my_python_path)

sc_matrix<-singlets@assays$RNA@counts
sc_lable<-singlets$seurat_clusters %>% as.factor()
sc_matrix2 <- createGiottoObject(raw_exprs = sc_matrix,instructions = instrs)
sc_matrix2norm <- normalizeGiotto(gobject = sc_matrix2)
sc_matrix2norm@cell_metadata$leiden_clus<-as.character(sc_lable)

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
head(colnames(Sig)); head(rownames(Sig))



AE <- AverageExpression(object =singlets, features = rownames(Sig), 
                        return.seurat = FALSE, group.by = "seurat_clusters", slot = "counts", verbose = TRUE)

corr <- cor(AE$RNA)
head(corr)
write_csv(cbind(rownames(as.data.frame(corr)), as.data.frame(corr)), "Heart_sc_correlation_fig1b.csv")


# Compute correlation P-value
p.mat <- corrplot::cor.mtest(mat = corr, conf.level = 0.95)

# Visualize
d <- ggcorrplot::ggcorrplot(
  corr = corr,
  p.mat = p.mat[[1]],
  hc.order = TRUE,
  sig.level = 1,
  type = "full",
  insig = "blank",
  lab = TRUE,
  outline.col = "white",
  method = "square",
  colors = c("#4477AA", "white", "#BB4444"),
  #colors = c("#6D9EC1", "white", "#E46726"),
  title = "Cell Type per Spot Correlation",
  legend.title = "Correlation\n(Pearson)") + ggplot2::labs(x = 'Cluster Number', y = 'Cluster Number') + 
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90, size = 18),
    axis.text.y = ggplot2::element_text(size = 18),
    axis.text = ggplot2::element_text(size = 18, vjust = 0.5))

d + scale_fill_gradient2(limit = c(0,1), low = "white", high =  "red", mid = "orange", midpoint = 0.5) + scale_x_continuous(breaks=seq(0,14,1)) +  scale_y_continuous(breaks=seq(0,14,1))



