#library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(tidyverse)
library(reshape2)
library(S4Vectors)
library(png)
library(RColorBrewer)

counts1 <- read.csv("counts.dat", quote = "", sep = "\t") #read filtered sc data

rownames(counts1) = counts1$cell
counts1 <- dplyr::select(counts1, -1)

counts1 <- t(counts1)

counts_meta1 <- read.csv("count_meta.csv", quote = "", sep = "\t") #read filtered meta data
rownames(counts_meta1) = counts_meta1$cell
counts_meta1 <- dplyr::select(counts_meta1, -1)


counts.sp1 <- as(as.matrix(counts1), "dgCMatrix") #Convert to sparse matrix to save space
counts.metasc1 <-CreateSeuratObject(counts = counts.sp1)

counts.metasc1 <- AddMetaData(object = counts.metasc1, metadata = counts_meta1$res.0.7, col.name = "RNA_snn_res.0.7_meta")


lower <- stringr::str_to_title(rownames(counts.metasc1))


counts.metasc2 <- counts.metasc1@assays$RNA@counts
rownames(counts.metasc2) <- lower
counts.metasc3 <-CreateSeuratObject(counts = counts.metasc2, assay = "RNA")

counts.metasc3@meta.data$seurat_clusters <- counts_meta1$origident
counts.metasc3 <- SetIdent(counts.metasc3, value = "seurat_clusters")
Idents(counts.metasc3) <- factor(x = Idents(counts.metasc3), levels = sort(levels(counts.metasc3)))


###scRNA SEQ Brain
require(data.table)
counts_meta <- read.table(paste0("sc_mta_data.tsv")) #read filtered meta data
rownames(counts_meta) = counts_meta$V1
counts_meta <- dplyr::select(counts_meta, -1)


counts.metasc <- readRDS("seurat_counts.rds")

counts2 <- GetAssayData(counts.metasc, assay = "RNA")

counts2 <- counts2[-c(1),]


counts.metasc <-CreateSeuratObject(counts = counts2)

counts.metasc <- AddMetaData(object = counts.metasc, metadata = counts_meta$V2, col.name = "RNA_snn_res.0.7_meta")

counts.metasc@meta.data$seurat_clusters <- counts_meta$V2
counts.metasc@meta.data$seurat_clusters

counts.metasc <- SetIdent(counts.metasc, value = "seurat_clusters")
Idents(counts.metasc) <- factor(x = Idents(counts.metasc), levels = sort(levels(counts.metasc)))


sub_deneme <- c("Astrocytes_14", "Astrocytes_40", 
                "Astrocytes_41", "Blood_73", "Ependymal_47", "Excluded_38", "Immune_32", 
                "Immune_34", "Immune_35",  "Neurons_25",
                "Neurons_26","Neurons_27", "Neurons_63",
                "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                "Vascular_67", "Vascular_69")

sub_remain_deneme <- c("Astrocytes_42", "Oligos_5","Neurons_48",
                       "Oligos_53","Vascular_68",
                       "Neurons_12", "Neurons_14",
                       "Neurons_21","Neurons_52","Neurons_51",
                       "Neurons_18","Neurons_15",
                       "Neurons_11","Neurons_24",
                       "Neurons_23")

sub_obj <- subset(counts.metasc, idents = sub_deneme)

sub_obj2 <- t(sub_obj@assays$RNA@counts)

write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_0_50/subset_2923cells_brain_sc.dat")
write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_10_40/subset_2923cells_brain_sc.dat")
write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_20_30/subset_2923cells_brain_sc.dat")
write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_30_20/subset_2923cells_brain_sc.dat")
write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_40_10/subset_2923cells_brain_sc.dat")
write_csv(cbind(rownames(sub_obj2), as.data.frame(sub_obj2)), "mixture_50_50_0/subset_2923cells_brain_sc.dat")


write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_0_50/subset_2923cells_brain_sc_meta.dat")
write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_10_40/subset_2923cells_brain_sc_meta.dat")
write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_20_30/subset_2923cells_brain_sc_meta.dat")
write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_30_20/subset_2923cells_brain_sc_meta.dat")
write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_40_10/subset_2923cells_brain_sc_meta.dat")
write_csv(cbind(colnames(sub_obj), as.data.frame(sub_obj$seurat_clusters)), "mixture_50_50_0/subset_2923cells_brain_sc_meta.dat")


#random selection of 50_0_50 of heart cell

counts_hc <- counts.metasc3@assays$RNA@counts

library(dplyr)

random_2923_cols <- counts_hc[, sample(ncol(counts_hc), 2923)]

selection <- colnames(random_2923_cols)

counts.metasc4 <- counts.metasc3[,selection]


inter <- intersect(rownames(counts.metasc4), rownames(sub_obj))

counts.metasc5 <- counts.metasc4[inter,]
sub_obj3 <- sub_obj[inter,]
table(counts.metasc5$seurat_clusters)

merged_object <- merge(x = counts.metasc5, y = sub_obj3, add.cell.ids = c("Human", "Mouse"))
write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_0_50/Brain_nc50_brainunc0_heart50_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_0_50/Brain_nc50_brainunc0_heart50_reference_labels.dat")
#random selection of 50_10_40 of heart cell

counts_hc <- counts.metasc3@assays$RNA@counts
random_2338_cols <- counts_hc[, sample(ncol(counts_hc), 2338)]
selection <- colnames(random_2338_cols)

counts.metasc4 <- counts.metasc3[,selection]

sub_obj_10 <- subset(counts.metasc, idents = sub_remain_deneme)

brain_585_cols <- sub_obj_10[, sample(ncol(sub_obj_10), 585)]

inter <- intersect(rownames(counts.metasc4), rownames(sub_obj))
counts.metasc5 <- counts.metasc4[inter,]
sub_obj3 <- sub_obj[inter,]
brain_585_cols_2 <- brain_585_cols[inter,]

merged_object <- merge(x = counts.metasc5, y = c(sub_obj3, brain_585_cols_2), add.cell.ids = c("Human", "Mouse","Mouse"))
table(merged_object$seurat_clusters)

write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_10_40/Brain_nc50_brainunc10_heart40_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_10_40/Brain_nc50_brainunc10_heart40_reference_labels.dat")

#random selection of 50_20_30 of heart cell

counts_bc <- counts.metasc3@assays$RNA@counts


random_1754_cols <- counts_bc[, sample(ncol(counts_bc), 1754)]
selection <- colnames(random_1754_cols)

counts.metasc4 <- counts.metasc3[,selection]
sub_obj_10 <- subset(counts.metasc, idents = sub_remain_deneme)

brain_1169_cols <- sub_obj_10[, sample(ncol(sub_obj_10), 1169)]

inter <- intersect(rownames(counts.metasc4), rownames(sub_obj))
counts.metasc5 <- counts.metasc4[inter,]
sub_obj3 <- sub_obj[inter,]
brain_1169_cols_2 <- brain_1169_cols[inter,]
dim(counts.metasc5); dim(sub_obj3); dim(brain_1169_cols_2)

merged_object <- merge(x = counts.metasc5, y = c(sub_obj3, brain_1169_cols_2), add.cell.ids = c("Human", "Mouse","Mouse"))

write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_20_30/Brain_nc50_brainunc20_heart30_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_20_30/Brain_nc50_brainunc20_heart30_reference_labels.dat")

#random selection of 50_30_20 of heart cell

counts_hc <- counts.metasc3@assays$RNA@counts
random_1169_cols <- counts_hc[, sample(ncol(counts_hc), 1169)]
head(colnames(random_1169_cols))
selection <- colnames(random_1169_cols)

counts.metasc4 <- counts.metasc3[,selection]


sub_obj_10 <- subset(counts.metasc, idents = sub_remain_deneme)

brain_1754_cols <- sub_obj_10[, sample(ncol(sub_obj_10), 1754)]

inter <- intersect(rownames(counts.metasc4), rownames(sub_obj))
counts.metasc5 <- counts.metasc4[inter,]
sub_obj3 <- sub_obj[inter,]
brain_1754_cols_2 <- brain_1754_cols[inter,]

merged_object <- merge(x = counts.metasc5, y = c(sub_obj3, brain_1754_cols_2), add.cell.ids = c("Human", "Mouse","Mouse"))
table(merged_object$seurat_clusters)

write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_30_20/Brain_nc50_brainunc30_heart20_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_30_20/Brain_nc50_brainunc30_heart20_reference_labels.dat")

#random selection of 50_40_10 of heart cell
idents2 <- c("0","1", "2", "3", "4", "5", "7")

counts_hc <- subset(counts.metasc3, idents = idents2)
table(counts_hc$seurat_clusters)

counts_hc2 <- counts_hc@assays$RNA@counts
random_585_cols <- counts_hc2[, sample(ncol(counts_hc2), 585)]
selection <- colnames(random_585_cols)

counts.metasc4 <- counts_hc[,selection]


sub_obj_10 <- subset(counts.metasc, idents = sub_remain_deneme)

brain_2338_cols <- sub_obj_10[, sample(ncol(sub_obj_10), 2338)]

inter <- intersect(rownames(counts.metasc4), rownames(sub_obj))
counts.metasc5 <- counts.metasc4[inter,]
sub_obj3 <- sub_obj[inter,]
brain_2338_cols_2 <- brain_2338_cols[inter,]

merged_object <- merge(x = counts.metasc5, y = c(sub_obj3, brain_2338_cols_2), add.cell.ids = c("Human", "Mouse","Mouse"))
table(merged_object$seurat_clusters)

write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_40_10/Brain_nc50_brainunc40_heart10_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_40_10/Brain_nc50_brainunc40_heart10_reference_labels.dat")

#random selection of 50_50_0 of heart cell

brain_2923_cols <- sub_obj_10[, sample(ncol(sub_obj_10), 2923)]

merged_object <- merge(x = brain_2923_cols, y = sub_obj, add.cell.ids = c("Mouse", "Mouse"))
table(merged_object$seurat_clusters)

write_csv(cbind(rownames(merged_object), as.data.frame(merged_object@assays$RNA@counts)), "mixture_50_50_0/Brain_nc50_brainunc50_heart0_reference.dat")

write_csv(cbind(colnames(merged_object), as.data.frame(merged_object$seurat_clusters)), "mixture_50_50_0/Brain_nc50_brainunc50_heart0_reference_labels.dat")
