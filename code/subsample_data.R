library(Seurat)
library(dplyr)
library(readr)

# Load your Visium data (adjust path to your actual data location)
data_dir <- "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data"
expression_matrix <- Read10X(data.dir = data_dir)

# Load coordinates
coor_ST <- read.csv(paste0(data_dir, "/tissue_positions_list.csv"), header = FALSE)
rownames(coor_ST) <- coor_ST$V1
coor_ST <- coor_ST[coor_ST$V2 == "1", ]  # Only spots under tissue

# Match spots
match_n <- match(colnames(expression_matrix), rownames(coor_ST))
coor_ST_reor <- coor_ST[match_n, ]

# SUBSAMPLE TO 800 SPOTS
set.seed(42)  # Reproducibility
n_spots <- min(800, ncol(expression_matrix))
selected_spots <- sample(1:ncol(expression_matrix), n_spots)

expression_matrix_sub <- expression_matrix[, selected_spots]
coor_ST_sub <- coor_ST_reor[selected_spots, ]

# Create Seurat object
st_obj <- CreateSeuratObject(counts = expression_matrix_sub, assay = "RNA")
st_obj <- AddMetaData(st_obj, metadata = coor_ST_sub$V5, col.name = "X_coor")
st_obj <- AddMetaData(st_obj, metadata = coor_ST_sub$V6, col.name = "Y_coor")

# Save
saveRDS(st_obj, "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/st_800spots.rds")
write.csv(coor_ST_sub, 
          "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/coor_ST_800.csv")

cat("Subsampled to", n_spots, "spots\n")