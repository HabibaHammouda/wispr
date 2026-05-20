install.packages(c("Seurat", "ggplot2", "dplyr", "Matrix", "tidyverse", "spatstat", "devtools"))

install.packages(c("hydroGOF"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Giotto")

# List of all the packages you tried to install
packages_to_check <- c(
  "BiocManager", "Giotto", "Seurat", "ggplot2", "dplyr", 
  "Matrix", "tidyverse", "reshape2", "S4Vectors", "png", 
  "RColorBrewer", "data.table", "reticulate", "hydroGOF"
)

# Check which packages are successfully installed
status <- sapply(packages_to_check, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

# Create a clean summary table
verification_results <- data.frame(
  Package = names(status),
  Installed_Successfully = ifelse(status, "✅ YES", "❌ NO - Failed to install")
)

print(verification_results, row.names = FALSE)

# 1. Define your list of target packages
packages <- c(
  "Seurat", "ggplot2", "dplyr", "Matrix", "tidyverse", 
  "reshape2", "png", "RColorBrewer", "data.table", "reticulate"
)

# 2. Check which packages are missing
missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

# 3. Print a summary report
if (length(missing_pkgs) == 0) {
  message("✨ Success! All packages are fully installed and ready to go.")
} else {
  warning("⚠️ The following packages failed to download or are missing:\n", 
          paste("- ", missing_pkgs, collapse = "\n"))
  message("\n💡 Tip: Try running install.packages(c(", 
          paste0('"', missing_pkgs, '"', collapse = ", "), ")) again.")
}