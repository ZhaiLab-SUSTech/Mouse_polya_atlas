# install_packages.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("tidyverse")
bioc_packages <- c("WGCNA")

cat("--- Installing CRAN packages ---\n")
install.packages(cran_packages, repos = "https://cloud.r-project.org/")

cat("\n--- Installing Bioconductor packages ---\n")
BiocManager::install(bioc_packages)

cat("\n--- All packages installed successfully! ---\n")