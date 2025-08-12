# ---
# Usage:       Rscript cutree_expression.R <input_csv_path> <output_dir>
# ---

library(WGCNA)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript cutree_expression.R <input_csv_path> <output_dir>", call. = FALSE)
}

input_csv_path <- args[1]
output_dir <- args[2]

min_cluster_size <- 12
deep_split_level <- 1
hclust_method <- "ward.D2"

if (!file.exists(input_csv_path)) {
  stop("Error: Input file not found at the provided path.")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dissTOM <- read.csv(input_csv_path, row.names = 1, check.names = FALSE)
geneTree <- hclust(as.dist(dissTOM), method = hclust_method)

pdf_path <- file.path(output_dir, "gene_dendrogram_expression_wardD2.pdf")
pdf(pdf_path, width = 8, height = 10)
plot(geneTree, main = "Gene Cluster Dendrogram (Expression)", xlab = "", sub = "")
dev.off()

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = as.matrix(dissTOM),
  deepSplit = deep_split_level,
  pamRespectsDendro = FALSE,
  minClusterSize = min_cluster_size
)

module_df <- data.frame(
  Gene = colnames(dissTOM),
  Module = dynamicMods
)

csv_output_path <- file.path(output_dir, "expression_modules.csv")
write.csv(module_df, file = csv_output_path, row.names = FALSE)
