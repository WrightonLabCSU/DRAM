# Load necessary libraries
library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\\t", " ", labels_to_color_raw)

# Create a vector of colors for the tip labels
tip_colors <- rep("black", length(tree$tip.label))
tip_colors[tree$tip.label %in% labels_to_color] <- "red"

# Create a new tree object without duplicate labels for the colored tips
unique_labels <- unique(c(tree$tip.label, labels_to_color))
tree$tip.label <- unique_labels

# Plot the tree with the correct colors and avoid overlapping text
pdf(output_pdf, width = 15, height = 15)  # Increase size to avoid compression
plot(tree, type = "unrooted", show.tip.label = FALSE)
tiplabels(pch = 16, col = tip_colors, cex = 1.5)
tiplabels(tree$tip.label, col = tip_colors, bg = tip_colors, cex = 0.7, adj = 1, offset = 0.5)
dev.off()
