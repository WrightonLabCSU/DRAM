# Load required libraries
library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_png <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color <- readLines(labels_file)

# Plot the unrooted tree
png(output_png, width = 1024, height = 1024)
plot(tree, type = "unrooted", cex = 0.6)

# Color the specified labels
tiplabels(pch = 19, tip = which(tree$tip.label %in% labels_to_color), col = "red", cex = 1)
tiplabels(tree$tip.label, tip = which(tree$tip.label %in% labels_to_color), frame = "none", col = "red", adj = c(1,0.5), cex = 0.6)

dev.off()
