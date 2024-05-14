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

# Print the labels to be colored for debugging
cat("Labels to be colored:\n")
print(labels_to_color)

# Ensure labels to color are present in the tree
valid_labels <- tree$tip.label %in% labels_to_color

# Print the valid labels for debugging
cat("Valid labels found in the tree:\n")
print(tree$tip.label[valid_labels])

# Plot the unrooted tree
png(output_png, width = 1024, height = 1024)
plot(tree, type = "unrooted", cex = 0.6)

# Color the specified labels
colored_tips <- which(tree$tip.label %in% labels_to_color)
tiplabels(pch = 19, tip = colored_tips, col = "red", cex = 1)

# Add text labels in red for the specified tips
for (i in colored_tips) {
  tiplabels(tree$tip.label[i], tip = i, frame = "none", col = "red", adj = c(1, 0.5), cex = 0.6)
}

dev.off()
