library(ape)

# Define the input files
args <- commandArgs(trailingOnly = TRUE)
newick_file <- args[1]
labels_file <- args[2]
output_pdf <- args[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color <- readLines(labels_file)

# Remove tab spaces from labels to color
labels_to_color <- sub("\t", " ", labels_to_color)

# Set up colors
tip_colors <- rep("black", length(tree$tip.label))
names(tip_colors) <- tree$tip.label

# Color the tips that match the labels
for (label in labels_to_color) {
  if (label %in% names(tip_colors)) {
    tip_colors[label] <- "red"
  }
}

# Increase the plot size
pdf(output_pdf, width = 15, height = 15)
plot(tree, type = "unrooted", show.tip.label = FALSE, cex = 0.8)

# Add colored labels
tiplabels(tree$tip.label, tip = 1:length(tree$tip.label), bg = tip_colors, frame = "none", col = tip_colors, cex = 0.8, adj = 1)

dev.off()
