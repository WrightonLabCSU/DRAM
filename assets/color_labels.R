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
labels_to_color <- gsub("\t", "_", labels_to_color_raw) # Replace tabs with underscores

# Prepare the color vector
tip_colors <- rep("black", length(tree$tip.label))

# Apply red color to the specified labels
for (i in seq_along(tree$tip.label)) {
  if (tree$tip.label[i] %in% labels_to_color) {
    tip_colors[i] <- "red"
  }
}

# Set plot size
pdf(output_pdf, width = 30, height = 30) # Increase the size of the plot

# Plot the tree without tip labels to adjust spacing
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Add tip labels with custom colors and larger font size
tiplabels(tree$tip.label, frame = "none", adj = c(1, 1), col = tip_colors, cex = 0.7)

# Close the PDF device
dev.off()
