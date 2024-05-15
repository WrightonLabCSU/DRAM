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

# Debug: Print the labels to be colored
cat("Labels to be colored:\n")
print(labels_to_color)

# Prepare the color vector
tip_colors <- rep("black", length(tree$tip.label))

# Debug: Print the labels found in the tree
cat("Valid labels found in the tree:\n")
print(tree$tip.label)

# Apply red color to the specified labels
for (i in seq_along(tree$tip.label)) {
  if (tree$tip.label[i] %in% labels_to_color) {
    tip_colors[i] <- "red"
  }
}

# Debug: Print the colored tips
cat("Tips to color (indices):\n")
print(which(tip_colors == "red"))

# Set plot size
pdf(output_pdf, width = 20, height = 20)

# Plot the tree without tip labels to adjust spacing
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Add tip labels with custom colors and larger font size
tiplabels(tree$tip.label, frame = "none", adj = c(1, 1), col = tip_colors, cex = 0.7)

# Close the PDF device
dev.off()
