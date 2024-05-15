# Load required libraries
library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)

# Process the labels to extract query IDs
labels_to_color <- sapply(labels_to_color_raw, function(x) {
  parts <- unlist(strsplit(x, "\t"))
  return(parts[2])  # Consider only the second part (query ID)
})

# Print the labels to be colored for debugging
cat("Processed Labels to be colored:\n")
print(labels_to_color)

# Print the tree tip labels for debugging
cat("Tree tip labels:\n")
print(tree$tip.label)

# Ensure labels to color are present in the tree
valid_labels <- tree$tip.label %in% labels_to_color

# Print the valid labels for debugging
cat("Valid labels found in the tree:\n")
print(tree$tip.label[valid_labels])

# Plot the unrooted tree with larger size
pdf(output_pdf, width = 30, height = 30)  # Export to PDF with larger size
plot(tree, type = "unrooted", cex = 0.8)

# Color the specified labels
colored_tips <- which(tree$tip.label %in% labels_to_color)

# Add text labels in red for the specified tips without duplicating
for (i in colored_tips) {
  # Apply jitter to avoid overlapping
  x_jitter <- runif(1, -0.02, 0.02)
  y_jitter <- runif(1, -0.02, 0.02)
  tiplabels(pch = 19, tip = i, col = "red", cex = 1)
  x <- tree$edge.length[i] + x_jitter
  y <- i + y_jitter
  text(x, y, tree$tip.label[i], col = "red", pos = 4, cex = 0.8)
}

dev.off()
