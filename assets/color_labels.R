library(ape)

# Define the input files
args <- commandArgs(trailingOnly = TRUE)
newick_file <- args[1]
labels_file <- args[2]
output_pdf <- args[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- sub(".*\t", "", labels_to_color_raw)

# Ensure unique labels
labels_to_color <- unique(labels_to_color)

# Open PDF device with increased size
pdf(output_pdf, width = 15, height = 15)

# Plot the tree without labels
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Add colored tip labels without duplicating
tiplabels(pch = 19, col = "black", adj = c(1, 1))

for (i in seq_along(tree$tip.label)) {
  if (tree$tip.label[i] %in% labels_to_color) {
    tiplabels(tree$tip.label[i], i, col = "red", adj = c(1, 1), cex = 0.6, font = 2)
  } else {
    tiplabels(tree$tip.label[i], i, col = "black", adj = c(1, 1), cex = 0.6, font = 1)
  }
}

# Close the PDF device
dev.off()
