# Load necessary libraries
library(ape)
library(phytools)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\t", "_", labels_to_color_raw)

# Define the color for the specified labels
color <- "red"

# Adjust the plotting parameters for a larger plot
pdf(output_pdf, width = 20, height = 20)

# Plot the tree with radiating tips
plot(tree, type = "unrooted", show.tip.label = FALSE, cex = 0.8, no.margin = TRUE, edge.width = 2)

# Add tip labels with radiating text
tiplabels(
  text = tree$tip.label,
  adj = c(0.5, 0.5),
  frame = "none",
  pos = 4,
  offset = 0.5,
  cex = 0.8,
  font = 2,
  col = ifelse(tree$tip.label %in% labels_to_color, color, "black")
)

# Save the plot to a file
dev.off()

cat("Plot saved to", output_pdf, "\n")
