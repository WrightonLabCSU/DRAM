# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(ggrepel)

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

# Create a data frame to associate labels with their colors
label_colors <- data.frame(label = tree$tip.label, color = "black")
label_colors$color[label_colors$label %in% labels_to_color] <- "red"

# Convert the tree to a ggtree object
p <- ggtree(tree, layout = "daylight") +
  geom_tiplab(aes(label = label), size = 2, align = TRUE, linesize = 0.5) +
  geom_tiplab(data = label_colors, aes(label = label, color = I(color)), size = 2, align = TRUE, linesize = 0.5) +
  geom_text_repel(data = label_colors[label_colors$color == "red", ], aes(label = label), color = "red", size = 2, max.overlaps = 20)

# Save the plot to a PDF file with increased size
ggsave(output_pdf, plot = p, width = 20, height = 20)
