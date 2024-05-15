#!/usr/bin/env Rscript

# Load necessary libraries
library(ape)
library(ggplot2)
library(ggtree)

# Define the input files
args <- commandArgs(trailingOnly = TRUE)
newick_file <- args[1]
labels_file <- args[2]
output_pdf <- args[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\t", " ", labels_to_color_raw)

# Check if the labels to color are in the tree
valid_labels <- intersect(tree$tip.label, labels_to_color)
cat("Labels to be colored:\n", labels_to_color, "\n")
cat("Valid labels found in the tree:\n", valid_labels, "\n")

# Create a vector for tip labels colors
tip_colors <- rep("black", length(tree$tip.label))
names(tip_colors) <- tree$tip.label

# Set the color for the labels we want to color
tip_colors[valid_labels] <- "red"

# Create a data frame for plotting
tip_labels <- data.frame(label = tree$tip.label, color = tip_colors)

# Plot the tree using ggtree
p <- ggtree(tree, layout = "unrooted") +
  geom_tiplab(aes(label = label, color = color), data = tip_labels, size = 2, offset = 0.05, align = TRUE, linetype = "dotted") +
  scale_color_identity() +
  theme_tree2()

# Save the plot
ggsave(output_pdf, plot = p, height = 30, width = 30, units = "in")

cat("Colored tree saved to", output_pdf, "\n")
