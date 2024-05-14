library(ape)
library(ggtree)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
nwk_file <- args[1]
labels_file <- args[2]
output_file <- args[3]

# Read the Newick file
tree <- read.tree(nwk_file)

# Read the labels to be colored
labels <- readLines(labels_file)

# Prepare a data frame for labels to be colored
colored_labels <- data.frame(label = labels, color = "red", stringsAsFactors = FALSE)

# Plot the tree
p <- ggtree(tree, layout = "unrooted") +
  geom_tiplab(aes(label=label), color="black", size=2, offset=0.02) +
  geom_tiplab(data=colored_labels, aes(label=label), color="red", size=2, offset=0.02) +
  theme_tree2()

# Save to PDF
ggsave(output_file, p, width = 16, height = 16, units = "in")
