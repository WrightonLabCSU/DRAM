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
labels_to_color <- sub(".*\t", "", labels_to_color_raw)  # Extract the second part after tab

# Plot the tree
pdf(output_pdf, width = 20, height = 20)  # Increase size for better visibility
plot(tree, type = "unrooted", show.tip.label = FALSE, cex = 0.6)

# Define colors
colors <- rep("black", length(tree$tip.label))
colors[tree$tip.label %in% labels_to_color] <- "red"

# Adjust the labels to radiate outwards
tiplabels(text = tree$tip.label, tip = 1:Ntip(tree), 
          col = colors, cex = 0.6, adj = 1, srt = 90)

# Use ape's tiplabels with angle adjustment
for (i in 1:Ntip(tree)) {
  angle <- ifelse(tree$edge.length[i] > 0, tree$edge.length[i], 1)
  angle <- ifelse(i <= Ntip(tree) / 2, angle * 180 / pi, (angle + pi) * 180 / pi)
  text(tree$tip.label[i], pos = angle, cex = 0.6, col = colors[i])
}

dev.off()
