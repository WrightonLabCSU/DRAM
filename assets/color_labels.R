library(ape)
library(phytools)

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

# Function to adjust label positions and avoid overlap
tiplabels_radiate <- function(tree, labels_to_color, colors, cex = 0.6) {
  n_tips <- Ntip(tree)
  angles <- 360 * (1:n_tips / n_tips)
  radians <- angles * pi / 180
  
  for (i in 1:n_tips) {
    x <- cos(radians[i])
    y <- sin(radians[i])
    label <- tree$tip.label[i]
    color <- colors[i]
    
    text(x, y, labels = label, col = color, cex = cex, srt = angles[i], adj = ifelse(angles[i] > 180, 1, 0))
  }
}

# Adjust the labels to radiate outwards
tiplabels_radiate(tree, labels_to_color, colors, cex = 0.6)

dev.off()
