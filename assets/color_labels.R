# Load necessary libraries
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
labels_to_color <- readLines(labels_file)

# Prepare the colors for the tip labels
tip_colors <- rep("black", length(tree$tip.label))
tip_colors[tree$tip.label %in% labels_to_color] <- "red"

# Open a PDF device for plotting
pdf(output_pdf, width = 20, height = 20)

# Plot the tree without labels
plot(tree, type = "unrooted", show.tip.label = FALSE, edge.width = 1)

# Add colored tip labels
for (i in 1:length(tree$tip.label)) {
    tip_label <- tree$tip.label[i]
    tip_angle <- 360 * i / length(tree$tip.label)
    angle <- ifelse(tip_angle > 180, tip_angle - 180, tip_angle)
    if (tip_label %in% labels_to_color) {
        text(x = tree$tip[i, 1], y = tree$tip[i, 2], labels = tip_label, col = "red", srt = angle, adj = 0.5)
    } else {
        text(x = tree$tip[i, 1], y = tree$tip[i, 2], labels = tip_label, col = "black", srt = angle, adj = 0.5)
    }
}

# Close the PDF device
dev.off()
