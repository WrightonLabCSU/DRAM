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
labels_to_color <- sapply(labels_to_color_raw, function(x) strsplit(x, "\t")[[1]][2])

# Function to color the labels
color_labels <- function(tree, labels, color) {
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips <- tree$tip.label[tree$edge[is_tip, 2]]
  col <- rep("black", length(ordered_tips))
  col[ordered_tips %in% labels] <- color
  return(col)
}

# Apply coloring to the labels
tip_colors <- color_labels(tree, labels_to_color, "red")

# Plot the tree with colored labels
pdf(output_pdf, width = 16, height = 16)
plot(tree, type = "unrooted", show.tip.label = TRUE, tip.color = tip_colors, cex = 0.5)
tiplabels(tree$tip.label[tree$tip.label %in% labels_to_color], tip = which(tree$tip.label %in% labels_to_color), col = "red", frame = "none", cex = 0.5)
dev.off()

cat("Tree with colored labels saved to", output_pdf)
