library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- sub("\t", " ", labels_to_color_raw)

# Define the labels present in the tree
valid_labels <- tree$tip.label

# Identify the indices of the tips to color
tips_to_color <- which(valid_labels %in% labels_to_color)

# Plot the tree with colored tips
pdf(output_pdf, width=20, height=20)
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Color the tips
tiplabels(pch = 16, col = "red", tip = tips_to_color)
tiplabels(text = valid_labels[tips_to_color], tip = tips_to_color, col = "red", adj = 0, cex = 0.5, offset = 0.5, frame = "none", srt = 45)

# Add non-colored tip labels
non_colored_tips <- setdiff(seq_along(valid_labels), tips_to_color)
tiplabels(text = valid_labels[non_colored_tips], tip = non_colored_tips, col = "black", adj = 0, cex = 0.5, offset = 0.5, frame = "none", srt = 45)

dev.off()
