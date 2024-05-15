# Load necessary libraries
library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\\t", " ", labels_to_color_raw)

# Create a vector of colors for the tip labels
tip_colors <- rep("black", length(tree$tip.label))
tip_colors[tree$tip.label %in% labels_to_color] <- "red"

# Print out information for debugging
cat("Tree tip labels:\n")
print(tree$tip.label)
cat("\nLabels to color:\n")
print(labels_to_color)
cat("\nTip colors:\n")
print(tip_colors)

# Increase the size of the PDF
pdf(output_pdf, width = 20, height = 20)  # Increase size to avoid compression
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Use tryCatch to handle errors during plotting
tryCatch({
    tiplabels(pch = 16, col = tip_colors, cex = 1.5)
    tiplabels(tree$tip.label, col = tip_colors, bg = tip_colors, cex = 0.7, adj = 1, offset = 0.5)
}, error = function(e) {
    cat("Error during plotting: ", e$message, "\n")
})

dev.off()
