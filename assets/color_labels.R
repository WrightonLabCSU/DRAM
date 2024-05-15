library(ape)

# Define the input files
args <- commandArgs(trailingOnly = TRUE)
newick_file <- args[1]
labels_file <- args[2]
output_pdf <- args[3]

# Debug: Print the input arguments
cat("Newick file:", newick_file, "\n")
cat("Labels file:", labels_file, "\n")
cat("Output PDF:", output_pdf, "\n")

# Read the Newick tree
cat("Reading the Newick tree...\n")
tree <- read.tree(newick_file)

# Debug: Print the number of tips in the tree
cat("Number of tips in the tree:", length(tree$tip.label), "\n")

# Read the labels to be colored
cat("Reading the labels to be colored...\n")
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\t", "", labels_to_color_raw)

# Debug: Print the labels to be colored
cat("Labels to be colored:\n")
print(labels_to_color)

# Find and match the labels
cat("Matching the labels...\n")
matched_labels <- labels_to_color[labels_to_color %in% tree$tip.label]

# Debug: Print the matched labels
cat("Matched labels:\n")
print(matched_labels)

# Set the colors for the labels
label_colors <- ifelse(tree$tip.label %in% matched_labels, "red", "black")

# Debug: Print the label colors
cat("Label colors set.\n")

# Plot the tree
cat("Plotting the tree...\n")
pdf(output_pdf, width = 20, height = 20)
plot(tree, type = "unrooted", show.tip.label = TRUE, cex = 0.6, label.offset = 0.1)
tiplabels(tree$tip.label, adj = 1, frame = "none", col = label_colors, cex = 0.6)
dev.off()

cat("Tree plotting complete. PDF saved as:", output_pdf, "\n")
