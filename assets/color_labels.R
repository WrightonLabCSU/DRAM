library(ape)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)

# Clean and prepare the labels for matching
labels_to_color <- gsub("\t", "_", labels_to_color_raw)

# Create a vector for label colors, default to black
label_colors <- rep("black", length(tree$tip.label))

# Color the specified labels in red
label_colors[tree$tip.label %in% labels_to_color] <- "red"

# Calculate the angles for the labels
n_tips <- length(tree$tip.label)
angles <- 360 * (1:n_tips) / n_tips
angles <- ifelse(angles > 180, angles - 180, angles)

# Plot the tree
pdf(output_pdf, width = 10, height = 10)
plot(tree, type = "unrooted", show.tip.label = FALSE)

# Add labels with colors and angles
for (i in 1:n_tips) {
  label <- tree$tip.label[i]
  color <- label_colors[i]
  angle <- angles[i]
  text(tree$tip[i, 1], tree$tip[i, 2], labels = label, pos = 4, srt = angle, cex = 0.6, col = color)
}

dev.off()
