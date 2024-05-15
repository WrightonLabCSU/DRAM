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

# Plot the tree and extract coordinates
pdf(output_pdf, width = 10, height = 10)
plot(tree, type = "unrooted", show.tip.label = FALSE)
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Add labels with colors and angles
for (i in 1:length(tree$tip.label)) {
  label <- tree$tip.label[i]
  color <- label_colors[i]
  angle <- last_plot$theta[i]
  if (!is.na(angle) && angle > 90 && angle < 270) {
    angle <- angle + 180
  }
  text(last_plot$xx[i], last_plot$yy[i], labels = label, pos = 4, cex = 0.6, col = color, srt = angle)
}

dev.off()
