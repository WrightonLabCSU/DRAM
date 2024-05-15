library(ape)

args <- commandArgs(trailingOnly = TRUE)
newick_file <- args[1]
labels_file <- args[2]
output_pdf <- args[3]

tree <- read.tree(newick_file)
labels_to_color_raw <- readLines(labels_file)

labels_to_color <- gsub("\t", " ", labels_to_color_raw)
labels_to_color <- trimws(labels_to_color)

tip_colors <- rep("black", length(tree$tip.label))
tip_labels <- tree$tip.label

for (i in seq_along(tree$tip.label)) {
  for (label in labels_to_color) {
    if (grepl(label, tree$tip.label[i], fixed = TRUE)) {
      tip_colors[i] <- "red"
      tip_labels[i] <- tree$tip.label[i]
      break
    }
  }
}

pdf(output_pdf, width = 16, height = 16)
plot(tree, type = "unrooted", show.tip.label = FALSE)

offset <- 0.02
cex <- 0.5
for (i in seq_along(tree$tip.label)) {
  tiplabels(tip_labels[i], adj = c(0.5, 0.5), frame = "none", col = tip_colors[i], cex = cex, offset = offset)
}

dev.off()
