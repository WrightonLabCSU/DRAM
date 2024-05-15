library(ggtree)
library(ggplot2)
library(ggrepel)

# Define the input files
newick_file <- commandArgs(trailingOnly = TRUE)[1]
labels_file <- commandArgs(trailingOnly = TRUE)[2]
output_pdf <- commandArgs(trailingOnly = TRUE)[3]

# Read the Newick tree
tree <- read.tree(newick_file)

# Read the labels to be colored
labels_to_color_raw <- readLines(labels_file)
labels_to_color <- gsub("\t", "", labels_to_color_raw)

# Create a data frame for the labels
label_data <- data.frame(label = tree$tip.label, color = "black", stringsAsFactors = FALSE)
label_data$color[label_data$label %in% labels_to_color] <- "red"

# Ensure no duplicates in the labels
label_data <- label_data[!duplicated(label_data$label),]

# Create a ggtree plot
p <- ggtree(tree, layout = "unrooted") +
  geom_tiplab(aes(label = label, color = color), data = label_data, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  geom_text_repel(aes(label = label, color = color), data = label_data, size = 2.5, show.legend = FALSE, nudge_x = 0.5, nudge_y = 0.5)

# Save the plot to a PDF file
ggsave(output_pdf, plot = p, width = 20, height = 20)

# Print message
cat("Colored tree saved to", output_pdf, "\n")
