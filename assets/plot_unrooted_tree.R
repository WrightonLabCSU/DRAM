# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(ggrepel)

# Read the Newick tree file
tree_file <- "colored_tree.nwk"
tree <- read.tree(tree_file)

# Read the labels to be colored
labels_to_color <- scan("extracted_query_ids.txt", what = "", sep = "\n")
labels_to_color <- gsub("\t", "_", labels_to_color)  # Replace tabs with underscores

# Create a data frame for labels
label_data <- data.frame(label = tree$tip.label, color = "black")
label_data$color[label_data$label %in% labels_to_color] <- "red"

# Plot the tree using ggtree
p <- ggtree(tree, layout = "unrooted") +
  geom_tiplab(aes(label = label, color = color), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  theme_tree2()

# Add labels using geom_text_repel to avoid overlapping
p <- p + geom_text_repel(
  data = label_data[label_data$color == "red", ],
  aes(label = label, color = color),
  size = 3,
  nudge_x = 0.5,
  nudge_y = 0.5,
  max.overlaps = Inf
)

# Save the plot as a PDF
ggsave("colored_tree.pdf", plot = p, width = 20, height = 20)
