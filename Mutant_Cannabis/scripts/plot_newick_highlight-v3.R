#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))

# ---------------------------
# Argument parsing
# ---------------------------
parser <- ArgumentParser(
  description = "Render a Newick tree with highlighted tips and zoomed-in subtree"
)

parser$add_argument("-n", "--newick", required = TRUE, help = "Input Newick file")
parser$add_argument("-o", "--output_prefix", default = "tree", help = "Output file prefix")
parser$add_argument("-H", "--highlight", nargs = "+", help = "Sample IDs to highlight")
parser$add_argument("--layout", default = "circular", choices = c("circular","rectangular"), help="Tree layout")
parser$add_argument("--size", type="double", default = 0.2, help="Branch linewidth")
parser$add_argument("--full_width", type="double", default=10, help="Full tree width in inches")
parser$add_argument("--full_height", type="double", default=10, help="Full tree height in inches")
parser$add_argument("--zoom_width", type="double", default=8, help="Zoomed tree width")
parser$add_argument("--zoom_height", type="double", default=8, help="Zoomed tree height")
parser$add_argument("--combine", action="store_true", help="Save combined full+zoom plot")

args <- parser$parse_args()

# ---------------------------
# Read tree
# ---------------------------
phylo <- read.tree(args$newick)
if (is.null(phylo)) stop("Failed to read Newick file")

# ---------------------------
# Highlight data
# ---------------------------
tip_data <- NULL
if (!is.null(args$highlight)) {
  parsed <- lapply(args$highlight, function(x) {
    if (grepl(":", x)) {
      parts <- strsplit(x, ":", fixed = TRUE)[[1]]
      data.frame(
        label = parts[1],
        #display = paste(parts[1], parts[2], sep = " - "),
        #display = paste(parts[1], parts[2], sep = " - "),
        display = paste(
          parts[1],
          gsub('^"|"$', "", parts[2]),
          sep = " - "
        ),

        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        label = x,
        display = x,
        stringsAsFactors = FALSE
      )
    }
  })

  tip_data <- do.call(rbind, parsed)
  #tip_data <- data.frame(label = args$highlight)
}

# ---------------------------
# Full tree plot
# ---------------------------
p_full <- ggtree(phylo, layout=args$layout, linewidth=args$size)
p_full <- p_full +
  theme(plot.margin = margin(5, 0, 5, 5)) #trying to get the two plots together


if (!is.null(tip_data)) {
  tip_coords <- p_full$data %>% filter(label %in% tip_data$label)
  p_full <- p_full + geom_tippoint(data=tip_coords, aes(x=x, y=y), color="red", size=2)
}

full_file <- paste0(args$output_prefix, "_full.png")
ggsave(filename = full_file, plot=p_full, width=args$full_width, height=args$full_height, dpi=300)
message("Saved full tree: ", full_file)

# ---------------------------
# Zoomed-in subtree
# ---------------------------
zoom_file <- NULL
if (!is.null(tip_data)) {
  #node <- getMRCA(phylo, args$highlight)
  node <- getMRCA(phylo, tip_data$label)

  
  if (!is.na(node)) {
    # Extract subtree around MRCA
    subtree <- extract.clade(phylo, node = node)
    
    p_zoom <- ggtree(subtree, layout = args$layout, linewidth = args$size)
    p_zoom <- p_zoom +
      theme(plot.margin = margin(5, 5, 5, 0)) #trying to get two plots close together

    # Highlight tips in the subtree
    #tip_coords_zoom <- p_zoom$data %>% filter(label %in% tip_data$label)
    tip_coords_zoom <- p_zoom$data %>%
      filter(label %in% tip_data$label) %>%
      left_join(tip_data, by = "label")
    p_zoom <- p_zoom +
      geom_tippoint(data = tip_coords_zoom, aes(x=x, y=y), color="red", size=2) +
      #geom_tiplab(data = tip_coords_zoom, aes(x=x, y=y, label=label),
      #            color="red", size=3, align=TRUE, offset=0.5)

      #geom_tiplab2(
      #  data = tip_coords_zoom,
      #  aes(label = label),
      #  color = "red",
      #  size = 3,
      #  align = TRUE,
      #  linetype = "dotted",
      #  linesize = 0.4,
      #  offset = 0.25,
      #  angle = 0
      #)

      geom_tiplab2(
        data = tip_coords_zoom,
        aes(label = display),
        color = "red",
        size = 3,
        align = TRUE,
        linetype = "dotted",
        linesize = 0.4,
        offset = 0.25,
        angle = 0
      )
    
    zoom_file <- paste0(args$output_prefix, "_zoom.png")
    ggsave(filename = zoom_file, plot=p_zoom,
           width=args$zoom_width, height=args$zoom_height, dpi=300)
    message("Saved zoomed tree: ", zoom_file)
    
    # Optional combined plot
    if (args$combine) {
      combined_file <- paste0(args$output_prefix, "_combined.png")
      combined_plot <- plot_grid(
        p_full,
        ggdraw() + draw_plot(p_zoom, y = 0.15, height = 0.85),
        ncol = 2,
        rel_widths = c(1, 1),
        labels = c("Full Tree", "Zoomed Subtree"),
        label_size = 12
      )

      combined_plot <- combined_plot + theme(plot.background = element_rect(fill = "white", color = NA)) #otherwise it's grey in between

      ggsave(
        filename = combined_file,
        plot = combined_plot,
        width = args$full_width + args$zoom_width,
        height = max(args$full_height, args$zoom_height) * 0.55,
        dpi = 300
      )

      #combined_plot <- plot_grid(p_full, p_zoom, ncol=2, labels=c("Full Tree","Zoomed Subtree"))
      #ggsave(filename = combined_file, plot=combined_plot,
      #       width=args$full_width + args$zoom_width,
      #       height=max(args$full_height, args$zoom_height), dpi=300)
      message("Saved combined plot: ", combined_file)
    }
  } else {
    warning("Cannot find MRCA of highlighted tips; skipping zoomed plot")
  }
}
