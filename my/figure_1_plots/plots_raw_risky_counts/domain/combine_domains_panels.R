#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

# ==============================
# INPUT FILES
# ==============================

file_euk  <- "eukaryotes_all_species.png"
file_pro  <- "Prokaryotes_all_species.png"
file_arc  <- "archea_all_species.png"
file_vir  <- "viruses_all_species.png"

# ==============================
# FUNCTION TO LOAD PNG AS GG OBJECT
# ==============================

load_png_plot <- function(path) {
  img <- png::readPNG(path)
  grid::rasterGrob(img, interpolate = TRUE)
}

g1 <- wrap_elements(load_png_plot(file_euk))
g2 <- wrap_elements(load_png_plot(file_pro))
g3 <- wrap_elements(load_png_plot(file_arc))
g4 <- wrap_elements(load_png_plot(file_vir))

# ==============================
# COMBINE WITH PANEL LABELS
# ==============================

final_plot <- (g1 + g2) /
              (g3 + g4) +
  plot_annotation(
    tag_levels = "A"
  )

# ==============================
# SAVE FINAL FIGURE
# ==============================

ggsave(
  "Figure1_Domain_Combined.png",
  final_plot,
  width = 16,
  height = 12,
  dpi = 300
)

cat("DONE. Saved as Figure1_Domain_Combined.png\n")
