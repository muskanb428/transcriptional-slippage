#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("Rscript Figure1_all_domains_heatmap.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]

domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, aa, obs_risky) %>%
    mutate(
      obs_risky = suppressWarnings(as.numeric(obs_risky)),
      domain = domain_name
    ) %>%
    filter(!is.na(obs_risky))
}

# ==============================
# LOAD ALL DATA
# ==============================

all_data <- map_dfr(domains, function(dm){

  runs_files <- list.files(
    file.path(results_root, dm),
    pattern="\\.runs\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  map_dfr(runs_files, ~safe_read(.x, dm))
})

# ==============================
# SUMMARISE
# ==============================

summary_df <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(total_risky = sum(obs_risky), .groups="drop") %>%
  mutate(value = log10(total_risky + 1))

# ==============================
# HEATMAP FUNCTION
# ==============================

make_heatmap <- function(df, title_text){

  ggplot(df, aes(x=aa, y=species, fill=value)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma") +
    theme_minimal(base_size=11) +
    labs(title=title_text,
         x="Amino Acid",
         y="Species",
         fill="log10\n(Risky+1)") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))
}

# ==============================
# PANELS
# ==============================

pA <- make_heatmap(summary_df, "All Domains Combined")
pB <- make_heatmap(filter(summary_df, domain=="eukaryotes"), "Eukaryotes")
pC <- make_heatmap(filter(summary_df, domain=="prokaryotes"), "Prokaryotes")
pD <- make_heatmap(filter(summary_df, domain=="archea"), "Archaea")
pE <- make_heatmap(filter(summary_df, domain=="viruses"), "Viruses")

# ==============================
# COMBINE
# ==============================

final_plot <- pA /
             (pB | pC) /
             (pD | pE) +
  plot_annotation(tag_levels = "A")

ggsave("Figure1_Heatmap_AllDomains.png",
       final_plot,
       width = 18,
       height = 20,
       dpi = 300)

cat("DONE. Saved as Figure1_Heatmap_AllDomains.png\n")
