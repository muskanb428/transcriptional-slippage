#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("Rscript Figure_SuppressionIndex_Panels.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, aa, obs_risky, random_expected_risky_count) %>%
    mutate(
      obs_risky = as.numeric(obs_risky),
      random_expected_risky_count = as.numeric(random_expected_risky_count),
      domain = domain_name
    ) %>%
    filter(!is.na(obs_risky))
}

# ==============================
# LOAD DATA
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
# FUNCTION TO BUILD DOMAIN PANEL
# ==============================

make_domain_panel <- function(dm){

  df <- all_data %>%
    filter(domain == dm) %>%
    group_by(species, aa) %>%
    summarise(
      obs = sum(obs_risky),
      exp = sum(random_expected_risky_count, na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(
      SI = (obs + 1) / (exp + 1),
      log2SI = log2(SI)
    )

  ggplot(df,
         aes(x=aa, y=species, fill=log2SI)) +
    geom_tile(color="white") +
    scale_fill_gradient2(
      low="blue",
      mid="white",
      high="red",
      midpoint=0,
      limits=c(-3,3),
      name="log2(SI)"
    ) +
    theme_minimal(base_size=10) +
    labs(title=dm,
         x="Amino Acid",
         y="Species") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}

# ==============================
# BUILD PANELS
# ==============================

pA <- make_domain_panel("eukaryotes")
pB <- make_domain_panel("prokaryotes")
pC <- make_domain_panel("archea")
pD <- make_domain_panel("viruses")

final_plot <- (pA | pB) /
              (pC | pD) +
  plot_annotation(tag_levels="A")

ggsave("Figure_SuppressionIndex_AllDomains_Panels.png",
       final_plot,
       width=18,
       height=20,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_AllDomains_Panels.png\n")
