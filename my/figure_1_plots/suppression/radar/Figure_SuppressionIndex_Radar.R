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
  cat("Rscript Figure_SuppressionIndex_Radar.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(aa, obs_risky, random_expected_risky_count) %>%
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
# DOMAIN-LEVEL SUPPRESSION INDEX
# ==============================

si_data <- all_data %>%
  group_by(domain, aa) %>%
  summarise(
    obs = sum(obs_risky),
    exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    SI = (obs + 1) / (exp + 1),
    log2SI = log2(SI)
  )

# ==============================
# RADAR FUNCTION
# ==============================

make_radar <- function(dm){

  df <- si_data %>%
    filter(domain == dm) %>%
    arrange(aa)

  # close polygon
  df <- rbind(df, df[1,])

  ggplot(df, aes(x=aa, y=log2SI, group=1)) +
    geom_polygon(fill="steelblue", alpha=0.4) +
    geom_line(color="steelblue", size=1) +
    coord_polar() +
    theme_minimal(base_size=12) +
    labs(title=dm,
         x=NULL,
         y="log2(SI)") +
    theme(axis.text.x=element_text(size=8))
}

pA <- make_radar("eukaryotes")
pB <- make_radar("prokaryotes")
pC <- make_radar("archea")
pD <- make_radar("viruses")

final_plot <- (pA | pB) /
              (pC | pD) +
  plot_annotation(tag_levels="A")

ggsave("Figure_SuppressionIndex_Radar_AllDomains.png",
       final_plot,
       width=14,
       height=12,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_Radar_AllDomains.png\n")
