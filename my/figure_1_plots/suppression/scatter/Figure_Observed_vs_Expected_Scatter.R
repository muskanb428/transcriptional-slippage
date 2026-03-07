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
  cat("Rscript Figure_Observed_vs_Expected_Scatter.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# -----------------------------
# Safe Reader
# -----------------------------
safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, obs_risky, random_expected_risky_count) %>%
    mutate(
      obs_risky = as.numeric(obs_risky),
      random_expected_risky_count = as.numeric(random_expected_risky_count),
      domain = domain_name
    ) %>%
    filter(!is.na(obs_risky))
}

# -----------------------------
# Load Data
# -----------------------------
all_data <- map_dfr(domains, function(dm){

  runs_files <- list.files(
    file.path(results_root, dm),
    pattern="\\.runs\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  map_dfr(runs_files, ~safe_read(.x, dm))
})

# -----------------------------
# Species-level totals
# -----------------------------
species_totals <- all_data %>%
  group_by(domain, species) %>%
  summarise(
    total_obs = sum(obs_risky),
    total_exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  )

# -----------------------------
# Plot function
# -----------------------------
make_plot <- function(dm){

  df <- species_totals %>% filter(domain == dm)

  ggplot(df, aes(x=total_exp, y=total_obs)) +
    geom_point(size=2, alpha=0.8) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal(base_size=12) +
    labs(title=dm,
         x="Expected Risky Runs (log10)",
         y="Observed Risky Runs (log10)")
}

pA <- make_plot("eukaryotes")
pB <- make_plot("prokaryotes")
pC <- make_plot("archea")
pD <- make_plot("viruses")

final_plot <- (pA | pB) /
              (pC | pD)

ggsave("Figure_Observed_vs_Expected_Scatter.png",
       final_plot,
       width=12,
       height=10,
       dpi=300)

cat("DONE. Saved as Figure_Observed_vs_Expected_Scatter.png\n")
