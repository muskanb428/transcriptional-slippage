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
  cat("Rscript Figure1_Volcano_Domains.R <results_root>\n")
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
# FUNCTION TO BUILD VOLCANO
# ==============================

make_volcano <- function(dm){

  df <- all_data %>%
    filter(domain == dm) %>%
    group_by(aa) %>%
    summarise(
      obs = sum(obs_risky),
      exp = sum(random_expected_risky_count, na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(
      log2OE = log2((obs + 1) / (exp + 1)),
      pvalue = ppois(obs - 1, lambda = exp, lower.tail = FALSE),
      neglog10p = -log10(pvalue + 1e-300)
    )

  ggplot(df, aes(x=log2OE, y=neglog10p)) +
    geom_point(size=3, alpha=0.8) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    theme_minimal(base_size=12) +
    labs(title=dm,
         x="log2(Observed / Expected)",
         y="-log10(p-value)")
}

pA <- make_volcano("eukaryotes")
pB <- make_volcano("prokaryotes")
pC <- make_volcano("archea")
pD <- make_volcano("viruses")

final_plot <- (pA | pB) /
              (pC | pD) +
  plot_annotation(tag_levels="A")

ggsave("Figure1_Volcano_Domains.png",
       final_plot,
       width=14,
       height=12,
       dpi=300)

cat("DONE. Saved as Figure1_Volcano_Domains.png\n")
