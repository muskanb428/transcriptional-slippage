#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("Rscript Figure_GC_vs_Risk.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# -----------------------------
# Safe Reader
# -----------------------------
safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species,
           obs_risky,
           dna_G_frac,
           dna_C_frac) %>%
    mutate(
      obs_risky = as.numeric(obs_risky),
      dna_G_frac = as.numeric(dna_G_frac),
      dna_C_frac = as.numeric(dna_C_frac),
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
# Species-Level Summary
# -----------------------------
species_gc_risk <- all_data %>%
  group_by(domain, species) %>%
  summarise(
    total_risky = sum(obs_risky),
    mean_G = mean(dna_G_frac, na.rm=TRUE),
    mean_C = mean(dna_C_frac, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    GC_content = mean_G + mean_C
  )

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(species_gc_risk,
            aes(x=GC_content,
                y=total_risky,
                color=domain)) +
  geom_point(size=2, alpha=0.8) +
  geom_smooth(method="lm", se=FALSE) +
  scale_y_log10() +
  theme_minimal(base_size=13) +
  labs(title="GC Content vs Homopolymer Risk",
       x="GC Content",
       y="Total Risky Runs (log10)",
       color="Domain")

ggsave("Figure_GC_vs_Risk.png",
       p,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Saved as Figure_GC_vs_Risk.png\n")
