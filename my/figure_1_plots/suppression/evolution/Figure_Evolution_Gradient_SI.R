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
  cat("Rscript Figure_Evolution_Gradient_SI.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]

# Ordered evolution axis
domain_order <- c("archea", "prokaryotes", "eukaryotes", "viruses")

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
all_data <- map_dfr(domain_order, function(dm){

  runs_files <- list.files(
    file.path(results_root, dm),
    pattern="\\.runs\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  map_dfr(runs_files, ~safe_read(.x, dm))
})

# -----------------------------
# Species-Level SI
# -----------------------------
species_si <- all_data %>%
  group_by(domain, species) %>%
  summarise(
    total_obs = sum(obs_risky),
    total_exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    SI = (total_obs + 1)/(total_exp + 1),
    log2SI = log2(SI)
  )

species_si$domain <- factor(species_si$domain,
                             levels = domain_order)

# Domain medians
domain_medians <- species_si %>%
  group_by(domain) %>%
  summarise(median_log2SI = median(log2SI),
            .groups="drop")

# -----------------------------
# Plot
# -----------------------------
p <- ggplot() +
  geom_jitter(data=species_si,
              aes(x=domain, y=log2SI),
              width=0.1,
              alpha=0.6,
              size=2,
              color="grey40") +
  geom_line(data=domain_medians,
            aes(x=domain, y=median_log2SI, group=1),
            size=1.2,
            color="steelblue") +
  geom_point(data=domain_medians,
             aes(x=domain, y=median_log2SI),
             size=4,
             color="steelblue") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_minimal(base_size=13) +
  labs(title="Evolutionary Gradient of Homopolymer Suppression",
       x="Domain (Evolutionary Order)",
       y="Median log2(Suppression Index)")

ggsave("Figure_Evolution_Gradient_SI.png",
       p,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Saved as Figure_Evolution_Gradient_SI.png\n")
