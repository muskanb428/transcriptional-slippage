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
  cat("Rscript Figure_SuppressionIndex_SpeciesDistribution.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

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
# SPECIES-LEVEL SUPPRESSION INDEX
# ==============================

species_si <- all_data %>%
  group_by(domain, species) %>%
  summarise(
    total_obs = sum(obs_risky),
    total_exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    SI = (total_obs + 1) / (total_exp + 1),
    log2SI = log2(SI)
  )

# ==============================
# PLOT
# ==============================

p <- ggplot(species_si,
            aes(x=domain, y=log2SI, fill=domain)) +
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  geom_jitter(width=0.15, size=2, alpha=0.8) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_minimal(base_size=13) +
  labs(title="Species-Level Suppression Index Across Domains",
       x="Domain",
       y="log2(Suppression Index)") +
  theme(legend.position="none")

ggsave("Figure_SuppressionIndex_SpeciesDistribution.png",
       p,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_SpeciesDistribution.png\n")
