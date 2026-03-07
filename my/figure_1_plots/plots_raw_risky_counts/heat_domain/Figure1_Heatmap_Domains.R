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
  cat("Rscript Figure1_Heatmap_Domains.R <results_root>\n")
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
# DOMAIN LEVEL SUMMARY
# ==============================

heatmap_data <- all_data %>%
  group_by(domain, aa) %>%
  summarise(
    obs = sum(obs_risky),
    exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    log2OE = log2((obs + 1) / (exp + 1))
  )

# Optional: Order amino acids alphabetically
heatmap_data$aa <- factor(heatmap_data$aa,
                          levels = sort(unique(heatmap_data$aa)))

# ==============================
# PLOT HEATMAP
# ==============================

p <- ggplot(heatmap_data,
            aes(x=aa, y=domain, fill=log2OE)) +
  geom_tile(color="white") +
  scale_fill_gradient2(
    low="blue",
    mid="white",
    high="red",
    midpoint=0,
    name="log2(O/E)"
  ) +
  theme_minimal(base_size=13) +
  labs(title="Risky Amino Acid Enrichment Across Domains",
       x="Amino Acid",
       y="Domain") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("Figure1_Heatmap_Domains.png",
       p,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Saved as Figure1_Heatmap_Domains.png\n")
