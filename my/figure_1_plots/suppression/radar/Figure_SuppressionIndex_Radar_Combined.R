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
  cat("Rscript Figure_SuppressionIndex_Radar_Combined.R <results_root>\n")
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
# DOMAIN-LEVEL SI
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

# Order amino acids consistently
aa_levels <- sort(unique(si_data$aa))
si_data$aa <- factor(si_data$aa, levels = aa_levels)

# Close polygons
si_closed <- si_data %>%
  group_by(domain) %>%
  arrange(aa) %>%
  do(rbind(., .[1,])) %>%
  ungroup()

# ==============================
# PLOT
# ==============================

p <- ggplot(si_closed,
            aes(x=aa, y=log2SI, group=domain, color=domain, fill=domain)) +
  geom_polygon(alpha=0.15) +
  geom_line(size=1) +
  coord_polar() +
  theme_minimal(base_size=13) +
  labs(title="Combined Suppression Index Radar",
       y="log2(SI)",
       x=NULL) +
  theme(axis.text.x = element_text(size=8))

ggsave("Figure_SuppressionIndex_Radar_Combined.png",
       p,
       width=8,
       height=8,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_Radar_Combined.png\n")
