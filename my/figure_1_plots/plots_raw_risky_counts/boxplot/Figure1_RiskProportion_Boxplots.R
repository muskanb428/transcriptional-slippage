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
  cat("Rscript Figure1_RiskProportion_Boxplots.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, obs_risky) %>%
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
# CALCULATE RISK PROPORTION PER SPECIES
# ==============================

species_risk <- all_data %>%
  group_by(domain, species) %>%
  summarise(
    total_runs = n(),
    total_risky = sum(obs_risky),
    risk_prop = total_risky / total_runs,
    .groups="drop"
  )

# ==============================
# PANEL A (ALL DOMAINS)
# ==============================

pA <- ggplot(species_risk,
             aes(x=domain, y=risk_prop, fill=domain)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) +
  geom_jitter(width=0.15, size=2, alpha=0.7) +
  theme_minimal(base_size=13) +
  labs(title="A  |  Risk Proportion Across Domains",
       x="Domain",
       y="Risk Proportion") +
  theme(legend.position="none")

# ==============================
# DOMAIN-SPECIFIC PANELS
# ==============================

make_domain_panel <- function(dm){

  df <- species_risk %>% filter(domain==dm)

  ggplot(df,
         aes(x=species, y=risk_prop)) +
    geom_col(fill="steelblue") +
    theme_minimal(base_size=11) +
    labs(title=dm,
         x="Species",
         y="Risk Proportion") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))
}

pB <- make_domain_panel("eukaryotes")
pC <- make_domain_panel("Prokaryotes")
pD <- make_domain_panel("archea")
pE <- make_domain_panel("viruses")

# ==============================
# COMBINE PANELS
# ==============================

final_plot <- pA /
             (pB | pC) /
             (pD | pE) +
  plot_annotation(tag_levels="A")

ggsave("Figure1_RiskProportion_AllDomains.png",
       final_plot,
       width=18,
       height=22,
       dpi=300)

cat("DONE. Saved as Figure1_RiskProportion_AllDomains.png\n")
