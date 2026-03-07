#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

# ==============================
# ARGUMENT
# ==============================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("Rscript plot_stacked_AA_composition.R <results_root>\n")
  cat("Example:\n")
  cat("Rscript plot_stacked_AA_composition.R ../results\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, aa, obs_risky) %>%
    mutate(
      obs_risky = suppressWarnings(as.numeric(obs_risky)),
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
# CALCULATE PROPORTIONS PER SPECIES
# ==============================

species_summary <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(total_risky=sum(obs_risky), .groups="drop") %>%
  group_by(domain, species) %>%
  mutate(prop = total_risky / sum(total_risky)) %>%
  ungroup()

# ==============================
# DOMAIN-LEVEL PROPORTIONS
# ==============================

domain_summary <- all_data %>%
  group_by(domain, aa) %>%
  summarise(total_risky=sum(obs_risky), .groups="drop") %>%
  group_by(domain) %>%
  mutate(prop = total_risky / sum(total_risky)) %>%
  ungroup()

# ==============================
# 1️⃣ SEPARATE DOMAIN PLOTS
# ==============================

for (dm in domains) {

  df <- species_summary %>% filter(domain==dm)

  p <- ggplot(df,
              aes(x=species,
                  y=prop,
                  fill=aa)) +
    geom_col() +
    theme_minimal(base_size=11) +
    labs(title=paste("AA Composition of Risk -", dm),
         x="Species",
         y="Proportion of Risky Runs",
         fill="Amino Acid") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))

  ggsave(paste0(dm, "_AA_composition_species.png"),
         p,
         width=14,
         height=6,
         dpi=300)
}

# ==============================
# 2️⃣ COMBINED DOMAIN PLOT
# ==============================

p_combined <- ggplot(domain_summary,
                     aes(x=domain,
                         y=prop,
                         fill=aa)) +
  geom_col() +
  theme_minimal(base_size=12) +
  labs(title="AA Composition of Risk Across Domains",
       x="Domain",
       y="Proportion of Risky Runs",
       fill="Amino Acid")

ggsave("Combined_Domain_AA_composition.png",
       p_combined,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Generated separate and combined stacked barplots.\n")
