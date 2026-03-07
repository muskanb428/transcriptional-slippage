#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("Rscript Figure1_Bubble_AllDomains.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# ==============================
# SAFE READER
# ==============================

safe_read <- function(f, domain_name){
  read_tsv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    select(species, aa, obs_risky, random_expected_risky_count) %>%
    mutate(
      obs_risky = suppressWarnings(as.numeric(obs_risky)),
      random_expected_risky_count = suppressWarnings(as.numeric(random_expected_risky_count)),
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
# DOMAIN LEVEL SUMMARY
# ==============================

domain_summary <- all_data %>%
  group_by(domain, aa) %>%
  summarise(
    obs = sum(obs_risky),
    exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    prop = obs / sum(obs),
    enrich = log2((obs + 1) / (exp + 1))
  )

# ==============================
# PANEL A (ALL DOMAINS)
# ==============================

pA <- ggplot(domain_summary,
             aes(x=aa,
                 y=domain,
                 size=prop,
                 color=enrich)) +
  geom_point(alpha=0.8) +
  scale_color_gradient2(low="blue",
                        mid="white",
                        high="red",
                        midpoint=0) +
  theme_minimal(base_size=12) +
  labs(title="A  |  Risky AA Landscape Across Domains",
       x="Amino Acid",
       y="Domain",
       size="Risk Proportion",
       color="log2(O/E)") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

# ==============================
# FUNCTION FOR DOMAIN PANEL
# ==============================

make_domain_panel <- function(dm){

  df <- all_data %>%
    filter(domain==dm) %>%
    group_by(species, aa) %>%
    summarise(
      obs = sum(obs_risky),
      exp = sum(random_expected_risky_count, na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(
      prop = obs / sum(obs),
      enrich = log2((obs + 1) / (exp + 1))
    )

  ggplot(df,
         aes(x=aa,
             y=species,
             size=prop,
             color=enrich)) +
    geom_point(alpha=0.8) +
    scale_color_gradient2(low="blue",
                          mid="white",
                          high="red",
                          midpoint=0) +
    theme_minimal(base_size=10) +
    labs(title=dm,
         x="Amino Acid",
         y="Species",
         size="Risk Prop",
         color="log2(O/E)") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))
}

pB <- make_domain_panel("eukaryotes")
pC <- make_domain_panel("prokaryotes")
pD <- make_domain_panel("archea")
pE <- make_domain_panel("viruses")

# ==============================
# COMBINE PANELS
# ==============================

final_plot <- pA /
             (pB | pC) /
             (pD | pE) +
  plot_annotation(tag_levels="A")

ggsave("Figure1_Bubble_AllDomains.png",
       final_plot,
       width=18,
       height=22,
       dpi=300)

cat("DONE. Saved as Figure1_Bubble_AllDomains.png\n")
