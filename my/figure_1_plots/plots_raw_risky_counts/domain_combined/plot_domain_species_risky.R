#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

results_root <- ifelse(length(args) >= 1, args[1], "../results")
output_dir   <- ifelse(length(args) >= 2, args[2], "domain_combined")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

runs_files <- list.files(
  results_root,
  pattern = "\\.runs\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

message("Found ", length(runs_files), " runs.tsv files")

safe_read <- function(f) {
  
  df <- read_tsv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  path_parts <- str_split(f, "/", simplify = TRUE)
  domain <- path_parts[2]
  
  df %>%
    select(species, aa, obs_risky) %>%
    mutate(
      obs_risky = suppressWarnings(as.numeric(obs_risky)),
      domain = domain
    ) %>%
    filter(!is.na(obs_risky))
}

all_data <- map_dfr(runs_files, safe_read)

species_summary <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(
    total_risky = sum(obs_risky),
    .groups = "drop"
  )

# ==============================
# ONE PLOT PER DOMAIN (ALL SPECIES INSIDE)
# ==============================

for (dm in unique(species_summary$domain)) {
  
  dm_data <- species_summary %>%
    filter(domain == dm)
  
  p <- ggplot(dm_data,
              aes(x = aa,
                  y = total_risky,
                  color = species,
                  group = species)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Observed Risky Runs per AA -", dm),
      x = "Amino Acid",
      y = "Total Risky Runs",
      color = "Species"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  ggsave(
    file.path(output_dir, paste0(dm, "_all_species.png")),
    p,
    width = 12,
    height = 6,
    dpi = 300
  )
}

message("DONE.")
