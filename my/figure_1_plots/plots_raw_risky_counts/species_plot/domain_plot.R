#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

# ==============================
# ARGUMENTS
# ==============================

args <- commandArgs(trailingOnly = TRUE)

results_root <- ifelse(length(args) >= 1, args[1], "results")
output_dir   <- ifelse(length(args) >= 2, args[2], "plots/domain_species_risky")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================
# FIND FILES
# ==============================

runs_files <- list.files(
  results_root,
  pattern = "\\.runs\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(runs_files) == 0) {
  stop("No .runs.tsv files found!")
}

message("Found ", length(runs_files), " runs.tsv files")

# ==============================
# SAFE LOADER
# ==============================

safe_read <- function(f) {
  
  # read everything as character to avoid type conflicts
  df <- read_tsv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  # extract domain
  path_parts <- str_split(f, "/", simplify = TRUE)
  domain <- path_parts[2]
  
  # keep only required columns
  df <- df %>%
    select(species, aa, obs_risky) %>%
    mutate(
      obs_risky = suppressWarnings(as.numeric(obs_risky)),
      domain = domain
    )
  
  df
}

# ==============================
# LOAD ALL
# ==============================

all_data <- map_dfr(runs_files, safe_read)

# remove NAs from obs_risky
all_data <- all_data %>%
  filter(!is.na(obs_risky))

# ==============================
# SPECIES SUMMARY
# ==============================

species_summary <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(
    total_risky = sum(obs_risky),
    .groups = "drop"
  )

# ==============================
# DOMAIN SUMMARY
# ==============================

domain_summary <- all_data %>%
  group_by(domain, aa) %>%
  summarise(
    total_risky = sum(obs_risky),
    .groups = "drop"
  )

# ==============================
# SPECIES PLOTS
# ==============================

message("Generating species plots...")

for (sp in unique(species_summary$species)) {
  
  sp_data <- species_summary %>%
    filter(species == sp)
  
  p <- ggplot(sp_data, aes(x = aa, y = total_risky)) +
    geom_col() +
    facet_wrap(~ domain, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Observed Risky Runs per AA -", sp),
      x = "Amino Acid",
      y = "Total Risky Runs"
    )
  
  ggsave(
    file.path(output_dir, paste0(sp, "_risky.png")),
    p,
    width = 8,
    height = 4,
    dpi = 300
  )
}

# ==============================
# DOMAIN PLOTS
# ==============================

message("Generating domain plots...")

for (dm in unique(domain_summary$domain)) {
  
  dm_data <- domain_summary %>%
    filter(domain == dm)
  
  p <- ggplot(dm_data, aes(x = aa, y = total_risky)) +
    geom_col() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Overall Observed Risky Runs -", dm),
      x = "Amino Acid",
      y = "Total Risky Runs"
    )
  
  ggsave(
    file.path(output_dir, paste0(dm, "_overall.png")),
    p,
    width = 8,
    height = 4,
    dpi = 300
  )
}

# ==============================
# SAVE TABLES
# ==============================

write_tsv(species_summary,
          file.path(output_dir, "species_risky_summary.tsv"))

write_tsv(domain_summary,
          file.path(output_dir, "domain_risky_summary.tsv"))

message("DONE successfully.")
