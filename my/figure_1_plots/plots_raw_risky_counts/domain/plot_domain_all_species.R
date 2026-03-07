#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
})

# ==============================
# ARGUMENTS
# ==============================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage:\n")
  cat("Rscript plot_domain_all_species.R <results_root> <domain>\n\n")
  cat("Example:\n")
  cat("Rscript plot_domain_all_species.R ../results eukaryotes\n")
  quit(save = "no", status = 1)
}

results_root <- args[1]
domain_name  <- args[2]

domain_path <- file.path(results_root, domain_name)

if (!dir.exists(domain_path)) {
  stop("Domain folder not found: ", domain_path)
}

# ==============================
# FIND runs.tsv FILES
# ==============================

runs_files <- list.files(
  domain_path,
  pattern = "\\.runs\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(runs_files), "runs.tsv files in", domain_name, "\n")

if (length(runs_files) == 0) {
  stop("No runs.tsv files found inside domain.")
}

# ==============================
# SAFE READ FUNCTION
# ==============================

safe_read <- function(f) {
  
  df <- read_tsv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  df %>%
    select(species, aa, obs_risky) %>%
    mutate(obs_risky = suppressWarnings(as.numeric(obs_risky))) %>%
    filter(!is.na(obs_risky))
}

# ==============================
# LOAD ALL SPECIES
# ==============================

all_data <- map_dfr(runs_files, safe_read)

# ==============================
# SUMMARISE
# ==============================

species_summary <- all_data %>%
  group_by(species, aa) %>%
  summarise(
    total_risky = sum(obs_risky),
    .groups = "drop"
  )

# ==============================
# PLOT
# ==============================

output_file <- paste0(domain_name, "_all_species.png")

p <- ggplot(species_summary,
            aes(x = aa,
                y = total_risky,
                color = species,
                group = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste("Observed Risky Runs per Amino Acid -", domain_name),
    x = "Amino Acid",
    y = "Total Risky Runs",
    color = "Species"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(output_file, p, width = 12, height = 6, dpi = 300)

cat("DONE. Output saved as:", output_file, "\n")
