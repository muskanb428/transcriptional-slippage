#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage:\n")
  cat("Rscript plot_domain_heatmap.R <results_root> <domain>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domain_name  <- args[2]

domain_path <- file.path(results_root, domain_name)

runs_files <- list.files(
  domain_path,
  pattern="\\.runs\\.tsv$",
  recursive=TRUE,
  full.names=TRUE
)

cat("Found", length(runs_files), "files\n")

safe_read <- function(f){
  read_tsv(f, col_types = cols(.default=col_character()), progress=FALSE) %>%
    select(species, aa, obs_risky) %>%
    mutate(obs_risky = suppressWarnings(as.numeric(obs_risky))) %>%
    filter(!is.na(obs_risky))
}

all_data <- map_dfr(runs_files, safe_read)

summary_df <- all_data %>%
  group_by(species, aa) %>%
  summarise(total_risky=sum(obs_risky), .groups="drop")

heatmap_df <- summary_df %>%
  pivot_wider(names_from=aa, values_from=total_risky, values_fill=0)

heatmap_long <- heatmap_df %>%
  pivot_longer(-species, names_to="aa", values_to="value")

p <- ggplot(heatmap_long,
            aes(x=aa, y=species, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma") +
  theme_minimal(base_size=12) +
  labs(title=paste("Risky Runs Heatmap -", domain_name),
       x="Amino Acid",
       y="Species",
       fill="Risky\nRuns") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

output_file <- paste0(domain_name, "_heatmap.png")

ggsave(output_file, p, width=10, height=8, dpi=300)

cat("Saved:", output_file, "\n")
