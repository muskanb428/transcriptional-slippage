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
  cat("Rscript Figure_SuppressionIndex_AA_Ranking.R <results_root>\n")
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
# SPECIES × AA SUPPRESSION INDEX
# ==============================

species_aa_si <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(
    obs = sum(obs_risky),
    exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    SI = (obs + 1) / (exp + 1),
    log2SI = log2(SI)
  )

# ==============================
# AA-LEVEL SUMMARY
# ==============================

aa_summary <- species_aa_si %>%
  group_by(aa) %>%
  summarise(
    mean_log2SI = mean(log2SI, na.rm=TRUE),
    sd_log2SI = sd(log2SI, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(mean_log2SI)

# ==============================
# PLOT
# ==============================

p <- ggplot(aa_summary,
            aes(x=reorder(aa, mean_log2SI),
                y=mean_log2SI)) +
  geom_col(fill="steelblue") +
  geom_errorbar(aes(ymin=mean_log2SI - sd_log2SI,
                    ymax=mean_log2SI + sd_log2SI),
                width=0.2) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_minimal(base_size=13) +
  labs(title="Amino Acid Suppression Ranking Across Life",
       x="Amino Acid",
       y="Mean log2(Suppression Index)") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("Figure_SuppressionIndex_AA_Ranking.png",
       p,
       width=10,
       height=6,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_AA_Ranking.png\n")
