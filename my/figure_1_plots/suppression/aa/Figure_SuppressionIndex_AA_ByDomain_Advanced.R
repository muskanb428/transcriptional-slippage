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
  cat("Rscript Figure_SuppressionIndex_AA_ByDomain_Advanced.R <results_root>\n")
  quit(save="no", status=1)
}

results_root <- args[1]
domains <- c("eukaryotes", "prokaryotes", "archea", "viruses")

# -----------------------------
# Amino Acid Classification
# -----------------------------
aa_class <- data.frame(
  aa = c("A","V","L","I","M","F","W","Y",
         "S","T","N","Q","C",
         "K","R","H",
         "D","E",
         "G","P"),
  class = c(rep("Hydrophobic",8),
            rep("Polar",5),
            rep("Basic",3),
            rep("Acidic",2),
            rep("Special",2))
)

# -----------------------------
# Safe Reader
# -----------------------------
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

# -----------------------------
# Load Data
# -----------------------------
all_data <- map_dfr(domains, function(dm){

  runs_files <- list.files(
    file.path(results_root, dm),
    pattern="\\.runs\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  map_dfr(runs_files, ~safe_read(.x, dm))
})

# -----------------------------
# Compute Species × AA SI
# -----------------------------
species_aa_si <- all_data %>%
  group_by(domain, species, aa) %>%
  summarise(
    obs = sum(obs_risky),
    exp = sum(random_expected_risky_count, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    SI = (obs + 1)/(exp + 1),
    log2SI = log2(SI)
  )

# -----------------------------
# Domain-Level Summary
# -----------------------------
domain_summary <- species_aa_si %>%
  group_by(domain, aa) %>%
  summarise(
    mean_log2SI = mean(log2SI, na.rm=TRUE),
    sd_log2SI = sd(log2SI, na.rm=TRUE),
    p_value = wilcox.test(log2SI, mu=0)$p.value,
    .groups="drop"
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "*", ""),
    neg_signif = ifelse(p_value < 0.05 & mean_log2SI < 0, "*", "")
  ) %>%
  left_join(aa_class, by="aa")

# -----------------------------
# Plot Function
# -----------------------------
make_plot <- function(dm){

  df <- domain_summary %>% filter(domain == dm)

  ggplot(df,
         aes(x=reorder(aa, mean_log2SI),
             y=mean_log2SI,
             fill=class)) +
    geom_col() +
    geom_errorbar(aes(ymin=mean_log2SI - sd_log2SI,
                      ymax=mean_log2SI + sd_log2SI),
                  width=0.2) +
    geom_text(aes(label=neg_signif,
                  y=mean_log2SI - 0.2),
              size=4) +
    geom_hline(yintercept=0, linetype="dashed") +
    theme_minimal(base_size=12) +
    labs(title=dm,
         x="Amino Acid",
         y="Mean log2(SI)") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5),
          legend.position="none")
}

pA <- make_plot("eukaryotes")
pB <- make_plot("prokaryotes")
pC <- make_plot("archea")
pD <- make_plot("viruses")

final_plot <- (pA | pB) /
              (pC | pD)

ggsave("Figure_SuppressionIndex_AA_ByDomain_Advanced.png",
       final_plot,
       width=14,
       height=12,
       dpi=300)

cat("DONE. Saved as Figure_SuppressionIndex_AA_ByDomain_Advanced.png\n")
