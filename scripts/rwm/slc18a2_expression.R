# File description --------------------------------------------------------

### Goal
# slc18a2_expression.R summarises DGEA results and TPM values for the SLC18A2 gene.

### Execution steps
# x rwm/make_feature_metadata.R
# x rwm/dgea.R
# ===>>> rwm/slc18a2_expression.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - gene_metadata.csv produced by make_feature_metadata.R
# - DGEA results produced by dgea.R
# - rnaseqtools package available from GitHub with 'pak::pak("benvallin/rnaseqtools")'
# - paramstream package available from GitHub with 'pak::pak("benvallin/paramstream")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (ben.vallin@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Cramb_2026/data/rwm/slc18a2_expression"

# Path to directory containing DGEA results
dgea_dir_path <- "publications/Cramb_2026/data/rwm/dgea"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Cramb_2026/data/rwm/make_feature_metadata/v46/gene_metadata.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(rnaseqtools)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE) 

# Define DGE analysis parameters ------------------------------------------

dge_params <- list(
  coef_var = "genotype",
  coef_ref_lvl = "CTRL",
  coef_test_lvl = "TRIP",
  model_formula = "~ genotype",
  padj_dge = 0.05,
  min_cnt_incl = 10,
  min_samp_incl = 3,
  protein_coding_only = TRUE)

# Stream analysis parameter values ----------------------------------------

input_param_log_id <- get_param_log_id(param_list = dge_params,
                                       in_dir_path = dgea_dir_path)

write_param_log(param_list = dge_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = dge_params,
                                        in_dir_path = out_dir_path)

# Load feature metadata ---------------------------------------------------

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Load DGEA results -------------------------------------------------------

sample_metadata <- read_rds(file.path(dgea_dir_path, paste0("sample_metadata_", input_param_log_id, ".rds")))

test_names <- c("deseq2_wald", "deseq2_lrt", 
                "limma_trend", "limma_voom", 
                "edger_glm", "edger_glmql")

for(i in test_names) {
  
  assign(x = paste0("lfc_", i),
         value = read_csv(file = file.path(dgea_dir_path, paste0("lfc_", i, "_", input_param_log_id, ".csv"))))
  
  assign(x = paste0("lfc_sign_", i),
         value = read_csv(file = file.path(dgea_dir_path, paste0("lfc_sign_", i, "_", input_param_log_id, ".csv"))))
  
}

tpm <- as_tibble(fread(file = file.path(dgea_dir_path, paste0("tpm_", input_param_log_id, ".csv")))) 

# Compute mean TPMs -------------------------------------------------------

mean_tpm <- compute_mean_tpm(input = tpm,
                             gene_id = "ensembl_gene_id_version",
                             sample_metadata = sample_metadata,
                             sample_id_var = "sample_id",
                             group_id_vars = "genotype",
                             cnt_nm = "tpm")

# Focus on SLC18A2 --------------------------------------------------------

# SLC18A2 DGEA results
get_gene_dgea_results <- function(lfc_data, gene_symbol = "SLC18A2") {
  
  lfc_data[lfc_data$gene_symbol == gene_symbol,
           c("ensembl_gene_id_version", "gene_symbol", 
             "log2FoldChange", "pvalue", "padj")]
  
}

lfc_slc18a2 <- lfc_deseq2_wald |> 
  get_gene_dgea_results() |> 
  mutate(dgea_method = "deseq2_wald") |> 
  bind_rows(lfc_deseq2_lrt |> 
              get_gene_dgea_results() |> 
              mutate(dgea_method = "deseq2_lrt")) |> 
  bind_rows(lfc_edger_glm |> 
              get_gene_dgea_results() |> 
              mutate(dgea_method = "edger_glm")) |> 
  bind_rows(lfc_edger_glmql |> 
              get_gene_dgea_results() |> 
              mutate(dgea_method = "edger_glmql")) |> 
  bind_rows(lfc_limma_trend |> 
              get_gene_dgea_results() |> 
              mutate(dgea_method = "limma_trend")) |> 
  bind_rows(lfc_limma_voom |> 
              get_gene_dgea_results() |> 
              mutate(dgea_method = "limma_voom")) |> 
  select(dgea_method, everything())

# Sample-level SLC18A2 TPMs 
tpm_slc18a2 <- tpm |> 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) |> 
  filter(gene_symbol == "SLC18A2") |> 
  pivot_longer(-c(ensembl_gene_id_version, gene_symbol),
               names_to = "sample_id",
               values_to = "TPM") |> 
  left_join(sample_metadata[, c("sample_id", "line_name", "line_name_generic", "genotype")],
            by = join_by(sample_id)) |> 
  select(line_name_generic, line_name, genotype, ensembl_gene_id_version, gene_symbol, TPM) |> 
  arrange(gene_symbol, line_name_generic)

# Genotype-level SLC18A2 mean TPMs
mean_tpm_scl18a2 <- mean_tpm |> 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) |> 
  filter(gene_symbol == "SLC18A2") |>  
  select(genotype, ensembl_gene_id_version, gene_symbol, everything(), -sample_data) |> 
  arrange(gene_symbol, genotype)

# Write results to file ---------------------------------------------------

write_csv(lfc_slc18a2, file.path(out_dir_path, "lfc.csv"))

write_csv(tpm_slc18a2, file.path(out_dir_path, "tpm.csv"))

write_csv(mean_tpm_scl18a2, file.path(out_dir_path, "mean_tpm.csv"))



