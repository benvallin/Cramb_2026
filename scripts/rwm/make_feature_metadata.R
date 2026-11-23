# File description --------------------------------------------------------

### Goal
# make_feature_metadata.R produces gene and transcript metadata from GENECODE annotation.

### Execution steps
# ===>>> rwm/make_feature_metadata.R
# - rwm/dgea.R
# - rwm/slc18a2_expression.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - rnaseqtools package available from GitHub with 'pak::pak("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (ben.vallin@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Cramb_2026/data/rwm/make_feature_metadata"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define annotation parameters --------------------------------------------

species <- "human"

release <- 46

# Update paths ------------------------------------------------------------

sr_dir <- ifelse(species == "mouse",
                 paste0("vM", release),
                 paste0("v", release))

out_dir_sub_path <- file.path(out_dir_path, sr_dir)

# Download GENCODE annotation ---------------------------------------------

download_gencode_annotation(species = species,
                            release = release,
                            out_dir_path = out_dir_sub_path)

# Import GENCODE annotation -----------------------------------------------

gencode_annotation <- import_gencode_annotation(file = list.files(path = out_dir_sub_path,
                                                                  pattern = "^.*\\.gtf\\.gz$",
                                                                  full.names = TRUE))

# Add transgene annotation ------------------------------------------------

transgenes <- c("mTagBFP2", "puromycinR")

transgene_annotation <- tibble(gene_id = transgenes,
                               transcript_id = transgenes,
                               gene_name = transgenes,
                               gene_type = "protein_coding")

gencode_annotation <- bind_rows(gencode_annotation, transgene_annotation)

write_csv(gencode_annotation, file.path(out_dir_sub_path, "gencode_annotation.csv"))

# Make feature metadata ---------------------------------------------------

gene_metadata <- make_gene_metadata(input = gencode_annotation)

transcript_metadata <- make_transcript_metadata(input = gencode_annotation)

# Write feature metadata to file ------------------------------------------

write_csv(gene_metadata, file.path(out_dir_sub_path, "gene_metadata.csv"))

write_csv(transcript_metadata, file.path(out_dir_sub_path, "transcript_metadata.csv"))


