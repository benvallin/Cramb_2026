# File description --------------------------------------------------------

### Goal
# dgea.R runs differential gene expression analysis (DGEA) using DESeq2, edgeR and limma methods.

### Execution steps
# x rwm/make_feature_metadata.R
# ===>>> rwm/dgea.R
# - rwm/slc18a2_expression.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - sample_metadata.csv available on GitHub repository 'benvallin/Cramb_2026', at metadata/rwm/sample_metadata.csv
# - gene_metadata.csv produced by make_feature_metadata.R
# - transcript_metadata.csv produced by make_feature_metadata.R
# - Salmon's quant.sf files available from the Gene Expression Omnibus (GEO) under accession number GSE353879
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
out_dir_path <- "publications/Cramb_2026/data/rwm/dgea"

# Path to sample metadata file (csv)
sample_metadata_path <- "publications/Cramb_2026/metadata/rwm/sample_metadata.csv"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Cramb_2026/data/rwm/make_feature_metadata/v46/gene_metadata.csv"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "publications/Cramb_2026/data/rwm/make_feature_metadata/v46/transcript_metadata.csv"

# Path to directory containing Salmon quant files
salmon_quant_dir_path <- "publications/Thomas-Wright_2026/seq_data_Thomas-Wright_2026"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(arrow)
library(tximport)
library(DESeq2)
library(edgeR)
library(limma)
library(rnaseqtools)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Set n_workers for parallel processing
n_workers <- 48

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

# Load feature metadata ---------------------------------------------------

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Load transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path) 

# Load sample metadata ----------------------------------------------------

sample_metadata <- read_csv(sample_metadata_path)

# Check DGE analysis parameters -------------------------------------------

# Check unique dge_params entries
n_dge_params <- length(names(dge_params))
n_unique_dge_params <- length(unique(names(dge_params)))

if(n_dge_params != n_unique_dge_params) {
  
  non_unique <- setNames(object = as.integer(table(names(dge_params))), nm = names(table(names(dge_params))))
  non_unique <- names(non_unique[non_unique > 1])
  
  stop("The following dge_params entries are not unique: ", paste0(non_unique, collapse = ", "))
  
}

# Define required dge_params entries
req_dge_params <- c("coef_var", "coef_ref_lvl", "coef_test_lvl", 
                    "model_formula", "padj_dge", "min_cnt_incl", "min_samp_incl")

# Define optional dge_params entries of type logical
opt_lgl_dge_params <- c("protein_coding_only", "exclude_mt_genes")

# Define optional dge_params entries of type numeric
opt_num_dge_params <- c()

# Define optional dge_params entries of type character
opt_chr_dge_params <- c()

# Define dge_params entries which should be names of logical columns in sample_metadata
met_lgl_dge_params <- c()

# Define reserved dge_params entries
res_dge_params <- unique(c(req_dge_params, opt_lgl_dge_params, opt_num_dge_params, 
                           opt_chr_dge_params, met_lgl_dge_params))

# Define reserved entries (required + optional) included in user's dge_params
for(i in c("req_dge_params", "opt_lgl_dge_params", "opt_num_dge_params", 
           "opt_chr_dge_params", "met_lgl_dge_params", "res_dge_params")) {
  
  assign(x = paste0("incl_", i),
         value = names(dge_params)[names(dge_params) %in% eval(as.symbol(i))])
  
}

# Check mandatory dge_params entries
missing <- req_dge_params[!req_dge_params %in% incl_req_dge_params]

if(length(missing) > 0L) {
  
  stop("The following mandatory entries are missing from dge_params: ", paste0(missing, collapse = ", "))
  
}

# Check optional dge_params entries of type logical
not_lgl <- incl_opt_lgl_dge_params[map_lgl(.x = incl_opt_lgl_dge_params, 
                                           .f = ~ !is.logical(dge_params[[.x]]))]

if(length(not_lgl) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected logical type: ", paste0(not_lgl, collapse = ", "))
  
}

# Check optional dge_params entries of type numeric
not_num <- incl_opt_num_dge_params[map_lgl(.x = incl_opt_num_dge_params, 
                                           .f = ~ !is.numeric(dge_params[[.x]]))]

if(length(not_num) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected numeric type: ", paste0(not_num, collapse = ", "))
  
}

# Check optional dge_params entries of type character
not_chr <- incl_opt_chr_dge_params[map_lgl(.x = incl_opt_chr_dge_params, 
                                           .f = ~ !is.character(dge_params[[.x]]))]

if(length(not_chr) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected character type: ", paste0(not_chr, collapse = ", "))
  
}

# Check dge_params entries which should be names of logical columns in sample_metadata
met_lgl_exp <- setNames(object = map_chr(.x = incl_met_lgl_dge_params, 
                                         .f = ~ as.character(dge_params[[.x]])), 
                        nm = incl_met_lgl_dge_params)

not_met_lgl <- met_lgl_exp[!map_lgl(.x = met_lgl_exp, 
                                    .f = ~ .x %in% colnames(sample_metadata) && is.logical(sample_metadata[[.x]]))]   

if(length(not_met_lgl) > 0L) {
  
  stop("The following dge_params entries are not names of logical columns in sample_metadata: ", 
       paste0(paste0(names(not_met_lgl), " - ", not_met_lgl), collapse = ", "))
  
}

# Define expected cell filter vars
# => dge_params entries that are not reserved
exp_cell_filt_vars <- setdiff(names(dge_params), incl_res_dge_params)

# Define actual cell filter vars
# => Expected cell filter vars which are column names of sample_metadata
cell_filt_vars <- exp_cell_filt_vars[exp_cell_filt_vars %in% colnames(sample_metadata)]

# Define problematic cell filter vars
# => Expected cell filter vars which are not column names of sample_metadata
prob_cell_filt_vars <- exp_cell_filt_vars[!exp_cell_filt_vars %in% cell_filt_vars]

if(length(prob_cell_filt_vars) > 0) {
  
  stop("The following dge_params entries are not column names of sample_metadata: ", paste0(prob_cell_filt_vars, collapse = ", "))
  
}

# Filter sample metadata --------------------------------------------------

# Subset to sample filter var values defined in dge_params
if(length(cell_filt_vars) > 0) {
  
  for(i in cell_filt_vars) {
    
    sample_metadata <- sample_metadata[sample_metadata[[i]] %in% dge_params[[i]],]
    
  }
  
}

# Set coef_var levels
sample_metadata[[dge_params$coef_var]] = fct_relevel(sample_metadata[[dge_params$coef_var]],
                                                     c(dge_params$coef_ref_lvl, dge_params$coef_test_lvl))

# Drop unused levels in categorical variables
sample_metadata <- sample_metadata %>%
  mutate(across(.cols = where(is.factor),
                .fns = droplevels))

# Construct design matrix -------------------------------------------------

# Create design matrix
design <- model.matrix(object = formula(dge_params$model_formula),
                       data = sample_metadata)

# Make the matrix full-rank
# => Remove columns with all 0 from design matrix
design <- make_full_rank(model_matrix = design)

# Construct coef_nm 
coef_nm <- paste0(dge_params$coef_var, dge_params$coef_test_lvl)

# Load Salmon's quantification results ------------------------------------

# Construct paths to Salmon's quant.sf files
salmon_quant_files <- setNames(object = map_chr(.x = sample_metadata$sample_id,
                                                .f = ~ file.path(salmon_quant_dir_path, paste0(.x, "_quant.sf"))),
                               nm = sample_metadata$sample_id)

# Construct txi object
txi <- tximport(files = salmon_quant_files,
                type = "salmon",
                txOut = FALSE,
                countsFromAbundance = "no",
                tx2gene = transcript_metadata)


# Construct dds
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = sample_metadata,
                                design = design)

# Filter genes ------------------------------------------------------------

# Filter out poorly expressed genes
if(all(c("min_cnt_incl", "min_samp_incl") %in% incl_res_dge_params)) {
  
  keep <- rowSums(counts(dds) >= dge_params$min_cnt_incl) >= dge_params$min_samp_incl
  dds <- dds[keep,]
  
}

# Subset to protein-coding genes only
if("protein_coding_only" %in% incl_res_dge_params && 
   dge_params$protein_coding_only) {
  
  protein_coding_genes <- gene_metadata %>%
    filter(gene_type == "protein_coding") %>%
    pull(ensembl_gene_id_version)
  
  dds <- dds[rownames(dds) %in% protein_coding_genes,]
  
}

# Exclude mitochondrial genes
if("exclude_mt_genes" %in% incl_res_dge_params && 
   dge_params$exclude_mt_genes) {
  
  mt_genes <- gene_metadata %>%
    filter(str_detect(gene_symbol, "^MT-")) %>%
    pull(ensembl_gene_id_version)
  
  dds <- dds[!rownames(dds) %in% mt_genes,]
  
}

# DESeq2 workflow ---------------------------------------------------------

# Perform DGE
# # => Using Wald test
dds_wald <- DESeq(object = dds,
                  test = "Wald",
                  parallel = TRUE,
                  BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))

# => Using LRT
dds_lrt <- DESeq(object = dds,
                 test = "LRT",
                 full = design,
                 reduced = design[, which(colnames(design) != coef_nm),
                                  drop = F],
                 parallel = TRUE,
                 BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))

# Tidy coef_nm 
coef_nm_deseq2 <- coef_nm

for(i in c("-", "\\s", "\\+", "/")) {
  
  coef_nm_deseq2 <- gsub(pattern = i,
                         replacement = "\\.",
                         x = coef_nm_deseq2)
  
}

# Extract results table for each test
for(i in c("wald", "lrt")) {
  
  # Construct res
  assign(x = paste0("res_", i),
         value = results(object = eval(sym(paste0("dds_", i))),
                         name = coef_nm_deseq2,
                         alpha = dge_params$padj_dge,
                         cooksCutoff = FALSE,
                         independentFiltering = FALSE))
  
  # Construct res_shrunken_lfc
  assign(x = paste0("res_", i, "_shrunken_lfc"),
         value = lfcShrink(dds = eval(sym(paste0("dds_", i))),
                           coef = coef_nm_deseq2,
                           type = "apeglm",
                           parallel = T,
                           BPPARAM = BiocParallel::MulticoreParam(workers = n_workers)))
  
  # Construct lfc
  assign(x = paste0("lfc_deseq2_", i),
         value = tidy_DESeqResults(input = eval(sym(paste0("res_", i))),
                                   gene_id = "ensembl_gene_id_version") %>%
           full_join(tidy_DESeqResults(input = eval(sym(paste0("res_", i, "_shrunken_lfc"))),
                                       gene_id = "ensembl_gene_id_version") %>%
                       select(ensembl_gene_id_version,
                              shrunken_log2FoldChange = log2FoldChange,
                              shrunken_lfcSE = lfcSE),
                     by = join_by(ensembl_gene_id_version)) %>%
           left_join(unique(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")]),
                     by = join_by(ensembl_gene_id_version)))
  
  # Construct lfc_sign
  assign(x = paste0("lfc_sign_deseq2_", i),
         value = eval(sym(paste0("lfc_deseq2_", i))) %>%
           filter(padj < dge_params$padj_dge))
  
}

# limma-voom workflow -----------------------------------------------------

y <- DGEList(counts = counts(object = dds, normalized = FALSE))

y <- normLibSizes(object = y)

# limma-trend framework
fit_limma_trend <- lmFit(object = edgeR::cpm(y = y, log = TRUE, prior.count = 3),
                         design = design)

fit_limma_trend <- eBayes(fit = fit_limma_trend, trend = TRUE)

# limma-voom framework
fit_limma_voom <- voomLmFit(counts = y,
                            design = design,
                            sample.weights = TRUE,
                            plot = FALSE)

fit_limma_voom <- eBayes(fit = fit_limma_voom)

for(i in c("trend", "voom")) {
  
  assign(x = paste0("lfc_limma_", i),
         value = topTable(fit = eval(sym(paste0("fit_limma_", i))), 
                          coef = coef_nm, 
                          number = Inf) %>% 
           rownames_to_column(var = "ensembl_gene_id_version") %>%
           mutate(log2FoldChange = logFC,
                  pvalue = P.Value,
                  padj = adj.P.Val,
                  sign_log2fc_times_minus_log10pvalue = sign(log2FoldChange) * -log(x = pvalue, base = 10)) %>%
           left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")],
                     by = join_by(ensembl_gene_id_version)))
  
  assign(x = paste0("lfc_sign_limma_", i),
         value = eval(sym(paste0("lfc_limma_", i))) %>%
           filter(padj < dge_params$padj_dge))
  
}

# edgeR workflow ----------------------------------------------------------

y <- DGEList(counts = counts(object = dds, normalized = FALSE))

y <- normLibSizes(object = y)

# GLM-QL framework
fit_edger_glmql <- glmQLFit(y = y, 
                            design = design, 
                            robust = TRUE)

res_edger_glmql <- glmQLFTest(glmfit = fit_edger_glmql, 
                              coef = coef_nm)

# GLM framework
y <- estimateDisp(y = y, 
                  design = design, 
                  robust = TRUE)

fit_edger_glm <- glmFit(y = y, 
                        design = design, 
                        robust = TRUE)

res_edger_glm <- glmLRT(glmfit = fit_edger_glm, 
                        coef = coef_nm)

for(i in c("glm", "glmql")) {
  
  assign(x = paste0("lfc_edger_", i),
         value = topTags(object = eval(sym(paste0("res_edger_", i))), 
                         n = Inf, 
                         sort.by = "PValue")$table %>% 
           rownames_to_column(var = "ensembl_gene_id_version") %>%
           mutate(log2FoldChange = logFC,
                  pvalue = PValue,
                  padj = FDR,
                  sign_log2fc_times_minus_log10pvalue = sign(log2FoldChange) * -log(x = pvalue, base = 10)) %>%
           left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")],
                     by = join_by(ensembl_gene_id_version)))
  
  assign(x = paste0("lfc_sign_edger_", i),
         value = eval(sym(paste0("lfc_edger_", i))) %>%
           filter(padj < dge_params$padj_dge))
  
}

# Stream output analysis parameter values ---------------------------------

write_param_log(param_list = dge_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = dge_params,
                                        in_dir_path = out_dir_path)

# Write results to files --------------------------------------------------

write_rds(x = sample_metadata,
          file = file.path(out_dir_path, paste0("sample_metadata_", output_param_log_id, ".rds")))

write_parquet(x = sample_metadata,
              sink = file.path(out_dir_path, paste0("sample_metadata_", output_param_log_id, ".parquet")))

fwrite(x = as_tibble(txi$counts, rownames = "ensembl_gene_id_version"),
       file = file.path(out_dir_path, paste0("cnt_", output_param_log_id, ".csv")))

fwrite(x = as_tibble(txi$abundance, rownames = "ensembl_gene_id_version"),
       file = file.path(out_dir_path, paste0("tpm_", output_param_log_id, ".csv")))

test_names <- c("deseq2_wald", "deseq2_lrt", 
                "limma_trend", "limma_voom", 
                "edger_glm", "edger_glmql")

walk(.x = test_names,
     .f = function(x) {
       
       write_csv(x = eval(sym(paste0("lfc_", x))),
                 file = file.path(out_dir_path, paste0("lfc_", x, "_", output_param_log_id, ".csv")))
       
       write_csv(x = eval(sym(paste0("lfc_sign_", x))),
                 file = file.path(out_dir_path, paste0("lfc_sign_", x, "_", output_param_log_id, ".csv")))
       
     })



