# File description --------------------------------------------------------

### Goal
# random_effects_metaanalysis.R computes a precision-weighted average of SLC18A2 log2 fold change across publicly available snRNA-seq datasets.

### Execution steps
# ===>>> public/random_effects_metaanalysis.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - DESeq2 DGEA results for Kamath 2022 (PMID: 35513515) / scripts & results available upon request
# - DESeq2 DGEA results for Smajic 2021 (PMID: 34919646) / scripts & results available upon request
# - DESeq2 DGEA results for Yaghmaeian Salmani 2024 (PMID: 38587883) / scripts & results available upon request
# - DESeq2 DGEA results for Fushiki 2026 (DOI: https://doi.org/10.1101/2024.12.20.629776) / scripts & results available upon request
# - paramstream package available from GitHub with 'pak::pak("benvallin/paramstream")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (ben.vallin@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Cramb_2026/data/random_effects_metaanalysis"

# Path to directory containing Kamath DGE results
dge_kamath_dir_path <- "data/public/Kamath_2022_pmid_35513515/dge_pseudobulk_deseq2"

# Path to directory containing Smajic DGE results
dge_smajic_dir_path <- "data/public/Smajic_2021_pmid_34919646/dge_pseudobulk_deseq2"

# Path to directory containing Yaghmaeian Salmani DGE results
dge_yaghmaeian_salmani_dir_path <- "data/public/Yaghmaeian_Salmani_2024_pmid_38587883/dge_pseudobulk_deseq2"

# Path to directory containing Fushiki DGE results
dge_fushiki_dir_path <- "data/public/Fushiki_2026/dge_pseudobulk_deseq2"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(paramstream)
library(ggpubr)
library(metafor)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define functions --------------------------------------------------------

# run_rma()
run_rma <- function(data, yi = "log2FoldChange", sei = "lfcSE", mods = NULL, comp_nm, exp_var_nm, slab) {
  
  if(is.null(slab)) {
    
    slab <- paste0("dataset", seq(nrow(data)))
    
  }
  
  if(is.null(mods)) {
    
    fit <- rma(yi = data[[yi]],
               sei = data[[sei]],
               data = data,
               slab = slab,
               method = "REML")
    
    coef <- broom::tidy(fit, conf.int = TRUE) %>% 
      mutate(k = fit$k,
             tau2 = fit$tau2,
             I2 = fit$I2,
             QEp = fit$QEp,
             comp_nm = comp_nm,
             exp_var_nm = exp_var_nm) %>% 
      select(comp_nm, exp_var_nm, everything())
    
    return(list(coefficients = coef,
                fit = fit,
                predictions = NULL))
    
  } else {
    
    fit <- rma(yi = data[[yi]],
               sei = data[[sei]],
               mods = as.formula(mods),
               data = data,
               slab = slab,
               method = "REML")
    
    coef <- broom::tidy(fit, conf.int = TRUE) %>% 
      mutate(k = fit$k,
             tau2_residual = fit$tau2,
             I2_residual = fit$I2,
             QMp = fit$QMp,
             QM = fit$QM,
             QEp = fit$QEp,
             comp_nm = comp_nm,
             exp_var_nm = exp_var_nm) %>% 
      select(comp_nm, exp_var_nm, everything())
    
    mod_vars <- all.vars(as.formula(mods))
    
    mod_grid <- data %>%
      select(all_of(mod_vars)) %>%
      mutate(across(.cols = everything(), .fns = factor)) %>%  
      distinct() %>%
      arrange(across(.cols = everything())) 
    
    newmods_matrix <- model.matrix(as.formula(mods), data = mod_grid)[, -1, drop = FALSE]
    
    pred <- as_tibble(predict(fit, newmods = newmods_matrix)) %>%
      bind_cols(mod_grid) %>%
      mutate(comp_nm = comp_nm,
             exp_var_nm = exp_var_nm) %>%
      select(comp_nm, exp_var_nm, all_of(mod_vars), everything())
    
    return(list(coefficients = coef,
                fit = fit,
                predictions = pred))
  }
  
}

# get_rma_estimates()
extract_estimates <- function(rma_res, comp_nm, exp_var_nm) {
  
  fit <- rma_res$fit
  coef <- rma_res$coefficients
  
  res_datasets <- tibble(comp_nm = comp_nm,
                         exp_var_nm = exp_var_nm,
                         dataset = fit$slab,
                         type = "study",
                         estimate = as.numeric(fit$yi),
                         se = sqrt(fit$vi),
                         conf.low = as.numeric(fit$yi) - 1.96 * sqrt(fit$vi),
                         conf.high = as.numeric(fit$yi) + 1.96 * sqrt(fit$vi))
  
  res_pooled <- coef %>% 
    mutate(dataset = "Pooled", type = "summary") %>% 
    select(comp_nm, exp_var_nm, dataset, type, estimate, se = std.error, conf.low, conf.high)
  
  bind_rows(res_datasets, res_pooled)
  
}

# plot_rma_estimates()
plot_rma_estimates <- function(input, 
                               test_cond = "test cond.", 
                               ref_cond = "ref cond.", 
                               title = NULL, 
                               xintercepts = NULL, 
                               fill_colors = c( "#0b8c1f", "#ff6d00"),
                               padj_lab = TRUE,
                               lab_nudge = 0.08) {
  
  if((!is.null(xintercepts) && !is.numeric(xintercepts)) ||
     (is.numeric(xintercepts) && (length(xintercepts) == 0L ||
                                  any(is.na(xintercepts))))) {
    stop("Invalid xintercepts argument.",
         call. = F)
  }
  
  datasets <- unique(input$dataset)
  
  dataset_levels <- c("Pooled", rev(setdiff(datasets, "Pooled")))
  
  df <- input %>% 
    mutate(dataset = fct_relevel(dataset, dataset_levels),
           padj_lab = if("padj" %in% colnames(input)) { formatC(.data$padj, format = "e", digits = 1) } else { "" },
           estimate_ci_lab = paste0(round(.data$estimate, 2), " [", round(.data$conf.low, 2), ", ", round(.data$conf.high, 2), "]"),
           lab_val = if(.env$padj_lab & "padj" %in% colnames(input)) { paste0(.data$estimate_ci_lab, " | ", .data$padj_lab) } else { .data$estimate_ci_lab },
           lab_pos = estimate + se + lab_nudge)
  
  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = .data$estimate,
                                              y = .data$dataset,
                                              fill = .data$type)) +
    # annotate(geom = "rect",
    #          xmin = -Inf,
    #          xmax = Inf,
    #          ymin = datasets[[1]],
    #          ymax = datasets[[length(datasets)-1]],
    #          alpha = 0.3, fill = "#5d76cb") +
    annotate(geom = "rect",
             xmin = -Inf,
             xmax = 0,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.3, fill = "#f0396a") +
    annotate(geom = "rect",
             xmin = 0,
             xmax = Inf,
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.3, fill = "#0077b6") +
    ggplot2::geom_hline(yintercept = df$dataset,
                        color = "grey",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed",
                        color = "#c9184a",
                        linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = xintercepts,
                        linetype = "dashed",
                        color = "orange",
                        linewidth = 0.5) +
    ggplot2::geom_linerange(mapping = ggplot2::aes(xmin = .data$estimate - .data$se,
                                                   xmax = .data$estimate + .data$se)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$type),
                        color = "black",
                        size = 5,
                        stroke = 1) +
    ggplot2::geom_text(mapping = ggplot2::aes(x = .data$lab_pos,
                                              label = .data$lab_val),
                       hjust = 0) +
    ggpubr::theme_pubr(legend = "bottom") +
    ggplot2::theme(plot.title = element_text(size = 14, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   strip.text.x = ggplot2::element_text(size = 14)) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_shape_manual(values = c(21, 23)) +
    ggplot2::labs(title = title,
                  x = bquote(log[2] * " fold change (" * .(test_cond) * " vs " * .(ref_cond) * ")"),
                  y = NULL,
                  fill = NULL,
                  shape = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 5)))
  
  if(!is.null(xintercepts)) {
    
    p <- p +
      ggplot2::scale_x_continuous(breaks = c(0, xintercepts))
    
  }
  
  p
  
}

# Define DGE analysis parameters ------------------------------------------

## Kamath
# PD vs CTRL DaNs
dge_params_kamath = list(
  broad_cell_type = "DaN",
  cell_type = NULL,
  disease_status = c("CTRL", "PD"),
  input_pseudobulks = "pb_broad_cell_type",
  coef_var = "disease_status",
  coef_ref_lvl = "CTRL",
  coef_test_lvl = "PD",
  model_formula = "~ donor_sex + disease_status",
  protein_coding_only = TRUE,
  min_cnt_incl = 0,
  padj_dge = 0.1
)

## Smajic
# PD vs CTRL DaNs
dge_params_smajic = list(
  cell_type = "DaNs",
  disease_status = c("CTRL", "PD"),
  input_pseudobulks = "pb_cell_type",
  coef_var = "disease_status",
  coef_ref_lvl = "CTRL",
  coef_test_lvl = "PD",
  model_formula = "~ donor_sex + disease_status",
  protein_coding_only = TRUE,
  min_cnt_incl = 0,
  padj_dge = 0.1
)

## Yaghmaeian Salmani
# Lesioned vs intact DaNs
dge_params_yaghmaeian_salmani = list(
  input_pseudobulks = "pb_cell_type",
  coef_var = "condition",
  coef_ref_lvl = "intact",
  coef_test_lvl = "lesion",
  model_formula = "~ animal_id + condition",
  protein_coding_only = TRUE,
  min_cnt_incl = 0,
  padj_dge = 0.1
)

## Fushiki
# MitoPark vs control DaNs
dge_params_fushiki = list(
  genotype = c("MitoPark littermate controls", "MitoPark mice"),
  age = c("8 weeks old",  "16 weeks old", "24 weeks old"),
  cell_type = c("DaN"),
  input_pseudobulks = "pb_cell_type",
  coef_var = "genotype",
  coef_ref_lvl = "MitoPark littermate controls",
  coef_test_lvl = "MitoPark mice",
  # model_formula = "~ sex + genotype",
  model_formula = "~ sex + age + genotype",
  protein_coding_only = TRUE,
  min_cnt_incl = 0,
  padj_dge = 0.1
)


# Stream analysis parameter values ----------------------------------------

dge_param_log_id_kamath <- get_param_log_id(param_list = dge_params_kamath,
                                            in_dir_path = dge_kamath_dir_path)

dge_param_log_id_smajic <- get_param_log_id(param_list = dge_params_smajic,
                                            in_dir_path = dge_smajic_dir_path)

dge_param_log_id_yaghmaeian_salmani <- get_param_log_id(param_list = dge_params_yaghmaeian_salmani,
                                                        in_dir_path = dge_yaghmaeian_salmani_dir_path)

dge_param_log_id_fushiki <- get_param_log_id(param_list = dge_params_fushiki,
                                             in_dir_path = dge_fushiki_dir_path)

# Load DGE results --------------------------------------------------------

dataset_summary <- tibble(
  dataset = as_factor(c("kamath", 
                        "smajic",
                        "yaghmaeian_salmani", 
                        "fushiki")),
  dataset_pretty = as_factor(c("Kamath et al., 2022",
                               "Smajic et al., 2021",
                               "Yaghmaeian Salmani et al., 2024", 
                               "Fushiki et al., 2026"))) %>% 
  mutate(species = case_when(dataset %in% c("kamath", "smajic", "tiklova") ~ "human",
                             dataset %in% c("yaghmaeian_salmani", "fushiki") ~"mouse",
                             TRUE ~ NA_character_),
         tech = case_when(dataset %in% c("kamath", "yaghmaeian_salmani", "fushiki") ~ "FANS-snRNA-seq",
                          dataset %in% c("smajic") ~ "snRNA-seq",
                          TRUE ~ NA_character_))

sample_metadata <- list()  
lfc_data <- list()

for(i in dataset_summary$dataset) {
  
  dge_param_log_id <- paste0("dge_param_log_id_", i)
  
  if(exists(x = dge_param_log_id)) {
    
    # Load sample metadata
    sample_metadata[[i]] <- read_rds(file.path(eval(as.symbol(paste0("dge_", i, "_dir_path"))), 
                                               paste0("pseudobulk_metadata_", eval(as.symbol(dge_param_log_id)), 
                                                      ".rds")))
    
    # Load DGE results
    lfc_data[[i]] <- read_csv(file.path(eval(as.symbol(paste0("dge_", i, "_dir_path"))),
                                        paste0("lfc_deseq2_wald_", eval(as.symbol(dge_param_log_id)), 
                                               ".csv"))) %>% 
      mutate(gene_symbol = toupper(gene_symbol),
             dataset = i,
             comparison = "diseased vs control DaNs") %>% 
      filter(gene_symbol == "SLC18A2") %>% 
      select(dataset, comparison, gene_symbol, log2FoldChange, lfcSE, padj)
    
  }
}

lfc_meta <- bind_rows(lfc_data) %>% 
  left_join(dataset_summary,
            by = join_by(dataset)) %>% 
  dplyr::arrange(dataset_pretty) 

# Meta-Analysis via Linear (Mixed-Effects) Models -------------------------

meta_datasets <- c("kamath", 
                   "smajic",
                   "yaghmaeian_salmani", 
                   "fushiki")

rma_res <- list()
rma_mod_res <- list()

for(i in unique(lfc_meta$gene_symbol)) {
  
  data <- lfc_meta %>% 
    filter(gene_symbol == i,
           dataset %in% meta_datasets)
  
  slab <- unique(data$dataset_pretty)
  
  rma_res[[i]] <- run_rma(data = data, 
                          comp_nm = "diseased vs control DaNs", 
                          exp_var_nm = i,
                          slab = slab)
  
  rma_mod_res[[i]] <- run_rma(data = data, 
                              mods = "~ species",
                              comp_nm = "diseased vs control DaNs", 
                              exp_var_nm = i,
                              slab = slab)
  
}

rma_coef <- lapply(X = names(rma_res),
                   FUN = function(x) rma_res[[x]]$coefficients) %>% 
  setNames(nm = names(rma_res)) %>% 
  bind_rows() %>% 
  mutate(padj = p.adjust(p = p.value, method = "BH"),
         log2FoldChange = estimate,
         lfcSE = std.error)

rma_coef_sign <- rma_coef %>% 
  filter(padj < 0.05)

df_padj <- bind_rows(lfc_meta %>% 
                       select(dataset = dataset_pretty, 
                              comp_nm = comparison, 
                              exp_var_nm = gene_symbol, 
                              padj),
                     rma_coef %>% 
                       mutate(dataset = "Pooled") %>% 
                       select(dataset, 
                              comp_nm, 
                              exp_var_nm, 
                              padj))

rma_estimates <- list()

for(i in unique(lfc_meta$gene_symbol)) {
  
  rma_estimates[[i]] <- extract_estimates(rma_res = rma_res[[i]],
                                          comp_nm = "diseased vs control DaNs", 
                                          exp_var_nm = i) %>% 
    left_join(df_padj,
              by = join_by(comp_nm, exp_var_nm, dataset))
  
}

rma_mod_pred <- lapply(X = names(rma_mod_res),
                       FUN = function(x) rma_mod_res[[x]]$predictions) %>% 
  setNames(nm = names(rma_mod_res)) %>% 
  bind_rows()

rma_mod_coef <- lapply(X = names(rma_mod_res),
                       FUN = function(x) rma_mod_res[[x]]$coefficients) %>% 
  setNames(nm = names(rma_mod_res)) %>% 
  bind_rows() %>% 
  group_by(comp_nm, term) %>% 
  mutate(padj = p.adjust(p = p.value, method = "BH")) %>% 
  ungroup() %>% 
  nest(data = c(term:k, padj)) %>% 
  group_by(comp_nm) %>% 
  mutate(QMp_adj = p.adjust(p = QMp, method = "BH")) %>% 
  ungroup() %>% 
  unnest(cols = c(data))

rma_mod_coef_sign <- rma_mod_coef %>% 
  filter(padj < 0.05)

rma_mod_coef_QMp_sign <- rma_mod_coef %>% 
  filter(QMp_adj < 0.05)

tmp <- lfc_meta |> 
  select(dataset = dataset_pretty,
         species, tech,
         comparison,
         everything(),
         - dataset)

write_csv(tmp, file.path(out_dir_path, "lfc_meta.csv"))

write_csv(rma_coef, file.path(out_dir_path, "rma_coef.csv"))

write_csv(rma_mod_coef_sign, file.path(out_dir_path, "rma_mod_coef_sign.csv"))

# Summary plots -----------------------------------------------------------

# SLC18A2 (VMAT2)
plot_rma_estimates(input = rma_estimates$SLC18A2, 
                   title = "SLC18A2",
                   ref_cond = "CTRL DaNs",
                   test_cond = "diseased",
                   xintercepts = c(seq(-2, -1, 1),  seq(1, 2, 1)))

ggsave(filename = file.path(out_dir_path, "lfc_SLC18A2_diseased_vs_ctrl_dans.pdf"),
       width = 8,
       height = 3+0.22*length(rma_estimates$SLC18A2$exp_var_nm),
       device = cairo_pdf)


