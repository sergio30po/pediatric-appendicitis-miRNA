# ==============================================================================
# 05_DE_validation_edgeR_DESeq2.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Purpose:
#   Validate the RPM-based non-parametric differential-expression layer using
#   two independent count-based frameworks:
#
#       edgeR with TMM normalization and quasi-likelihood GLMs
#       DESeq2 with native median-of-ratios normalization
#
#   All miRNAs frozen upstream by script 01 are tested. No historical
#   candidate list is used to restrict or prioritize features.
#
# Clinical variables inherited from script 01:
#   Appendicitis_status = NIA / Appendicitis
#   Phenotype           = NIA / Mild / Severe
#
# Models:
#   Primary:
#       unadjusted
#
#   Sensitivity:
#       adjusted for standardized age + sex
#
#   Batch is not included in the frozen script-05 model set.
#
# Prespecified analyses:
#   Appendicitis_vs_NIA
#   Phenotype_global
#   NIA_vs_Mild
#   NIA_vs_Severe
#   Mild_vs_Severe
#
# Multiple testing:
#   BH-FDR is calculated independently for each
#   method x model x analysis.
#
#   Holm is retained as a conservative audit.
#   Raw p < 0.05 is exploratory.
#
# Feature-set harmonization:
#   - no additional edgeR::filterByExpr() filtering;
#   - DESeq2 independent filtering disabled;
#   - RPM, edgeR and DESeq2 therefore operate on the same upstream-frozen
#     filtered feature set.
#
# Data-driven candidate tracking:
#   The primary candidate recurrence summary is based on the unadjusted
#   edgeR + DESeq2 analyses only.
#
#   Age+sex-adjusted results are summarized separately as sensitivity evidence
#   and do not inflate the primary recurrence count.
#
# Interpretation:
#   Candidate PCA/heatmaps are post-selection visualizations only.
#   They are not independent evidence of group separation.
#
# ==============================================================================

source(here::here("R", "00_analysis_config.R"))
source(here::here("R", "00_utils.R"))

required_packages <- c(
  "edgeR",
  "DESeq2",
  "dplyr",
  "ggplot2",
  "ComplexHeatmap",
  "circlize"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing package(s) required by script 05: ",
    paste(missing_packages, collapse = ", "),
    ". Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}

context <- get_analysis_context()
results_dir <- context$results_dir
figures_dir <- context$figures_dir

validation_results_dir <- file.path(
  results_dir,
  "de_val"
)

validation_figures_dir <- file.path(
  figures_dir,
  "de_val"
)

volcano_dir <- file.path(
  validation_figures_dir,
  "vol"
)

heatmap_dir <- file.path(
  validation_figures_dir,
  "hm"
)

candidate_pca_dir <- file.path(
  validation_figures_dir,
  "pca"
)

for (directory in c(
  validation_results_dir,
  validation_figures_dir,
  volcano_dir,
  heatmap_dir,
  candidate_pca_dir
)) {
  ensure_dir(directory)
}

message(
  "\n============================================================"
)
message(
  "05_DE_validation_edgeR_DESeq2.R"
)
message(
  "Analysis ID: ",
  context$analysis_id
)
message(
  "Role: ",
  context$role
)
message(
  "Expected n: ",
  context$expected_n
)
message(
  "Count-based validation: edgeR-TMM + DESeq2"
)
message(
  "Primary recurrence: unadjusted models only"
)
message(
  "============================================================\n"
)

model_code <- function(model_name) {
  codes <- c(
    unadjusted = "unadj",
    adjusted_age_sex = "age_sex"
  )

  code <- unname(codes[model_name])

  if (length(code) != 1L || is.na(code)) {
    stop(
      "Unknown model name: ",
      model_name,
      call. = FALSE
    )
  }

  code
}

method_code <- function(method_name) {
  codes <- c(
    edgeR_TMM = "edger",
    DESeq2 = "deseq2"
  )

  code <- unname(codes[method_name])

  if (length(code) != 1L || is.na(code)) {
    stop(
      "Unknown method name: ",
      method_name,
      call. = FALSE
    )
  }

  code
}

contrast_code <- function(contrast_name) {
  codes <- c(
    Appendicitis_vs_NIA = "AvN",
    Phenotype_global = "Pg",
    NIA_vs_Mild = "NvM",
    NIA_vs_Severe = "NvS",
    Mild_vs_Severe = "MvS"
  )

  code <- unname(codes[contrast_name])

  if (length(code) != 1L || is.na(code)) {
    stop(
      "Unknown contrast name: ",
      contrast_name,
      call. = FALSE
    )
  }

  code
}

output_path_dictionary <- rbind(
  data.frame(
    type = "model",
    code = c("unadj", "age_sex"),
    label = c(
      "unadjusted",
      "adjusted_age_sex"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    type = "method",
    code = c("edger", "deseq2"),
    label = c("edgeR_TMM", "DESeq2"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    type = "contrast",
    code = c("AvN", "Pg", "NvM", "NvS", "MvS"),
    label = c(
      "Appendicitis_vs_NIA",
      "Phenotype_global",
      "NIA_vs_Mild",
      "NIA_vs_Severe",
      "Mild_vs_Severe"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    type = "candidate",
    code = c("raw", "fdr", "holm"),
    label = c(
      "raw p < 0.05",
      "BH-FDR < 0.05",
      "Holm < 0.05"
    ),
    stringsAsFactors = FALSE
  )
)

utils::write.csv(
  output_path_dictionary,
  file.path(
    validation_results_dir,
    "output_path_dictionary.csv"
  ),
  row.names = FALSE
)

required <- c(
  counts = file.path(results_dir, "raw_ftd_counts.rds"),
  metadata = file.path(results_dir, "metadata_analysis_samples.rds")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop("Run script 01 first. Missing:\n", paste("-", missing, collapse = "\n"))

counts <- as.matrix(readRDS(required[["counts"]]))
metadata <- as.data.frame(readRDS(required[["metadata"]]))
storage.mode(counts) <- "numeric"

if (any(abs(counts - round(counts)) > 1e-8)) stop("Count matrix contains non-integer values.")
counts <- round(counts)
storage.mode(counts) <- "integer"

if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("Count matrix must contain miRNA row names and sample column names.")
}
if (anyDuplicated(rownames(counts)) > 0) {
  stop("Duplicated miRNA identifiers detected in the count matrix.")
}
if (anyDuplicated(colnames(counts)) > 0) {
  stop("Duplicated sample identifiers detected in the count matrix.")
}
if (anyDuplicated(metadata$Sample) > 0) {
  stop("Duplicated sample identifiers detected in metadata.")
}
if (!setequal(colnames(counts), metadata$Sample)) {
  stop("Count matrix and metadata contain different samples.")
}
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
  stop("Count matrix contains missing, non-finite or negative values.")
}

metadata <- metadata[match(colnames(counts), metadata$Sample), , drop = FALSE]
stopifnot(identical(metadata$Sample, colnames(counts)))

message(
  "Loaded and aligned: ",
  nrow(counts),
  " miRNAs x ",
  ncol(counts),
  " samples."
)


required_metadata_columns <- c(
  "Sample",
  "Phenotype",
  "Phenotype_full",
  "Appendicitis_status",
  "Sex",
  "Age"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Missing metadata columns generated by script 01: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata$Phenotype <- factor(
  as.character(metadata$Phenotype),
  levels = c("NIA", "Mild", "Severe")
)

metadata$Phenotype_full <- factor(
  as.character(metadata$Phenotype_full),
  levels = c(
    "Non-inflamed appendix",
    "Mild appendicitis",
    "Severe appendicitis"
  )
)

metadata$Appendicitis_status <- factor(
  as.character(metadata$Appendicitis_status),
  levels = c("NIA", "Appendicitis")
)

metadata$Sex <- droplevels(
  factor(as.character(metadata$Sex))
)

metadata$Age <- as.numeric(metadata$Age)

if (
  anyNA(metadata$Phenotype) ||
  anyNA(metadata$Phenotype_full) ||
  anyNA(metadata$Appendicitis_status) ||
  anyNA(metadata$Age) ||
  anyNA(metadata$Sex)
) {
  stop(
    "Phenotype, Phenotype_full, Appendicitis_status, Age or Sex contains ",
    "missing/unexpected values in the selected cohort.",
    call. = FALSE
  )
}

expected_full <- c(
  NIA = "Non-inflamed appendix",
  Mild = "Mild appendicitis",
  Severe = "Severe appendicitis"
)

if (!identical(
  as.character(metadata$Phenotype_full),
  unname(expected_full[as.character(metadata$Phenotype)])
)) {
  stop(
    "Phenotype_full is inconsistent with Phenotype. Re-run script 01.",
    call. = FALSE
  )
}

expected_status <- ifelse(
  as.character(metadata$Phenotype) == "NIA",
  "NIA",
  "Appendicitis"
)

if (!identical(
  as.character(metadata$Appendicitis_status),
  expected_status
)) {
  stop(
    "Appendicitis_status is inconsistent with Phenotype. Re-run script 01.",
    call. = FALSE
  )
}

if (
  nrow(counts) == 0L ||
  ncol(counts) != nrow(metadata)
) {
  stop(
    "Count matrix is empty or inconsistent with metadata.",
    call. = FALSE
  )
}

message(
  "Appendicitis status: ",
  paste(
    names(
      table(
        metadata$Appendicitis_status
      )
    ),
    as.integer(
      table(
        metadata$Appendicitis_status
      )
    ),
    collapse = "; "
  )
)

message(
  "Clinical groups: ",
  paste(
    names(
      table(
        metadata$Phenotype
      )
    ),
    as.integer(
      table(
        metadata$Phenotype
      )
    ),
    collapse = "; "
  )
)

# Standardize age for adjusted GLMs. This is a linear reparameterization of the
# same continuous covariate and improves numerical conditioning in DESeq2.
age_mean <- mean(metadata$Age)
age_sd <- stats::sd(metadata$Age)
if (!is.finite(age_sd) || age_sd <= 0) stop("Age has zero or invalid variance; adjusted models cannot be fitted.")
metadata$Age_z <- (metadata$Age - age_mean) / age_sd

age_scaling <- data.frame(
  variable = "Age",
  mean = age_mean,
  sd = age_sd,
  transformed_variable = "Age_z",
  stringsAsFactors = FALSE
)
utils::write.csv(
  age_scaling,
  file.path(validation_results_dir, "age_scaling.csv"),
  row.names = FALSE
)

analysis_manifest <- data.frame(
  analysis_id = context$analysis_id,
  n_samples = ncol(counts),
  n_miRNAs_tested = nrow(counts),
  edgeR_version = as.character(utils::packageVersion("edgeR")),
  DESeq2_version = as.character(utils::packageVersion("DESeq2")),
  edgeR_normalization = "TMM",
  DESeq2_normalization = "native median-of-ratios size factors",
  additional_edgeR_filterByExpr = FALSE,
  DESeq2_independentFiltering = FALSE,
  primary_model = "unadjusted",
  adjusted_model_role = "age + sex sensitivity analysis",
  clinical_group_variable = "Phenotype",
  binary_status_variable = "Appendicitis_status",
  NIA_definition = "Non-inflamed appendix; not a healthy-control cohort",
  primary_multiple_testing = "Benjamini-Hochberg FDR",
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_manifest,
  file.path(validation_results_dir, "manifest.csv"),
  row.names = FALSE
)

utils::write.csv(
  data.frame(miRNA = rownames(counts), stringsAsFactors = FALSE),
  file.path(validation_results_dir, "tested_miRNAs.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 1. Model specifications and audit
# ------------------------------------------------------------------------------

model_specs <- list(
  unadjusted = list(
    condition_terms = c("Appendicitis_status"),
    severity_terms = c("Phenotype"),
    description = "Primary unadjusted model"
  )
)

if (ADJUST_FOR_AGE_SEX) {
  model_specs$adjusted_age_sex <- list(
    condition_terms = c("Age_z", "Sex", "Appendicitis_status"),
    severity_terms = c("Age_z", "Sex", "Phenotype"),
    description = "Sensitivity model adjusted for age and sex"
  )
}

if (!is.na(BATCH_COLUMN)) {
  message(
    "BATCH_COLUMN is configured as '",
    BATCH_COLUMN,
    "', but batch adjustment is not part of the frozen script-05 model set. ",
    "Only unadjusted and age+sex models are run here."
  )
}

formula_from_terms <- function(terms) {
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

is_estimable <- function(formula, data) {
  mm <- tryCatch(
    stats::model.matrix(formula, data = data),
    error = function(e) e
  )
  if (inherits(mm, "error")) {
    return(list(
      estimable = FALSE,
      rank = NA_integer_,
      n_columns = NA_integer_,
      columns = "",
      error = conditionMessage(mm)
    ))
  }

  list(
    estimable = qr(mm)$rank == ncol(mm),
    rank = qr(mm)$rank,
    n_columns = ncol(mm),
    columns = paste(colnames(mm), collapse = ";"),
    error = ""
  )
}

audit_rows <- list()
valid_models <- list()
for (model_name in names(model_specs)) {
  spec <- model_specs[[model_name]]
  cond_formula <- formula_from_terms(spec$condition_terms)
  sev_formula <- formula_from_terms(spec$severity_terms)
  ca <- is_estimable(cond_formula, metadata)
  sa <- is_estimable(sev_formula, metadata)

  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    model = model_name,
    description = spec$description,
    status_formula = paste(deparse(cond_formula), collapse = ""),
    severity_formula = paste(deparse(sev_formula), collapse = ""),
    status_estimable = ca$estimable,
    severity_estimable = sa$estimable,
    status_rank = ca$rank,
    status_columns = ca$n_columns,
    severity_rank = sa$rank,
    severity_columns = sa$n_columns,
    status_model_error = ca$error,
    severity_model_error = sa$error,
    stringsAsFactors = FALSE
  )

  if (ca$estimable && sa$estimable) {
    valid_models[[model_name]] <- spec
  } else {
    warning("Skipping non-estimable model: ", model_name)
  }
}
model_audit <- dplyr::bind_rows(audit_rows)
utils::write.csv(model_audit, file.path(validation_results_dir, "model_audit.csv"), row.names = FALSE)

if (length(valid_models) == 0) stop("No estimable models available.")

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

contrast_effect_definition <- function(contrast) {
  definitions <- c(
    Appendicitis_vs_NIA = "Positive log2FC = higher in Appendicitis than NIA",
    NIA_vs_Mild = "Positive log2FC = higher in Mild than NIA",
    NIA_vs_Severe = "Positive log2FC = higher in Severe than NIA",
    Mild_vs_Severe = "Positive log2FC = higher in Severe than Mild",
    Phenotype_global = "Global omnibus test; no single log2FC direction"
  )
  unname(definitions[contrast])
}

contrast_direction <- function(log2FC, contrast) {
  if (contrast == "Phenotype_global") {
    return(rep(NA_character_, length(log2FC)))
  }

  positive_label <- switch(
    contrast,
    Appendicitis_vs_NIA = "Higher_in_Appendicitis",
    NIA_vs_Mild = "Higher_in_Mild",
    NIA_vs_Severe = "Higher_in_Severe",
    Mild_vs_Severe = "Higher_in_Severe"
  )

  negative_label <- switch(
    contrast,
    Appendicitis_vs_NIA = "Higher_in_NIA",
    NIA_vs_Mild = "Higher_in_NIA",
    NIA_vs_Severe = "Higher_in_NIA",
    Mild_vs_Severe = "Higher_in_Mild"
  )

  dplyr::case_when(
    is.na(log2FC) ~ NA_character_,
    log2FC > 0 ~ positive_label,
    log2FC < 0 ~ negative_label,
    TRUE ~ "No_change"
  )
}

standardize_table <- function(tab, method, model, contrast, global = FALSE) {
  if (!"miRNA" %in% colnames(tab)) stop("Internal error: miRNA column missing.")
  if (!"p_value" %in% colnames(tab)) stop("Internal error: p_value column missing.")

  tab$FDR <- stats::p.adjust(tab$p_value, method = "BH")
  tab$Holm <- stats::p.adjust(tab$p_value, method = "holm")
  tab$raw_significant <- !is.na(tab$p_value) & tab$p_value < ALPHA_RAW
  tab$fdr_significant <- !is.na(tab$FDR) & tab$FDR < ALPHA_FDR
  tab$holm_significant <- !is.na(tab$Holm) & tab$Holm < ALPHA_FDR
  tab$method <- method
  tab$model <- model
  tab$contrast <- contrast
  tab$global_test <- global
  if (!"log2FC" %in% colnames(tab)) tab$log2FC <- NA_real_
  tab$effect_definition <- contrast_effect_definition(contrast)
  tab$direction <- contrast_direction(tab$log2FC, contrast)
  tab
}

write_table <- function(
    tab,
    model_name,
    method_dir,
    filename = NULL
) {
  d <- file.path(
    validation_results_dir,
    model_code(model_name),
    method_code(method_dir),
    "tab"
  )

  ensure_dir(d)

  contrast_name <- unique(
    as.character(
      tab$contrast
    )
  )

  if (length(contrast_name) != 1L) {
    stop(
      "Result table must contain exactly one contrast.",
      call. = FALSE
    )
  }

  utils::write.csv(
    tab,
    file.path(
      d,
      paste0(
        contrast_code(contrast_name),
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

write_candidate_sets <- function(
    tab,
    model_name,
    method_dir,
    contrast
) {
  base_dir <- file.path(
    validation_results_dir,
    model_code(model_name),
    method_code(method_dir),
    "cand"
  )

  raw_dir <- file.path(
    base_dir,
    "raw"
  )

  fdr_dir <- file.path(
    base_dir,
    "fdr"
  )

  holm_dir <- file.path(
    base_dir,
    "holm"
  )

  for (d in c(
    raw_dir,
    fdr_dir,
    holm_dir
  )) {
    ensure_dir(d)
  }

  short_name <- paste0(
    contrast_code(contrast),
    ".csv"
  )

  utils::write.csv(
    tab[
      tab$raw_significant,
      ,
      drop = FALSE
    ],
    file.path(
      raw_dir,
      short_name
    ),
    row.names = FALSE
  )

  utils::write.csv(
    tab[
      tab$fdr_significant,
      ,
      drop = FALSE
    ],
    file.path(
      fdr_dir,
      short_name
    ),
    row.names = FALSE
  )

  utils::write.csv(
    tab[
      tab$holm_significant,
      ,
      drop = FALSE
    ],
    file.path(
      holm_dir,
      short_name
    ),
    row.names = FALSE
  )

  invisible(NULL)
}

plot_volcano <- function(tab, model_name, method_dir, contrast) {
  if (
    contrast == "Phenotype_global" ||
      all(is.na(tab$log2FC))
  ) {
    return(invisible(NULL))
  }

  plotting <- tab[
    is.finite(tab$log2FC) & !is.na(tab$p_value),
    ,
    drop = FALSE
  ]

  if (nrow(plotting) == 0) {
    return(invisible(NULL))
  }

  plotting$minus_log10_p <- -log10(
    pmax(plotting$p_value, .Machine$double.xmin)
  )

  plotting$significance_class <- dplyr::case_when(
    plotting$fdr_significant ~ "BH-FDR < 0.05",
    plotting$raw_significant ~ "Raw p < 0.05 only",
    TRUE ~ "Not significant"
  )

  labels <- plotting |>
    dplyr::filter(fdr_significant) |>
    dplyr::arrange(FDR, p_value) |>
    utils::head(12)

  p <- ggplot2::ggplot(
    plotting,
    ggplot2::aes(
      x = log2FC,
      y = minus_log10_p
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(ALPHA_RAW),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = significance_class,
        alpha = significance_class
      ),
      size = 2.2
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "BH-FDR < 0.05" = 16,
        "Raw p < 0.05 only" = 1,
        "Not significant" = 1
      ),
      drop = FALSE
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "BH-FDR < 0.05" = 1,
        "Raw p < 0.05 only" = 0.75,
        "Not significant" = 0.30
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = paste(method_dir, "—", contrast),
      subtitle = paste0(
        model_name,
        "; BH-FDR across ",
        nrow(tab),
        " tested miRNAs; no log2FC cutoff"
      ),
      x = "log2 fold-change",
      y = expression(-log[10](p)),
      shape = NULL,
      alpha = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (nrow(labels) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = labels,
        ggplot2::aes(label = miRNA),
        size = 3,
        check_overlap = TRUE,
        vjust = -0.6
      )
  }

  d <- file.path(
    volcano_dir,
    model_code(model_name),
    method_code(method_dir)
  )

  ensure_dir(d)

  short_name <- contrast_code(
    contrast
  )

  save_ggplot_pdf(
    plot = p,
    filename = file.path(
      d,
      paste0(
        short_name,
        ".pdf"
      )
    ),
    width = 7,
    height = 5.5
  )

  save_ggplot_tiff(
    plot = p,
    filename = file.path(
      d,
      paste0(
        short_name,
        ".tiff"
      )
    ),
    width = 7,
    height = 5.5,
    dpi = 600
  )

  invisible(p)
}


# ------------------------------------------------------------------------------
# Candidate visualization helpers
# ------------------------------------------------------------------------------

# Candidate PCA and candidate heatmaps are post-selection visualizations.
# They are generated separately for:
#   - exploratory raw p < 0.05 candidates
#   - BH-FDR < 0.05 candidates
#
# They must not be interpreted as independent evidence of group separation.

candidate_group_spec <- function(
    contrast
) {

  if (contrast == "Appendicitis_vs_NIA") {
    return(
      list(
        variable = "Appendicitis_status",
        levels = c("NIA", "Appendicitis"),
        labels = c("NIA", "Appendicitis")
      )
    )
  }

  if (contrast == "Phenotype_global") {
    return(
      list(
        variable = "Phenotype",
        levels = c("NIA", "Mild", "Severe"),
        labels = c("NIA", "Mild", "Severe")
      )
    )
  }

  phenotype_specs <- list(
    NIA_vs_Mild = list(
      levels = c("NIA", "Mild"),
      labels = c("NIA", "Mild")
    ),
    NIA_vs_Severe = list(
      levels = c("NIA", "Severe"),
      labels = c("NIA", "Severe")
    ),
    Mild_vs_Severe = list(
      levels = c("Mild", "Severe"),
      labels = c("Mild", "Severe")
    )
  )

  if (!contrast %in% names(phenotype_specs)) {
    stop(
      "Unknown contrast for visualization: ",
      contrast,
      call. = FALSE
    )
  }

  list(
    variable = "Phenotype",
    levels = phenotype_specs[[contrast]]$levels,
    labels = phenotype_specs[[contrast]]$labels
  )
}


candidate_output_info <- function(
    candidate_type = c(
      "raw",
      "fdr"
    )
) {

  candidate_type <- match.arg(
    candidate_type
  )

  if (candidate_type == "raw") {
    return(
      list(
        label = "raw p < 0.05",
        tag = "raw"
      )
    )
  }

  list(
    label = "BH-FDR < 0.05",
    tag = "fdr"
  )
}


candidate_mirnas_from_table <- function(
    tab,
    candidate_type = c(
      "raw",
      "fdr"
    )
) {

  candidate_type <- match.arg(
    candidate_type
  )

  if (candidate_type == "raw") {
    return(
      unique(
        tab$miRNA[
          tab$raw_significant
        ]
      )
    )
  }

  unique(
    tab$miRNA[
      tab$fdr_significant
    ]
  )
}


save_method_candidate_heatmap <- function(
    tab,
    expression_matrix,
    expression_label,
    model_name,
    method_dir,
    contrast,
    candidate_type = c(
      "raw",
      "fdr"
    )
) {

  candidate_type <- match.arg(
    candidate_type
  )

  info <- candidate_output_info(
    candidate_type
  )

  selected_miRNAs <- intersect(
    candidate_mirnas_from_table(
      tab,
      candidate_type
    ),
    rownames(
      expression_matrix
    )
  )

  if (length(
    selected_miRNAs
  ) == 0) {
    message(
      "Candidate heatmap skipped: ",
      method_dir,
      " / ",
      model_name,
      " / ",
      contrast,
      " [",
      info$label,
      "] — no candidates."
    )

    return(
      invisible(
        NULL
      )
    )
  }

  spec <- candidate_group_spec(
    contrast
  )

  keep <- as.character(
    metadata[[
      spec$variable
    ]]
  ) %in%
    spec$levels

  plot_metadata <- metadata[
    keep,
    ,
    drop = FALSE
  ]

  plot_metadata$Plot_group <- factor(
    as.character(
      plot_metadata[[
        spec$variable
      ]]
    ),
    levels = spec$levels,
    labels = spec$labels
  )

  plot_metadata <- plot_metadata |>
    dplyr::arrange(
      Plot_group,
      Sample
    )

  ordered_samples <- plot_metadata$Sample

  matrix_selected <- expression_matrix[
    selected_miRNAs,
    ordered_samples,
    drop = FALSE
  ]

  z_matrix <- t(
    scale(
      t(
        matrix_selected
      )
    )
  )

  z_matrix[
    !is.finite(
      z_matrix
    )
  ] <- 0

  annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = plot_metadata$Plot_group,
    annotation_name_side = "left"
  )

  heatmap_object <- ComplexHeatmap::Heatmap(
    z_matrix,
    name = "Row Z-score",
    col = circlize::colorRamp2(
      c(
        -2,
        0,
        2
      ),
      c(
        "navy",
        "white",
        "firebrick3"
      )
    ),
    top_annotation = annotation,
    column_split = plot_metadata$Plot_group,
    cluster_rows = length(
      selected_miRNAs
    ) > 1L,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(
      fontsize = 8
    ),
    column_title = paste0(
      method_dir,
      " — ",
      model_name,
      " — ",
      contrast,
      " — ",
      info$label
    ),
    use_raster = FALSE
  )

  output_dir <- file.path(
    heatmap_dir,
    model_code(model_name),
    method_code(method_dir),
    info$tag
  )

  ensure_dir(output_dir)

  figure_height <- max(
    4.5,
    min(
      11,
      0.28 * length(
        selected_miRNAs
      ) + 3
    )
  )

  short_name <- contrast_code(
    contrast
  )

  pdf_file <- file.path(
    output_dir,
    paste0(
      short_name,
      ".pdf"
    )
  )

  tiff_file <- file.path(
    output_dir,
    paste0(
      short_name,
      ".tiff"
    )
  )

  draw_heatmap <- function() {
    ComplexHeatmap::draw(
      heatmap_object,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  }

  draw_to_pdf(
    filename = pdf_file,
    width = 7,
    height = figure_height,
    draw = draw_heatmap
  )

  draw_to_tiff(
    filename = tiff_file,
    width = 7,
    height = figure_height,
    units = "in",
    res = 600,
    compression = "lzw",
    draw = draw_heatmap
  )

  manifest <- data.frame(
    model = model_name,
    method = method_dir,
    contrast = contrast,
    candidate_definition = info$label,
    expression_matrix = expression_label,
    miRNA = selected_miRNAs,
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    manifest,
    file.path(
      output_dir,
      paste0(
        short_name,
        "_mirnas.csv"
      )
    ),
    row.names = FALSE
  )

  invisible(
    heatmap_object
  )
}


save_method_candidate_pca <- function(
    tab,
    expression_matrix,
    expression_label,
    model_name,
    method_dir,
    contrast,
    candidate_type = c(
      "raw",
      "fdr"
    ),
    ellipse_level = 0.95
) {

  candidate_type <- match.arg(
    candidate_type
  )

  info <- candidate_output_info(
    candidate_type
  )

  selected_miRNAs <- intersect(
    candidate_mirnas_from_table(
      tab,
      candidate_type
    ),
    rownames(
      expression_matrix
    )
  )

  if (length(
    selected_miRNAs
  ) < 2L) {
    message(
      "Candidate PCA skipped: ",
      method_dir,
      " / ",
      model_name,
      " / ",
      contrast,
      " [",
      info$label,
      "] — fewer than 2 candidates."
    )

    return(
      invisible(
        NULL
      )
    )
  }

  spec <- candidate_group_spec(
    contrast
  )

  keep <- as.character(
    metadata[[
      spec$variable
    ]]
  ) %in%
    spec$levels

  plot_metadata <- metadata[
    keep,
    ,
    drop = FALSE
  ]

  plot_metadata$Plot_group <- factor(
    as.character(
      plot_metadata[[
        spec$variable
      ]]
    ),
    levels = spec$levels,
    labels = spec$labels
  )

  pca_matrix <- expression_matrix[
    selected_miRNAs,
    plot_metadata$Sample,
    drop = FALSE
  ]

  feature_sd <- apply(
    pca_matrix,
    1,
    stats::sd
  )

  pca_matrix <- pca_matrix[
    is.finite(
      feature_sd
    ) &
      feature_sd > 0,
    ,
    drop = FALSE
  ]

  if (nrow(
    pca_matrix
  ) < 2L) {
    message(
      "Candidate PCA skipped: ",
      method_dir,
      " / ",
      model_name,
      " / ",
      contrast,
      " [",
      info$label,
      "] — fewer than 2 non-zero-variance candidates."
    )

    return(
      invisible(
        NULL
      )
    )
  }

  pca <- stats::prcomp(
    t(
      pca_matrix
    ),
    center = TRUE,
    scale. = TRUE
  )

  if (ncol(
    pca$x
  ) < 2L) {
    return(
      invisible(
        NULL
      )
    )
  }

  variance <- (
    pca$sdev^2
  ) / sum(
    pca$sdev^2
  )

  pca_df <- data.frame(
    Sample = rownames(
      pca$x
    ),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  if (ncol(
    pca$x
  ) >= 3L) {
    pca_df$PC3 <- pca$x[, 3]
  }

  pca_df <- dplyr::left_join(
    pca_df,
    plot_metadata,
    by = "Sample"
  )

  group_counts <- table(
    pca_df$Plot_group
  )

  eligible_groups <- names(
    group_counts[
      group_counts >= 3
    ]
  )

  ellipse_data <- pca_df[
    as.character(
      pca_df$Plot_group
    ) %in%
      eligible_groups,
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = Plot_group
    )
  )

  if (nrow(
    ellipse_data
  ) > 0) {
    p <- p +
      ggplot2::stat_ellipse(
        data = ellipse_data,
        ggplot2::aes(
          fill = Plot_group
        ),
        geom = "polygon",
        type = "norm",
        level = ellipse_level,
        alpha = 0.12,
        color = NA,
        show.legend = FALSE
      ) +
      ggplot2::stat_ellipse(
        data = ellipse_data,
        type = "norm",
        level = ellipse_level,
        linewidth = 0.6,
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.95
    ) +
    ggplot2::labs(
      title = paste0(
        method_dir,
        " candidate PCA — ",
        contrast
      ),
      subtitle = paste0(
        nrow(
          pca_matrix
        ),
        " miRNAs selected at ",
        info$label,
        "; ",
        model_name,
        "; post-selection visualization"
      ),
      x = sprintf(
        "PC1 (%.1f%%)",
        100 * variance[1]
      ),
      y = sprintf(
        "PC2 (%.1f%%)",
        100 * variance[2]
      ),
      color = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )

  output_dir <- file.path(
    candidate_pca_dir,
    model_code(model_name),
    method_code(method_dir),
    info$tag
  )

  ensure_dir(output_dir)

  prefix <- contrast_code(
    contrast
  )

  save_ggplot_pdf(
    plot = p,
    filename = file.path(
      output_dir,
      paste0(
        prefix,
        "_pc12.pdf"
      )
    ),
    width = 7,
    height = 5
  )

  save_ggplot_tiff(
    plot = p,
    filename = file.path(
      output_dir,
      paste0(
        prefix,
        "_pc12.tiff"
      )
    ),
    width = 7,
    height = 5,
    dpi = 600
  )

  utils::write.csv(
    pca_df,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_coords.csv"
      )
    ),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(
      PC = paste0(
        "PC",
        seq_along(
          variance
        )
      ),
      variance_explained_pct = 100 * variance,
      cumulative_variance_explained_pct =
        100 * cumsum(
          variance
        ),
      stringsAsFactors = FALSE
    ),
    file.path(
      output_dir,
      paste0(
        prefix,
        "_var.csv"
      )
    ),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(
      model = model_name,
      method = method_dir,
      contrast = contrast,
      candidate_definition = info$label,
      expression_matrix = expression_label,
      miRNA = rownames(
        pca_matrix
      ),
      stringsAsFactors = FALSE
    ),
    file.path(
      output_dir,
      paste0(
        prefix,
        "_mirnas.csv"
      )
    ),
    row.names = FALSE
  )

  invisible(
    p
  )
}


save_method_candidate_visuals <- function(
    tab,
    expression_matrix,
    expression_label,
    model_name,
    method_dir,
    contrast
) {

  for (candidate_type in c(
    "raw",
    "fdr"
  )) {

    save_method_candidate_heatmap(
      tab = tab,
      expression_matrix = expression_matrix,
      expression_label = expression_label,
      model_name = model_name,
      method_dir = method_dir,
      contrast = contrast,
      candidate_type = candidate_type
    )

    save_method_candidate_pca(
      tab = tab,
      expression_matrix = expression_matrix,
      expression_label = expression_label,
      model_name = model_name,
      method_dir = method_dir,
      contrast = contrast,
      candidate_type = candidate_type
    )
  }

  invisible(
    NULL
  )
}


# ------------------------------------------------------------------------------
# 3. edgeR
# ------------------------------------------------------------------------------

run_edger_model <- function(model_name, spec) {
  message("edgeR — ", model_name)

  # Feature filtering was frozen upstream in script 01.
  # Do not call edgeR::filterByExpr() here.
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::normLibSizes(dge, method = "TMM")

  edgeR_expression <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 1
  )

  # Appendicitis status
  cond_formula <- formula_from_terms(spec$condition_terms)
  design_c <- stats::model.matrix(cond_formula, data = metadata)
  dge_c <- edgeR::estimateDisp(dge, design_c, robust = TRUE)
  fit_c <- edgeR::glmQLFit(dge_c, design_c, robust = TRUE)

  coef_c <- grep("^Appendicitis_status", colnames(design_c))
  if (length(coef_c) != 1) stop("Unexpected Appendicitis_status coefficients in edgeR design.")
  qlf_c <- edgeR::glmQLFTest(fit_c, coef = coef_c)
  tab_c <- edgeR::topTags(qlf_c, n = Inf, sort.by = "none")$table
  tab_c <- data.frame(
    miRNA = rownames(tab_c),
    log2FC = tab_c$logFC,
    average_logCPM = tab_c$logCPM,
    statistic = tab_c$F,
    p_value = tab_c$PValue,
    stringsAsFactors = FALSE
  )
  tab_c <- standardize_table(tab_c, "edgeR_TMM", model_name, "Appendicitis_vs_NIA")
  write_table(tab_c, model_name, "edgeR_TMM", "edgeR_TMM_Appendicitis_vs_NIA_all_results.csv")
  write_candidate_sets(tab_c, model_name, "edgeR_TMM", "Appendicitis_vs_NIA")
  plot_volcano(tab_c, model_name, "edgeR_TMM", "Appendicitis_vs_NIA")
  save_method_candidate_visuals(
    tab = tab_c,
    expression_matrix = edgeR_expression,
    expression_label = "TMM log2-CPM",
    model_name = model_name,
    method_dir = "edgeR_TMM",
    contrast = "Appendicitis_vs_NIA"
  )

  # Three-group phenotype
  sev_formula <- formula_from_terms(spec$severity_terms)
  design_s <- stats::model.matrix(sev_formula, data = metadata)
  dge_s <- edgeR::estimateDisp(dge, design_s, robust = TRUE)
  fit_s <- edgeR::glmQLFit(dge_s, design_s, robust = TRUE)
  sev_coef <- grep("^Phenotype", colnames(design_s))
  if (length(sev_coef) != 2) stop("Expected two Phenotype coefficients in edgeR three-group design.")

  # global
  qlf_global <- edgeR::glmQLFTest(fit_s, coef = sev_coef)
  tab_g <- edgeR::topTags(qlf_global, n = Inf, sort.by = "none")$table
  tab_g <- data.frame(
    miRNA = rownames(tab_g),
    statistic = tab_g$F,
    p_value = tab_g$PValue,
    stringsAsFactors = FALSE
  )
  tab_g <- standardize_table(tab_g, "edgeR_TMM", model_name, "Phenotype_global", global = TRUE)
  write_table(tab_g, model_name, "edgeR_TMM", "edgeR_TMM_Phenotype_global_all_results.csv")
  write_candidate_sets(tab_g, model_name, "edgeR_TMM", "Phenotype_global")
  save_method_candidate_visuals(
    tab = tab_g,
    expression_matrix = edgeR_expression,
    expression_label = "TMM log2-CPM",
    model_name = model_name,
    method_dir = "edgeR_TMM",
    contrast = "Phenotype_global"
  )

  # pairwise helper
  contrast_vector <- function(num_level, den_level) {
    v <- rep(0, ncol(design_s))
    names(v) <- colnames(design_s)

    coefficient_name <- function(level) {
      if (level == "NIA") return(NULL)
      paste0("Phenotype", level)
    }
    num <- coefficient_name(num_level)
    den <- coefficient_name(den_level)
    if (!is.null(num)) v[num] <- v[num] + 1
    if (!is.null(den)) v[den] <- v[den] - 1
    v
  }

  pairwise <- list(
    NIA_vs_Mild = c("Mild", "NIA"),
    NIA_vs_Severe = c("Severe", "NIA"),
    Mild_vs_Severe = c("Severe", "Mild")
  )

  for (nm in names(pairwise)) {
    lv <- pairwise[[nm]]
    qlf <- edgeR::glmQLFTest(fit_s, contrast = contrast_vector(lv[1], lv[2]))
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    tab <- data.frame(
      miRNA = rownames(tt),
      log2FC = tt$logFC,
      average_logCPM = tt$logCPM,
      statistic = tt$F,
      p_value = tt$PValue,
      stringsAsFactors = FALSE
    )
    tab <- standardize_table(tab, "edgeR_TMM", model_name, nm)
    write_table(tab, model_name, "edgeR_TMM", paste0("edgeR_TMM_", nm, "_all_results.csv"))
    write_candidate_sets(tab, model_name, "edgeR_TMM", nm)
    plot_volcano(tab, model_name, "edgeR_TMM", nm)
    save_method_candidate_visuals(
      tab = tab,
      expression_matrix = edgeR_expression,
      expression_label = "TMM log2-CPM",
      model_name = model_name,
      method_dir = "edgeR_TMM",
      contrast = nm
    )
  }

  # Save normalized factors for audit
  norm <- data.frame(
    Sample = colnames(counts),
    library_size = dge$samples$lib.size,
    TMM_norm_factor = dge$samples$norm.factors,
    effective_library_size = (
      dge$samples$lib.size *
        dge$samples$norm.factors
    ),
    stringsAsFactors = FALSE
  )
  d <- file.path(
    validation_results_dir,
    model_code(model_name),
    method_code("edgeR_TMM")
  )
  ensure_dir(d)
  utils::write.csv(
    norm,
    file.path(
      d,
      "norm.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. DESeq2
# ------------------------------------------------------------------------------

run_deseq2_model <- function(model_name, spec) {
  message("DESeq2 — ", model_name)

  coldata <- metadata
  rownames(coldata) <- coldata$Sample

  # Appendicitis status
  cond_formula <- formula_from_terms(spec$condition_terms)
  dds_c <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = cond_formula
  )
  dds_c <- DESeq2::DESeq(dds_c, quiet = TRUE)

  DESeq2_expression <- log2(
    DESeq2::counts(
      dds_c,
      normalized = TRUE
    ) + 1
  )

  res_c <- DESeq2::results(dds_c, contrast = c("Appendicitis_status", "Appendicitis", "NIA"), independentFiltering = FALSE)
  tab_c <- data.frame(
    miRNA = rownames(res_c),
    log2FC = res_c$log2FoldChange,
    standard_error = res_c$lfcSE,
    statistic = res_c$stat,
    p_value = res_c$pvalue,
    baseMean = res_c$baseMean,
    stringsAsFactors = FALSE
  )
  tab_c <- standardize_table(tab_c, "DESeq2", model_name, "Appendicitis_vs_NIA")
  write_table(tab_c, model_name, "DESeq2", "DESeq2_Appendicitis_vs_NIA_all_results.csv")
  write_candidate_sets(tab_c, model_name, "DESeq2", "Appendicitis_vs_NIA")
  plot_volcano(tab_c, model_name, "DESeq2", "Appendicitis_vs_NIA")
  save_method_candidate_visuals(
    tab = tab_c,
    expression_matrix = DESeq2_expression,
    expression_label = "log2(DESeq2 normalized counts + 1)",
    model_name = model_name,
    method_dir = "DESeq2",
    contrast = "Appendicitis_vs_NIA"
  )

  # Three-group phenotype pairwise
  sev_formula <- formula_from_terms(spec$severity_terms)
  dds_s <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = sev_formula
  )
  dds_s <- DESeq2::DESeq(dds_s, quiet = TRUE)

  pairwise <- list(
    NIA_vs_Mild = c("Mild", "NIA"),
    NIA_vs_Severe = c("Severe", "NIA"),
    Mild_vs_Severe = c("Severe", "Mild")
  )
  for (nm in names(pairwise)) {
    lv <- pairwise[[nm]]
    rr <- DESeq2::results(dds_s, contrast = c("Phenotype", lv[1], lv[2]), independentFiltering = FALSE)
    tab <- data.frame(
      miRNA = rownames(rr),
      log2FC = rr$log2FoldChange,
      standard_error = rr$lfcSE,
      statistic = rr$stat,
      p_value = rr$pvalue,
      baseMean = rr$baseMean,
      stringsAsFactors = FALSE
    )
    tab <- standardize_table(tab, "DESeq2", model_name, nm)
    write_table(tab, model_name, "DESeq2", paste0("DESeq2_", nm, "_all_results.csv"))
    write_candidate_sets(tab, model_name, "DESeq2", nm)
    plot_volcano(tab, model_name, "DESeq2", nm)
    save_method_candidate_visuals(
      tab = tab,
      expression_matrix = DESeq2_expression,
      expression_label = "log2(DESeq2 normalized counts + 1)",
      model_name = model_name,
      method_dir = "DESeq2",
      contrast = nm
    )
  }

  # Global phenotype LRT. Keep covariates in the reduced model.
  reduced_terms <- setdiff(spec$severity_terms, "Phenotype")
  reduced_formula <- if (length(reduced_terms) == 0) {
    stats::as.formula("~ 1")
  } else {
    formula_from_terms(reduced_terms)
  }

  dds_lrt <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = sev_formula
  )
  dds_lrt <- DESeq2::DESeq(dds_lrt, test = "LRT", reduced = reduced_formula, quiet = TRUE)
  rr <- DESeq2::results(dds_lrt, independentFiltering = FALSE)
  tab_g <- data.frame(
    miRNA = rownames(rr),
    statistic = rr$stat,
    p_value = rr$pvalue,
    baseMean = rr$baseMean,
    stringsAsFactors = FALSE
  )
  tab_g <- standardize_table(tab_g, "DESeq2", model_name, "Phenotype_global", global = TRUE)
  write_table(tab_g, model_name, "DESeq2", "DESeq2_Phenotype_global_all_results.csv")
  write_candidate_sets(tab_g, model_name, "DESeq2", "Phenotype_global")
  save_method_candidate_visuals(
    tab = tab_g,
    expression_matrix = DESeq2_expression,
    expression_label = "log2(DESeq2 normalized counts + 1)",
    model_name = model_name,
    method_dir = "DESeq2",
    contrast = "Phenotype_global"
  )

  # Size factors
  sf <- data.frame(
    Sample = colnames(dds_c),
    DESeq2_size_factor = DESeq2::sizeFactors(dds_c),
    stringsAsFactors = FALSE
  )
  d <- file.path(
    validation_results_dir,
    model_code(model_name),
    method_code("DESeq2")
  )
  ensure_dir(d)
  utils::write.csv(
    sf,
    file.path(
      d,
      "norm.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Execute all estimable models
# ------------------------------------------------------------------------------

for (model_name in names(valid_models)) {
  run_edger_model(model_name, valid_models[[model_name]])
  run_deseq2_model(model_name, valid_models[[model_name]])
}

# Summary
summary_rows <- list()
for (model_name in names(valid_models)) {
  for (method in c("edgeR_TMM", "DESeq2")) {
    table_dir <- file.path(
      validation_results_dir,
      model_code(model_name),
      method_code(method),
      "tab"
    )
    files <- list.files(
      table_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    )
    for (f in files) {
      x <- utils::read.csv(f, stringsAsFactors = FALSE)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        model = model_name,
        method = method,
        contrast = if ("contrast" %in% colnames(x)) x$contrast[1] else basename(f),
        file = basename(f),
        n_miRNAs_expected = nrow(counts),
        n_rows_returned = nrow(x),
        n_p_values_available = sum(!is.na(x$p_value)),
        n_p_values_NA = sum(is.na(x$p_value)),
        n_raw_p_lt_0.05 = sum(x$p_value < ALPHA_RAW, na.rm = TRUE),
        n_FDR_lt_0.05 = sum(x$FDR < ALPHA_FDR, na.rm = TRUE),
        n_Holm_lt_0.05 = sum(x$Holm < ALPHA_FDR, na.rm = TRUE),
        minimum_p_value = if (all(is.na(x$p_value))) NA_real_ else min(x$p_value, na.rm = TRUE),
        minimum_FDR = if (all(is.na(x$FDR))) NA_real_ else min(x$FDR, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}
summary_table <- dplyr::bind_rows(summary_rows)
utils::write.csv(
  summary_table,
  file.path(validation_results_dir, "summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 6. Data-driven candidate evidence summaries
# ------------------------------------------------------------------------------

# Re-read the complete result tables into one standardized long table.
#
# This avoids defining candidates from file names or external lists and creates
# a single audit object spanning method x model x contrast.

result_tables <- list()

result_index <- 0L

for (model_name in names(valid_models)) {

  for (method in c(
    "edgeR_TMM",
    "DESeq2"
  )) {

    table_dir <- file.path(
      validation_results_dir,
      model_code(model_name),
      method_code(method),
      "tab"
    )

    files <- list.files(
      table_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    )

    for (f in files) {

      result_index <- result_index + 1L

      x <- utils::read.csv(
        f,
        stringsAsFactors = FALSE
      )

      required_result_columns <- c(
        "miRNA",
        "p_value",
        "FDR",
        "Holm",
        "raw_significant",
        "fdr_significant",
        "holm_significant",
        "method",
        "model",
        "contrast",
        "log2FC",
        "direction"
      )

      missing_result_columns <- setdiff(
        required_result_columns,
        colnames(x)
      )

      if (
        length(
          missing_result_columns
        ) > 0L
      ) {

        stop(
          "Result table ",
          basename(f),
          " lacks required column(s): ",
          paste(
            missing_result_columns,
            collapse = ", "
          ),
          call. = FALSE
        )
      }

      x$source_file <- basename(f)

      result_tables[[result_index]] <- x
    }
  }
}


all_count_results <- dplyr::bind_rows(
  result_tables
)


expected_contrasts <- c(
  "Appendicitis_vs_NIA",
  "Phenotype_global",
  "NIA_vs_Mild",
  "NIA_vs_Severe",
  "Mild_vs_Severe"
)


expected_methods <- c(
  "edgeR_TMM",
  "DESeq2"
)


if (!setequal(
  unique(
    all_count_results$contrast
  ),
  expected_contrasts
)) {

  stop(
    "Unexpected contrast set in combined count-based DE results.",
    call. = FALSE
  )
}


if (!setequal(
  unique(
    all_count_results$method
  ),
  expected_methods
)) {

  stop(
    "Unexpected method set in combined count-based DE results.",
    call. = FALSE
  )
}


expected_total_rows <-
  nrow(counts) *
  length(
    expected_contrasts
  ) *
  length(
    expected_methods
  ) *
  length(
    valid_models
  )


if (
  nrow(
    all_count_results
  ) != expected_total_rows
) {

  stop(
    "Combined count-based DE result table has ",
    nrow(
      all_count_results
    ),
    " rows; expected ",
    expected_total_rows,
    ".",
    call. = FALSE
  )
}


utils::write.csv(
  all_count_results,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_all_results_long.csv"
    )
  ),
  row.names = FALSE
)


safe_min <- function(
    x
) {

  x <- x[
    !is.na(
      x
    ) &
      is.finite(
        x
      )
  ]

  if (
    length(
      x
    ) == 0L
  ) {

    return(
      NA_real_
    )
  }

  min(
    x
  )
}


collapse_unique <- function(
    x
) {

  x <- unique(
    as.character(
      x[
        !is.na(
          x
        ) &
          nzchar(
            as.character(
              x
            )
          )
      ]
    )
  )

  if (
    length(
      x
    ) == 0L
  ) {

    return(
      ""
    )
  }

  paste(
    sort(
      x
    ),
    collapse = "; "
  )
}


summarize_candidate_evidence <- function(
    data,
    summary_label
) {

  out <- data |>
    dplyr::group_by(
      miRNA
    ) |>
    dplyr::summarise(

      evidence_scope = summary_label,

      n_results_tested = dplyr::n(),

      n_p_values_available = sum(
        !is.na(
          p_value
        )
      ),

      minimum_raw_p = safe_min(
        p_value
      ),

      minimum_FDR = safe_min(
        FDR
      ),

      minimum_Holm = safe_min(
        Holm
      ),

      n_nominal_results = sum(
        raw_significant,
        na.rm = TRUE
      ),

      n_FDR_results = sum(
        fdr_significant,
        na.rm = TRUE
      ),

      n_Holm_results = sum(
        holm_significant,
        na.rm = TRUE
      ),

      n_nominal_methods = dplyr::n_distinct(
        method[
          raw_significant %in%
            TRUE
        ]
      ),

      n_FDR_methods = dplyr::n_distinct(
        method[
          fdr_significant %in%
            TRUE
        ]
      ),

      n_nominal_contrasts = dplyr::n_distinct(
        contrast[
          raw_significant %in%
            TRUE
        ]
      ),

      n_FDR_contrasts = dplyr::n_distinct(
        contrast[
          fdr_significant %in%
            TRUE
        ]
      ),

      nominal_methods = collapse_unique(
        method[
          raw_significant %in%
            TRUE
        ]
      ),

      FDR_methods = collapse_unique(
        method[
          fdr_significant %in%
            TRUE
        ]
      ),

      nominal_contrasts = collapse_unique(
        contrast[
          raw_significant %in%
            TRUE
        ]
      ),

      FDR_contrasts = collapse_unique(
        contrast[
          fdr_significant %in%
            TRUE
        ]
      ),

      .groups = "drop"
    )

  out
}


primary_results <- all_count_results |>
  dplyr::filter(
    model == "unadjusted"
  )


if (
  nrow(
    primary_results
  ) !=
    nrow(
      counts
    ) *
      length(
        expected_contrasts
      ) *
      length(
        expected_methods
      )
) {

  stop(
    "Primary unadjusted count-based result set has an unexpected size.",
    call. = FALSE
  )
}


primary_candidate_summary <- summarize_candidate_evidence(
  data = primary_results,
  summary_label = "primary_unadjusted"
)


# Cross-method agreement is evaluated contrast-by-contrast in the primary
# unadjusted analyses.
#
# For pairwise contrasts, direction concordance is also recorded.

primary_cross_method_by_contrast <- primary_results |>
  dplyr::group_by(
    miRNA,
    contrast
  ) |>
  dplyr::summarise(

    n_methods_tested = dplyr::n_distinct(
      method
    ),

    n_methods_nominal = dplyr::n_distinct(
      method[
        raw_significant %in%
          TRUE
      ]
    ),

    n_methods_FDR = dplyr::n_distinct(
      method[
        fdr_significant %in%
          TRUE
      ]
    ),

    edgeR_p = p_value[
      method == "edgeR_TMM"
    ][1L],

    DESeq2_p = p_value[
      method == "DESeq2"
    ][1L],

    edgeR_FDR = FDR[
      method == "edgeR_TMM"
    ][1L],

    DESeq2_FDR = FDR[
      method == "DESeq2"
    ][1L],

    edgeR_log2FC = log2FC[
      method == "edgeR_TMM"
    ][1L],

    DESeq2_log2FC = log2FC[
      method == "DESeq2"
    ][1L],

    .groups = "drop"
  )


primary_cross_method_by_contrast$nominal_both_methods <-
  primary_cross_method_by_contrast$n_methods_nominal == 2L


primary_cross_method_by_contrast$FDR_both_methods <-
  primary_cross_method_by_contrast$n_methods_FDR == 2L


primary_cross_method_by_contrast$direction_concordant <- with(
  primary_cross_method_by_contrast,
  ifelse(
    contrast == "Phenotype_global",
    NA,
    ifelse(
      is.na(
        edgeR_log2FC
      ) |
        is.na(
          DESeq2_log2FC
        ),
      NA,
      sign(
        edgeR_log2FC
      ) ==
        sign(
          DESeq2_log2FC
        )
    )
  )
)


utils::write.csv(
  primary_cross_method_by_contrast,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_primary_cross_method_by_contrast.csv"
    )
  ),
  row.names = FALSE
)


cross_method_mirna_summary <- primary_cross_method_by_contrast |>
  dplyr::group_by(
    miRNA
  ) |>
  dplyr::summarise(

    n_contrasts_nominal_both_methods = sum(
      nominal_both_methods,
      na.rm = TRUE
    ),

    n_contrasts_FDR_both_methods = sum(
      FDR_both_methods,
      na.rm = TRUE
    ),

    n_pairwise_nominal_both_methods_direction_concordant = sum(
      contrast !=
        "Phenotype_global" &
        nominal_both_methods &
        direction_concordant %in%
          TRUE,
      na.rm = TRUE
    ),

    cross_method_nominal_contrasts = collapse_unique(
      contrast[
        nominal_both_methods %in%
          TRUE
      ]
    ),

    cross_method_FDR_contrasts = collapse_unique(
      contrast[
        FDR_both_methods %in%
          TRUE
      ]
    ),

    .groups = "drop"
  )


primary_candidate_summary <- dplyr::left_join(
  primary_candidate_summary,
  cross_method_mirna_summary,
  by = "miRNA"
)


primary_candidate_summary$signal_tier <- dplyr::case_when(

  primary_candidate_summary$n_FDR_results > 0L ~
    "FDR_supported_primary",

  primary_candidate_summary$n_contrasts_nominal_both_methods > 0L ~
    "cross_method_nominal",

  primary_candidate_summary$n_nominal_results >= 2L ~
    "recurrent_primary_nominal",

  primary_candidate_summary$n_nominal_results == 1L ~
    "single_primary_nominal",

  TRUE ~
    "no_primary_nominal_signal"
)


primary_tier_order <- c(
  "FDR_supported_primary",
  "cross_method_nominal",
  "recurrent_primary_nominal",
  "single_primary_nominal",
  "no_primary_nominal_signal"
)


primary_candidate_summary$signal_tier <- factor(
  primary_candidate_summary$signal_tier,
  levels = primary_tier_order
)


primary_candidate_summary <- primary_candidate_summary |>
  dplyr::arrange(
    signal_tier,
    dplyr::desc(
      n_contrasts_FDR_both_methods
    ),
    dplyr::desc(
      n_FDR_results
    ),
    dplyr::desc(
      n_contrasts_nominal_both_methods
    ),
    dplyr::desc(
      n_nominal_results
    ),
    minimum_FDR,
    minimum_raw_p
  )


primary_candidate_summary$signal_tier <- as.character(
  primary_candidate_summary$signal_tier
)


if (
  nrow(
    primary_candidate_summary
  ) !=
    nrow(
      counts
    ) ||
  !setequal(
    primary_candidate_summary$miRNA,
    rownames(
      counts
    )
  )
) {

  stop(
    "Primary candidate evidence summary does not match the complete ",
    "upstream-frozen feature set.",
    call. = FALSE
  )
}


primary_any_nominal <- primary_candidate_summary |>
  dplyr::filter(
    n_nominal_results >= 1L
  )


primary_recurrent_nominal <- primary_candidate_summary |>
  dplyr::filter(
    n_nominal_results >= 2L
  )


primary_cross_method_nominal <- primary_candidate_summary |>
  dplyr::filter(
    n_contrasts_nominal_both_methods >= 1L
  )


primary_FDR_supported <- primary_candidate_summary |>
  dplyr::filter(
    n_FDR_results >= 1L
  )


utils::write.csv(
  primary_candidate_summary,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_candidate_primary_unadjusted_summary.csv"
    )
  ),
  row.names = FALSE
)


utils::write.csv(
  primary_any_nominal,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_candidate_primary_any_nominal.csv"
    )
  ),
  row.names = FALSE
)


utils::write.csv(
  primary_recurrent_nominal,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_candidate_primary_recurrent_nominal.csv"
    )
  ),
  row.names = FALSE
)


utils::write.csv(
  primary_cross_method_nominal,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_candidate_primary_cross_method_nominal.csv"
    )
  ),
  row.names = FALSE
)


utils::write.csv(
  primary_FDR_supported,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "DE_count_candidate_primary_FDR_supported.csv"
    )
  ),
  row.names = FALSE
)


# Age + sex sensitivity evidence is intentionally summarized separately.

if (
  "adjusted_age_sex" %in%
    unique(
      all_count_results$model
    )
) {

  adjusted_results <- all_count_results |>
    dplyr::filter(
      model == "adjusted_age_sex"
    )

  adjusted_candidate_summary <- summarize_candidate_evidence(
    data = adjusted_results,
    summary_label = "adjusted_age_sex"
  )

  adjusted_candidate_summary <- adjusted_candidate_summary |>
    dplyr::arrange(
      dplyr::desc(
        n_FDR_results
      ),
      dplyr::desc(
        n_nominal_results
      ),
      minimum_FDR,
      minimum_raw_p
    )

  utils::write.csv(
    adjusted_candidate_summary,
    prepare_output_file(
      file.path(
        validation_results_dir,
        "DE_count_candidate_adjusted_age_sex_summary.csv"
      )
    ),
    row.names = FALSE
  )

} else {

  adjusted_candidate_summary <- data.frame()
}


# ------------------------------------------------------------------------------
# 7. Final QC / audit
# ------------------------------------------------------------------------------

final_qc <- data.frame(

  analysis_id = context$analysis_id,

  n_samples = ncol(
    counts
  ),

  n_NIA = sum(
    metadata$Phenotype ==
      "NIA"
  ),

  n_Mild = sum(
    metadata$Phenotype ==
      "Mild"
  ),

  n_Severe = sum(
    metadata$Phenotype ==
      "Severe"
  ),

  n_miRNAs = nrow(
    counts
  ),

  n_valid_models = length(
    valid_models
  ),

  valid_models = paste(
    names(
      valid_models
    ),
    collapse = "; "
  ),

  n_methods = length(
    expected_methods
  ),

  n_contrasts = length(
    expected_contrasts
  ),

  combined_result_rows = nrow(
    all_count_results
  ),

  expected_combined_result_rows =
    expected_total_rows,

  complete_result_grid_valid =
    nrow(
      all_count_results
    ) ==
      expected_total_rows,

  primary_unadjusted_rows = nrow(
    primary_results
  ),

  n_primary_any_nominal = nrow(
    primary_any_nominal
  ),

  n_primary_recurrent_nominal = nrow(
    primary_recurrent_nominal
  ),

  n_primary_cross_method_nominal = nrow(
    primary_cross_method_nominal
  ),

  n_primary_FDR_supported = nrow(
    primary_FDR_supported
  ),

  edgeR_additional_filtering = FALSE,

  DESeq2_independent_filtering = FALSE,

  historical_candidate_list_used = FALSE,

  samples_removed_by_script = 0L,

  stringsAsFactors = FALSE
)


utils::write.csv(
  final_qc,
  prepare_output_file(
    file.path(
      validation_results_dir,
      "05_final_QC_summary.csv"
    )
  ),
  row.names = FALSE
)


message(
  ""
)

message(
  "Primary unadjusted count-based candidate summary:"
)

message(
  "- nominal p < 0.05 in >=1 method x contrast result: ",
  nrow(
    primary_any_nominal
  )
)

message(
  "- recurrent nominal in >=2 primary results: ",
  nrow(
    primary_recurrent_nominal
  )
)

message(
  "- nominal in both edgeR and DESeq2 for >=1 same contrast: ",
  nrow(
    primary_cross_method_nominal
  )
)

message(
  "- FDR-supported in >=1 primary result: ",
  nrow(
    primary_FDR_supported
  )
)


if (
  nrow(
    primary_any_nominal
  ) > 0L
) {

  print(
    utils::head(
      primary_any_nominal[
        ,
        c(
          "miRNA",
          "signal_tier",
          "n_nominal_results",
          "n_FDR_results",
          "n_contrasts_nominal_both_methods",
          "minimum_raw_p",
          "minimum_FDR",
          "nominal_methods",
          "nominal_contrasts",
          "cross_method_nominal_contrasts"
        ),
        drop = FALSE
      ],
      30L
    )
  )
}


write_session_info(results_dir)

message("")
message("05 complete — count-based models generated for ", context$analysis_id, ".")
message(
  "Tested feature set: ",
  nrow(counts),
  " miRNAs; ",
  ncol(counts),
  " samples."
)
message("No additional edgeR feature filtering was performed.")
message("DESeq2 independent filtering was disabled.")
message("Primary model: unadjusted; age + sex is a separate sensitivity model when enabled.")
message("BH-FDR is calculated separately for each method x model x contrast.")
message("Primary recurrence is based only on unadjusted edgeR + DESeq2 results.")
message("No historical candidate list was used for testing or candidate definition.")
message("Method-specific candidate sets, volcano plots, heatmaps and candidate PCA were exported.")
message("")
message("Main outputs:")
message("- ", file.path(validation_results_dir, "manifest.csv"))
message("- ", file.path(validation_results_dir, "model_audit.csv"))
message("- ", file.path(validation_results_dir, "summary.csv"))
message("- ", file.path(validation_results_dir, "DE_count_all_results_long.csv"))
message("- ", file.path(validation_results_dir, "DE_count_candidate_primary_unadjusted_summary.csv"))
message("- ", file.path(validation_results_dir, "DE_count_primary_cross_method_by_contrast.csv"))
message("- ", file.path(validation_results_dir, "DE_count_candidate_primary_cross_method_nominal.csv"))
message("- ", file.path(validation_results_dir, "DE_count_candidate_primary_FDR_supported.csv"))
if (
  "adjusted_age_sex" %in%
    names(
      valid_models
    )
) {
  message("- ", file.path(validation_results_dir, "DE_count_candidate_adjusted_age_sex_summary.csv"))
}
message("- ", file.path(validation_results_dir, "05_final_QC_summary.csv"))
message("- ", validation_results_dir)
message("- ", volcano_dir)
message("- ", heatmap_dir)
message("- ", candidate_pca_dir)
