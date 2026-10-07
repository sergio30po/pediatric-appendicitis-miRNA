# ==============================================================================
# 00_utils.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Script:
#   00_utils.R
#
# Purpose:
#   Provide shared utility functions used throughout the complete miRNA-seq
#   analysis pipeline.
#
# Role in pipeline:
#   Shared infrastructure only.
#
#   This script centralizes:
#     - output-directory and file-path handling;
#     - graphics-device handling;
#     - standardized figure export;
#     - project input loading;
#     - base input validation;
#     - dynamic reconstruction of analysis-cohort membership;
#     - cohort selection and alignment;
#     - cohort traceability;
#     - session-information export;
#     - common statistical helpers;
#     - legacy integrated-candidate audit access.
#
#   No differential-expression analysis, candidate selection or biological
#   interpretation is performed here.
#
# Inputs:
#   Expected project files:
#
#     data/counts.tsv
#     data/metadata_master.xlsx
#     data/Metrics.xlsx
#
#   Required objects from R/00_analysis_config.R:
#
#     project_root
#     ANALYSIS_CATALOG
#     PRIMARY_ANALYSIS_ID
#
# Outputs:
#   This script does not generate files when sourced.
#
#   Utility functions defined here can generate:
#
#     analysis_manifest.csv
#     cohort_samples.csv
#     sessionInfo.txt
#     PDF/TIFF/PNG figures
#
# Methodological decisions:
#   - Raw count matrices must contain non-negative integer-like values.
#   - Sample IDs must be unique and aligned across counts, metadata and metrics.
#   - Condition and Apend_type must be internally consistent:
#
#       Condition C -> Apend_type 0
#       Condition P -> Apend_type 1 or 2
#
#   - Cohort definitions and exclusions are obtained exclusively from
#     R/00_analysis_config.R.
#   - Cohort-membership variables stored historically in metadata_master.xlsx
#     are never trusted as the source of truth.
#   - Cohort-membership columns are reconstructed dynamically from
#     ANALYSIS_CATALOG every time the pipeline is run.
#
# Required packages:
#   readxl
#   ggplot2
#
# Reproducibility:
#   Downstream scripts should reuse these functions rather than implement
#   separate versions of input loading, cohort selection, output handling or
#   p-value adjustment.
#
# Usage:
#
#   source("R/00_analysis_config.R")
#   source("R/00_utils.R")
#
# ==============================================================================


# ==============================================================================
# 1. DEPENDENCIES
# ==============================================================================

required_utility_packages <- c(
  "readxl",
  "ggplot2"
)


missing_utility_packages <- required_utility_packages[
  !vapply(
    required_utility_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(missing_utility_packages) > 0L) {
  
  stop(
    "Missing package(s) required by R/00_utils.R:\n",
    paste(
      "-",
      missing_utility_packages,
      collapse = "\n"
    ),
    "\nRun R/00_packages_setup.R first.",
    call. = FALSE
  )
}


# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# This utility file does not define scientific configuration.
#
# project_root and all cohort/analysis parameters must come from:
#
#   R/00_analysis_config.R


if (!exists(
  "project_root",
  inherits = TRUE
)) {
  
  stop(
    "project_root is not defined. ",
    "Source R/00_analysis_config.R before R/00_utils.R.",
    call. = FALSE
  )
}


# ==============================================================================
# 3. INPUTS
# ==============================================================================

# ==============================================================================
# 3.1 Output-path utilities
# ==============================================================================

ensure_dir <- function(path) {
  
  if (
    length(path) != 1L ||
    is.na(path) ||
    !nzchar(path)
  ) {
    
    stop(
      "'path' must be one non-empty path.",
      call. = FALSE
    )
  }
  
  
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  if (!dir.exists(path)) {
    
    stop(
      "Could not create directory: ",
      path,
      call. = FALSE
    )
  }
  
  
  invisible(path)
}


ensure_parent_dir <- function(path) {
  
  ensure_dir(
    dirname(path)
  )
  
  invisible(path)
}


warn_if_long_path <- function(
    path,
    threshold = 240L
) {
  
  if (
    .Platform$OS.type == "windows" &&
    nchar(
      path,
      type = "chars"
    ) >= threshold
  ) {
    
    warning(
      "Long Windows output path (",
      nchar(
        path,
        type = "chars"
      ),
      " characters): ",
      path,
      "\nShorten downstream directory/file names. ",
      "Do not use machine-specific drive mappings as a pipeline requirement.",
      call. = FALSE
    )
  }
  
  
  invisible(path)
}


prepare_output_file <- function(path) {
  
  ensure_parent_dir(path)
  
  warn_if_long_path(path)
  
  path
}


# ==============================================================================
# 3.2 Project input reader
# ==============================================================================

read_project_inputs <- function() {
  
  data_dir <- file.path(
    project_root,
    "data"
  )
  
  
  files <- c(
    
    counts = file.path(
      data_dir,
      "counts.tsv"
    ),
    
    metadata_master = file.path(
      data_dir,
      "metadata_master.xlsx"
    ),
    
    metrics = file.path(
      data_dir,
      "Metrics.xlsx"
    )
  )
  
  
  missing_files <- files[
    !file.exists(files)
  ]
  
  
  if (length(missing_files) > 0L) {
    
    stop(
      "Missing required project input file(s):\n",
      paste(
        "-",
        missing_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Counts
  # ---------------------------------------------------------------------------
  
  counts <- utils::read.delim(
    files[["counts"]],
    row.names = 1,
    check.names = FALSE
  )
  
  
  counts <- as.matrix(
    counts
  )
  
  
  storage.mode(
    counts
  ) <- "numeric"
  
  
  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------
  
  metadata_master <- readxl::read_excel(
    files[["metadata_master"]],
    sheet = "metadata_master"
  )
  
  
  metadata_master <- as.data.frame(
    metadata_master,
    stringsAsFactors = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Sequencing metrics
  # ---------------------------------------------------------------------------
  
  metrics <- readxl::read_excel(
    files[["metrics"]]
  )
  
  
  metrics <- as.data.frame(
    metrics,
    stringsAsFactors = FALSE
  )
  
  
  list(
    counts = counts,
    metadata_master = metadata_master,
    metrics = metrics,
    input_files = files
  )
}


# ==============================================================================
# 4. VALIDATION
# ==============================================================================

# ==============================================================================
# 4.1 Base-input validation
# ==============================================================================

validate_base_inputs <- function(
    counts,
    metadata_master,
    metrics
) {
  
  # ---------------------------------------------------------------------------
  # Counts structure
  # ---------------------------------------------------------------------------
  
  if (!is.matrix(
    counts
  )) {
    
    stop(
      "'counts' must be a matrix.",
      call. = FALSE
    )
  }
  
  
  if (
    is.null(
      rownames(counts)
    ) ||
    is.null(
      colnames(counts)
    )
  ) {
    
    stop(
      "counts.tsv must contain miRNA row names and sample column names.",
      call. = FALSE
    )
  }
  
  
  if (
    anyNA(
      rownames(counts)
    ) ||
    anyNA(
      colnames(counts)
    ) ||
    any(
      !nzchar(
        rownames(counts)
      )
    ) ||
    any(
      !nzchar(
        colnames(counts)
      )
    )
  ) {
    
    stop(
      "counts.tsv contains missing or empty miRNA/sample identifiers.",
      call. = FALSE
    )
  }
  
  
  if (
    anyDuplicated(
      rownames(counts)
    ) > 0L
  ) {
    
    stop(
      "Duplicated miRNA IDs in counts.tsv.",
      call. = FALSE
    )
  }
  
  
  if (
    anyDuplicated(
      colnames(counts)
    ) > 0L
  ) {
    
    stop(
      "Duplicated sample IDs in counts.tsv.",
      call. = FALSE
    )
  }
  
  
  if (
    anyNA(counts) ||
    any(!is.finite(counts)) ||
    any(counts < 0)
  ) {
    
    stop(
      "counts.tsv contains missing, non-finite or negative values.",
      call. = FALSE
    )
  }
  
  
  non_integer_like <- abs(
    counts -
      round(counts)
  ) > sqrt(
    .Machine$double.eps
  )
  
  
  if (any(
    non_integer_like
  )) {
    
    stop(
      "counts.tsv contains non-integer-like values. ",
      "Raw counts are required for the count-based DE methods.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Metadata structure
  # ---------------------------------------------------------------------------
  
  required_metadata <- c(
    "Sample",
    "Condition",
    "Sex",
    "Age",
    "Apend_type"
  )
  
  
  missing_metadata_columns <- setdiff(
    required_metadata,
    colnames(
      metadata_master
    )
  )
  
  
  if (
    length(
      missing_metadata_columns
    ) > 0L
  ) {
    
    stop(
      "metadata_master.xlsx lacks required column(s): ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  
  metadata_samples <- as.character(
    metadata_master$Sample
  )
  
  
  if (
    anyNA(
      metadata_samples
    ) ||
    any(
      !nzchar(
        metadata_samples
      )
    )
  ) {
    
    stop(
      "metadata_master.xlsx contains missing/empty Sample IDs.",
      call. = FALSE
    )
  }
  
  
  if (
    anyDuplicated(
      metadata_samples
    ) > 0L
  ) {
    
    stop(
      "Duplicated Sample IDs in metadata_master.xlsx.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Metrics structure
  # ---------------------------------------------------------------------------
  
  required_metrics <- c(
    "Sample",
    "Reads"
  )
  
  
  missing_metric_columns <- setdiff(
    required_metrics,
    colnames(
      metrics
    )
  )
  
  
  if (
    length(
      missing_metric_columns
    ) > 0L
  ) {
    
    stop(
      "Metrics.xlsx lacks required column(s): ",
      paste(
        missing_metric_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  
  metric_samples <- as.character(
    metrics$Sample
  )
  
  
  if (
    anyNA(
      metric_samples
    ) ||
    any(
      !nzchar(
        metric_samples
      )
    )
  ) {
    
    stop(
      "Metrics.xlsx contains missing/empty Sample IDs.",
      call. = FALSE
    )
  }
  
  
  if (
    anyDuplicated(
      metric_samples
    ) > 0L
  ) {
    
    stop(
      "Duplicated Sample IDs in Metrics.xlsx.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Cross-file sample coverage
  # ---------------------------------------------------------------------------
  
  missing_metadata_samples <- setdiff(
    colnames(counts),
    metadata_samples
  )
  
  
  missing_metric_samples <- setdiff(
    colnames(counts),
    metric_samples
  )
  
  
  if (
    length(
      missing_metadata_samples
    ) > 0L
  ) {
    
    stop(
      "Samples present in counts.tsv but absent from metadata_master.xlsx: ",
      paste(
        missing_metadata_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  if (
    length(
      missing_metric_samples
    ) > 0L
  ) {
    
    stop(
      "Samples present in counts.tsv but absent from Metrics.xlsx: ",
      paste(
        missing_metric_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  invisible(TRUE)
}


# ==============================================================================
# 4.2 Selected-cohort phenotype validation
# ==============================================================================

validate_selected_metadata <- function(
    metadata
) {
  
  required_columns <- c(
    "Sample",
    "Condition",
    "Sex",
    "Age",
    "Apend_type"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(
      metadata
    )
  )
  
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "Selected metadata lacks required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Sample
  # ---------------------------------------------------------------------------
  
  if (
    anyNA(
      metadata$Sample
    ) ||
    any(
      !nzchar(
        as.character(
          metadata$Sample
        )
      )
    ) ||
    anyDuplicated(
      metadata$Sample
    ) > 0L
  ) {
    
    stop(
      "Invalid Sample identifiers in selected cohort.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Condition
  # ---------------------------------------------------------------------------
  
  if (
    anyNA(
      metadata$Condition
    )
  ) {
    
    stop(
      "Unexpected or missing Condition in selected cohort.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Apend_type
  # ---------------------------------------------------------------------------
  
  if (
    anyNA(
      metadata$Apend_type
    )
  ) {
    
    stop(
      "Unexpected or missing Apend_type in selected cohort.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Condition / Apend_type consistency
  # ---------------------------------------------------------------------------
  
  condition_chr <- as.character(
    metadata$Condition
  )
  
  
  apend_chr <- as.character(
    metadata$Apend_type
  )
  
  
  invalid_control <- (
    condition_chr == "C" &
      apend_chr != "0"
  )
  
  
  invalid_patient <- (
    condition_chr == "P" &
      !apend_chr %in%
      c(
        "1",
        "2"
      )
  )
  
  
  if (
    any(
      invalid_control
    ) ||
    any(
      invalid_patient
    )
  ) {
    
    bad_samples <- metadata$Sample[
      invalid_control |
        invalid_patient
    ]
    
    
    stop(
      "Condition/Apend_type inconsistency detected for sample(s): ",
      paste(
        bad_samples,
        collapse = ", "
      ),
      ". Expected C -> 0 and P -> 1/2.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Sex
  # ---------------------------------------------------------------------------
  
  if (
    anyNA(
      metadata$Sex
    )
  ) {
    
    stop(
      "Sex contains missing values in selected cohort.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Age
  # ---------------------------------------------------------------------------
  
  if (
    anyNA(
      metadata$Age
    ) ||
    any(
      !is.finite(
        metadata$Age
      )
    )
  ) {
    
    stop(
      "Age contains missing or non-finite values in selected cohort.",
      call. = FALSE
    )
  }
  
  
  invisible(TRUE)
}


# ==============================================================================
# 4.3 Rebuild cohort-membership metadata
# ==============================================================================

rebuild_cohort_membership_metadata <- function(
    metadata_master,
    counts
) {
  
  if (!exists(
    "ANALYSIS_CATALOG",
    inherits = TRUE
  )) {
    
    stop(
      "ANALYSIS_CATALOG is not defined. ",
      "Source R/00_analysis_config.R first.",
      call. = FALSE
    )
  }
  
  
  if (!exists(
    "PRIMARY_ANALYSIS_ID",
    inherits = TRUE
  )) {
    
    stop(
      "PRIMARY_ANALYSIS_ID is not defined. ",
      "Source R/00_analysis_config.R first.",
      call. = FALSE
    )
  }
  
  
  if (
    !"Sample" %in%
    colnames(
      metadata_master
    )
  ) {
    
    stop(
      "metadata_master must contain a Sample column.",
      call. = FALSE
    )
  }
  
  
  metadata_master$Sample <- as.character(
    metadata_master$Sample
  )
  
  
  final_qc_samples <- colnames(
    counts
  )
  
  
  if (
    length(
      final_qc_samples
    ) == 0L
  ) {
    
    stop(
      "Count matrix contains no samples.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Remove/rebuild configured cohort membership columns
  # ---------------------------------------------------------------------------
  
  cohort_membership_summary <- list()
  
  
  for (
    analysis_id in names(
      ANALYSIS_CATALOG
    )
  ) {
    
    cfg <- ANALYSIS_CATALOG[[
      analysis_id
    ]]
    
    
    if (
      is.null(
        cfg$exclude_samples
      )
    ) {
      
      cfg$exclude_samples <- character(0)
    }
    
    
    cohort_samples <- setdiff(
      final_qc_samples,
      cfg$exclude_samples
    )
    
    
    membership_column <- paste0(
      "cohort_",
      analysis_id
    )
    
    
    metadata_master[[
      membership_column
    ]] <- metadata_master$Sample %in%
      cohort_samples
    
    
    # Count only samples actually belonging to the final count matrix.
    observed_n <- sum(
      metadata_master$Sample %in%
        cohort_samples,
      na.rm = TRUE
    )
    
    
    if (
      !is.na(cfg$expected_n) &&
      observed_n != cfg$expected_n
    ) {
      
      stop(
        "Reconstructed cohort '",
        analysis_id,
        "' contains ",
        observed_n,
        " samples but ANALYSIS_CATALOG expects ",
        cfg$expected_n,
        ".",
        call. = FALSE
      )
    }
    
    
    cohort_membership_summary[[
      analysis_id
    ]] <- data.frame(
      
      analysis_id = analysis_id,
      
      role = cfg$role,
      
      expected_n = cfg$expected_n,
      
      reconstructed_n = observed_n,
      
      excluded_samples = paste(
        cfg$exclude_samples,
        collapse = ";"
      ),
      
      stringsAsFactors = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Rebuild primary-analysis membership
  # ---------------------------------------------------------------------------
  
  primary_cfg <- ANALYSIS_CATALOG[[
    PRIMARY_ANALYSIS_ID
  ]]
  
  
  primary_exclusions <- primary_cfg$exclude_samples
  
  
  if (
    is.null(
      primary_exclusions
    )
  ) {
    
    primary_exclusions <- character(0)
  }
  
  
  primary_samples <- setdiff(
    final_qc_samples,
    primary_exclusions
  )
  
  
  metadata_master$primary_analysis <-
    metadata_master$Sample %in%
    primary_samples
  
  
  # ---------------------------------------------------------------------------
  # Rebuild primary-analysis exclusion reason
  # ---------------------------------------------------------------------------
  
  metadata_master$primary_exclusion_reason <- ""
  
  
  not_in_final_counts <- !metadata_master$Sample %in%
    final_qc_samples
  
  
  metadata_master$primary_exclusion_reason[
    not_in_final_counts
  ] <- "Not present in final post-sequencing QC count matrix"
  
  
  excluded_by_primary_config <-
    metadata_master$Sample %in%
    primary_exclusions
  
  
  metadata_master$primary_exclusion_reason[
    excluded_by_primary_config
  ] <- paste0(
    "Excluded by primary cohort configuration: ",
    PRIMARY_ANALYSIS_ID
  )
  
  
  # ---------------------------------------------------------------------------
  # Validate primary cohort
  # ---------------------------------------------------------------------------
  
  observed_primary_n <- sum(
    metadata_master$primary_analysis,
    na.rm = TRUE
  )
  
  
  if (
    !is.na(primary_cfg$expected_n) &&
    observed_primary_n !=
      primary_cfg$expected_n
  ) {
    
    stop(
      "Reconstructed primary cohort contains ",
      observed_primary_n,
      " samples but expected ",
      primary_cfg$expected_n,
      ".",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Attach audit table
  # ---------------------------------------------------------------------------
  
  cohort_membership_summary <- do.call(
    rbind,
    cohort_membership_summary
  )
  
  
  rownames(
    cohort_membership_summary
  ) <- NULL
  
  
  attr(
    metadata_master,
    "cohort_membership_summary"
  ) <- cohort_membership_summary
  
  
  metadata_master
}


# ==============================================================================
# 5. ANALYSIS / PROCESSING
# ==============================================================================

# ==============================================================================
# 5.1 Cohort selection
# ==============================================================================

select_analysis_cohort <- function(
    counts,
    metadata_master,
    metrics,
    context
) {
  
  validate_base_inputs(
    counts = counts,
    metadata_master = metadata_master,
    metrics = metrics
  )
  
  
  if (
    is.null(
      context$analysis_id
    ) ||
    is.null(
      context$exclude_samples
    ) ||
    is.null(
      context$expected_n
    )
  ) {
    
    stop(
      "'context' lacks analysis_id, exclude_samples or expected_n.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Rebuild all cohort-membership variables from ANALYSIS_CATALOG
  # ---------------------------------------------------------------------------
  
  metadata_master <- rebuild_cohort_membership_metadata(
    metadata_master = metadata_master,
    counts = counts
  )
  
  
  cohort_membership_summary <- attr(
    metadata_master,
    "cohort_membership_summary"
  )
  
  
  # ---------------------------------------------------------------------------
  # Standardize sample IDs
  # ---------------------------------------------------------------------------
  
  metadata_master$Sample <- as.character(
    metadata_master$Sample
  )
  
  
  metrics$Sample <- as.character(
    metrics$Sample
  )
  
  
  # ---------------------------------------------------------------------------
  # Select current analysis cohort
  # ---------------------------------------------------------------------------
  
  keep_samples <- setdiff(
    colnames(
      counts
    ),
    context$exclude_samples
  )
  
  
  if (
    !is.na(context$expected_n) &&
    length(
      keep_samples
    ) !=
      context$expected_n
  ) {
    
    stop(
      "Cohort '",
      context$analysis_id,
      "' has ",
      length(
        keep_samples
      ),
      " samples but expected ",
      context$expected_n,
      ". Check R/00_analysis_config.R and the project input files.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Align metadata
  # ---------------------------------------------------------------------------
  
  metadata <- metadata_master[
    match(
      keep_samples,
      metadata_master$Sample
    ),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # Align metrics
  # ---------------------------------------------------------------------------
  
  metrics_selected <- metrics[
    match(
      keep_samples,
      metrics$Sample
    ),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # Align counts
  # ---------------------------------------------------------------------------
  
  counts_selected <- counts[
    ,
    keep_samples,
    drop = FALSE
  ]
  
  
  if (
    anyNA(
      metadata$Sample
    ) ||
    anyNA(
      metrics_selected$Sample
    )
  ) {
    
    stop(
      "Cohort alignment failed.",
      call. = FALSE
    )
  }
  
  
  if (!identical(
    metadata$Sample,
    colnames(
      counts_selected
    )
  )) {
    
    stop(
      "Metadata/count alignment failed.",
      call. = FALSE
    )
  }
  
  
  if (!identical(
    metrics_selected$Sample,
    colnames(
      counts_selected
    )
  )) {
    
    stop(
      "Metrics/count alignment failed.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Standardize phenotype variables
  # ---------------------------------------------------------------------------
  
  metadata$Condition <- factor(
    as.character(
      metadata$Condition
    ),
    levels = c(
      "C",
      "P"
    )
  )
  
  
  metadata$Apend_type <- factor(
    as.character(
      metadata$Apend_type
    ),
    levels = c(
      "0",
      "1",
      "2"
    )
  )
  
  
  metadata$Sex <- factor(
    as.character(
      metadata$Sex
    )
  )
  
  
  metadata$Age <- suppressWarnings(
    as.numeric(
      metadata$Age
    )
  )
  
  
  metrics_selected$Reads <- suppressWarnings(
    as.numeric(
      metrics_selected$Reads
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Validate phenotype and technical variables
  # ---------------------------------------------------------------------------
  
  validate_selected_metadata(
    metadata
  )
  
  
  if (
    anyNA(
      metrics_selected$Reads
    ) ||
    any(
      !is.finite(
        metrics_selected$Reads
      )
    ) ||
    any(
      metrics_selected$Reads <= 0
    )
  ) {
    
    stop(
      "Invalid Metrics$Reads values in selected cohort.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Validate current cohort membership
  # ---------------------------------------------------------------------------
  
  current_membership_column <- paste0(
    "cohort_",
    context$analysis_id
  )
  
  
  if (
    !current_membership_column %in%
    colnames(
      metadata
    )
  ) {
    
    stop(
      "Reconstructed cohort-membership column is missing: ",
      current_membership_column,
      call. = FALSE
    )
  }
  
  
  if (
    !all(
      metadata[[
        current_membership_column
      ]]
    )
  ) {
    
    stop(
      "At least one selected sample is not marked as belonging to cohort '",
      context$analysis_id,
      "'.",
      call. = FALSE
    )
  }
  
  
  # Primary cohort must agree exactly with primary_analysis.
  
  if (
    identical(
      context$analysis_id,
      PRIMARY_ANALYSIS_ID
    ) &&
    !all(
      metadata$primary_analysis
    )
  ) {
    
    stop(
      "Selected primary cohort disagrees with reconstructed ",
      "primary_analysis membership.",
      call. = FALSE
    )
  }
  
  
  list(
    
    counts = counts_selected,
    
    metadata = metadata,
    
    metrics = metrics_selected,
    
    metadata_master = metadata_master,
    
    cohort_membership_summary =
      cohort_membership_summary
  )
}


# ==============================================================================
# 6. OUTPUTS
# ==============================================================================

# ==============================================================================
# 6.1 Graphics devices
# ==============================================================================

open_pdf_device <- function(
    filename,
    width,
    height
) {
  
  filename <- prepare_output_file(
    filename
  )
  
  
  cairo_error <- NULL
  
  
  if (
    isTRUE(
      capabilities(
        "cairo"
      )
    )
  ) {
    
    opened <- tryCatch(
      {
        
        grDevices::cairo_pdf(
          filename = filename,
          width = width,
          height = height
        )
        
        TRUE
      },
      error = function(e) {
        
        cairo_error <<- conditionMessage(
          e
        )
        
        FALSE
      }
    )
    
    
    if (
      opened
    ) {
      
      return(
        invisible(
          list(
            file = filename,
            device = "cairo_pdf"
          )
        )
      )
    }
  }
  
  
  pdf_error <- NULL
  
  
  opened <- tryCatch(
    {
      
      grDevices::pdf(
        file = filename,
        width = width,
        height = height,
        useDingbats = FALSE
      )
      
      TRUE
    },
    error = function(e) {
      
      pdf_error <<- conditionMessage(
        e
      )
      
      FALSE
    }
  )
  
  
  if (
    !opened
  ) {
    
    stop(
      "Unable to open PDF device for:\n",
      filename,
      if (
        !is.null(
          cairo_error
        )
      ) {
        paste0(
          "\nCairo error: ",
          cairo_error
        )
      } else {
        ""
      },
      if (
        !is.null(
          pdf_error
        )
      ) {
        paste0(
          "\nBase PDF error: ",
          pdf_error
        )
      } else {
        ""
      },
      call. = FALSE
    )
  }
  
  
  invisible(
    list(
      file = filename,
      device = "pdf"
    )
  )
}


open_tiff_device <- function(
    filename,
    width,
    height,
    units = "in",
    res = 600,
    compression = "lzw"
) {
  
  filename <- prepare_output_file(
    filename
  )
  
  
  grDevices::tiff(
    filename = filename,
    width = width,
    height = height,
    units = units,
    res = res,
    compression = compression
  )
  
  
  invisible(
    filename
  )
}


close_graphics_device <- function() {
  
  if (
    grDevices::dev.cur() > 1L
  ) {
    
    grDevices::dev.off()
  }
  
  
  invisible(
    NULL
  )
}


draw_to_pdf <- function(
    filename,
    width,
    height,
    draw
) {
  
  if (
    !is.function(
    draw
    )
  ) {
    
    stop(
      "'draw' must be a function.",
      call. = FALSE
    )
  }
  
  
  open_pdf_device(
    filename = filename,
    width = width,
    height = height
  )
  
  
  on.exit(
    close_graphics_device(),
    add = TRUE
  )
  
  
  draw()
  
  
  invisible(
    filename
  )
}


draw_to_tiff <- function(
    filename,
    width,
    height,
    draw,
    units = "in",
    res = 600,
    compression = "lzw"
) {
  
  if (
    !is.function(
    draw
    )
  ) {
    
    stop(
      "'draw' must be a function.",
      call. = FALSE
    )
  }
  
  
  open_tiff_device(
    filename = filename,
    width = width,
    height = height,
    units = units,
    res = res,
    compression = compression
  )
  
  
  on.exit(
    close_graphics_device(),
    add = TRUE
  )
  
  
  draw()
  
  
  invisible(
    filename
  )
}


# ==============================================================================
# 6.2 ggplot export
# ==============================================================================

save_ggplot_pdf <- function(
    plot,
    filename,
    width,
    height
) {
  
  draw_to_pdf(
    filename = filename,
    width = width,
    height = height,
    draw = function() {
      print(
        plot
      )
    }
  )
}


save_ggplot_tiff <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600,
    compression = "lzw"
) {
  
  draw_to_tiff(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    compression = compression,
    draw = function() {
      print(
        plot
      )
    }
  )
}


save_ggplot_png <- function(
    plot,
    filename,
    width,
    height,
    dpi = 300
) {
  
  filename <- prepare_output_file(
    filename
  )
  
  
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  
  
  invisible(
    filename
  )
}


# ==============================================================================
# 6.3 Analysis manifest
# ==============================================================================

write_analysis_manifest <- function(
    context,
    metadata,
    results_dir
) {
  
  ensure_dir(
    results_dir
  )
  
  
  validate_selected_metadata(
    metadata
  )
  
  
  manifest <- data.frame(
    
    analysis_id =
      context$analysis_id,
    
    role =
      context$role,
    
    description =
      context$description,
    
    n =
      nrow(
        metadata
      ),
    
    n_condition_C =
      sum(
        metadata$Condition == "C"
      ),
    
    n_condition_P =
      sum(
        metadata$Condition == "P"
      ),
    
    n_apend_type_0 =
      sum(
        metadata$Apend_type == "0"
      ),
    
    n_apend_type_1_mild =
      sum(
        metadata$Apend_type == "1"
      ),
    
    n_apend_type_2_severe =
      sum(
        metadata$Apend_type == "2"
      ),
    
    excluded_samples =
      paste(
        context$exclude_samples,
        collapse = ";"
      ),
    
    run_timestamp =
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S %z"
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  output_file <- prepare_output_file(
    file.path(
      results_dir,
      "analysis_manifest.csv"
    )
  )
  
  
  utils::write.csv(
    manifest,
    output_file,
    row.names = FALSE
  )
  
  
  manifest
}


# ==============================================================================
# 6.4 Exact cohort-sample manifest
# ==============================================================================

write_cohort_sample_manifest <- function(
    context,
    metadata,
    results_dir
) {
  
  ensure_dir(
    results_dir
  )
  
  
  validate_selected_metadata(
    metadata
  )
  
  
  preferred_columns <- c(
    "Sample",
    "Condition",
    "Apend_type",
    "Sex",
    "Age",
    "primary_analysis",
    "primary_exclusion_reason",
    paste0(
      "cohort_",
      names(
        ANALYSIS_CATALOG
      )
    )
  )
  
  
  cohort_manifest <- metadata[
    ,
    intersect(
      preferred_columns,
      colnames(
        metadata
      )
    ),
    drop = FALSE
  ]
  
  
  cohort_manifest$analysis_id <-
    context$analysis_id
  
  
  cohort_manifest$analysis_role <-
    context$role
  
  
  cohort_manifest <- cohort_manifest[
    ,
    c(
      "analysis_id",
      "analysis_role",
      setdiff(
        colnames(
          cohort_manifest
        ),
        c(
          "analysis_id",
          "analysis_role"
        )
      )
    ),
    drop = FALSE
  ]
  
  
  output_file <- prepare_output_file(
    file.path(
      results_dir,
      "cohort_samples.csv"
    )
  )
  
  
  utils::write.csv(
    cohort_manifest,
    output_file,
    row.names = FALSE
  )
  
  
  cohort_manifest
}


# ==============================================================================
# 6.5 Session information
# ==============================================================================

write_session_info <- function(
    results_dir
) {
  
  ensure_dir(
    results_dir
  )
  
  
  output_file <- prepare_output_file(
    file.path(
      results_dir,
      "sessionInfo.txt"
    )
  )
  
  
  capture.output(
    sessionInfo(),
    file = output_file
  )
  
  
  invisible(
    output_file
  )
}


# ==============================================================================
# 7. QC / AUDIT
# ==============================================================================

# ==============================================================================
# 7.1 Legacy integrated-candidate audit reader
# ==============================================================================

# IMPORTANT:
#
# This file is an output of script 06 and can be useful for documenting
# concordance across RPM, edgeR and DESeq2.
#
# Candidate definition and prioritization belong to the downstream integration
# workflow, not to this shared utility file.
#
# The function name is retained for backwards compatibility.

read_candidate_priority <- function(
    results_dir
) {
  
  candidate_file <- file.path(
    results_dir,
    "integrated",
    "downstream_candidate_priority.csv"
  )
  
  
  if (
    !file.exists(
      candidate_file
    )
  ) {
    
    stop(
      "Candidate-priority audit file not found. ",
      "Run R/06_integrated_method_comparison.R first.",
      call. = FALSE
    )
  }
  
  
  candidate_priority <- utils::read.csv(
    candidate_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  required_columns <- c(
    "miRNA",
    "downstream_category"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(
      candidate_priority
    )
  )
  
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "downstream_candidate_priority.csv lacks required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  
  candidate_priority
}


# ==============================================================================
# 8. REPRODUCIBILITY / COMMON HELPERS
# ==============================================================================

# ==============================================================================
# 8.1 Safe log2 transformation
# ==============================================================================

safe_log2 <- function(
    x
) {
  
  if (
    any(
      x < 0,
      na.rm = TRUE
    )
  ) {
    
    stop(
      "safe_log2() received negative values.",
      call. = FALSE
    )
  }
  
  
  log2(
    x + 1
  )
}


# ==============================================================================
# 8.2 Benjamini-Hochberg adjustment
# ==============================================================================

bh_adjust <- function(
    p
) {
  
  stats::p.adjust(
    p,
    method = "BH"
  )
}


# ==============================================================================
# 8.3 Holm adjustment
# ==============================================================================

holm_adjust <- function(
    p
) {
  
  stats::p.adjust(
    p,
    method = "holm"
  )
}


# ==============================================================================
# 8.4 Filename sanitization
# ==============================================================================

sanitize_filename <- function(
    x
) {
  
  original_x <- x
  
  
  x <- iconv(
    x,
    from = "",
    to = "ASCII//TRANSLIT",
    sub = "_"
  )
  
  
  x[
    is.na(
      x
    )
  ] <- original_x[
    is.na(
      x
    )
  ]
  
  
  x <- gsub(
    "[^A-Za-z0-9._-]",
    "_",
    x
  )
  
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  
  x[
    is.na(
      x
    ) |
      !nzchar(
        x
      )
  ] <- "unnamed"
  
  
  x
}


# ==============================================================================
# 8.5 Effect direction
# ==============================================================================

effect_direction <- function(
    x
) {
  
  ifelse(
    is.na(
      x
    ),
    NA_character_,
    ifelse(
      x > 0,
      "up",
      ifelse(
        x < 0,
        "down",
        "flat"
      )
    )
  )
}
