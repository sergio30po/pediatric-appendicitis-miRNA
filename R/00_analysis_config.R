# ==============================================================================
# 00_analysis_config.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Purpose:
#   Central configuration file for the complete miRNA-seq analysis pipeline.
#   This script defines the frozen analytical cohort, primary analysis, expression
#   filtering thresholds, statistical significance thresholds, covariate
#   adjustment settings, technical QC parameters, and common output paths.
#
# Role in pipeline:
#   Configuration only. No statistical analysis is performed in this script.
#   All downstream scripts must obtain cohort definitions and shared analysis
#   parameters from this file instead of redefining them locally.
#
# Primary analysis cohort:
#   Cohort membership is defined from study input data after prespecified QC.
#   No study-specific sample identifiers or observed cohort sizes are hard-coded
#   in this public methodological configuration.
#
# Main statistical framework:
#   - Expression filtering:
#       RPM >= 100 in >= 50% of samples according to the group-based filtering
#       strategy implemented downstream.
#
#   - Nominal significance:
#       p < 0.05
#
#   - Multiple-testing significance:
#       FDR < 0.05
#
#   - Main covariate-adjusted model:
#       age + sex
#
#   - Batch:
#       disabled unless a validated experimental batch variable is available.
#
#   - Hemolysis:
#       exploratory sequencing-based miRNA proxy.
#
# Inputs:
#   - Repository structure
#   - Environment variables:
#       MIRNA_PROJECT_ROOT
#       MIRNA_ANALYSIS_ID
#       MIRNA_BATCH_COLUMN
#
# Outputs:
#   This script does not generate analytical result files.
#   It provides configuration objects and helper functions used downstream.
#
# Required packages:
#   here
#
# Reproducibility:
#   This file is the single source of truth for global analysis parameters.
#   Sample exclusions, primary-cohort definitions, filtering thresholds and
#   top-level output paths must not be independently modified in downstream
#   scripts.
#
# ==============================================================================


# ==============================================================================
# 1. DEPENDENCIES
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}


# ==============================================================================
# 2. PROJECT ROOT
# ==============================================================================

# Optional override for automated environments or CI.
#
# In normal interactive use, the repository root is discovered using
# here::here(), allowing the project to be cloned or moved without editing
# absolute paths.

project_root_override <- Sys.getenv(
  "MIRNA_PROJECT_ROOT",
  unset = ""
)


if (nzchar(project_root_override)) {
  
  project_root <- normalizePath(
    project_root_override,
    winslash = "/",
    mustWork = TRUE
  )
  
} else {
  
  project_root <- normalizePath(
    here::here(),
    winslash = "/",
    mustWork = TRUE
  )
}


# Validate minimum expected repository structure.

required_project_dirs <- c(
  file.path(project_root, "R"),
  file.path(project_root, "data")
)


missing_project_dirs <- required_project_dirs[
  !dir.exists(required_project_dirs)
]


if (length(missing_project_dirs) > 0L) {
  
  stop(
    "Project root is not valid. Missing required directory/directories:\n",
    paste(
      "-",
      missing_project_dirs,
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. ANALYSIS COHORTS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Cohort strategy
# ------------------------------------------------------------------------------
#
# Public methodological configuration:
#
#   primary_analysis
#   role = primary
#
# Cohort membership is reconstructed from the supplied post-QC input data.
#
# ------------------------------------------------------------------------------


ANALYSIS_CATALOG <- list(
  primary_analysis = list(
    role = "primary",
    description = "Primary analysis cohort reconstructed from post-QC input data.",
    exclude_samples = character(0),
    expected_n = NA_integer_
  )
)


# ==============================================================================
# 4. PRIMARY ANALYSIS
# ==============================================================================

# Final manuscript primary cohort.
#
# All downstream scripts using get_analysis_id() without an explicit
# MIRNA_ANALYSIS_ID environment override will run this cohort.

PRIMARY_ANALYSIS_ID <- "primary_analysis"


# ==============================================================================
# 5. EXPRESSION FILTERING
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1 Main filtering rule
# ------------------------------------------------------------------------------

# Minimum abundance threshold.

FILTER_MIN_RPM <- 100


# Required prevalence.

FILTER_PREVALENCE <- 0.50


# Clinical grouping variable used by the filtering procedure.

FILTER_GROUP_VARIABLE <- "Apend_type"


# ------------------------------------------------------------------------------
# 5.2 Filtering robustness grid
# ------------------------------------------------------------------------------

# Prespecified RPM thresholds for filtering robustness analyses.

FILTER_MIN_RPM_GRID <- c(
  10,
  50,
  100,
  200,
  500
)


# Prespecified prevalence thresholds for filtering robustness analyses.

FILTER_PREVALENCE_GRID <- c(
  0.25,
  0.50,
  0.75
)


# ==============================================================================
# 6. STATISTICAL SIGNIFICANCE THRESHOLDS
# ==============================================================================

# Nominal significance threshold.

ALPHA_RAW <- 0.05


# Benjamini-Hochberg false discovery rate threshold.

ALPHA_FDR <- 0.05


# ==============================================================================
# 7. COVARIATE ADJUSTMENT
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Main adjustment model
# ------------------------------------------------------------------------------

# Age and sex are available in the current metadata and constitute the
# prespecified basic adjusted model.

ADJUST_FOR_AGE_SEX <- TRUE


# ------------------------------------------------------------------------------
# 7.2 Technical batch
# ------------------------------------------------------------------------------

# Batch adjustment remains disabled unless a validated experimental batch
# variable is explicitly available in metadata_master.xlsx.
#
# Batch must not be inferred from sample identifiers, file order, sequencing
# order or other undocumented surrogate variables.

BATCH_COLUMN <- Sys.getenv(
  "MIRNA_BATCH_COLUMN",
  unset = ""
)


if (!nzchar(BATCH_COLUMN)) {
  BATCH_COLUMN <- NA_character_
}


# ==============================================================================
# 8. TECHNICAL QC PARAMETERS
# ==============================================================================

# ------------------------------------------------------------------------------
# 8.1 Hemolysis proxy
# ------------------------------------------------------------------------------

# Exploratory sequencing-based hemolysis proxy.
#
# This is not a qPCR delta-Cq analysis and qPCR-specific clinical thresholds
# must not be transferred directly to these sequencing measurements.

HEMOLYSIS_MIRNA_RBC <- "hsa-miR-451a"

HEMOLYSIS_MIRNA_REFERENCE <- "hsa-miR-23a-3p"


# ==============================================================================
# 9. ANALYSIS-ID RESOLVER
# ==============================================================================

get_analysis_id <- function() {
  
  analysis_id <- Sys.getenv(
    "MIRNA_ANALYSIS_ID",
    unset = PRIMARY_ANALYSIS_ID
  )
  
  
  if (!analysis_id %in% names(ANALYSIS_CATALOG)) {
    
    stop(
      "Unknown MIRNA_ANALYSIS_ID='",
      analysis_id,
      "'. Valid values: ",
      paste(
        names(ANALYSIS_CATALOG),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  analysis_id
}


# ==============================================================================
# 10. ANALYSIS CONTEXT
# ==============================================================================

get_analysis_context <- function(
    analysis_id = get_analysis_id()
) {
  
  cfg <- ANALYSIS_CATALOG[[analysis_id]]
  
  
  # ---------------------------------------------------------------------------
  # 10.1 Output directories
  # ---------------------------------------------------------------------------
  
  results_dir <- file.path(
    project_root,
    "results",
    analysis_id
  )
  
  
  figures_dir <- file.path(
    project_root,
    "figures",
    analysis_id
  )
  
  
  dir.create(
    results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    figures_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 10.2 Directory validation
  # ---------------------------------------------------------------------------
  
  if (!dir.exists(results_dir)) {
    
    stop(
      "Could not create results directory: ",
      results_dir,
      call. = FALSE
    )
  }
  
  
  if (!dir.exists(figures_dir)) {
    
    stop(
      "Could not create figures directory: ",
      figures_dir,
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 10.3 Return analysis context
  # ---------------------------------------------------------------------------
  
  list(
    analysis_id = analysis_id,
    role = cfg$role,
    description = cfg$description,
    exclude_samples = cfg$exclude_samples,
    expected_n = cfg$expected_n,
    project_root = project_root,
    results_dir = results_dir,
    figures_dir = figures_dir
  )
}


# ==============================================================================
# 11. CONFIGURATION VALIDATION
# ==============================================================================

# Validate primary analysis identifier.

if (!PRIMARY_ANALYSIS_ID %in% names(ANALYSIS_CATALOG)) {
  
  stop(
    "PRIMARY_ANALYSIS_ID is not present in ANALYSIS_CATALOG.",
    call. = FALSE
  )
}


# Validate analytical-cohort expected sample counts.

expected_n_values <- vapply(
  ANALYSIS_CATALOG,
  function(x) x$expected_n,
  integer(1)
)


if (any(!is.na(expected_n_values) & expected_n_values <= 0L)) {
  stop(
    "Configured expected_n values must be positive when specified.",
    call. = FALSE
  )
}


# Validate filtering parameters.

if (
  !is.numeric(FILTER_MIN_RPM) ||
  length(FILTER_MIN_RPM) != 1L ||
  !is.finite(FILTER_MIN_RPM) ||
  FILTER_MIN_RPM < 0
) {
  
  stop(
    "FILTER_MIN_RPM must be a single non-negative numeric value.",
    call. = FALSE
  )
}


if (
  !is.numeric(FILTER_PREVALENCE) ||
  length(FILTER_PREVALENCE) != 1L ||
  !is.finite(FILTER_PREVALENCE) ||
  FILTER_PREVALENCE <= 0 ||
  FILTER_PREVALENCE > 1
) {
  
  stop(
    "FILTER_PREVALENCE must be in the interval (0, 1].",
    call. = FALSE
  )
}


# Validate significance thresholds.

if (
  ALPHA_RAW <= 0 ||
  ALPHA_RAW >= 1 ||
  ALPHA_FDR <= 0 ||
  ALPHA_FDR >= 1
) {
  
  stop(
    "ALPHA_RAW and ALPHA_FDR must be between 0 and 1.",
    call. = FALSE
  )
}
