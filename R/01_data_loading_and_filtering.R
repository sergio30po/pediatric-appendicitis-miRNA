# ==============================================================================
# 01_data_loading_and_filtering.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Script:
#   01_data_loading_and_filtering.R
#
# Purpose:
#   Load and validate the project input data, construct the selected analysis
#   cohort, define the canonical clinical-group labels, normalize raw miRNA
#   counts to reads per million (RPM), apply the prespecified
#   abundance/prevalence filter, and generate the frozen input objects used by
#   all downstream analyses.
#
# Role in pipeline:
#   First analytical step.
#
#   This script:
#     - loads raw counts, clinical metadata and sequencing metrics;
#     - applies the cohort definition from 00_analysis_config.R;
#     - rebuilds cohort-membership variables dynamically from ANALYSIS_CATALOG;
#     - validates sample identities and phenotype coding;
#     - creates canonical clinical-group labels used by all downstream scripts;
#     - audits cohort membership and phenotype coding;
#     - verifies the RPM normalization denominator;
#     - normalizes the COMPLETE count matrix to RPM;
#     - applies the prespecified group-aware expression filter;
#     - evaluates filtering sensitivity across the predefined threshold grid;
#     - preserves matched filtered raw counts for edgeR and DESeq2;
#     - generates standardized matrices for visualization;
#     - records complete sample and cohort traceability.
#
# Primary analysis cohort:
#   The cohort is reconstructed from the supplied post-QC input data using the
#   phenotype coding and validation rules below.
#
# Canonical derived clinical variables:
#
#   Phenotype:
#       NIA
#       Mild
#       Severe
#
#   Phenotype_full:
#       Non-inflamed appendix
#       Mild appendicitis
#       Severe appendicitis
#
#   Appendicitis_status:
#       NIA
#       Appendicitis
#
# IMPORTANT:
#   NIA = non-inflamed appendix.
#
#   The NIA group is NOT a healthy-control cohort. It represents symptomatic
#   surgical patients whose appendix showed no evidence of acute appendiceal
#   inflammation.
#
#   Apend_type remains the original coded variable (0/1/2) and is never
#   overwritten.
#
#   Condition remains the original source variable (C/P) for traceability but
#   should not be used downstream as a display label for clinical groups.
#
# Cohort handling:
#   Study-specific sample identifiers are not hard-coded in this public script.
#
# Normalization:
#   RPM = raw count / complete miRNA library size * 1,000,000
#
#   The denominator is calculated from the COMPLETE raw count matrix before
#   filtering. Filtering must never be performed before RPM normalization.
#
# Main expression filter:
#   Minimum abundance:
#       RPM >= 100
#
#   Minimum prevalence:
#       >= 50% of samples
#
#   Scope:
#       any Apend_type group
#
#   Minimum groups:
#       1
#
#   The prevalence requirement is evaluated from the observed size of each
#   phenotype group. A miRNA is retained when it meets the criterion in at least
#   one phenotype group.
#
# Inputs:
#   data/counts.tsv
#   data/metadata_master.xlsx
#   data/Metrics.xlsx
#
# Configuration:
#   R/00_analysis_config.R
#
# Utilities:
#   R/00_utils.R
#
# Main outputs:
#   results/<analysis_id>/
#
#     analysis_manifest.csv
#     cohort_samples.csv
#     cohort_group_summary.csv
#     cohort_membership_summary.csv
#     special_sample_audit.csv
#     sample_input_traceability.csv
#
#     RPM_denominator_check.csv
#     RPM_normalization_manifest.csv
#
#     RPM_filter_group_requirements.csv
#     RPM_filter_sensitivity_summary.csv
#     RPM_filter_diagnostics.csv
#     filter_manifest.csv
#
#     miRNA_ftd.csv
#     raw_ftd_counts.csv
#     scaled_data.csv
#
#     counts_analysis_samples.rds
#     metadata_analysis_samples.rds
#     Metrics_analysis_samples.rds
#     RPM_counts_all_miRNAs.rds
#     filtered_RPM.rds
#     raw_ftd_counts.rds
#     scaled_data.rds
#     RPM_object_unscaled.rds
#     RPM_object_scaled.rds
#
#     01_final_QC_summary.csv
#     sessionInfo.txt
#
# Methodological decisions:
#   - Raw counts are never replaced by RPM values in count-based models.
#   - edgeR and DESeq2 downstream analyses must use raw_ftd_counts.
#   - RPM-based analyses use filtered_RPM / miRNA_ftd.
#   - Filtering is independent of downstream statistical significance.
#   - Cohort-membership columns from metadata_master.xlsx are not trusted as
#     analytical truth; they are rebuilt from ANALYSIS_CATALOG.
#   - Clinical display labels are derived once here and inherited downstream.
#
# Reproducibility:
#   All sample inclusion/exclusion decisions, cohort definitions, clinical
#   labels and filtering parameters are recorded explicitly.
#
# ==============================================================================


# ==============================================================================
# 1. DEPENDENCIES
# ==============================================================================

if (!requireNamespace(
  "here",
  quietly = TRUE
)) {
  
  stop(
    "Package 'here' is required. Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}


source(
  here::here(
    "R",
    "00_analysis_config.R"
  )
)


source(
  here::here(
    "R",
    "00_utils.R"
  )
)


if (!requireNamespace(
  "miRPM",
  quietly = TRUE
)) {
  
  stop(
    "Package 'miRPM' is required. Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}


if (
  utils::compareVersion(
    as.character(
      utils::packageVersion("miRPM")
    ),
    "0.2.0"
  ) < 0
) {
  
  stop(
    "miRPM >= 0.2.0 is required. ",
    "Run R/00_packages_setup.R to update it.",
    call. = FALSE
  )
}


# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

context <- get_analysis_context()


results_dir <- context$results_dir

figures_dir <- context$figures_dir


# ------------------------------------------------------------------------------
# Canonical clinical-group labels
# ------------------------------------------------------------------------------

# These labels are created once in script 01 and are inherited by all
# downstream scripts through metadata_analysis_samples.rds.

APEND_TYPE_SHORT_LABELS <- c(
  "0" = "NIA",
  "1" = "Mild",
  "2" = "Severe"
)


APEND_TYPE_FULL_LABELS <- c(
  "0" = "Non-inflamed appendix",
  "1" = "Mild appendicitis",
  "2" = "Severe appendicitis"
)


APEND_TYPE_ORDER <- c(
  "0",
  "1",
  "2"
)


PHENOTYPE_ORDER <- unname(
  APEND_TYPE_SHORT_LABELS[
    APEND_TYPE_ORDER
  ]
)


PHENOTYPE_FULL_ORDER <- unname(
  APEND_TYPE_FULL_LABELS[
    APEND_TYPE_ORDER
  ]
)


APPENDICITIS_STATUS_ORDER <- c(
  "NIA",
  "Appendicitis"
)


message(
  "\n============================================================"
)

message(
  "01_data_loading_and_filtering.R"
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
  "Description: ",
  context$description
)

message(
  "Expected n: ",
  context$expected_n
)

message(
  "============================================================\n"
)


# ==============================================================================
# 3. INPUTS
# ==============================================================================

# ==============================================================================
# 3.1 Load project inputs
# ==============================================================================

inputs <- read_project_inputs()


# ==============================================================================
# 3.2 Select analysis cohort
# ==============================================================================

cohort <- select_analysis_cohort(
  counts = inputs$counts,
  metadata_master = inputs$metadata_master,
  metrics = inputs$metrics,
  context = context
)


counts <- cohort$counts

metadata <- cohort$metadata

Metrics <- cohort$metrics


# ==============================================================================
# 3.3 Define canonical clinical-group labels
# ==============================================================================

add_clinical_group_labels <- function(
    metadata_df
) {
  
  if (
    !"Apend_type" %in%
    colnames(
      metadata_df
    )
  ) {
    
    stop(
      "Cannot derive clinical-group labels: Apend_type is missing.",
      call. = FALSE
    )
  }
  
  
  apend_chr <- as.character(
    metadata_df$Apend_type
  )
  
  
  unknown_apend_types <- setdiff(
    unique(
      apend_chr[
        !is.na(
          apend_chr
        )
      ]
    ),
    APEND_TYPE_ORDER
  )
  
  
  if (
    length(
      unknown_apend_types
    ) > 0L
  ) {
    
    stop(
      "Unknown Apend_type value(s) while creating clinical labels: ",
      paste(
        unknown_apend_types,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  short_labels <- unname(
    APEND_TYPE_SHORT_LABELS[
      apend_chr
    ]
  )
  
  
  full_labels <- unname(
    APEND_TYPE_FULL_LABELS[
      apend_chr
    ]
  )
  
  
  if (
    anyNA(
      short_labels
    ) ||
    anyNA(
      full_labels
    )
  ) {
    
    stop(
      "Failed to assign clinical-group labels to all samples.",
      call. = FALSE
    )
  }
  
  
  metadata_df$Phenotype <- factor(
    short_labels,
    levels = PHENOTYPE_ORDER
  )
  
  
  metadata_df$Phenotype_full <- factor(
    full_labels,
    levels = PHENOTYPE_FULL_ORDER
  )
  
  
  metadata_df$Appendicitis_status <- factor(
    ifelse(
      apend_chr == "0",
      "NIA",
      "Appendicitis"
    ),
    levels = APPENDICITIS_STATUS_ORDER
  )
  
  
  metadata_df
}


metadata <- add_clinical_group_labels(
  metadata
)


# Also reconstruct the labels in the full traceability metadata object so that
# exported audit files do not contain ambiguous clinical-group terminology.

cohort$metadata_master <- add_clinical_group_labels(
  cohort$metadata_master
)


message(
  "Canonical clinical labels created: ",
  paste(
    PHENOTYPE_ORDER,
    collapse = " / "
  )
)


# ==============================================================================
# 4. VALIDATION
# ==============================================================================

# ==============================================================================
# 4.1 Basic dimensions
# ==============================================================================

if (
  ncol(counts) != nrow(metadata) ||
  ncol(counts) != nrow(Metrics)
) {
  
  stop(
    "Selected count, metadata and metrics objects have inconsistent ",
    "sample dimensions.",
    call. = FALSE
  )
}


if (!identical(
  colnames(counts),
  metadata$Sample
)) {
  
  stop(
    "Count matrix and metadata are not identically sample-aligned.",
    call. = FALSE
  )
}


if (!identical(
  colnames(counts),
  Metrics$Sample
)) {
  
  stop(
    "Count matrix and sequencing metrics are not identically sample-aligned.",
    call. = FALSE
  )
}


# ==============================================================================
# 4.2 Clinical-label consistency
# ==============================================================================

expected_short <- unname(
  APEND_TYPE_SHORT_LABELS[
    as.character(
      metadata$Apend_type
    )
  ]
)


expected_full <- unname(
  APEND_TYPE_FULL_LABELS[
    as.character(
      metadata$Apend_type
    )
  ]
)


if (!identical(
  as.character(
    metadata$Phenotype
  ),
  expected_short
)) {
  
  stop(
    "Phenotype labels are inconsistent with Apend_type.",
    call. = FALSE
  )
}


if (!identical(
  as.character(
    metadata$Phenotype_full
  ),
  expected_full
)) {
  
  stop(
    "Phenotype_full labels are inconsistent with Apend_type.",
    call. = FALSE
  )
}


expected_status <- ifelse(
  as.character(
    metadata$Apend_type
  ) == "0",
  "NIA",
  "Appendicitis"
)


if (!identical(
  as.character(
    metadata$Appendicitis_status
  ),
  expected_status
)) {
  
  stop(
    "Appendicitis_status is inconsistent with Apend_type.",
    call. = FALSE
  )
}


# ==============================================================================
# 4.3 Exact cohort manifest
# ==============================================================================

analysis_manifest <- write_analysis_manifest(
  context = context,
  metadata = metadata,
  results_dir = results_dir
)


cohort_sample_manifest <- write_cohort_sample_manifest(
  context = context,
  metadata = metadata,
  results_dir = results_dir
)


# ==============================================================================
# 4.4 Phenotype-group audit
# ==============================================================================

observed_apend_counts <- table(
  factor(
    metadata$Apend_type,
    levels = APEND_TYPE_ORDER
  )
)


cohort_group_summary <- data.frame(
  
  analysis_id = context$analysis_id,
  
  Apend_type = APEND_TYPE_ORDER,
  
  Phenotype = unname(
    APEND_TYPE_SHORT_LABELS[
      APEND_TYPE_ORDER
    ]
  ),
  
  Phenotype_full = unname(
    APEND_TYPE_FULL_LABELS[
      APEND_TYPE_ORDER
    ]
  ),
  
  n = as.integer(
    observed_apend_counts
  ),
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  cohort_group_summary,
  prepare_output_file(
    file.path(
      results_dir,
      "cohort_group_summary.csv"
    )
  ),
  row.names = FALSE
)


message(
  "Observed phenotype groups: ",
  paste(
    paste0(
      cohort_group_summary$Phenotype,
      "=",
      cohort_group_summary$n
    ),
    collapse = "; "
  )
)


# ==============================================================================
# 4.5 Primary-analysis cohort audit
# ==============================================================================

if (
  identical(context$analysis_id, PRIMARY_ANALYSIS_ID) &&
  identical(context$role, "primary")
) {
  message(
    "Primary analysis cohort validated: ",
    nrow(metadata),
    " samples."
  )
}


# ==============================================================================
# 4.6 Global sample traceability
# ==============================================================================

# IMPORTANT:
# Use the metadata reconstructed by select_analysis_cohort(), not the original
# historical cohort flags stored in metadata_master.xlsx.

traceability <- cohort$metadata_master


traceability$Sample <- as.character(
  traceability$Sample
)


traceability$present_in_counts <- traceability$Sample %in%
  colnames(
    inputs$counts
  )


traceability$present_in_metrics <- traceability$Sample %in%
  as.character(
    inputs$metrics$Sample
  )


traceability$excluded_by_analysis_config <- traceability$Sample %in%
  context$exclude_samples


traceability$included_in_current_analysis <- traceability$Sample %in%
  metadata$Sample


traceability$current_analysis_id <- context$analysis_id


traceability$current_analysis_status <- ifelse(
  
  traceability$included_in_current_analysis,
  
  "INCLUDED",
  
  ifelse(
    
    traceability$excluded_by_analysis_config,
    
    "EXCLUDED_BY_ANALYSIS_CONFIG",
    
    ifelse(
      
      !traceability$present_in_counts,
      
      "NOT_PRESENT_IN_COUNT_MATRIX",
      
      "NOT_INCLUDED_OTHER_REASON"
    )
  )
)


utils::write.csv(
  traceability,
  prepare_output_file(
    file.path(
      results_dir,
      "sample_input_traceability.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 4.8 Cohort-membership audit
# ==============================================================================

utils::write.csv(
  cohort$cohort_membership_summary,
  prepare_output_file(
    file.path(
      results_dir,
      "cohort_membership_summary.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 5. ANALYSIS / PROCESSING
# ==============================================================================

# ==============================================================================
# 5.1 Verify RPM denominator
# ==============================================================================

count_library_size <- colSums(
  counts
)


if (
  anyNA(
    count_library_size
  ) ||
  any(
    !is.finite(
      count_library_size
    )
  ) ||
  any(
    count_library_size <= 0
  )
) {
  
  stop(
    "At least one sample has an invalid count-derived library size.",
    call. = FALSE
  )
}


names(
  count_library_size
) <- colnames(
  counts
)


metric_library_size <- Metrics$Reads


names(
  metric_library_size
) <- Metrics$Sample


metric_library_size <- metric_library_size[
  colnames(
    counts
  )
]


library_size_delta <- as.numeric(
  metric_library_size -
    count_library_size
)


library_size_relative_delta <- library_size_delta /
  count_library_size


denominator_check <- data.frame(
  
  Sample = colnames(
    counts
  ),
  
  Metrics_Reads = as.numeric(
    metric_library_size
  ),
  
  miRNA_count_sum = as.numeric(
    count_library_size
  ),
  
  difference = library_size_delta,
  
  relative_difference = library_size_relative_delta,
  
  exact_match = abs(
    library_size_delta
  ) < 1e-8,
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  denominator_check,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_denominator_check.csv"
    )
  ),
  row.names = FALSE
)


if (!all(
  denominator_check$exact_match
)) {
  
  warning(
    "Metrics$Reads is not identical to the complete miRNA count-matrix ",
    "column sum for all samples. RPM normalization will use the explicitly ",
    "calculated complete count-matrix library sizes. Inspect ",
    "RPM_denominator_check.csv.",
    call. = FALSE
  )
  
} else {
  
  message(
    "Metrics$Reads exactly matches the complete quantified miRNA count ",
    "library size for every selected sample."
  )
}


# ==============================================================================
# 5.2 RPM normalization
# ==============================================================================

RPM_counts <- miRPM::normalize_rpm(
  
  count_matrix = counts,
  
  library_sizes = count_library_size
)


if (
  anyNA(
    RPM_counts
  ) ||
  any(
    !is.finite(
      RPM_counts
    )
  ) ||
  any(
    RPM_counts < 0
  )
) {
  
  stop(
    "RPM normalization produced invalid values.",
    call. = FALSE
  )
}


if (!identical(
  dim(
    RPM_counts
  ),
  dim(
    counts
  )
)) {
  
  stop(
    "RPM normalization changed matrix dimensions unexpectedly.",
    call. = FALSE
  )
}


if (!identical(
  dimnames(
    RPM_counts
  ),
  dimnames(
    counts
  )
)) {
  
  stop(
    "RPM normalization changed matrix dimnames unexpectedly.",
    call. = FALSE
  )
}


# ==============================================================================
# 5.3 Validate RPM sums
# ==============================================================================

rpm_column_sums <- colSums(
  RPM_counts
)


rpm_sum_difference <- rpm_column_sums -
  1e6


if (
  any(
    abs(
      rpm_sum_difference
    ) > 1e-4
  )
) {
  
  stop(
    "RPM normalization validation failed: at least one sample does not ",
    "sum approximately to 1,000,000 RPM.",
    call. = FALSE
  )
}


normalization_manifest <- data.frame(
  
  Sample = colnames(
    counts
  ),
  
  raw_library_size = as.numeric(
    count_library_size
  ),
  
  Metrics_Reads = as.numeric(
    metric_library_size
  ),
  
  metrics_minus_count_sum = library_size_delta,
  
  Metrics_matches_count_sum =
    denominator_check$exact_match,
  
  RPM_sum = as.numeric(
    rpm_column_sums
  ),
  
  RPM_sum_minus_1e6 = as.numeric(
    rpm_sum_difference
  ),
  
  normalization_denominator =
    "complete_count_matrix_colSum",
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  normalization_manifest,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_normalization_manifest.csv"
    )
  ),
  row.names = FALSE
)


message(
  "RPM normalization validated: all samples sum to approximately 1,000,000."
)


# ==============================================================================
# 5.4 Filtering function
# ==============================================================================

run_miRPM_filter <- function(
    rpm,
    metadata,
    min_rpm,
    prevalence,
    group_var
) {
  
  miRPM::filter_mirnas(
    
    rpm_matrix = rpm,
    
    metadata = metadata,
    
    min_rpm = min_rpm,
    
    min_prevalence = prevalence,
    
    group_column = group_var,
    
    prevalence_scope = "any_group",
    
    min_groups = 1,
    
    return_diagnostics = TRUE
  )
}


# ==============================================================================
# 5.5 Filter group requirements
# ==============================================================================

filter_group_values <- as.character(
  metadata[[
    FILTER_GROUP_VARIABLE
  ]]
)


filter_group_sizes <- table(
  filter_group_values
)


filter_group_requirements <- data.frame(
  
  group = names(
    filter_group_sizes
  ),
  
  phenotype = unname(
    APEND_TYPE_SHORT_LABELS[
      names(
        filter_group_sizes
      )
    ]
  ),
  
  phenotype_full = unname(
    APEND_TYPE_FULL_LABELS[
      names(
        filter_group_sizes
      )
    ]
  ),
  
  n_samples = as.integer(
    filter_group_sizes
  ),
  
  minimum_RPM = FILTER_MIN_RPM,
  
  minimum_prevalence = FILTER_PREVALENCE,
  
  required_samples = ceiling(
    as.integer(
      filter_group_sizes
    ) *
      FILTER_PREVALENCE
  ),
  
  prevalence_scope = "any_group",
  
  minimum_groups_required = 1L,
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  filter_group_requirements,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_filter_group_requirements.csv"
    )
  ),
  row.names = FALSE
)


message(
  "Main RPM filter requirements: ",
  paste(
    paste0(
      filter_group_requirements$phenotype,
      ": ",
      filter_group_requirements$required_samples,
      "/",
      filter_group_requirements$n_samples,
      " samples >= ",
      FILTER_MIN_RPM,
      " RPM"
    ),
    collapse = "; "
  )
)


# ==============================================================================
# 5.6 Filter sensitivity grid
# ==============================================================================

sensitivity_rows <- list()


row_index <- 0L


for (
  min_rpm in FILTER_MIN_RPM_GRID
) {
  
  for (
    prevalence in FILTER_PREVALENCE_GRID
  ) {
    
    row_index <- row_index + 1L
    
    
    tmp_filter <- run_miRPM_filter(
      
      rpm = RPM_counts,
      
      metadata = metadata,
      
      min_rpm = min_rpm,
      
      prevalence = prevalence,
      
      group_var = FILTER_GROUP_VARIABLE
    )
    
    
    sensitivity_rows[[
      row_index
    ]] <- data.frame(
      
      min_rpm = min_rpm,
      
      prevalence = prevalence,
      
      retained_miRNAs =
        tmp_filter$summary$retained_mirnas,
      
      removed_miRNAs =
        tmp_filter$summary$removed_mirnas,
      
      retained_fraction =
        tmp_filter$summary$retained_fraction,
      
      retained_pct = round(
        100 *
          tmp_filter$summary$retained_fraction,
        2
      ),
      
      only_one_group = sum(
        tmp_filter$diagnostics$groups_passing == 1L
      ),
      
      two_or_more_groups = sum(
        tmp_filter$diagnostics$groups_passing >= 2L
      ),
      
      stringsAsFactors = FALSE
    )
  }
}


filter_sensitivity <- do.call(
  rbind,
  sensitivity_rows
)


filter_sensitivity <- filter_sensitivity[
  order(
    filter_sensitivity$min_rpm,
    filter_sensitivity$prevalence
  ),
  ,
  drop = FALSE
]


rownames(
  filter_sensitivity
) <- NULL


utils::write.csv(
  filter_sensitivity,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_filter_sensitivity_summary.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 5.7 Prespecified manuscript filter
# ==============================================================================

selected_filter <- run_miRPM_filter(
  
  rpm = RPM_counts,
  
  metadata = metadata,
  
  min_rpm = FILTER_MIN_RPM,
  
  prevalence = FILTER_PREVALENCE,
  
  group_var = FILTER_GROUP_VARIABLE
)


miRNA_ftd <- selected_filter$filtered_matrix


if (
  nrow(
    miRNA_ftd
  ) == 0L
) {
  
  stop(
    "The prespecified RPM filter retained zero miRNAs.",
    call. = FALSE
  )
}


# ==============================================================================
# 5.8 Preserve matched raw counts
# ==============================================================================

raw_ftd_counts <- counts[
  rownames(
    miRNA_ftd
  ),
  colnames(
    miRNA_ftd
  ),
  drop = FALSE
]


if (!identical(
  colnames(
    miRNA_ftd
  ),
  colnames(
    raw_ftd_counts
  )
)) {
  
  stop(
    "Filtered RPM and raw-count matrices are not sample-aligned.",
    call. = FALSE
  )
}


if (!identical(
  rownames(
    miRNA_ftd
  ),
  rownames(
    raw_ftd_counts
  )
)) {
  
  stop(
    "Filtered RPM and raw-count matrices are not miRNA-aligned.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. OUTPUTS
# ==============================================================================

# ==============================================================================
# 6.1 Filter diagnostics
# ==============================================================================

utils::write.csv(
  selected_filter$diagnostics,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_filter_diagnostics.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 6.2 Filter manifest
# ==============================================================================

filter_manifest <- data.frame(
  
  analysis_id = context$analysis_id,
  
  filter_group_variable =
    selected_filter$parameters$group_column,
  
  minimum_RPM =
    selected_filter$parameters$min_rpm,
  
  prevalence =
    selected_filter$parameters$min_prevalence,
  
  prevalence_scope =
    selected_filter$parameters$prevalence_scope,
  
  min_groups =
    selected_filter$parameters$min_groups,
  
  initial_miRNAs =
    selected_filter$summary$input_mirnas,
  
  retained_miRNAs =
    selected_filter$summary$retained_mirnas,
  
  removed_miRNAs =
    selected_filter$summary$removed_mirnas,
  
  retained_fraction =
    selected_filter$summary$retained_fraction,
  
  retained_percent =
    100 *
    selected_filter$summary$retained_fraction,
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  filter_manifest,
  prepare_output_file(
    file.path(
      results_dir,
      "filter_manifest.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 6.3 Visualization scaling
# ==============================================================================

scaled_data <- t(
  scale(
    t(
      miRNA_ftd
    )
  )
)


scaled_data[
  !is.finite(
    scaled_data
  )
] <- 0


storage.mode(
  scaled_data
) <- "numeric"


# ==============================================================================
# 6.4 Human-readable matrices
# ==============================================================================

utils::write.csv(
  miRNA_ftd,
  prepare_output_file(
    file.path(
      results_dir,
      "miRNA_ftd.csv"
    )
  ),
  row.names = TRUE
)


utils::write.csv(
  raw_ftd_counts,
  prepare_output_file(
    file.path(
      results_dir,
      "raw_ftd_counts.csv"
    )
  ),
  row.names = TRUE
)


utils::write.csv(
  scaled_data,
  prepare_output_file(
    file.path(
      results_dir,
      "scaled_data.csv"
    )
  ),
  row.names = TRUE
)


# ==============================================================================
# 6.5 RDS analysis objects
# ==============================================================================

saveRDS(
  counts,
  prepare_output_file(
    file.path(
      results_dir,
      "counts_analysis_samples.rds"
    )
  )
)


saveRDS(
  metadata,
  prepare_output_file(
    file.path(
      results_dir,
      "metadata_analysis_samples.rds"
    )
  )
)


saveRDS(
  Metrics,
  prepare_output_file(
    file.path(
      results_dir,
      "Metrics_analysis_samples.rds"
    )
  )
)


saveRDS(
  RPM_counts,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_counts_all_miRNAs.rds"
    )
  )
)


saveRDS(
  miRNA_ftd,
  prepare_output_file(
    file.path(
      results_dir,
      "filtered_RPM.rds"
    )
  )
)


saveRDS(
  raw_ftd_counts,
  prepare_output_file(
    file.path(
      results_dir,
      "raw_ftd_counts.rds"
    )
  )
)


saveRDS(
  scaled_data,
  prepare_output_file(
    file.path(
      results_dir,
      "scaled_data.rds"
    )
  )
)


# ==============================================================================
# 6.6 Lightweight analysis objects
# ==============================================================================

RPM_object_unscaled <- list(
  
  assay = miRNA_ftd,
  
  colData = metadata,
  
  rowData = data.frame(
    miRNA = rownames(
      miRNA_ftd
    ),
    row.names = rownames(
      miRNA_ftd
    ),
    stringsAsFactors = FALSE
  )
)


RPM_object_scaled <- list(
  
  assay = scaled_data,
  
  colData = metadata,
  
  rowData = data.frame(
    miRNA = rownames(
      scaled_data
    ),
    row.names = rownames(
      scaled_data
    ),
    stringsAsFactors = FALSE
  )
)


saveRDS(
  RPM_object_unscaled,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_object_unscaled.rds"
    )
  )
)


saveRDS(
  RPM_object_scaled,
  prepare_output_file(
    file.path(
      results_dir,
      "RPM_object_scaled.rds"
    )
  )
)


# ==============================================================================
# 7. QC / AUDIT
# ==============================================================================

# ==============================================================================
# 7.1 Final dimensions and group labels
# ==============================================================================

final_qc <- data.frame(
  
  analysis_id = context$analysis_id,
  
  n_samples = ncol(
    counts
  ),
  
  n_NIA = sum(
    metadata$Phenotype == "NIA"
  ),
  
  n_Mild = sum(
    metadata$Phenotype == "Mild"
  ),
  
  n_Severe = sum(
    metadata$Phenotype == "Severe"
  ),
  
  n_Appendicitis = sum(
    metadata$Appendicitis_status == "Appendicitis"
  ),
  
  total_input_miRNAs = nrow(
    RPM_counts
  ),
  
  retained_miRNAs = nrow(
    miRNA_ftd
  ),
  
  removed_miRNAs =
    nrow(
      RPM_counts
    ) -
    nrow(
      miRNA_ftd
    ),
  
  retained_fraction =
    nrow(
      miRNA_ftd
    ) /
    nrow(
      RPM_counts
    ),
  
  all_RPM_sums_valid = all(
    abs(
      colSums(
        RPM_counts
      ) -
        1e6
    ) <= 1e-4
  ),
  
  raw_RPM_sample_alignment = identical(
    colnames(
      raw_ftd_counts
    ),
    colnames(
      miRNA_ftd
    )
  ),
  
  raw_RPM_feature_alignment = identical(
    rownames(
      raw_ftd_counts
    ),
    rownames(
      miRNA_ftd
    )
  ),
  
  phenotype_labels_valid = identical(
    as.character(
      metadata$Phenotype
    ),
    unname(
      APEND_TYPE_SHORT_LABELS[
        as.character(
          metadata$Apend_type
        )
      ]
    )
  ),
  
  primary_membership_valid = if (
    "primary_analysis" %in%
    colnames(
      metadata
    )
  ) {
    all(
      metadata$primary_analysis
    )
  } else {
    NA
  },
  
  stringsAsFactors = FALSE
)


utils::write.csv(
  final_qc,
  prepare_output_file(
    file.path(
      results_dir,
      "01_final_QC_summary.csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. REPRODUCIBILITY
# ==============================================================================

write_session_info(
  results_dir
)


message(
  "\n============================================================"
)

message(
  "01 complete"
)

message(
  "Analysis: ",
  context$analysis_id
)

message(
  "Samples: ",
  ncol(
    counts
  )
)

message(
  "NIA / Mild / Severe: ",
  sum(
    metadata$Phenotype == "NIA"
  ),
  " / ",
  sum(
    metadata$Phenotype == "Mild"
  ),
  " / ",
  sum(
    metadata$Phenotype == "Severe"
  )
)

message(
  "Appendicitis / NIA: ",
  sum(
    metadata$Appendicitis_status == "Appendicitis"
  ),
  " / ",
  sum(
    metadata$Appendicitis_status == "NIA"
  )
)

message(
  "miRNAs retained: ",
  nrow(
    miRNA_ftd
  ),
  " / ",
  nrow(
    RPM_counts
  ),
  " (",
  round(
    100 *
      nrow(
        miRNA_ftd
      ) /
      nrow(
        RPM_counts
      ),
    2
  ),
  "%)"
)

message(
  "Filter: RPM >= ",
  FILTER_MIN_RPM,
  " in >= ",
  100 *
    FILTER_PREVALENCE,
  "% of at least one ",
  FILTER_GROUP_VARIABLE,
  " group."
)

message(
  "Clinical labels frozen as NIA / Mild / Severe."
)

message(
  "NIA = non-inflamed appendix; not a healthy-control cohort."
)

message(
  "Cohort-membership variables rebuilt from ANALYSIS_CATALOG."
)

message(
  "============================================================\n"
)
