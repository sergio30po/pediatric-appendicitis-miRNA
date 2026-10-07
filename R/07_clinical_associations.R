# ==============================================================================
# 07_clinical_associations.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Purpose:
#   Evaluate clinical associations of the integrated candidate miRNAs generated
#   by script 06, including inflammatory markers, perioperative outcomes,
#   detailed complications, hemolysis sensitivity and exploratory ROC analyses.
#
# Candidate hierarchy:
#   Defined upstream by script 06 and loaded from:
#     integrated/downstream_candidate_priority.csv
#
#   Script 07 does not redefine or re-rank the integrated candidate hierarchy.
#   It performs downstream clinical analyses on the carried-forward set only.
#
# Clinical terminology:
#   NIA = non-inflamed appendix.
#   NIA is not a healthy-control cohort.
#
# Complications:
#   The free-text clinical field `Complicaciones` is imported from
#   data/metadata_clinical_deidentified.xlsx. It is NOT required to be present
#   in metadata_master or metadata_analysis_samples.rds.
#
#   The raw text is preserved and reproducibly classified into:
#
#     Postoperative complications:
#       - intra-abdominal abscess
#       - bowel obstruction
#       - readmission for infection
#       - wound dehiscence
#
#     Intraoperative event:
#       - conversion/reconversion to open surgery
#
#     Broad perioperative adverse event:
#       postoperative complication OR conversion/reconversion
#
#   Entries describing alternative diagnoses/pathology or baseline disease
#   characteristics are retained for audit but are NOT automatically classified
#   as postoperative complications. This includes:
#       - gastroenteritis / Sapovirus / Giardia
#       - Yersinia-associated pathology
#       - enterobiasis
#       - negative appendectomy / "blanca"
#       - Meckel diverticulitis
#       - plastron/abscess already describing baseline severe disease
#
# IMPORTANT:
#   Complications and hospital stay are downstream clinical outcomes. They are
#   tested descriptively against candidate miRNA expression but are NEVER added
#   as covariates to the severity differential-expression models.
#
#   Blank `Complicaciones` cells are interpreted as no recorded complication/
#   event in the available deidentified clinical dataset.
#
# No candidate is promoted on the basis of this script.
# ==============================================================================

source(here::here("R", "00_analysis_config.R"))
source(here::here("R", "00_utils.R"))

context <- get_analysis_context()
results_dir <- context$results_dir
figures_dir <- context$figures_dir

out_dir <- file.path(
  results_dir,
  "clinical_associations_ROC"
)

fig_dir <- file.path(
  figures_dir,
  "clinical_associations_ROC"
)

ensure_dir(out_dir)
ensure_dir(fig_dir)

writeLines(
  c(
    "NIA = non-inflamed appendix; it is not a healthy-control cohort.",
    "Clinical follow-up uses the integrated candidate hierarchy from script 06.",
    "Primary_robust = robust FDR-supported integrated candidates.",
    "Secondary_integrated = RPM + edgeR + DESeq2 nominally concordant candidates.",
    "Inflammatory covariates are evaluated one at a time to limit overfitting.",
    "Appendicitis-only Mild-vs-Severe models are included as sensitivity analyses.",
    "The free-text Complicaciones field is imported from metadata_clinical_deidentified.xlsx.",
    "Postoperative complications and intraoperative conversion are classified separately.",
    "Complications are downstream outcomes and are never used as severity adjustment covariates.",
    "No sample is excluded on the basis of the sequencing hemolysis proxy."
  ),
  file.path(out_dir, "analysis_notes.txt")
)

required_packages <- c(
  "pROC",
  "edgeR",
  "DESeq2",
  "dplyr",
  "ggplot2",
  "readxl",
  "openxlsx"
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
    "Missing package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}

rpm_file <- file.path(results_dir, "filtered_RPM.rds")
counts_file <- file.path(results_dir, "raw_ftd_counts.rds")
meta_file <- file.path(results_dir, "metadata_analysis_samples.rds")
clinical_file <- file.path(project_root, "data", "metadata_clinical_deidentified.xlsx")
hemolysis_file <- file.path(results_dir, "hemolysis_proxy.csv")

for (f in c(rpm_file, counts_file, meta_file, clinical_file, hemolysis_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

rpm <- as.matrix(readRDS(rpm_file))
counts <- as.matrix(readRDS(counts_file))
storage.mode(counts) <- "numeric"
if (any(abs(counts - round(counts)) > 1e-8)) {
  stop("raw_ftd_counts.rds contains non-integer values.")
}
counts <- round(counts)
storage.mode(counts) <- "integer"

analysis_meta <- as.data.frame(readRDS(meta_file))
clinical <- as.data.frame(readxl::read_excel(
  clinical_file,
  na = c("", "NA", "N/A", "#NULL!", "NULL")
))

analysis_meta <- analysis_meta[match(colnames(rpm), analysis_meta$Sample), , drop = FALSE]
clinical <- clinical[match(colnames(rpm), clinical$Sample), , drop = FALSE]
counts <- counts[, match(colnames(rpm), colnames(counts)), drop = FALSE]

if (anyNA(clinical$Sample)) stop("Clinical metadata missing for one or more analysis samples.")
if (!identical(colnames(counts), colnames(rpm))) stop("Counts and RPM sample order differs.")

# Avoid duplicate columns from clinical metadata.
clinical <- clinical[, setdiff(colnames(clinical), intersect(setdiff(colnames(clinical), "Sample"), colnames(analysis_meta))), drop = FALSE]
metadata <- cbind(
  analysis_meta,
  clinical[, setdiff(colnames(clinical), "Sample"), drop = FALSE]
)
rownames(metadata) <- metadata$Sample

candidate_priority_file <- file.path(
  results_dir,
  "integrated",
  "downstream_candidate_priority.csv"
)


candidate_hierarchy_file <- file.path(
  results_dir,
  "integrated",
  "candidate_hierarchy_summary.csv"
)


for (f in c(
  candidate_priority_file,
  candidate_hierarchy_file
)) {
  if (!file.exists(f)) {
    stop(
      "Missing script-06 candidate hierarchy output: ",
      f,
      call. = FALSE
    )
  }
}


integrated_priority <- utils::read.csv(
  candidate_priority_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


candidate_hierarchy_summary <- utils::read.csv(
  candidate_hierarchy_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_priority_columns <- c(
  "miRNA",
  "evidence_class",
  "downstream_category",
  "best_evidence_class",
  "downstream_rank"
)


if (!all(
  required_priority_columns %in%
    colnames(
      integrated_priority
    )
)) {
  stop(
    "Script-06 downstream candidate hierarchy lacks required columns: ",
    paste(
      setdiff(
        required_priority_columns,
        colnames(
          integrated_priority
        )
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}


if (
  anyDuplicated(
    integrated_priority$miRNA
  ) > 0L
) {
  stop(
    "Duplicated miRNAs detected in the script-06 downstream hierarchy.",
    call. = FALSE
  )
}


if (
  !all(
    integrated_priority$downstream_category %in%
      c(
        "Primary_robust",
        "Secondary_integrated"
      )
  )
) {
  stop(
    "Unexpected downstream category in script-06 candidate hierarchy.",
    call. = FALSE
  )
}


candidate_miRNAs <- intersect(
  integrated_priority$miRNA,
  rownames(rpm)
)


if (
  length(
    candidate_miRNAs
  ) != nrow(
    integrated_priority
  )
) {
  stop(
    "One or more script-06 carried-forward miRNAs are absent from filtered_RPM.rds.",
    call. = FALSE
  )
}


candidate_panel <- integrated_priority |>
  dplyr::filter(
    miRNA %in%
      candidate_miRNAs
  ) |>
  dplyr::arrange(
    downstream_rank
  )


utils::write.csv(
  candidate_panel,
  file.path(
    out_dir,
    "candidate_panel.csv"
  ),
  row.names = FALSE
)


primary_miRNAs <- candidate_panel |>
  dplyr::filter(
    downstream_category ==
      "Primary_robust"
  ) |>
  dplyr::pull(
    miRNA
  ) |>
  unique()


secondary_miRNAs <- candidate_panel |>
  dplyr::filter(
    downstream_category ==
      "Secondary_integrated"
  ) |>
  dplyr::pull(
    miRNA
  ) |>
  unique()


targeted_miRNAs <- candidate_panel$miRNA


utils::write.csv(
  candidate_hierarchy_summary,
  file.path(
    out_dir,
    "candidate_hierarchy_summary.csv"
  ),
  row.names = FALSE
)

message(
  "Clinical candidate hierarchy from script 06: ",
  length(primary_miRNAs),
  " Primary_robust + ",
  length(secondary_miRNAs),
  " Secondary_integrated = ",
  length(targeted_miRNAs),
  " miRNAs."
)


# ------------------------------------------------------------------------------
# 1. Clinical variable definitions, detailed complication parsing and cleaning
# ------------------------------------------------------------------------------

to_numeric_clean <- function(x) {

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- trimws(
    as.character(x)
  )

  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "#NULL!",
      "NULL"
    )
  ] <- NA_character_

  x <- gsub(
    ",",
    ".",
    x,
    fixed = TRUE
  )

  suppressWarnings(
    as.numeric(x)
  )
}


recode_binary <- function(
    x,
    variable_name
) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "#NULL!",
      "NULL"
    )
  ] <- NA_character_


  if (
    variable_name == "Sexo"
  ) {

    # Clinical file: H = male, M = female.
    return(
      dplyr::case_when(
        x == "H" ~ 1,
        x == "M" ~ 0,
        TRUE ~ NA_real_
      )
    )
  }


  dplyr::case_when(
    x %in% c(
      "1",
      "SI",
      "SÍ",
      "YES",
      "Y",
      "POSITIVO",
      "PRESENTE"
    ) ~ 1,

    x %in% c(
      "0",
      "NO",
      "N",
      "NEGATIVO",
      "AUSENTE",
      "NO ECO"
    ) ~ 0,

    TRUE ~ NA_real_
  )
}


normalize_clinical_text <- function(x) {

  x_original <- as.character(x)

  x_original[
    is.na(
      x_original
    )
  ] <- ""

  x_original <- gsub(
    "<[^>]+>",
    " ",
    x_original
  )

  x_original <- gsub(
    "[\r\n\t]+",
    " ",
    x_original
  )

  x_ascii <- suppressWarnings(
    iconv(
      x_original,
      from = "",
      to = "ASCII//TRANSLIT"
    )
  )

  bad <- is.na(
    x_ascii
  )

  x_ascii[
    bad
  ] <- x_original[
    bad
  ]

  x_ascii <- toupper(
    x_ascii
  )

  x_ascii <- gsub(
    "\\s+",
    " ",
    x_ascii
  )

  trimws(
    x_ascii
  )
}


# ---------------------------------------------------------------------------
# 1.1 Detailed free-text complication field
# ---------------------------------------------------------------------------

if (
  "Complicaciones" %in%
    colnames(
      metadata
    )
) {

  metadata$Complicaciones_raw <- as.character(
    metadata$Complicaciones
  )

} else {

  metadata$Complicaciones_raw <- NA_character_

  warning(
    "Column `Complicaciones` was not found in metadata_clinical_deidentified.xlsx. ",
    "Detailed complication classification will be empty."
  )
}


metadata$Complicaciones_normalized <- normalize_clinical_text(
  metadata$Complicaciones_raw
)


metadata$Complicaciones_note_present <- as.integer(
  nzchar(
    metadata$Complicaciones_normalized
  )
)


# Postoperative complications.
#
# IMPORTANT:
# A PLASTRON entry containing ABSCESO is treated as baseline complex disease,
# not automatically as a postoperative intra-abdominal abscess.

metadata$Complication_intraabdominal_abscess <- as.integer(
  grepl(
    "ABSCESO",
    metadata$Complicaciones_normalized
  ) &
    !grepl(
      "PLASTRON",
      metadata$Complicaciones_normalized
    ) &
    !grepl(
      "SIN\\s+ABSCESO",
      metadata$Complicaciones_normalized
    )
)


metadata$Complication_bowel_obstruction <- as.integer(
  grepl(
    "OBSTRUCCION",
    metadata$Complicaciones_normalized
  )
)


metadata$Complication_readmission_infection <- as.integer(
  grepl(
    "REINGRESO",
    metadata$Complicaciones_normalized
  ) &
    grepl(
      "INFECCION",
      metadata$Complicaciones_normalized
    )
)


metadata$Complication_wound_dehiscence <- as.integer(
  grepl(
    "DEHISCENCIA",
    metadata$Complicaciones_normalized
  )
)


metadata$Postoperative_complication_any <- as.integer(
  metadata$Complication_intraabdominal_abscess == 1L |
    metadata$Complication_bowel_obstruction == 1L |
    metadata$Complication_readmission_infection == 1L |
    metadata$Complication_wound_dehiscence == 1L
)


# Intraoperative conversion is recorded separately and is not labelled as a
# postoperative complication.

metadata$Intraoperative_conversion <- as.integer(
  grepl(
    "RECONVERSION",
    metadata$Complicaciones_normalized
  ) |
    grepl(
      "CONVERSION.*ABIERT",
      metadata$Complicaciones_normalized
    )
)


metadata$Perioperative_adverse_event_any <- as.integer(
  metadata$Postoperative_complication_any == 1L |
    metadata$Intraoperative_conversion == 1L
)


# Audit-only categories: retained to prevent misclassification of free text.

metadata$Baseline_plastron_abscess_note <- as.integer(
  grepl(
    "PLASTRON",
    metadata$Complicaciones_normalized
  ) &
    grepl(
      "ABSCESO",
      metadata$Complicaciones_normalized
    )
)


metadata$Alternative_infectious_or_pathology_finding <- as.integer(
  grepl(
    paste(
      c(
        "GASTROENTERITIS",
        "SAPOVIRUS",
        "GIARDIA",
        "YERSIN",
        "ENTEROBIAS"
      ),
      collapse = "|"
    ),
    metadata$Complicaciones_normalized
  )
)


metadata$Negative_appendectomy_note <- as.integer(
  grepl(
    "BLANCA",
    metadata$Complicaciones_normalized
  )
)


metadata$Meckel_diverticulitis_note <- as.integer(
  grepl(
    "MECKEL",
    metadata$Complicaciones_normalized
  )
)


metadata$Complicaciones_note_classified <- as.integer(
  metadata$Postoperative_complication_any == 1L |
    metadata$Intraoperative_conversion == 1L |
    metadata$Baseline_plastron_abscess_note == 1L |
    metadata$Alternative_infectious_or_pathology_finding == 1L |
    metadata$Negative_appendectomy_note == 1L |
    metadata$Meckel_diverticulitis_note == 1L
)


metadata$Complicaciones_unclassified_note <- as.integer(
  metadata$Complicaciones_note_present == 1L &
    metadata$Complicaciones_note_classified == 0L
)


# Backward-compatible binary field, now explicitly defined.
metadata$Complicaciones_binary <-
  metadata$Postoperative_complication_any


complication_variables <- c(
  "Postoperative_complication_any",
  "Perioperative_adverse_event_any",
  "Complication_intraabdominal_abscess",
  "Complication_bowel_obstruction",
  "Complication_readmission_infection",
  "Complication_wound_dehiscence",
  "Intraoperative_conversion"
)


complication_labels <- c(
  Postoperative_complication_any =
    "Any postoperative complication",

  Perioperative_adverse_event_any =
    "Any perioperative adverse event",

  Complication_intraabdominal_abscess =
    "Intra-abdominal abscess",

  Complication_bowel_obstruction =
    "Bowel obstruction",

  Complication_readmission_infection =
    "Readmission for infection",

  Complication_wound_dehiscence =
    "Wound dehiscence",

  Intraoperative_conversion =
    "Conversion/reconversion to open surgery"
)


complication_text_audit <- metadata |>
  dplyr::transmute(
    Sample = Sample,
    Phenotype = as.character(
      Phenotype
    ),
    Appendicitis_status = as.character(
      Appendicitis_status
    ),
    Complicaciones_raw =
      Complicaciones_raw,
    Complicaciones_normalized =
      Complicaciones_normalized,
    note_present =
      Complicaciones_note_present,
    Postoperative_complication_any =
      Postoperative_complication_any,
    Perioperative_adverse_event_any =
      Perioperative_adverse_event_any,
    intraabdominal_abscess =
      Complication_intraabdominal_abscess,
    bowel_obstruction =
      Complication_bowel_obstruction,
    readmission_infection =
      Complication_readmission_infection,
    wound_dehiscence =
      Complication_wound_dehiscence,
    intraoperative_conversion =
      Intraoperative_conversion,
    baseline_plastron_abscess_note =
      Baseline_plastron_abscess_note,
    alternative_infectious_or_pathology_finding =
      Alternative_infectious_or_pathology_finding,
    negative_appendectomy_note =
      Negative_appendectomy_note,
    Meckel_diverticulitis_note =
      Meckel_diverticulitis_note,
    note_classified =
      Complicaciones_note_classified,
    unclassified_note =
      Complicaciones_unclassified_note
  )


utils::write.csv(
  complication_text_audit,
  file.path(
    out_dir,
    "complication_text_classification_audit.csv"
  ),
  row.names = FALSE
)


message(
  "Complication free-text notes detected: ",
  sum(
    metadata$Complicaciones_note_present,
    na.rm = TRUE
  )
)

message(
  "Derived postoperative complications: ",
  sum(
    metadata$Postoperative_complication_any,
    na.rm = TRUE
  )
)

message(
  "Derived intraoperative conversions: ",
  sum(
    metadata$Intraoperative_conversion,
    na.rm = TRUE
  )
)

message(
  "Unclassified non-empty complication notes: ",
  sum(
    metadata$Complicaciones_unclassified_note,
    na.rm = TRUE
  )
)

if (
  any(
    metadata$Complicaciones_unclassified_note == 1L,
    na.rm = TRUE
  )
) {
  warning(
    "One or more non-empty `Complicaciones` entries were not recognized by the ",
    "prespecified text classifier. Review complication_text_classification_audit.csv."
  )
}


# ---------------------------------------------------------------------------
# 1.2 Standard clinical variables
# ---------------------------------------------------------------------------

continuous_variables <- intersect(
  c(
    "Edad",
    "TiempoEvol",
    "Leucocitos",
    "Plaquetas",
    "Neutrófilos",
    "PCR",
    "PCT",
    "GrosorAp",
    "TiempoIQ",
    "TiempoIQAB",
    "TiempoHos"
  ),
  colnames(
    metadata
  )
)


binary_variables <- intersect(
  c(
    "Sexo",
    "Fiebre",
    "Vomitos",
    "Diarrea",
    "ECOComp",
    "LiqlibreEco",
    "ApDes",
    complication_variables
  ),
  colnames(
    metadata
  )
)


continuous_labels <- c(
  Edad = "Age",
  TiempoEvol = "Symptom duration",
  Leucocitos = "Leukocyte count",
  Plaquetas = "Platelet count",
  Neutrófilos = "Neutrophil count",
  PCR = "C-reactive protein",
  PCT = "Procalcitonin",
  GrosorAp = "Appendiceal diameter",
  TiempoIQ = "Surgery duration",
  TiempoIQAB = "Time from surgery to antibiotics",
  TiempoHos = "Hospital stay"
)


binary_labels <- c(
  Sexo = "Sex",
  Fiebre = "Fever",
  Vomitos = "Vomiting",
  Diarrea = "Diarrhea",
  ECOComp = "Complex ultrasound",
  LiqlibreEco = "Free fluid on ultrasound",
  ApDes = "Appendiceal destruction",
  complication_labels
)


for (
  v in continuous_variables
) {

  metadata[[
    v
  ]] <- to_numeric_clean(
    metadata[[
      v
    ]]
  )
}


nonderived_binary_variables <- setdiff(
  binary_variables,
  complication_variables
)


for (
  v in nonderived_binary_variables
) {

  metadata[[
    v
  ]] <- recode_binary(
    metadata[[
      v
    ]],
    v
  )
}


# ---------------------------------------------------------------------------
# 1.3 Complication frequencies and relationship with clinical phenotype
# ---------------------------------------------------------------------------

complication_frequency_rows <- list()


for (
  stratum_name in c(
    "All",
    "NIA",
    "Appendicitis",
    "Mild",
    "Severe"
  )
) {

  idx <- switch(
    stratum_name,
    All =
      rep(
        TRUE,
        nrow(
          metadata
        )
      ),
    NIA =
      metadata$Phenotype ==
        "NIA",
    Appendicitis =
      metadata$Appendicitis_status ==
        "Appendicitis",
    Mild =
      metadata$Phenotype ==
        "Mild",
    Severe =
      metadata$Phenotype ==
        "Severe"
  )


  for (
    v in complication_variables
  ) {

    values <- metadata[[
      v
    ]][
      idx
    ]

    complication_frequency_rows[[
      length(
        complication_frequency_rows
      ) +
        1L
    ]] <- data.frame(
      stratum = stratum_name,
      variable = v,
      variable_label = unname(
        complication_labels[
          v
        ]
      ),
      n_total = length(
        values
      ),
      n_event = sum(
        values == 1L,
        na.rm = TRUE
      ),
      event_pct = if (
        length(
          values
        ) > 0L
      ) {
        100 *
          mean(
            values == 1L,
            na.rm = TRUE
          )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
}


complication_frequency_summary <- dplyr::bind_rows(
  complication_frequency_rows
)


utils::write.csv(
  complication_frequency_summary,
  file.path(
    out_dir,
    "complication_frequency_summary.csv"
  ),
  row.names = FALSE
)


complication_phenotype_test_rows <- list()


for (
  v in complication_variables
) {

  # Global NIA / Mild / Severe Fisher test.
  tab_global <- table(
    metadata$Phenotype,
    factor(
      metadata[[
        v
      ]],
      levels = c(
        0,
        1
      )
    )
  )


  if (
    sum(
      tab_global[
        ,
        "1",
        drop = TRUE
      ]
    ) > 0L
  ) {

    ft <- stats::fisher.test(
      tab_global
    )

    complication_phenotype_test_rows[[
      length(
        complication_phenotype_test_rows
      ) +
        1L
    ]] <- data.frame(
      variable = v,
      variable_label = unname(
        complication_labels[
          v
        ]
      ),
      comparison = "Phenotype_global",
      test = "Fisher_exact",
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }


  # Appendicitis-only Mild vs Severe.
  idx_patient <- metadata$Appendicitis_status ==
    "Appendicitis"

  tab_patient <- table(
    droplevels(
      metadata$Phenotype[
        idx_patient
      ]
    ),
    factor(
      metadata[[
        v
      ]][
        idx_patient
      ],
      levels = c(
        0,
        1
      )
    )
  )


  if (
    nrow(
      tab_patient
    ) == 2L &&
    sum(
      tab_patient[
        ,
        "1",
        drop = TRUE
      ]
    ) > 0L
  ) {

    ft_patient <- stats::fisher.test(
      tab_patient
    )

    complication_phenotype_test_rows[[
      length(
        complication_phenotype_test_rows
      ) +
        1L
    ]] <- data.frame(
      variable = v,
      variable_label = unname(
        complication_labels[
          v
        ]
      ),
      comparison = "Appendicitis_Mild_vs_Severe",
      test = "Fisher_exact",
      p_value = ft_patient$p.value,
      stringsAsFactors = FALSE
    )
  }
}


complication_phenotype_tests <- dplyr::bind_rows(
  complication_phenotype_test_rows
)


if (
  nrow(
    complication_phenotype_tests
  ) > 0L
) {

  complication_phenotype_tests$FDR <- ave(
    complication_phenotype_tests$p_value,
    complication_phenotype_tests$comparison,
    FUN = function(p) {
      stats::p.adjust(
        p,
        method = "BH"
      )
    }
  )
}


utils::write.csv(
  complication_phenotype_tests,
  file.path(
    out_dir,
    "complication_by_phenotype_tests.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 2. Missingness audit
# ------------------------------------------------------------------------------

analysis_variables <- c(continuous_variables, binary_variables)
cohort_index <- list(
  All = rep(TRUE, nrow(metadata)),
  NIA = metadata$Phenotype == "NIA",
  Appendicitis = metadata$Appendicitis_status == "Appendicitis",
  Mild = metadata$Phenotype == "Mild",
  Severe = metadata$Phenotype == "Severe"
)

missing_rows <- list()
for (cohort_name in names(cohort_index)) {
  idx <- cohort_index[[cohort_name]]
  for (v in analysis_variables) {
    missing_rows[[length(missing_rows) + 1]] <- data.frame(
      cohort = cohort_name,
      variable = v,
      n_total = sum(idx),
      n_available = sum(!is.na(metadata[[v]][idx])),
      n_missing = sum(is.na(metadata[[v]][idx])),
      missing_pct = round(100 * mean(is.na(metadata[[v]][idx])), 2),
      stringsAsFactors = FALSE
    )
  }
}
missingness <- dplyr::bind_rows(missing_rows)
utils::write.csv(missingness, file.path(out_dir, "missingness_by_variable_and_group.csv"), row.names = FALSE)

sample_missing <- data.frame(
  Sample = metadata$Sample,
  n_missing = rowSums(is.na(metadata[, analysis_variables, drop = FALSE])),
  missing_pct = round(100 * rowMeans(is.na(metadata[, analysis_variables, drop = FALSE])), 2)
)
utils::write.csv(sample_missing, file.path(out_dir, "missingness_by_sample.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Clinical marker distributions by appendicitis severity
# ------------------------------------------------------------------------------

metadata$Severity_label <- factor(
  as.character(metadata$Phenotype),
  levels = c(
    "NIA",
    "Mild",
    "Severe"
  )
)

marker_summary_rows <- list()
for (v in continuous_variables) {
  for (g in levels(metadata$Severity_label)) {
    x <- metadata[[v]][metadata$Severity_label == g]
    x <- x[is.finite(x)]
    marker_summary_rows[[length(marker_summary_rows) + 1]] <- data.frame(
      variable = v,
      variable_label = unname(continuous_labels[v]),
      group = g,
      n_available = length(x),
      median = if (length(x)) stats::median(x) else NA_real_,
      Q1 = if (length(x)) unname(stats::quantile(x, 0.25)) else NA_real_,
      Q3 = if (length(x)) unname(stats::quantile(x, 0.75)) else NA_real_,
      mean = if (length(x)) mean(x) else NA_real_,
      sd = if (length(x) > 1) stats::sd(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
clinical_marker_summary <- dplyr::bind_rows(marker_summary_rows)
utils::write.csv(
  clinical_marker_summary,
  file.path(out_dir, "clinical_marker_group_summary.csv"),
  row.names = FALSE
)

clinical_marker_test_rows <- list()
for (v in continuous_variables) {
  ok <- !is.na(metadata$Severity_label) & is.finite(metadata[[v]])
  if (sum(ok) >= 10 && length(unique(metadata$Severity_label[ok])) >= 2) {
    kt <- stats::kruskal.test(metadata[[v]][ok] ~ metadata$Severity_label[ok])
    clinical_marker_test_rows[[length(clinical_marker_test_rows) + 1]] <- data.frame(
      variable = v,
      variable_label = unname(continuous_labels[v]),
      comparison = "Phenotype_global",
      n = sum(ok),
      statistic = unname(kt$statistic),
      p_value = kt$p.value,
      stringsAsFactors = FALSE
    )
  }

  pair_defs <- list(
    NIA_vs_Mild = c("NIA", "Mild"),
    NIA_vs_Severe = c("NIA", "Severe"),
    Mild_vs_Severe = c("Mild", "Severe")
  )
  for (nm in names(pair_defs)) {
    lev <- pair_defs[[nm]]
    x0 <- metadata[[v]][metadata$Severity_label == lev[1]]
    x1 <- metadata[[v]][metadata$Severity_label == lev[2]]
    x0 <- x0[is.finite(x0)]
    x1 <- x1[is.finite(x1)]
    if (length(x0) >= 4 && length(x1) >= 4) {
      wt <- stats::wilcox.test(x1, x0, exact = FALSE, correct = FALSE)
      clinical_marker_test_rows[[length(clinical_marker_test_rows) + 1]] <- data.frame(
        variable = v,
        variable_label = unname(continuous_labels[v]),
        comparison = nm,
        n = length(x0) + length(x1),
        statistic = unname(wt$statistic),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
clinical_marker_tests <- dplyr::bind_rows(clinical_marker_test_rows)
if (nrow(clinical_marker_tests) > 0) {
  clinical_marker_tests$FDR <- ave(
    clinical_marker_tests$p_value,
    clinical_marker_tests$comparison,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
}
utils::write.csv(
  clinical_marker_tests,
  file.path(out_dir, "clinical_marker_group_tests.csv"),
  row.names = FALSE
)

clinical_marker_fig_dir <- file.path(
  fig_dir,
  "clinical_markers"
)

ensure_dir(clinical_marker_fig_dir)
for (v in continuous_variables) {
  dat <- metadata[, c("Severity_label", v), drop = FALSE]
  dat <- dat[!is.na(dat$Severity_label) & is.finite(dat[[v]]), , drop = FALSE]
  if (nrow(dat) < 6) next
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = Severity_label, y = .data[[v]])
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.8) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(
      title = unname(continuous_labels[v]),
      x = "Group",
      y = unname(continuous_labels[v])
    )
  base_name <- file.path(
    clinical_marker_fig_dir,
    paste0(sanitize_filename(v), "_by_severity")
  )

  save_ggplot_pdf(
    plot = p,
    filename = paste0(
      base_name,
      ".pdf"
    ),
    width = 6.5,
    height = 5
  )

  save_ggplot_tiff(
    plot = p,
    filename = paste0(
      base_name,
      ".tiff"
    ),
    width = 6.5,
    height = 5,
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# 4. Continuous associations (Spearman)
# ------------------------------------------------------------------------------

strata <- list(
  All = rep(TRUE, nrow(metadata)),
  Appendicitis = metadata$Appendicitis_status == "Appendicitis",
  Mild = metadata$Phenotype == "Mild",
  Severe = metadata$Phenotype == "Severe"
)

continuous_results <- list()
for (mir in candidate_miRNAs) {
  expr <- as.numeric(rpm[mir, ])
  for (v in continuous_variables) {
    for (stratum in names(strata)) {
      idx <- strata[[stratum]]
      ok <- idx & is.finite(expr) & !is.na(metadata[[v]])
      if (sum(ok) < 5) next
      tst <- suppressWarnings(stats::cor.test(expr[ok], metadata[[v]][ok], method = "spearman", exact = FALSE))
      continuous_results[[length(continuous_results) + 1]] <- data.frame(
        miRNA = mir,
        variable = v,
        variable_label = unname(continuous_labels[v]),
        stratum = stratum,
        n = sum(ok),
        spearman_rho = unname(tst$estimate),
        p_value = tst$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
continuous_results <- dplyr::bind_rows(continuous_results)
if (nrow(continuous_results) > 0) {
  continuous_results$FDR <- ave(
    continuous_results$p_value,
    continuous_results$stratum,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
  continuous_results$fdr_significant <- continuous_results$FDR < ALPHA_FDR
}
utils::write.csv(
  continuous_results,
  file.path(out_dir, "continuous_clinical_associations.csv"),
  row.names = FALSE
)

primary_continuous_results <- continuous_results[
  continuous_results$miRNA %in% primary_miRNAs,
  ,
  drop = FALSE
]
utils::write.csv(
  primary_continuous_results,
  file.path(out_dir, "primary_miRNA_continuous_associations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Binary associations (Mann-Whitney; no automatic candidate promotion)
# ------------------------------------------------------------------------------

cliffs_delta <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  diffs <- outer(x, y, "-")
  (sum(diffs > 0) - sum(diffs < 0)) / (length(x) * length(y))
}

binary_results <- list()
for (mir in candidate_miRNAs) {
  expr <- as.numeric(rpm[mir, ])
  for (v in binary_variables) {
    for (stratum in c("All", "Appendicitis")) {
      idx <- if (
        stratum == "All"
      ) {
        rep(TRUE, nrow(metadata))
      } else {
        metadata$Appendicitis_status == "Appendicitis"
      }
      ok <- idx & is.finite(expr) & !is.na(metadata[[v]])
      if (sum(ok) < 6 || length(unique(metadata[[v]][ok])) < 2) next

      g0 <- expr[ok & metadata[[v]] == 0]
      g1 <- expr[ok & metadata[[v]] == 1]
      if (length(g0) < 3 || length(g1) < 3) next

      tst <- stats::wilcox.test(g1, g0, exact = FALSE, correct = FALSE)
      binary_results[[length(binary_results) + 1]] <- data.frame(
        miRNA = mir,
        variable = v,
        variable_label = unname(binary_labels[v]),
        stratum = stratum,
        n_group_0 = length(g0),
        n_group_1 = length(g1),
        median_group_0 = stats::median(g0),
        median_group_1 = stats::median(g1),
        cliffs_delta_group1_vs_group0 = cliffs_delta(g1, g0),
        p_value = tst$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
binary_results <- dplyr::bind_rows(binary_results)
if (nrow(binary_results) > 0) {
  binary_results$FDR <- ave(
    binary_results$p_value,
    binary_results$stratum,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
  binary_results$fdr_significant <- binary_results$FDR < ALPHA_FDR
}
utils::write.csv(
  binary_results,
  file.path(out_dir, "binary_clinical_associations.csv"),
  row.names = FALSE
)


complication_miRNA_associations <- binary_results |>
  dplyr::filter(
    variable %in%
      complication_variables
  )


if (
  nrow(
    complication_miRNA_associations
  ) > 0L
) {

  complication_miRNA_associations <- complication_miRNA_associations |>
    dplyr::group_by(
      stratum,
      variable
    ) |>
    dplyr::mutate(
      complication_specific_FDR = stats::p.adjust(
        p_value,
        method = "BH"
      )
    ) |>
    dplyr::ungroup()
}


utils::write.csv(
  complication_miRNA_associations,
  file.path(
    out_dir,
    "candidate_miRNA_complication_associations.csv"
  ),
  row.names = FALSE
)

primary_binary_results <- binary_results[
  binary_results$miRNA %in% primary_miRNAs,
  ,
  drop = FALSE
]
utils::write.csv(
  primary_binary_results,
  file.path(out_dir, "primary_miRNA_binary_associations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Hemolysis control
# ------------------------------------------------------------------------------

hemolysis_raw <- utils::read.csv(
  hemolysis_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hemolysis_available <- all(
  c("Sample", "hemolysis_logRPM_difference") %in% colnames(hemolysis_raw)
)

hemolysis_group_summary <- data.frame()
hemolysis_group_tests <- data.frame()
hemolysis_miRNA_associations <- data.frame()

if (hemolysis_available) {
  hemo <- hemolysis_raw[
    match(metadata$Sample, hemolysis_raw$Sample),
    c("Sample", "hemolysis_logRPM_difference"),
    drop = FALSE
  ]
  if (anyNA(hemo$Sample)) stop("Hemolysis proxy missing for one or more analysis samples.")

  metadata$hemolysis_logRPM_difference <- as.numeric(hemo$hemolysis_logRPM_difference)

  # Descriptive distribution of the sequencing-based hemolysis proxy by severity.
  hemo_summary_rows <- list()
  for (g in levels(metadata$Severity_label)) {
    x <- metadata$hemolysis_logRPM_difference[metadata$Severity_label == g]
    x <- x[is.finite(x)]
    hemo_summary_rows[[length(hemo_summary_rows) + 1]] <- data.frame(
      group = g,
      n_available = length(x),
      median = if (length(x)) stats::median(x) else NA_real_,
      Q1 = if (length(x)) unname(stats::quantile(x, 0.25)) else NA_real_,
      Q3 = if (length(x)) unname(stats::quantile(x, 0.75)) else NA_real_,
      mean = if (length(x)) mean(x) else NA_real_,
      sd = if (length(x) > 1) stats::sd(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  hemolysis_group_summary <- dplyr::bind_rows(hemo_summary_rows)

  # Group-level tests: the proxy itself should not simply reproduce disease severity.
  hemo_tests <- list()
  ok <- is.finite(metadata$hemolysis_logRPM_difference) & !is.na(metadata$Severity_label)
  if (sum(ok) >= 10) {
    kt <- stats::kruskal.test(
      metadata$hemolysis_logRPM_difference[ok] ~ metadata$Severity_label[ok]
    )
    hemo_tests[[length(hemo_tests) + 1]] <- data.frame(
      comparison = "Phenotype_global",
      n = sum(ok),
      statistic = unname(kt$statistic),
      p_value = kt$p.value,
      stringsAsFactors = FALSE
    )
  }

  hemo_pair_defs <- list(
    Appendicitis_vs_NIA = list(
      group = as.character(metadata$Appendicitis_status),
      negative = "NIA",
      positive = "Appendicitis"
    ),
    NIA_vs_Mild = list(
      group = as.character(metadata$Phenotype),
      negative = "NIA",
      positive = "Mild"
    ),
    NIA_vs_Severe = list(
      group = as.character(metadata$Phenotype),
      negative = "NIA",
      positive = "Severe"
    ),
    Mild_vs_Severe = list(
      group = as.character(metadata$Phenotype),
      negative = "Mild",
      positive = "Severe"
    )
  )
  for (nm in names(hemo_pair_defs)) {
    cfg <- hemo_pair_defs[[nm]]
    x0 <- metadata$hemolysis_logRPM_difference[cfg$group == cfg$negative]
    x1 <- metadata$hemolysis_logRPM_difference[cfg$group == cfg$positive]
    x0 <- x0[is.finite(x0)]
    x1 <- x1[is.finite(x1)]
    if (length(x0) >= 4 && length(x1) >= 4) {
      wt <- stats::wilcox.test(x1, x0, exact = FALSE, correct = FALSE)
      hemo_tests[[length(hemo_tests) + 1]] <- data.frame(
        comparison = nm,
        n = length(x0) + length(x1),
        statistic = unname(wt$statistic),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  hemolysis_group_tests <- dplyr::bind_rows(hemo_tests)
  if (nrow(hemolysis_group_tests) > 0) {
    hemolysis_group_tests$FDR <- stats::p.adjust(hemolysis_group_tests$p_value, method = "BH")
  }

  # Candidate expression vs hemolysis proxy.
  hemo_assoc <- list()
  hemo_strata <- list(
    All = rep(TRUE, nrow(metadata)),
    Appendicitis = metadata$Appendicitis_status == "Appendicitis"
  )
  for (mir in candidate_miRNAs) {
    expr <- as.numeric(rpm[mir, ])
    for (stratum in names(hemo_strata)) {
      idx <- hemo_strata[[stratum]] &
        is.finite(expr) &
        is.finite(metadata$hemolysis_logRPM_difference)
      if (sum(idx) < 8) next
      tst <- suppressWarnings(stats::cor.test(
        expr[idx],
        metadata$hemolysis_logRPM_difference[idx],
        method = "spearman",
        exact = FALSE
      ))
      hemo_assoc[[length(hemo_assoc) + 1]] <- data.frame(
        miRNA = mir,
        downstream_category = candidate_panel$downstream_category[
          match(mir, candidate_panel$miRNA)
        ],
        stratum = stratum,
        n = sum(idx),
        spearman_rho = unname(tst$estimate),
        p_value = tst$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  hemolysis_miRNA_associations <- dplyr::bind_rows(hemo_assoc)
  if (nrow(hemolysis_miRNA_associations) > 0) {
    hemolysis_miRNA_associations$FDR <- ave(
      hemolysis_miRNA_associations$p_value,
      hemolysis_miRNA_associations$stratum,
      FUN = function(p) stats::p.adjust(p, method = "BH")
    )
  }
} else {
  warning("Hemolysis proxy is unavailable; targeted hemolysis checks are skipped.")
}

utils::write.csv(
  hemolysis_group_summary,
  file.path(out_dir, "hemolysis_proxy_group_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  hemolysis_group_tests,
  file.path(out_dir, "hemolysis_proxy_group_tests.csv"),
  row.names = FALSE
)
utils::write.csv(
  hemolysis_miRNA_associations,
  file.path(out_dir, "candidate_miRNA_vs_hemolysis.csv"),
  row.names = FALSE
)

if (nrow(hemolysis_miRNA_associations) > 0) {
  utils::write.csv(
    hemolysis_miRNA_associations[
      hemolysis_miRNA_associations$miRNA %in% primary_miRNAs,
      ,
      drop = FALSE
    ],
    file.path(out_dir, "primary_miRNA_vs_hemolysis.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 6b. Correlation structure among inflammatory/technical covariates
# ------------------------------------------------------------------------------

covariate_corr_variables <- intersect(
  c("Neutrófilos", "Leucocitos", "PCR", "TiempoEvol", "hemolysis_logRPM_difference"),
  colnames(metadata)
)

covariate_corr_rows <- list()
if (length(covariate_corr_variables) >= 2) {
  corr_strata <- list(
    All = rep(TRUE, nrow(metadata)),
    Appendicitis = metadata$Appendicitis_status == "Appendicitis"
  )
  pairs <- utils::combn(covariate_corr_variables, 2, simplify = FALSE)
  for (stratum in names(corr_strata)) {
    idx0 <- corr_strata[[stratum]]
    for (pair in pairs) {
      x <- suppressWarnings(as.numeric(metadata[[pair[1]]]))
      y <- suppressWarnings(as.numeric(metadata[[pair[2]]]))
      ok <- idx0 & is.finite(x) & is.finite(y)
      if (sum(ok) < 8) next
      tst <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      covariate_corr_rows[[length(covariate_corr_rows) + 1]] <- data.frame(
        stratum = stratum,
        variable_1 = pair[1],
        variable_2 = pair[2],
        n = sum(ok),
        spearman_rho = unname(tst$estimate),
        p_value = tst$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
covariate_correlations <- dplyr::bind_rows(covariate_corr_rows)
if (nrow(covariate_correlations) > 0) {
  covariate_correlations$FDR <- ave(
    covariate_correlations$p_value,
    covariate_correlations$stratum,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
}
utils::write.csv(
  covariate_correlations,
  file.path(out_dir, "inflammatory_technical_covariate_correlations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Targeted count-model adjustment for inflammation and hemolysis
# ------------------------------------------------------------------------------

# These are pre-specified secondary sensitivity analyses on the dynamically
# selected integrated candidates (Primary_robust + Secondary_integrated).
# Each marker is added ONE AT A TIME to
# age + sex + severity. For every marker, the reference and extended model use
# the exact same complete-case sample set, so attenuation is not confounded by
# a change in sample composition.
#
# Hospital stay, surgical duration, complications and other likely downstream
# consequences of severity are deliberately not included as adjustment
# covariates. They remain descriptive clinical associations above.

adjustment_variables <- intersect(
  c("Neutrófilos", "Leucocitos", "PCR", "PCT", "TiempoEvol"),
  colnames(metadata)
)
if (hemolysis_available) {
  adjustment_variables <- c(adjustment_variables, "hemolysis_logRPM_difference")
}
adjustment_variables <- unique(adjustment_variables)

adjustment_labels <- c(
  `Neutrófilos` = "Neutrophil count",
  Leucocitos = "Leukocyte count",
  PCR = "C-reactive protein",
  PCT = "Procalcitonin",
  TiempoEvol = "Symptom duration",
  hemolysis_logRPM_difference = "Sequencing hemolysis proxy"
)

targeted_model_rows <- list()
targeted_audit_rows <- list()

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

fdr_from_p <- function(p) stats::p.adjust(p, method = "BH")

edger_extract_severity <- function(count_mat, dat, formula, model_role, covariate_name) {
  dge <- edgeR::DGEList(counts = count_mat)
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  design <- stats::model.matrix(formula, data = dat)
  if (qr(design)$rank != ncol(design)) {
    stop("Non-estimable edgeR design.")
  }
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)

  sev_coef <- grep("^Phenotype", colnames(design))
  if (length(sev_coef) != 2) stop("Expected two severity coefficients in edgeR.")

  out <- list()

  qlf_g <- edgeR::glmQLFTest(fit, coef = sev_coef)
  tt_g <- edgeR::topTags(qlf_g, n = Inf, sort.by = "none")$table
  pg <- tt_g$PValue
  names(pg) <- rownames(tt_g)
  fdrg <- fdr_from_p(pg)
  out[[length(out) + 1]] <- data.frame(
    miRNA = names(pg),
    contrast = "Phenotype_global",
    numerator = NA_character_,
    denominator = NA_character_,
    log2FC = NA_real_,
    p_value = unname(pg),
    FDR = unname(fdrg),
    stringsAsFactors = FALSE
  )

  contrast_vector <- function(num_level, den_level) {
    v <- rep(0, ncol(design))
    names(v) <- colnames(design)

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
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vector(lv[1], lv[2]))
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    p <- tt$PValue
    fdr <- fdr_from_p(p)
    out[[length(out) + 1]] <- data.frame(
      miRNA = rownames(tt),
      contrast = nm,
      numerator = lv[1],
      denominator = lv[2],
      log2FC = tt$logFC,
      p_value = p,
      FDR = fdr,
      stringsAsFactors = FALSE
    )
  }

  ans <- dplyr::bind_rows(out)
  ans$method <- "edgeR_TMM"
  ans$model_role <- model_role
  ans$adjustment_variable <- covariate_name
  ans
}

deseq_extract_severity <- function(count_mat, dat, formula, reduced_formula, model_role, covariate_name) {
  coldata <- dat
  rownames(coldata) <- coldata$Sample

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  out <- list()
  pairwise <- list(
    NIA_vs_Mild = c("Mild", "NIA"),
    NIA_vs_Severe = c("Severe", "NIA"),
    Mild_vs_Severe = c("Severe", "Mild")
  )

  for (nm in names(pairwise)) {
    lv <- pairwise[[nm]]
    rr <- DESeq2::results(
      dds,
      contrast = c("Phenotype", lv[1], lv[2]),
      independentFiltering = FALSE
    )
    p <- rr$pvalue
    fdr <- fdr_from_p(p)
    out[[length(out) + 1]] <- data.frame(
      miRNA = rownames(rr),
      contrast = nm,
      numerator = lv[1],
      denominator = lv[2],
      log2FC = rr$log2FoldChange,
      p_value = p,
      FDR = fdr,
      stringsAsFactors = FALSE
    )
  }

  dds_lrt <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = formula
  )
  dds_lrt <- DESeq2::DESeq(
    dds_lrt,
    test = "LRT",
    reduced = reduced_formula,
    quiet = TRUE
  )
  rr_g <- DESeq2::results(dds_lrt, independentFiltering = FALSE)
  pg <- rr_g$pvalue
  fdrg <- fdr_from_p(pg)
  out[[length(out) + 1]] <- data.frame(
    miRNA = rownames(rr_g),
    contrast = "Phenotype_global",
    numerator = NA_character_,
    denominator = NA_character_,
    log2FC = NA_real_,
    p_value = pg,
    FDR = fdrg,
    stringsAsFactors = FALSE
  )

  ans <- dplyr::bind_rows(out)
  ans$method <- "DESeq2"
  ans$model_role <- model_role
  ans$adjustment_variable <- covariate_name
  ans
}

if (length(targeted_miRNAs) > 0 && length(adjustment_variables) > 0) {
  # Model uses Age from analysis metadata, not the duplicated clinical Edad field.
  metadata$Age_model <- suppressWarnings(as.numeric(metadata$Age))
  metadata$Sex_model <- droplevels(factor(as.character(metadata$Sex)))
  metadata$Phenotype_model <- factor(
    as.character(metadata$Phenotype),
    levels = c(
      "NIA",
      "Mild",
      "Severe"
    )
  )

  for (v in adjustment_variables) {
    cov_value <- suppressWarnings(as.numeric(metadata[[v]]))
    ok <- is.finite(metadata$Age_model) &
      !is.na(metadata$Sex_model) &
      !is.na(metadata$Phenotype_model) &
      is.finite(cov_value)

    dat <- metadata[ok, , drop = FALSE]
    count_sub <- counts[, ok, drop = FALSE]
    dat$Age_z_targeted <- safe_z(dat$Age_model)
    dat$Covariate_z_targeted <- safe_z(cov_value[ok])
    dat$Sex_model <- droplevels(factor(dat$Sex_model))
    dat$Phenotype_model <- factor(
      as.character(dat$Phenotype_model),
      levels = c(
        "NIA",
        "Mild",
        "Severe"
      )
    )

    group_counts <- table(dat$Phenotype_model)
    n_NIA <- if ("NIA" %in% names(group_counts)) unname(group_counts["NIA"]) else 0
    n_mild <- if ("Mild" %in% names(group_counts)) unname(group_counts["Mild"]) else 0
    n_severe <- if ("Severe" %in% names(group_counts)) unname(group_counts["Severe"]) else 0

    usable <- nrow(dat) >= 18 &&
      min(c(n_NIA, n_mild, n_severe)) >= 4 &&
      all(is.finite(dat$Age_z_targeted)) &&
      all(is.finite(dat$Covariate_z_targeted)) &&
      nlevels(dat$Sex_model) >= 2

    targeted_audit_rows[[length(targeted_audit_rows) + 1]] <- data.frame(
      adjustment_variable = v,
      adjustment_label = unname(adjustment_labels[v]),
      n_complete = nrow(dat),
      n_NIA = n_NIA,
      n_mild = n_mild,
      n_severe = n_severe,
      usable = usable,
      note = if (usable) {
        "Reference and extended models use identical complete-case samples."
      } else {
        "Skipped: insufficient complete cases/group size, zero variance, or non-estimable sex factor."
      },
      stringsAsFactors = FALSE
    )

    if (!usable) next

    # Make model-column names simple and explicit for both packages.
    model_dat <- data.frame(
      Sample = dat$Sample,
      Age_z = dat$Age_z_targeted,
      Sex = dat$Sex_model,
      Phenotype = factor(
        as.character(dat$Phenotype_model),
        levels = c(
          "NIA",
          "Mild",
          "Severe"
        )
      ),
      Covariate_z = dat$Covariate_z_targeted,
      stringsAsFactors = FALSE
    )

    ref_formula <- ~ Age_z + Sex + Phenotype
    ext_formula <- ~ Age_z + Sex + Covariate_z + Phenotype
    ref_reduced <- ~ Age_z + Sex
    ext_reduced <- ~ Age_z + Sex + Covariate_z

    # Explicit estimability check before fitting.
    ref_mm <- stats::model.matrix(ref_formula, data = model_dat)
    ext_mm <- stats::model.matrix(ext_formula, data = model_dat)
    if (qr(ref_mm)$rank != ncol(ref_mm) || qr(ext_mm)$rank != ncol(ext_mm)) {
      targeted_audit_rows[[length(targeted_audit_rows)]]$usable <- FALSE
      targeted_audit_rows[[length(targeted_audit_rows)]]$note <- "Skipped: model matrix is rank deficient."
      next
    }

    edger_ref <- edger_extract_severity(
      count_sub, model_dat, ref_formula, "reference_age_sex_same_samples", v
    )
    edger_ext <- edger_extract_severity(
      count_sub, model_dat, ext_formula, "extended_age_sex_plus_covariate", v
    )
    deseq_ref <- deseq_extract_severity(
      count_sub, model_dat, ref_formula, ref_reduced,
      "reference_age_sex_same_samples", v
    )
    deseq_ext <- deseq_extract_severity(
      count_sub, model_dat, ext_formula, ext_reduced,
      "extended_age_sex_plus_covariate", v
    )

    all_models <- dplyr::bind_rows(edger_ref, edger_ext, deseq_ref, deseq_ext)
    all_models <- all_models[all_models$miRNA %in% targeted_miRNAs, , drop = FALSE]
    all_models$downstream_category <- candidate_panel$downstream_category[
      match(all_models$miRNA, candidate_panel$miRNA)
    ]
    all_models$n_complete <- nrow(model_dat)
    all_models$n_NIA <- n_NIA
    all_models$n_mild <- n_mild
    all_models$n_severe <- n_severe
    targeted_model_rows[[length(targeted_model_rows) + 1]] <- all_models
  }
}

targeted_model_long <- dplyr::bind_rows(targeted_model_rows)
targeted_model_audit <- dplyr::bind_rows(targeted_audit_rows)

utils::write.csv(
  targeted_model_audit,
  file.path(out_dir, "integrated_candidate_covariate_model_audit.csv"),
  row.names = FALSE
)

targeted_model_comparison <- data.frame()
if (nrow(targeted_model_long) > 0) {
  ref <- targeted_model_long[
    targeted_model_long$model_role == "reference_age_sex_same_samples",
    ,
    drop = FALSE
  ]
  ext <- targeted_model_long[
    targeted_model_long$model_role == "extended_age_sex_plus_covariate",
    ,
    drop = FALSE
  ]

  key <- c(
    "miRNA", "contrast", "numerator", "denominator", "method",
    "adjustment_variable", "n_complete", "n_NIA", "n_mild", "n_severe"
  )

  ref <- ref[, c(key, "log2FC", "p_value", "FDR"), drop = FALSE]
  ext <- ext[, c(key, "log2FC", "p_value", "FDR"), drop = FALSE]

  targeted_model_comparison <- dplyr::left_join(
    ref,
    ext,
    by = key,
    suffix = c("_reference", "_adjusted")
  )

  targeted_model_comparison$delta_log2FC_adjusted_minus_reference <-
    targeted_model_comparison$log2FC_adjusted -
    targeted_model_comparison$log2FC_reference

  targeted_model_comparison$effect_attenuation_pct <- ifelse(
    is.finite(targeted_model_comparison$log2FC_reference) &
      abs(targeted_model_comparison$log2FC_reference) > 1e-12 &
      is.finite(targeted_model_comparison$log2FC_adjusted),
    100 * (
      1 -
        abs(targeted_model_comparison$log2FC_adjusted) /
        abs(targeted_model_comparison$log2FC_reference)
    ),
    NA_real_
  )

  targeted_model_comparison$direction_stable <- ifelse(
    is.finite(targeted_model_comparison$log2FC_reference) &
      is.finite(targeted_model_comparison$log2FC_adjusted),
    sign(targeted_model_comparison$log2FC_reference) ==
      sign(targeted_model_comparison$log2FC_adjusted),
    NA
  )

  targeted_model_comparison$retains_FDR_lt_0.05_after_adjustment <-
    !is.na(targeted_model_comparison$FDR_adjusted) &
    targeted_model_comparison$FDR_adjusted < ALPHA_FDR

  targeted_model_comparison$interpretation <- dplyr::case_when(
    targeted_model_comparison$contrast == "Phenotype_global" &
      targeted_model_comparison$retains_FDR_lt_0.05_after_adjustment ~
      "Global severity association remains FDR-significant after adding the covariate.",
    targeted_model_comparison$contrast == "Phenotype_global" ~
      "Global severity association does not remain FDR-significant after adding the covariate.",
    targeted_model_comparison$retains_FDR_lt_0.05_after_adjustment &
      is.finite(targeted_model_comparison$effect_attenuation_pct) &
      targeted_model_comparison$effect_attenuation_pct < 25 ~
      "Association remains FDR-significant with limited effect attenuation.",
    targeted_model_comparison$retains_FDR_lt_0.05_after_adjustment ~
      "Association remains FDR-significant but effect attenuation should be reviewed.",
    TRUE ~
      "Association does not remain FDR-significant after adding the covariate; review effect attenuation and sample size."
  )
}

utils::write.csv(
  targeted_model_long,
  file.path(out_dir, "integrated_candidates_targeted_models_long.csv"),
  row.names = FALSE
)
utils::write.csv(
  targeted_model_comparison,
  file.path(out_dir, "integrated_candidates_inflammation_hemolysis_adjustment.csv"),
  row.names = FALSE
)


primary_targeted_model_comparison <- targeted_model_comparison[
  targeted_model_comparison$miRNA %in% primary_miRNAs,
  ,
  drop = FALSE
]

utils::write.csv(
  primary_targeted_model_comparison,
  file.path(out_dir, "primary_robust_candidate_inflammation_hemolysis_adjustment.csv"),
  row.names = FALSE
)

# Targeted scatterplots for integrated candidates and the main inflammatory/technical covariates.
targeted_fig_dir <- file.path(
  fig_dir,
  "covar"
)

ensure_dir(targeted_fig_dir)

plot_covariates <- intersect(
  c("Neutrófilos", "Leucocitos", "PCR", "PCT", "TiempoEvol", "hemolysis_logRPM_difference"),
  colnames(metadata)
)

for (mir in targeted_miRNAs) {
  for (v in plot_covariates) {
    dat <- data.frame(
      expression = as.numeric(rpm[mir, ]),
      covariate = suppressWarnings(as.numeric(metadata[[v]])),
      group = metadata$Severity_label,
      stringsAsFactors = FALSE
    )
    dat <- dat[
      is.finite(dat$expression) &
        is.finite(dat$covariate) &
        !is.na(dat$group),
      ,
      drop = FALSE
    ]
    if (nrow(dat) < 8) next

    rho <- suppressWarnings(stats::cor.test(
      dat$expression,
      dat$covariate,
      method = "spearman",
      exact = FALSE
    ))

    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = covariate, y = log2(expression + 1), shape = group)
    ) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::geom_smooth(method = "lm", se = TRUE) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::labs(
        title = paste(mir, "vs", unname(adjustment_labels[v])),
        subtitle = sprintf(
          "Spearman rho = %.2f; p = %.3g. Exploratory association.",
          unname(rho$estimate), rho$p.value
        ),
        x = unname(adjustment_labels[v]),
        y = "log2(RPM + 1)",
        shape = "Group"
      )

    base_name <- file.path(
      targeted_fig_dir,
      paste0(
        sanitize_filename(mir),
        "_vs_",
        sanitize_filename(v)
      )
    )

    save_ggplot_pdf(
      plot = p,
      filename = paste0(
        base_name,
        ".pdf"
      ),
      width = 7,
      height = 5.5
    )

    save_ggplot_tiff(
      plot = p,
      filename = paste0(
        base_name,
        ".tiff"
      ),
      width = 7,
      height = 5.5,
      dpi = 600
    )
  }
}

# ------------------------------------------------------------------------------
# 8. Appendicitis-only Mild-vs-Severe sensitivity
# ------------------------------------------------------------------------------
#
# Rationale:
# Group 0 consists of symptomatic surgical NIA cases, not healthy NIA cases. CRP and
# other inflammatory variables can therefore be heterogeneous in NIA cases for
# reasons unrelated to appendiceal inflammation. To answer the clinically focused
# question "does the miRNA remain associated with severity among appendicitis
# patients?", the following secondary models use only Mild (1) and Severe (2)
# patients. Each covariate is added ONE AT A TIME to age + sex + severity.
# FDR is recomputed across all tested miRNAs before retaining the integrated
# candidate panel. Upstream candidate tier is preserved.

patient_adjustment_variables <- intersect(
  c("Neutrófilos", "Leucocitos", "PCR", "PCT", "TiempoEvol", "hemolysis_logRPM_difference"),
  colnames(metadata)
)
patient_adjustment_variables <- unique(patient_adjustment_variables)

patient_model_rows <- list()
patient_audit_rows <- list()

edger_extract_patient_severity <- function(count_mat, dat, formula, model_role, covariate_name) {
  dge <- edgeR::DGEList(counts = count_mat)
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  design <- stats::model.matrix(formula, data = dat)
  if (qr(design)$rank != ncol(design)) stop("Non-estimable appendicitis-only edgeR design.")
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)

  sev_coef <- grep("^Severity_patient", colnames(design))
  if (length(sev_coef) != 1) stop("Expected one Mild-vs-Severe coefficient in appendicitis-only edgeR model.")

  qlf <- edgeR::glmQLFTest(fit, coef = sev_coef)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  p <- tt$PValue
  fdr <- stats::p.adjust(p, method = "BH")

  data.frame(
    miRNA = rownames(tt),
    contrast = "Appendicitis_Mild_vs_Severe",
    numerator = "Severe",
    denominator = "Mild",
    log2FC = tt$logFC,
    p_value = p,
    FDR = fdr,
    method = "edgeR_TMM",
    model_role = model_role,
    adjustment_variable = covariate_name,
    stringsAsFactors = FALSE
  )
}

deseq_extract_patient_severity <- function(count_mat, dat, formula, model_role, covariate_name) {
  coldata <- dat
  rownames(coldata) <- coldata$Sample
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  rr <- DESeq2::results(
    dds,
    contrast = c("Severity_patient", "Severe", "Mild"),
    independentFiltering = FALSE
  )
  p <- rr$pvalue
  fdr <- stats::p.adjust(p, method = "BH")

  data.frame(
    miRNA = rownames(rr),
    contrast = "Appendicitis_Mild_vs_Severe",
    numerator = "Severe",
    denominator = "Mild",
    log2FC = rr$log2FoldChange,
    p_value = p,
    FDR = fdr,
    method = "DESeq2",
    model_role = model_role,
    adjustment_variable = covariate_name,
    stringsAsFactors = FALSE
  )
}

if (length(targeted_miRNAs) > 0 && length(patient_adjustment_variables) > 0) {
  patient_base_idx <- metadata$Appendicitis_status == "Appendicitis"

  for (v in patient_adjustment_variables) {
    cov_value <- suppressWarnings(as.numeric(metadata[[v]]))
    age_value <- suppressWarnings(as.numeric(metadata$Age))
    sex_value <- factor(as.character(metadata$Sex))

    ok <- patient_base_idx &
      is.finite(age_value) &
      !is.na(sex_value) &
      is.finite(cov_value)

    dat <- metadata[ok, , drop = FALSE]
    count_sub <- counts[, ok, drop = FALSE]

    dat$Age_z_patient <- safe_z(age_value[ok])
    dat$Covariate_z_patient <- safe_z(cov_value[ok])
    dat$Sex_patient <- droplevels(factor(as.character(metadata$Sex[ok])))
    dat$Severity_patient <- factor(
      as.character(
        metadata$Phenotype[
          ok
        ]
      ),
      levels = c(
        "Mild",
        "Severe"
      )
    )

    group_counts <- table(dat$Severity_patient)
    n_mild <- if ("Mild" %in% names(group_counts)) unname(group_counts["Mild"]) else 0
    n_severe <- if ("Severe" %in% names(group_counts)) unname(group_counts["Severe"]) else 0

    usable <- nrow(dat) >= 20 &&
      min(c(n_mild, n_severe)) >= 5 &&
      all(is.finite(dat$Age_z_patient)) &&
      all(is.finite(dat$Covariate_z_patient)) &&
      nlevels(dat$Sex_patient) >= 2

    patient_audit_rows[[length(patient_audit_rows) + 1]] <- data.frame(
      adjustment_variable = v,
      adjustment_label = unname(adjustment_labels[v]),
      n_complete = nrow(dat),
      n_mild = n_mild,
      n_severe = n_severe,
      usable = usable,
      note = if (usable) {
        "Appendicitis-only reference and extended models use identical complete-case samples."
      } else {
        "Skipped: insufficient complete cases/group size, zero variance, or non-estimable sex factor."
      },
      stringsAsFactors = FALSE
    )

    if (!usable) next

    model_dat <- data.frame(
      Sample = dat$Sample,
      Age_z = dat$Age_z_patient,
      Sex = dat$Sex_patient,
      Severity_patient = dat$Severity_patient,
      Covariate_z = dat$Covariate_z_patient,
      stringsAsFactors = FALSE
    )

    ref_formula_patient <- ~ Age_z + Sex + Severity_patient
    ext_formula_patient <- ~ Age_z + Sex + Covariate_z + Severity_patient

    ref_mm <- stats::model.matrix(ref_formula_patient, data = model_dat)
    ext_mm <- stats::model.matrix(ext_formula_patient, data = model_dat)
    if (qr(ref_mm)$rank != ncol(ref_mm) || qr(ext_mm)$rank != ncol(ext_mm)) {
      patient_audit_rows[[length(patient_audit_rows)]]$usable <- FALSE
      patient_audit_rows[[length(patient_audit_rows)]]$note <- "Skipped: appendicitis-only model matrix is rank deficient."
      next
    }

    edger_ref <- edger_extract_patient_severity(
      count_sub, model_dat, ref_formula_patient,
      "reference_age_sex_same_appendicitis", v
    )
    edger_ext <- edger_extract_patient_severity(
      count_sub, model_dat, ext_formula_patient,
      "extended_age_sex_plus_covariate_same_appendicitis", v
    )
    deseq_ref <- deseq_extract_patient_severity(
      count_sub, model_dat, ref_formula_patient,
      "reference_age_sex_same_appendicitis", v
    )
    deseq_ext <- deseq_extract_patient_severity(
      count_sub, model_dat, ext_formula_patient,
      "extended_age_sex_plus_covariate_same_appendicitis", v
    )

    all_patient_models <- dplyr::bind_rows(edger_ref, edger_ext, deseq_ref, deseq_ext)
    all_patient_models <- all_patient_models[
      all_patient_models$miRNA %in% targeted_miRNAs,
      ,
      drop = FALSE
    ]
    all_patient_models$downstream_category <- candidate_panel$downstream_category[
      match(all_patient_models$miRNA, candidate_panel$miRNA)
    ]
    all_patient_models$n_complete <- nrow(model_dat)
    all_patient_models$n_mild <- n_mild
    all_patient_models$n_severe <- n_severe
    patient_model_rows[[length(patient_model_rows) + 1]] <- all_patient_models
  }
}

patient_models_long <- dplyr::bind_rows(patient_model_rows)
patient_model_audit <- dplyr::bind_rows(patient_audit_rows)

utils::write.csv(
  patient_model_audit,
  file.path(out_dir, "appendicitis_only_covariate_model_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  patient_models_long,
  file.path(out_dir, "integrated_candidates_appendicitis_only_models_long.csv"),
  row.names = FALSE
)

patient_model_comparison <- data.frame()
if (nrow(patient_models_long) > 0) {
  ref <- patient_models_long[
    patient_models_long$model_role == "reference_age_sex_same_appendicitis",
    ,
    drop = FALSE
  ]
  ext <- patient_models_long[
    patient_models_long$model_role == "extended_age_sex_plus_covariate_same_appendicitis",
    ,
    drop = FALSE
  ]

  key <- c(
    "miRNA", "contrast", "numerator", "denominator", "method",
    "adjustment_variable", "n_complete", "n_mild", "n_severe"
  )
  ref <- ref[, c(key, "log2FC", "p_value", "FDR"), drop = FALSE]
  ext <- ext[, c(key, "log2FC", "p_value", "FDR"), drop = FALSE]

  patient_model_comparison <- dplyr::left_join(
    ref, ext, by = key, suffix = c("_reference", "_adjusted")
  )
  patient_model_comparison$delta_log2FC_adjusted_minus_reference <-
    patient_model_comparison$log2FC_adjusted - patient_model_comparison$log2FC_reference
  patient_model_comparison$effect_attenuation_pct <- ifelse(
    is.finite(patient_model_comparison$log2FC_reference) &
      abs(patient_model_comparison$log2FC_reference) > 1e-12 &
      is.finite(patient_model_comparison$log2FC_adjusted),
    100 * (1 - abs(patient_model_comparison$log2FC_adjusted) /
      abs(patient_model_comparison$log2FC_reference)),
    NA_real_
  )
  patient_model_comparison$direction_stable <- ifelse(
    is.finite(patient_model_comparison$log2FC_reference) &
      is.finite(patient_model_comparison$log2FC_adjusted),
    sign(patient_model_comparison$log2FC_reference) ==
      sign(patient_model_comparison$log2FC_adjusted),
    NA
  )
  patient_model_comparison$retains_FDR_lt_0.05_after_adjustment <-
    !is.na(patient_model_comparison$FDR_adjusted) &
    patient_model_comparison$FDR_adjusted < ALPHA_FDR
  patient_model_comparison$interpretation <- dplyr::case_when(
    patient_model_comparison$retains_FDR_lt_0.05_after_adjustment &
      is.finite(patient_model_comparison$effect_attenuation_pct) &
      patient_model_comparison$effect_attenuation_pct < 25 ~
      "Appendicitis-only Mild-vs-Severe association remains FDR-significant with limited attenuation.",
    patient_model_comparison$retains_FDR_lt_0.05_after_adjustment ~
      "Appendicitis-only Mild-vs-Severe association remains FDR-significant; review attenuation.",
    TRUE ~
      "Appendicitis-only Mild-vs-Severe association does not remain FDR-significant after adding the covariate."
  )
}

utils::write.csv(
  patient_model_comparison,
  file.path(out_dir, "integrated_candidates_appendicitis_only_adjustment.csv"),
  row.names = FALSE
)


primary_patient_model_comparison <- patient_model_comparison[
  patient_model_comparison$miRNA %in% primary_miRNAs,
  ,
  drop = FALSE
]

utils::write.csv(
  primary_patient_model_comparison,
  file.path(out_dir, "primary_robust_candidate_appendicitis_only_adjustment.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Exploratory ROC/AUC — no Youden cutoff, sensitivity or specificity
# ------------------------------------------------------------------------------

roc_contrasts <- list(
  Appendicitis_vs_NIA = list(
    group = as.character(metadata$Appendicitis_status),
    negative = "NIA",
    positive = "Appendicitis"
  ),
  NIA_vs_Mild = list(
    group = as.character(metadata$Phenotype),
    negative = "NIA",
    positive = "Mild"
  ),
  NIA_vs_Severe = list(
    group = as.character(metadata$Phenotype),
    negative = "NIA",
    positive = "Severe"
  ),
  Mild_vs_Severe = list(
    group = as.character(metadata$Phenotype),
    negative = "Mild",
    positive = "Severe"
  )
)

roc_rows <- list()
for (mir in candidate_miRNAs) {
  expr <- as.numeric(rpm[mir, ])
  for (contrast in names(roc_contrasts)) {
    rc <- roc_contrasts[[contrast]]
    idx <- rc$group %in% c(rc$negative, rc$positive) & is.finite(expr)
    y <- rc$group[idx]
    x <- expr[idx]
    if (sum(y == rc$negative) < 5 || sum(y == rc$positive) < 5) next

    roc_obj <- pROC::roc(
      response = factor(y, levels = c(rc$negative, rc$positive)),
      predictor = x,
      levels = c(rc$negative, rc$positive),
      direction = "<",
      quiet = TRUE
    )
    ci <- pROC::ci.auc(roc_obj)

    roc_rows[[length(roc_rows) + 1]] <- data.frame(
      miRNA = mir,
      contrast = contrast,
      negative_group = rc$negative,
      positive_group = rc$positive,
      n_negative = sum(y == rc$negative),
      n_positive = sum(y == rc$positive),
      AUC = as.numeric(pROC::auc(roc_obj)),
      AUC_CI_low = as.numeric(ci[1]),
      AUC_CI_high = as.numeric(ci[3]),
      pROC_direction = roc_obj$direction,
      AUC_orientation = "fixed_direction_not_optimized",
      interpretation = paste0(
        "Exploratory in-sample directional AUC with fixed pROC direction='<'; ",
        "not external validation; no diagnostic cutoff estimated."
      ),
      stringsAsFactors = FALSE
    )
  }
}
roc_results <- dplyr::bind_rows(roc_rows)
utils::write.csv(roc_results, file.path(out_dir, "ROC_AUC_results.csv"), row.names = FALSE)

# Compact summary
summary <- candidate_panel |>
  dplyr::select(miRNA, downstream_category, best_evidence_class) |>
  dplyr::left_join(
    if (nrow(roc_results) > 0) {
      roc_results |>
        dplyr::group_by(miRNA) |>
        dplyr::summarise(max_exploratory_AUC = max(AUC, na.rm = TRUE), .groups = "drop")
    } else data.frame(miRNA = character(), max_exploratory_AUC = numeric()),
    by = "miRNA"
  )
utils::write.csv(summary, file.path(out_dir, "clinical_association_summary.csv"), row.names = FALSE)

# Workbook
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Candidate_panel")
openxlsx::writeData(wb, "Candidate_panel", candidate_panel)
openxlsx::addWorksheet(wb, "Candidate_hierarchy")
openxlsx::writeData(wb, "Candidate_hierarchy", candidate_hierarchy_summary)
openxlsx::addWorksheet(wb, "Missingness")
openxlsx::writeData(wb, "Missingness", missingness)
openxlsx::addWorksheet(wb, "Clinical_marker_summary")
openxlsx::writeData(wb, "Clinical_marker_summary", clinical_marker_summary)
openxlsx::addWorksheet(wb, "Clinical_marker_tests")
openxlsx::writeData(wb, "Clinical_marker_tests", clinical_marker_tests)
openxlsx::addWorksheet(wb, "Continuous")
openxlsx::writeData(wb, "Continuous", continuous_results)
openxlsx::addWorksheet(wb, "Primary_continuous")
openxlsx::writeData(wb, "Primary_continuous", primary_continuous_results)
openxlsx::addWorksheet(wb, "Binary")
openxlsx::writeData(wb, "Binary", binary_results)
openxlsx::addWorksheet(wb, "Complication_audit")
openxlsx::writeData(wb, "Complication_audit", complication_text_audit)
openxlsx::addWorksheet(wb, "Complication_frequency")
openxlsx::writeData(wb, "Complication_frequency", complication_frequency_summary)
openxlsx::addWorksheet(wb, "Complication_tests")
openxlsx::writeData(wb, "Complication_tests", complication_phenotype_tests)
openxlsx::addWorksheet(wb, "Complication_miRNA")
openxlsx::writeData(wb, "Complication_miRNA", complication_miRNA_associations)
openxlsx::addWorksheet(wb, "Primary_binary")
openxlsx::writeData(wb, "Primary_binary", primary_binary_results)
openxlsx::addWorksheet(wb, "Hemolysis_summary")
openxlsx::writeData(wb, "Hemolysis_summary", hemolysis_group_summary)
openxlsx::addWorksheet(wb, "Hemolysis_group")
openxlsx::writeData(wb, "Hemolysis_group", hemolysis_group_tests)
openxlsx::addWorksheet(wb, "Hemolysis_miRNA")
openxlsx::writeData(wb, "Hemolysis_miRNA", hemolysis_miRNA_associations)
openxlsx::addWorksheet(wb, "Targeted_model_audit")
openxlsx::writeData(wb, "Targeted_model_audit", targeted_model_audit)
openxlsx::addWorksheet(wb, "Inflammation_adjustment")
openxlsx::writeData(wb, "Inflammation_adjustment", targeted_model_comparison)
openxlsx::addWorksheet(wb, "Covariate_correlations")
openxlsx::writeData(wb, "Covariate_correlations", covariate_correlations)
openxlsx::addWorksheet(wb, "Appendicitis_only_audit")
openxlsx::writeData(wb, "Appendicitis_only_audit", patient_model_audit)
openxlsx::addWorksheet(wb, "Appendicitis_only_adjustment")
openxlsx::writeData(wb, "Appendicitis_only_adjustment", patient_model_comparison)
openxlsx::addWorksheet(wb, "ROC")
openxlsx::writeData(wb, "ROC", roc_results)
openxlsx::saveWorkbook(wb, file.path(out_dir, "clinical_associations_and_ROC.xlsx"), overwrite = TRUE)


final_qc <- data.frame(
  analysis_id = context$analysis_id,
  n_samples = nrow(metadata),
  n_NIA = sum(metadata$Phenotype == "NIA"),
  n_mild = sum(metadata$Phenotype == "Mild"),
  n_severe = sum(metadata$Phenotype == "Severe"),
  n_primary_robust_miRNAs = length(primary_miRNAs),
  n_secondary_integrated_miRNAs = length(secondary_miRNAs),
  n_clinical_candidate_miRNAs = length(targeted_miRNAs),
  complication_column_present = "Complicaciones" %in% colnames(metadata),
  n_complication_free_text_notes = sum(
    metadata$Complicaciones_note_present,
    na.rm = TRUE
  ),
  n_postoperative_complications = sum(
    metadata$Postoperative_complication_any,
    na.rm = TRUE
  ),
  n_intraoperative_conversions = sum(
    metadata$Intraoperative_conversion,
    na.rm = TRUE
  ),
  n_perioperative_adverse_events = sum(
    metadata$Perioperative_adverse_event_any,
    na.rm = TRUE
  ),
  n_unclassified_complication_notes = sum(
    metadata$Complicaciones_unclassified_note,
    na.rm = TRUE
  ),
  complications_used_as_adjustment_covariate = FALSE,
  samples_removed_by_script = 0L,
  stringsAsFactors = FALSE
)

utils::write.csv(
  final_qc,
  file.path(
    out_dir,
    "07_final_QC_summary.csv"
  ),
  row.names = FALSE
)

write_session_info(results_dir)
message(
  "07 complete — clinical associations, detailed complications, hemolysis control, ",
  "full-cohort and appendicitis-only covariate sensitivities, covariate audit and ",
  "exploratory ROC generated for the integrated candidate hierarchy from script 06."
)

message(
  "Primary_robust candidates: ",
  length(primary_miRNAs),
  if (length(primary_miRNAs) > 0) paste0(" — ", paste(primary_miRNAs, collapse = "; ")) else ""
)

message(
  "Secondary_integrated candidates: ",
  length(secondary_miRNAs),
  if (length(secondary_miRNAs) > 0) paste0(" — ", paste(secondary_miRNAs, collapse = "; ")) else ""
)

message(
  "The candidate hierarchy was inherited unchanged from script 06; miRNAs outside ",
  "the carried-forward set were not included in clinical or ROC analyses."
)

message(
  "Postoperative complications / intraoperative conversions / perioperative events: ",
  sum(metadata$Postoperative_complication_any, na.rm = TRUE),
  " / ",
  sum(metadata$Intraoperative_conversion, na.rm = TRUE),
  " / ",
  sum(metadata$Perioperative_adverse_event_any, na.rm = TRUE)
)

message(
  "Complications were analysed as downstream outcomes and were not used as ",
  "severity adjustment covariates."
)
