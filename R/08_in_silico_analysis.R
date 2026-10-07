# 08_in_silico_analysis.R
# Cohort-aware revision; candidate input comes directly from the downstream hierarchy frozen by script 06.
# ==============================================================================
# Final in-silico analysis using reliable local target databases
# ==============================================================================

#' This script avoids fragile online per-miRNA queries.
#'
#' Recommended target source:
#' - miRTarBase full download table: experimentally validated MTIs.
#'
#' Final analysis scope:
#' - integrated clinical candidates defined from script 06 evidence:
#'     * Primary_robust:
#'         FDR_concordant_robust / FDR_cross_supported_robust
#'     * Secondary_integrated:
#'         All_three_nominal_concordant
#' - enrichment is run once per miRNA;
#' - no historical miRNA panel is supplied to this script;
#' - no keyword-based clinical/pathway concordance;
#' - functional results are exploratory and never modify candidate priority;
#' - enrichment uses a target-database-derived universe matched to the evidence
#'   class used for each miRNA (strong validated or all validated targets).
#'
#' Optional sources:
#' - TarBase full table.
#' - TargetScan / miRDB prediction tables, used only as secondary evidence.
#'
#' Place downloaded target files in:
#' data/in_silico_targets/
#'
#' Accepted file types:
#' - .csv, .tsv, .txt, .xlsx, .xls
#'
#' If no reliable local target file is found, the script stops.
#' This is intentional: no targets will be invented.
# ==============================================================================


# 0. Configuration ------------------------------------------------------------

source(here::here("R", "00_analysis_config.R"))
source(here::here("R", "00_utils.R"))

context <- get_analysis_context()
results_dir <- context$results_dir
figures_dir <- context$figures_dir

integration_dir <- file.path(
  results_dir,
  "integrated"
)

clinical_dir <- file.path(
  results_dir,
  "clinical_associations_ROC"
)

target_input_dir <- file.path(
  project_root,
  "data",
  "in_silico_targets"
)

output_dir <- file.path(
  results_dir,
  "in_silico"
)

target_output_dir <- file.path(
  output_dir,
  "targets"
)

enrichment_dir <- file.path(
  output_dir,
  "enrichment"
)

figure_dir <- file.path(
  figures_dir,
  "in_silico"
)

target_figure_dir <- file.path(
  figure_dir,
  "target_summary"
)

for (directory in c(
  target_input_dir,
  output_dir,
  target_output_dir,
  enrichment_dir,
  figure_dir,
  target_figure_dir
)) {
  ensure_dir(directory)
}

# Set to TRUE only if you have high-confidence prediction files.
include_predicted_targets <- FALSE

# Use only stronger experimental evidence for enrichment when possible.
# This avoids very broad CLIP-only target lists dominating the functional analysis.
use_strong_validated_targets_for_enrichment <- TRUE

# Minimum number of retained target genes required to run enrichment.
minimum_genes_for_enrichment <- 10


# 1. Dependencies -------------------------------------------------------------

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "enrichplot",
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "readxl",
  "stringr",
  "purrr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages required by script 08: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ". Run R/00_packages_setup.R first.",
    call. = FALSE
  )
}


# 2. Candidate input -----------------------------------------------------------
#
# Candidate selection and prioritization are fixed upstream in script 06.
#
# Script 08 consumes:
#
#   integrated/downstream_candidate_priority.csv
#
# and does NOT recreate evidence classes or candidate tiers.
#
# Candidate categories are inherited from the evidence-integration rules in
# script 06; no study-specific candidate count is hard-coded here.
#
# Functional results are exploratory and never promote, demote or reorder
# candidates.

priority_file <- file.path(
  integration_dir,
  "downstream_candidate_priority.csv"
)

clinical_summary_file <- file.path(
  clinical_dir,
  "clinical_association_summary.csv"
)

if (!file.exists(priority_file)) {
  stop(
    "Missing script-06 downstream candidate hierarchy: ",
    priority_file,
    "\nRun R/06_integrated_method_comparison.R first.",
    call. = FALSE
  )
}

candidate_table <- utils::read.csv(
  priority_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_candidate_columns <- c(
  "miRNA",
  "evidence_class",
  "downstream_category",
  "best_evidence_class",
  "downstream_rank"
)

if (!all(
  required_candidate_columns %in%
    colnames(candidate_table)
)) {
  stop(
    "downstream_candidate_priority.csv must contain: ",
    paste(
      required_candidate_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  anyDuplicated(
    candidate_table$miRNA
  ) > 0L
) {
  stop(
    "Duplicated miRNAs detected in downstream_candidate_priority.csv.",
    call. = FALSE
  )
}

if (
  !all(
    candidate_table$downstream_category %in%
      c(
        "Primary_robust",
        "Secondary_integrated"
      )
  )
) {
  stop(
    "Unexpected downstream category in downstream_candidate_priority.csv.",
    call. = FALSE
  )
}

candidate_table <- candidate_table |>
  dplyr::filter(
    !is.na(miRNA),
    miRNA != ""
  ) |>
  dplyr::arrange(
    downstream_rank
  ) |>
  dplyr::mutate(
    final_priority = downstream_category
  )

candidate_miRNAs <- candidate_table$miRNA

primary_miRNAs <- candidate_table |>
  dplyr::filter(
    downstream_category == "Primary_robust"
  ) |>
  dplyr::pull(miRNA) |>
  unique()

secondary_miRNAs <- candidate_table |>
  dplyr::filter(
    downstream_category == "Secondary_integrated"
  ) |>
  dplyr::pull(miRNA) |>
  unique()

utils::write.csv(
  candidate_table,
  file.path(
    output_dir,
    "candidate_input_integrated.csv"
  ),
  row.names = FALSE
)

clinical_summary <- if (
  file.exists(clinical_summary_file)
) {

  utils::read.csv(
    clinical_summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

} else {

  warning(
    "Clinical summary from script 07 was not found. ",
    "Functional analysis will proceed without exploratory AUC context."
  )

  data.frame()
}

message(
  "Integrated candidate hierarchy inherited from script 06: ",
  length(primary_miRNAs),
  " Primary_robust + ",
  length(secondary_miRNAs),
  " Secondary_integrated = ",
  length(candidate_miRNAs),
  " miRNAs."
)

message(
  "Integrated candidate miRNAs: ",
  paste(
    candidate_miRNAs,
    collapse = ", "
  )
)


# 3. Create target-input template ---------------------------------------------

# Write a fresh candidate-specific template into the analysis output.
# Do not reuse/overwrite the legacy shared template in data/in_silico_targets/.
template_file <- file.path(
  output_dir,
  "target_input_template_current_candidates.csv"
)

template <- data.frame(
  miRNA = candidate_miRNAs,
  gene_symbol = NA_character_,
  source_database = NA_character_,
  evidence_type = NA_character_,
  species_miRNA = "Homo sapiens",
  species_target = "Homo sapiens",
  stringsAsFactors = FALSE
)

utils::write.csv(
  template,
  template_file,
  row.names = FALSE
)


# 4. Helper functions ----------------------------------------------------------

sanitize_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_\\-]+",
    "_",
    x
  )
}

standardize_colnames <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

read_target_file <- function(file_path) {

  extension <- tolower(
    tools::file_ext(file_path)
  )

  message(
    "Reading target file: ",
    basename(file_path)
  )

  if (extension %in% c("xlsx", "xls")) {

    sheets <- readxl::excel_sheets(file_path)

    tables <- lapply(
      sheets,
      function(sheet_name) {
        out <- tryCatch(
          {
            readxl::read_excel(
              file_path,
              sheet = sheet_name
            ) |>
              as.data.frame()
          },
          error = function(e) data.frame()
        )

        if (nrow(out) > 0) {
          out$source_sheet <- sheet_name
        }

        out
      }
    )

    table <- dplyr::bind_rows(tables)

  } else if (extension == "csv") {

    table <- utils::read.csv(
      file_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else {

    table <- tryCatch(
      {
        utils::read.delim(
          file_path,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      },
      error = function(e) {
        utils::read.table(
          file_path,
          sep = ",",
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )
  }

  if (nrow(table) == 0) {
    return(data.frame())
  }

  table$source_file <- basename(file_path)
  table
}

pick_column <- function(
    colnames_vector,
    patterns
) {

  clean_names <- standardize_colnames(
    colnames_vector
  )

  for (pattern in patterns) {
    match_index <- grep(
      pattern,
      clean_names,
      ignore.case = TRUE
    )

    if (length(match_index) > 0) {
      return(
        colnames_vector[match_index[1]]
      )
    }
  }

  NA_character_
}

is_human_species <- function(x) {

  x <- tolower(
    trimws(
      as.character(x)
    )
  )

  is.na(x) |
    x == "" |
    grepl(
      "homo sapiens|human|^hsa$",
      x
    )
}


standardize_target_table <- function(table) {

  if (nrow(table) == 0) {
    return(data.frame())
  }

  original_names <- colnames(table)
  clean_names <- standardize_colnames(original_names)

  names(table) <- clean_names

  miRNA_col <- pick_column(
    colnames(table),
    c(
      "^mirna$",
      "mirna_id",
      "mirna_name",
      "mature_mirna",
      "mirna"
    )
  )

  gene_col <- pick_column(
    colnames(table),
    c(
      "^gene_symbol$",
      "target_symbol",
      "target_gene",
      "target_symbol",
      "^symbol$",
      "gene"
    )
  )

  source_col <- pick_column(
    colnames(table),
    c(
      "source_database",
      "^database$",
      "database_name",
      "source_file"
    )
  )

  evidence_col <- pick_column(
    colnames(table),
    c(
      "evidence_type",
      "experiments",
      "experiment",
      "support_type",
      "validation",
      "method"
    )
  )

  species_miRNA_col <- pick_column(
    colnames(table),
    c(
      "species_mirna",
      "mirna_species",
      "species"
    )
  )

  species_target_col <- pick_column(
    colnames(table),
    c(
      "species_target",
      "target_species",
      "target_organism"
    )
  )

  if (
    is.na(miRNA_col) ||
    is.na(gene_col)
  ) {
    warning(
      "Could not identify miRNA/gene columns in file: ",
      paste(
        unique(table$source_file),
        collapse = "; "
      )
    )
    return(data.frame())
  }

  out <- data.frame(
    miRNA = as.character(table[[miRNA_col]]),
    gene_symbol = as.character(table[[gene_col]]),
    source_database = if (!is.na(source_col)) {
      as.character(table[[source_col]])
    } else {
      as.character(table$source_file)
    },
    evidence_type = if (!is.na(evidence_col)) {
      as.character(table[[evidence_col]])
    } else {
      NA_character_
    },
    species_miRNA = if (!is.na(species_miRNA_col)) {
      as.character(table[[species_miRNA_col]])
    } else {
      NA_character_
    },
    species_target = if (!is.na(species_target_col)) {
      as.character(table[[species_target_col]])
    } else {
      NA_character_
    },
    source_file = as.character(table$source_file),
    stringsAsFactors = FALSE
  )

  out |>
    dplyr::mutate(
      miRNA = trimws(miRNA),
      gene_symbol = trimws(gene_symbol),
      gene_symbol = toupper(gene_symbol),
      source_database = dplyr::coalesce(
        source_database,
        source_file
      ),
      source_database_lower = tolower(
        paste(
          source_database,
          source_file
        )
      ),
      evidence_lower = tolower(
        dplyr::coalesce(
          evidence_type,
          ""
        )
      ),
      source_class = dplyr::case_when(
        grepl("mirtarbase|hsa_mti|hsa mti|tarbase|mirecords", source_database_lower) ~ "validated",
        grepl("reporter|luciferase|western|qpcr|rt_qpcr|rt-qpcr|clip|ngs|microarray", evidence_lower) ~ "validated",
        grepl("targetscan|mirdb|mirwalk|microRNA_Target_Sites|target_sites", source_database_lower) ~ "predicted",
        TRUE ~ "unknown"
      ),
      strong_validation = source_class == "validated" &
        grepl(
          "reporter|luciferase|western|qpcr|rt_qpcr|rt-qpcr|functional|strong",
          evidence_lower,
          ignore.case = TRUE
        )
    ) |>
    dplyr::filter(
      !is.na(miRNA),
      miRNA != "",
      !is.na(gene_symbol),
      gene_symbol != ""
    )
}


# 5. Load local target files ---------------------------------------------------

target_files <- list.files(
  target_input_dir,
  pattern = "\\.(csv|tsv|txt|xlsx|xls)$",
  full.names = TRUE,
  ignore.case = TRUE
)

target_files <- target_files[
  !grepl(
    "target_input_template",
    basename(target_files),
    ignore.case = TRUE
  )
]

if (length(target_files) == 0) {
  stop(
    "No local target files were found in:\n",
    target_input_dir,
    "\n\nRecommended: download the full miRTarBase MTI table and place it in that folder.\n",
    "A template has been created here:\n",
    template_file
  )
}

raw_target_tables <- lapply(
  target_files,
  read_target_file
)

standardized_targets_all <- dplyr::bind_rows(
  lapply(
    raw_target_tables,
    standardize_target_table
  )
)

if (nrow(standardized_targets_all) == 0) {
  stop(
    "Target files were found, but no valid miRNA-target pairs could be parsed.\n",
    "Check column names. Required concepts: miRNA and gene symbol.\n",
    "Template file: ",
    template_file
  )
}



# Keep human records when species metadata are available.
standardized_targets_all <- standardized_targets_all |>
  dplyr::filter(
    is_human_species(species_miRNA),
    is_human_species(species_target)
  )

if (nrow(standardized_targets_all) == 0L) {
  stop(
    "No human miRNA-target records remained after species filtering.",
    call. = FALSE
  )
}

standardized_targets <- standardized_targets_all |>
  dplyr::filter(
    miRNA %in% candidate_miRNAs
  )

if (nrow(standardized_targets) == 0L) {
  stop(
    "Target databases were parsed successfully, but none of the integrated ",
    "candidate miRNAs had matching target records.",
    call. = FALSE
  )
}


# Reproducibility manifest for local target sources.
target_file_info <- file.info(
  target_files
)

target_input_manifest <- data.frame(
  source_file = basename(target_files),
  file_size_bytes = as.numeric(target_file_info$size),
  modified_time = as.character(target_file_info$mtime),
  md5 = unname(
    tools::md5sum(
      target_files
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  target_input_manifest,
  file.path(
    output_dir,
    "target_input_manifest.csv"
  ),
  row.names = FALSE
)


# 6. Build target summary ------------------------------------------------------

gene_target_summary_all <- standardized_targets_all |>
  dplyr::group_by(
    miRNA,
    gene_symbol
  ) |>
  dplyr::summarise(
    validated = any(
      source_class == "validated",
      na.rm = TRUE
    ),
    strong_validation = any(
      strong_validation,
      na.rm = TRUE
    ),
    predicted = any(
      source_class == "predicted",
      na.rm = TRUE
    ),
    n_sources = dplyr::n_distinct(
      source_database[
        !is.na(source_database)
      ]
    ),
    source_databases = paste(
      sort(
        unique(
          source_database[
            !is.na(source_database)
          ]
        )
      ),
      collapse = "; "
    ),
    evidence_types = paste(
      sort(
        unique(
          evidence_type[
            !is.na(evidence_type)
          ]
        )
      ),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    retained_for_enrichment =
      dplyr::case_when(
        use_strong_validated_targets_for_enrichment & strong_validation ~ TRUE,
        !use_strong_validated_targets_for_enrichment & validated ~ TRUE,
        include_predicted_targets & predicted & n_sources >= 2 ~ TRUE,
        TRUE ~ FALSE
      )
  ) |>
  dplyr::rename(
    input_miRNA = miRNA
  )

gene_target_summary <- gene_target_summary_all |>
  dplyr::filter(
    input_miRNA %in% candidate_miRNAs
  )

all_target_records <- standardized_targets |>
  dplyr::rename(
    input_miRNA = miRNA
  )

utils::write.csv(
  all_target_records,
  file = file.path(
    target_output_dir,
    "all_local_target_records.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  gene_target_summary,
  file = file.path(
    target_output_dir,
    "target_gene_summary.csv"
  ),
  row.names = FALSE
)


target_database_background_summary <- data.frame(
  metric = c(
    "human_miRNA_target_records",
    "unique_human_miRNAs",
    "unique_human_target_genes",
    "validated_target_genes",
    "strong_validated_target_genes"
  ),
  value = c(
    nrow(standardized_targets_all),
    dplyr::n_distinct(standardized_targets_all$miRNA),
    dplyr::n_distinct(standardized_targets_all$gene_symbol),
    dplyr::n_distinct(
      gene_target_summary_all$gene_symbol[
        gene_target_summary_all$validated
      ]
    ),
    dplyr::n_distinct(
      gene_target_summary_all$gene_symbol[
        gene_target_summary_all$strong_validation
      ]
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  target_database_background_summary,
  file.path(
    target_output_dir,
    "target_database_background_summary.csv"
  ),
  row.names = FALSE
)


# 7. Target figures ------------------------------------------------------------

target_count_summary <- gene_target_summary |>
  dplyr::group_by(
    input_miRNA
  ) |>
  dplyr::summarise(
    n_total_targets = dplyr::n_distinct(gene_symbol),
    n_retained_targets = sum(
      retained_for_enrichment,
      na.rm = TRUE
    ),
    n_validated_targets = sum(
      validated,
      na.rm = TRUE
    ),
    n_strong_validated_targets = sum(
      strong_validation,
      na.rm = TRUE
    ),
    n_predicted_targets = sum(
      predicted,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::right_join(
    data.frame(
      input_miRNA = candidate_miRNAs,
      stringsAsFactors = FALSE
    ),
    by = "input_miRNA"
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.numeric),
      ~ dplyr::coalesce(.x, 0)
    ),
    input_miRNA = factor(
      input_miRNA,
      levels = rev(candidate_miRNAs)
    )
  )

utils::write.csv(
  target_count_summary,
  file = file.path(
    target_output_dir,
    "target_count_summary.csv"
  ),
  row.names = FALSE
)

target_count_plot <- ggplot2::ggplot(
  target_count_summary,
  ggplot2::aes(
    x = input_miRNA,
    y = n_retained_targets
  )
) +
  ggplot2::geom_col(
    fill = "grey45"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Retained target genes per miRNA",
    subtitle = "Validated targets retained; predicted targets optional and disabled by default.",
    x = NULL,
    y = "Retained target genes"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )

ggplot2::ggsave(
  filename = prepare_output_file(
    file.path(
      target_figure_dir,
      "retained_target_genes_per_miRNA.png"
    )
  ),
  plot = target_count_plot,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

save_ggplot_pdf(
  plot = target_count_plot,
  filename = file.path(
    target_figure_dir,
    "retained_target_genes_per_miRNA.pdf"
  ),
  width = 7,
  height = 5
)


# 8. Analysis sets -------------------------------------------------------------
#
# Enrichment is run once per integrated miRNA.
#
# No pooled "Primary", "Secondary" or combined-candidate target union is used here.
# Pooling miRNAs would weight the functional analysis toward miRNAs with larger
# validated target sets and would answer a different biological question.

analysis_sets <- as.list(
  candidate_miRNAs
)

names(
  analysis_sets
) <- candidate_miRNAs


# 9. Enrichment helpers --------------------------------------------------------

convert_symbols_to_entrez <- function(gene_symbols) {

  gene_symbols <- unique(
    gene_symbols[
      !is.na(gene_symbols) &
        gene_symbols != ""
    ]
  )

  if (length(gene_symbols) == 0) {
    return(
      data.frame(
        SYMBOL = character(),
        ENTREZID = character()
      )
    )
  }

  mapped <- withCallingHandlers(
    suppressMessages(
      clusterProfiler::bitr(
        gene_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    ),
    warning = function(w) {
      msg <- conditionMessage(w)

      expected_mapping_warning <- grepl(
        "of input gene IDs are fail to map",
        msg,
        fixed = TRUE
      )

      if (expected_mapping_warning) {
        message(
          "Gene ID mapping: ",
          msg
        )
        invokeRestart("muffleWarning")
      }
    }
  )

  mapped |>
    dplyr::distinct(
      SYMBOL,
      ENTREZID
    )
}


# Build enrichment universes from the same local human target database(s).
#
# This avoids using the whole annotated genome as the default background when
# the target list itself comes from an experimentally ascertained miRNA-target
# resource.

background_gene_sets <- list(

  strong_validated =
    unique(
      gene_target_summary_all$gene_symbol[
        gene_target_summary_all$strong_validation
      ]
    ),

  validated =
    unique(
      gene_target_summary_all$gene_symbol[
        gene_target_summary_all$validated
      ]
    )
)


background_entrez_sets <- lapply(
  background_gene_sets,
  function(symbols) {

    mapped <- convert_symbols_to_entrez(
      symbols
    )

    unique(
      mapped$ENTREZID
    )
  }
)


background_mapping_audit <- data.frame(
  target_set = names(background_gene_sets),
  n_gene_symbols = vapply(
    background_gene_sets,
    length,
    integer(1)
  ),
  n_entrez_ids = vapply(
    background_entrez_sets,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)


utils::write.csv(
  background_mapping_audit,
  file.path(
    enrichment_dir,
    "enrichment_background_audit.csv"
  ),
  row.names = FALSE
)


run_enrichment_for_set <- function(
    set_name,
    miRNAs
) {

  message(
    "Running enrichment for set: ",
    set_name
  )

  genes <- gene_target_summary |>
    dplyr::filter(
      input_miRNA %in% miRNAs,
      retained_for_enrichment
    ) |>
    dplyr::pull(
      gene_symbol
    ) |>
    unique()

  target_set_used <- if (isTRUE(use_strong_validated_targets_for_enrichment)) {
    "strong_validated"
  } else {
    "validated"
  }

  if (
    length(genes) < minimum_genes_for_enrichment &&
    isTRUE(use_strong_validated_targets_for_enrichment)
  ) {

    fallback_genes <- gene_target_summary |>
      dplyr::filter(
        input_miRNA %in% miRNAs,
        validated
      ) |>
      dplyr::pull(
        gene_symbol
      ) |>
      unique()

    if (length(fallback_genes) >= minimum_genes_for_enrichment) {
      genes <- fallback_genes
      target_set_used <- "validated_fallback"
    }
  }

  if (length(genes) < minimum_genes_for_enrichment) {

    message(
      "Fewer than ",
      minimum_genes_for_enrichment,
      " retained genes for ",
      set_name,
      ". Enrichment skipped."
    )

    return(
      list(
        target_genes = genes,
        target_set_used = target_set_used,
        universe_key = if (
          target_set_used == "strong_validated"
        ) {
          "strong_validated"
        } else {
          "validated"
        },
        universe_entrez_n = NA_integer_,
        GO = data.frame(),
        KEGG = data.frame(),
        Reactome = data.frame(),
        GO_object = NULL,
        KEGG_object = NULL,
        Reactome_object = NULL
      )
    )
  }

  conversion <- convert_symbols_to_entrez(
    genes
  )

  entrez_ids <- unique(
    conversion$ENTREZID
  )

  universe_key <- if (
    target_set_used == "strong_validated"
  ) {
    "strong_validated"
  } else {
    "validated"
  }

  universe_entrez <- background_entrez_sets[[
    universe_key
  ]]

  if (
    length(universe_entrez) <
      minimum_genes_for_enrichment
  ) {
    message(
      "The matched enrichment universe for ",
      set_name,
      " contains fewer than ",
      minimum_genes_for_enrichment,
      " mapped genes. Enrichment skipped."
    )

    return(
      list(
        target_genes = genes,
        target_set_used = target_set_used,
        universe_key = universe_key,
        universe_entrez_n = length(universe_entrez),
        GO = data.frame(),
        KEGG = data.frame(),
        Reactome = data.frame(),
        GO_object = NULL,
        KEGG_object = NULL,
        Reactome_object = NULL
      )
    )
  }

  if (length(entrez_ids) < minimum_genes_for_enrichment) {

    message(
      "Fewer than ",
      minimum_genes_for_enrichment,
      " mapped Entrez IDs for ",
      set_name,
      ". Enrichment skipped."
    )

    return(
      list(
        target_genes = genes,
        target_set_used = target_set_used,
        universe_key = universe_key,
        universe_entrez_n = length(universe_entrez),
        GO = data.frame(),
        KEGG = data.frame(),
        Reactome = data.frame(),
        GO_object = NULL,
        KEGG_object = NULL,
        Reactome_object = NULL
      )
    )
  }

  go_result <- tryCatch(
    {
      clusterProfiler::enrichGO(
        gene = entrez_ids,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
    },
    error = function(e) NULL
  )

  kegg_result <- tryCatch(
    {
      clusterProfiler::enrichKEGG(
        gene = entrez_ids,
        universe = universe_entrez,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
    },
    error = function(e) NULL
  )

  reactome_result <- tryCatch(
    {
      ReactomePA::enrichPathway(
        gene = entrez_ids,
        universe = universe_entrez,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
    },
    error = function(e) NULL
  )

  list(
    target_genes = genes,
    target_set_used = target_set_used,
    universe_key = universe_key,
    universe_entrez_n = length(universe_entrez),
    GO = if (is.null(go_result)) data.frame() else as.data.frame(go_result),
    KEGG = if (is.null(kegg_result)) data.frame() else as.data.frame(kegg_result),
    Reactome = if (is.null(reactome_result)) data.frame() else as.data.frame(reactome_result),
    GO_object = go_result,
    KEGG_object = kegg_result,
    Reactome_object = reactome_result
  )
}


# 10. Run enrichment -----------------------------------------------------------

enrichment_results <- lapply(
  names(analysis_sets),
  function(set_name) {
    run_enrichment_for_set(
      set_name = set_name,
      miRNAs = analysis_sets[[set_name]]
    )
  }
)

names(enrichment_results) <- names(analysis_sets)

all_enrichment_tables <- list()

for (set_name in names(enrichment_results)) {

  set_result <- enrichment_results[[set_name]]

  set_output_dir <- file.path(
    enrichment_dir,
    sanitize_filename(set_name)
  )

  set_figure_dir <- file.path(
    figure_dir,
    sanitize_filename(set_name)
  )

  for (directory in c(
    set_output_dir,
    set_figure_dir
  )) {
    ensure_dir(directory)
  }

  utils::write.csv(
    data.frame(
      gene_symbol = set_result$target_genes
    ),
    file = file.path(
      set_output_dir,
      "retained_target_genes.csv"
    ),
    row.names = FALSE
  )

  for (database_name in c(
    "GO",
    "KEGG",
    "Reactome"
  )) {

    enrichment_table <- set_result[[database_name]]

    if (nrow(enrichment_table) == 0) {
      next
    }

    enrichment_table$analysis_set <- set_name
    enrichment_table$database <- database_name
    enrichment_table$target_set_used <- set_result$target_set_used
    enrichment_table$universe_key <- set_result$universe_key
    enrichment_table$universe_entrez_n <- set_result$universe_entrez_n

    enrichment_table$evidence_level <- dplyr::case_when(
      enrichment_table$p.adjust < 0.05 ~ "FDR < 0.05",
      enrichment_table$p.adjust < 0.10 ~ "FDR < 0.10",
      enrichment_table$pvalue < 0.05 ~ "Raw p < 0.05",
      TRUE ~ "Not significant"
    )

    utils::write.csv(
      enrichment_table,
      file = file.path(
        set_output_dir,
        paste0(
          database_name,
          "_enrichment.csv"
        )
      ),
      row.names = FALSE
    )

    all_enrichment_tables[[
      length(all_enrichment_tables) + 1
    ]] <- enrichment_table

    enrichment_object <- set_result[[
      paste0(
        database_name,
        "_object"
      )
    ]]

    plot_table <- enrichment_table |>
      dplyr::filter(
        p.adjust < 0.10 |
          pvalue < 0.05
      ) |>
      dplyr::arrange(
        p.adjust,
        pvalue
      ) |>
      dplyr::slice_head(
        n = 15
      )

    if (nrow(plot_table) > 0) {

      plot_table <- plot_table |>
        dplyr::mutate(
          Description = factor(
            Description,
            levels = rev(unique(Description))
          ),
          minus_log10_FDR = -log10(
            p.adjust + .Machine$double.eps
          )
        )

      dotplot_object <- ggplot2::ggplot(
        plot_table,
        ggplot2::aes(
          x = minus_log10_FDR,
          y = Description
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            size = Count,
            colour = evidence_level
          )
        ) +
        ggplot2::labs(
          title = paste(
            database_name,
            "enrichment -",
            set_name
          ),
          x = "-log10(FDR)",
          y = NULL,
          size = "Gene count",
          colour = "Evidence"
        ) +
        ggplot2::theme_classic(
          base_size = 11
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )

      ggplot2::ggsave(
        filename = prepare_output_file(
          file.path(
            set_figure_dir,
            paste0(
              database_name,
              "_dotplot.png"
            )
          )
        ),
        plot = dotplot_object,
        width = 9,
        height = 7,
        dpi = 600,
        bg = "white"
      )

      save_ggplot_pdf(
        plot = dotplot_object,
        filename = file.path(
          set_figure_dir,
          paste0(
            database_name,
            "_dotplot.pdf"
          )
        ),
        width = 9,
        height = 7
      )
    }
  }
}

all_enrichment <- dplyr::bind_rows(
  all_enrichment_tables
)

if (nrow(all_enrichment) == 0) {
  all_enrichment <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    analysis_set = character(),
    database = character(),
    target_set_used = character(),
    universe_key = character(),
    universe_entrez_n = integer(),
    evidence_level = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  all_enrichment,
  file = file.path(
    enrichment_dir,
    "all_enrichment_results.csv"
  ),
  row.names = FALSE
)



enrichment_run_audit <- dplyr::bind_rows(
  lapply(
    names(enrichment_results),
    function(set_name) {

      x <- enrichment_results[[
        set_name
      ]]

      data.frame(
        miRNA = set_name,
        target_set_used = x$target_set_used,
        universe_key = x$universe_key,
        universe_entrez_n = x$universe_entrez_n,
        n_retained_target_genes = length(
          unique(
            x$target_genes
          )
        ),
        n_GO_terms = nrow(
          x$GO
        ),
        n_KEGG_terms = nrow(
          x$KEGG
        ),
        n_Reactome_terms = nrow(
          x$Reactome
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  enrichment_run_audit,
  file.path(
    enrichment_dir,
    "enrichment_run_audit.csv"
  ),
  row.names = FALSE
)


# 11. Final functional interpretation ------------------------------------------
#
# No keyword-based clinical/pathway concordance is performed. Functional
# enrichment is descriptive and exploratory and is not used for candidate
# selection, promotion, or clinical validation.

summarize_database <- function(
    enrichment_table,
    database_name,
    candidate_miRNA,
    top_n = 5
) {

  x <- enrichment_table |>
    dplyr::filter(
      analysis_set == candidate_miRNA,
      database == database_name,
      is.finite(p.adjust)
    ) |>
    dplyr::arrange(
      p.adjust,
      pvalue
    )

  n_fdr <- sum(
    x$p.adjust < 0.05,
    na.rm = TRUE
  )

  top_terms <- x |>
    dplyr::filter(
      p.adjust < 0.05
    ) |>
    dplyr::slice_head(
      n = top_n
    ) |>
    dplyr::pull(
      Description
    )

  list(
    n_FDR_lt_0.05 = n_fdr,
    top_terms = if (
      length(top_terms) > 0
    ) {
      paste(
        top_terms,
        collapse = "; "
      )
    } else {
      "None"
    }
  )
}


functional_summary_rows <- lapply(
  candidate_miRNAs,
  function(mir) {

    target_row <- target_count_summary[
      as.character(
        target_count_summary$input_miRNA
      ) == mir,
      ,
      drop = FALSE
    ]

    go_summary <- summarize_database(
      all_enrichment,
      "GO",
      mir
    )

    kegg_summary <- summarize_database(
      all_enrichment,
      "KEGG",
      mir
    )

    reactome_summary <- summarize_database(
      all_enrichment,
      "Reactome",
      mir
    )

    candidate_row <- candidate_table[
      candidate_table$miRNA == mir,
      ,
      drop = FALSE
    ]

    clinical_auc <- NA_real_

    if (
      nrow(clinical_summary) > 0 &&
        "miRNA" %in%
          colnames(
            clinical_summary
          ) &&
        "max_exploratory_AUC" %in%
          colnames(
            clinical_summary
          )
    ) {

      tmp_auc <- clinical_summary$max_exploratory_AUC[
        clinical_summary$miRNA == mir
      ]

      if (length(tmp_auc) > 0) {
        clinical_auc <- tmp_auc[1]
      }
    }

    data.frame(
      miRNA = mir,
      final_priority = candidate_row$final_priority[1],
      best_evidence_class = if (
        "best_evidence_class" %in%
          colnames(
            candidate_row
          )
      ) {
        candidate_row$best_evidence_class[1]
      } else {
        NA_character_
      },
      max_exploratory_AUC = clinical_auc,
      n_total_validated_targets = if (
        nrow(target_row) > 0
      ) {
        target_row$n_validated_targets[1]
      } else {
        NA_integer_
      },
      n_strong_validated_targets = if (
        nrow(target_row) > 0
      ) {
        target_row$n_strong_validated_targets[1]
      } else {
        NA_integer_
      },
      n_retained_targets_initial_strong_filter = if (
        nrow(target_row) > 0
      ) {
        target_row$n_retained_targets[1]
      } else {
        NA_integer_
      },
      target_set_used = if (
        mir %in% enrichment_run_audit$miRNA
      ) {
        enrichment_run_audit$target_set_used[
          match(
            mir,
            enrichment_run_audit$miRNA
          )
        ]
      } else {
        NA_character_
      },
      n_target_genes_used_for_enrichment = if (
        mir %in% enrichment_run_audit$miRNA
      ) {
        enrichment_run_audit$n_retained_target_genes[
          match(
            mir,
            enrichment_run_audit$miRNA
          )
        ]
      } else {
        0L
      },
      enrichment_universe = if (
        mir %in% enrichment_run_audit$miRNA
      ) {
        enrichment_run_audit$universe_key[
          match(
            mir,
            enrichment_run_audit$miRNA
          )
        ]
      } else {
        NA_character_
      },
      enrichment_universe_entrez_n = if (
        mir %in% enrichment_run_audit$miRNA
      ) {
        enrichment_run_audit$universe_entrez_n[
          match(
            mir,
            enrichment_run_audit$miRNA
          )
        ]
      } else {
        NA_integer_
      },
      enrichment_performed = if (
        mir %in% enrichment_run_audit$miRNA
      ) {
        enrichment_run_audit$n_retained_target_genes[
          match(
            mir,
            enrichment_run_audit$miRNA
          )
        ] >= minimum_genes_for_enrichment
      } else {
        FALSE
      },
      n_GO_FDR_lt_0.05 = go_summary$n_FDR_lt_0.05,
      top_GO_terms = go_summary$top_terms,
      n_KEGG_FDR_lt_0.05 = kegg_summary$n_FDR_lt_0.05,
      top_KEGG_terms = kegg_summary$top_terms,
      n_Reactome_FDR_lt_0.05 =
        reactome_summary$n_FDR_lt_0.05,
      top_Reactome_terms =
        reactome_summary$top_terms,
      interpretation = paste0(
        "Functional target and enrichment results are exploratory. ",
        "The reported target_set_used and n_target_genes_used_for_enrichment ",
        "reflect the actual enrichment input after any validated-target fallback. ",
        "These results do not alter the candidate hierarchy fixed in script 06 ",
        "or constitute independent clinical validation."
      ),
      stringsAsFactors = FALSE
    )
  }
)

interpretation_table <- dplyr::bind_rows(
  functional_summary_rows
)

utils::write.csv(
  interpretation_table,
  file = file.path(
    output_dir,
    "final_in_silico_interpretation.csv"
  ),
  row.names = FALSE
)


# 12. Workbook -----------------------------------------------------------------



workbook <- openxlsx::createWorkbook()

add_sheet <- function(
    workbook,
    sheet_name,
    data
) {

  safe_sheet <- substr(
    gsub(
      "[\\[\\]\\*\\?/\\\\:]",
      "_",
      sheet_name
    ),
    1,
    31
  )

  openxlsx::addWorksheet(
    workbook,
    safe_sheet
  )

  openxlsx::writeData(
    workbook,
    safe_sheet,
    data
  )

  openxlsx::freezePane(
    workbook,
    safe_sheet,
    firstRow = TRUE
  )

  if (ncol(data) > 0) {
    openxlsx::setColWidths(
      workbook,
      safe_sheet,
      cols = seq_len(ncol(data)),
      widths = "auto"
    )
  }
}

add_sheet(
  workbook,
  "Candidate_input",
  candidate_table
)

add_sheet(
  workbook,
  "Target_input_manifest",
  target_input_manifest
)

add_sheet(
  workbook,
  "Background_summary",
  target_database_background_summary
)

add_sheet(
  workbook,
  "Enrichment_audit",
  enrichment_run_audit
)

add_sheet(
  workbook,
  "Target_counts",
  target_count_summary
)

add_sheet(
  workbook,
  "Target_summary",
  gene_target_summary
)

add_sheet(
  workbook,
  "All_enrichment",
  all_enrichment
)

add_sheet(
  workbook,
  "Final_interpretation",
  interpretation_table
)

openxlsx::saveWorkbook(
  workbook,
  file = file.path(
    output_dir,
    "in_silico_analysis.xlsx"
  ),
  overwrite = TRUE
)



# 13. Final QC / audit ----------------------------------------------------------

candidate_target_coverage <- target_count_summary |>
  dplyr::mutate(
    input_miRNA = as.character(
      input_miRNA
    )
  )

final_qc <- data.frame(
  analysis_id = context$analysis_id,
  n_integrated_candidates = length(candidate_miRNAs),
  n_primary_robust = length(primary_miRNAs),
  n_secondary_integrated = length(secondary_miRNAs),
  n_target_source_files = length(target_files),
  n_human_target_records = nrow(standardized_targets_all),
  n_candidates_with_any_target = sum(
    candidate_target_coverage$n_total_targets > 0,
    na.rm = TRUE
  ),
  n_candidates_with_retained_target = sum(
    candidate_target_coverage$n_retained_targets > 0,
    na.rm = TRUE
  ),
  n_candidates_with_enrichment_gene_set = sum(
    enrichment_run_audit$n_retained_target_genes >=
      minimum_genes_for_enrichment,
    na.rm = TRUE
  ),
  target_universe = "same local human target database; evidence-matched",
  predicted_targets_enabled = include_predicted_targets,
  historical_candidate_list_used = FALSE,
  functional_results_modify_candidate_priority = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  final_qc,
  file.path(
    output_dir,
    "08_final_QC_summary.csv"
  ),
  row.names = FALSE
)


# 14. Completion summary -------------------------------------------------------

message("")
message("08 complete — in-silico target/enrichment analysis completed for the integrated candidates.")
message("")
message("Target source files:")
message(
  paste(
    "-",
    basename(target_files),
    collapse = "\n"
  )
)
message("")
message(
  "Candidate hierarchy: ",
  length(primary_miRNAs),
  " Primary_robust + ",
  length(secondary_miRNAs),
  " Secondary_integrated = ",
  length(candidate_miRNAs)
)
message("")
message("Targets retained:")
message("- Strong validated targets are used for enrichment when available.")
message("- If too few strong targets are available, the script falls back to all validated targets.")
message("- Predicted targets are ignored unless include_predicted_targets = TRUE.")
message(
  "- Enrichment universe: genes represented in the same local human target ",
  "database and matched to the target evidence level used."
)
message(
  "- Functional enrichment does not alter candidate priority and is exploratory."
)
message("")
message("Main outputs:")
message("- ", file.path(output_dir, "candidate_input_integrated.csv"))
message("- ", file.path(output_dir, "target_input_manifest.csv"))
message("- ", file.path(target_output_dir, "target_gene_summary.csv"))
message("- ", file.path(target_output_dir, "target_count_summary.csv"))
message("- ", file.path(target_output_dir, "target_database_background_summary.csv"))
message("- ", file.path(enrichment_dir, "all_enrichment_results.csv"))
message("- ", file.path(enrichment_dir, "enrichment_run_audit.csv"))
message("- ", file.path(enrichment_dir, "enrichment_background_audit.csv"))
message("  Note: enrichment tables include evidence_level; prioritize FDR < 0.05, then FDR < 0.10.")
message("- ", file.path(output_dir, "final_in_silico_interpretation.csv"))
message("- ", file.path(output_dir, "in_silico_analysis.xlsx"))
message("- ", file.path(output_dir, "08_final_QC_summary.csv"))
message("- ", file.path(target_figure_dir, "retained_target_genes_per_miRNA.png"))
