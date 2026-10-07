# ==============================================================================
# 00_packages_setup.R
# ==============================================================================
#
# Project:
#   Pediatric appendicitis plasma miRNA-seq
#
# Script:
#   00_packages_setup.R
#
# Purpose:
#   Build, validate and document the R software environment required by the
#   complete miRNA-seq analysis pipeline.
#
# Role in pipeline:
#   Environment bootstrap and reproducibility audit.
#
#   This script:
#     - installs missing CRAN packages;
#     - installs missing Bioconductor packages;
#     - installs required GitHub packages;
#     - repairs packages that are installed but cannot be loaded;
#     - updates packages that do not meet pipeline-specific minimum versions;
#     - validates that all required packages are loadable;
#     - records package versions and installation provenance;
#     - records the R and Bioconductor environment;
#     - saves sessionInfo() for reproducibility.
#
#   No biological or statistical analysis is performed in this script.
#
# Inputs:
#   - Current R installation
#   - CRAN
#   - Bioconductor
#   - GitHub repositories:
#       cran/PoissonSeq
#       LXQin/DANA
#       sergio30po/miRPM
#
# Outputs:
#   results/package_versions_setup.csv
#   results/package_environment_setup.csv
#   results/sessionInfo_package_setup.txt
#
# Methodological decisions:
#   - CRAN packages are installed from https://cloud.r-project.org.
#   - Bioconductor packages are installed using BiocManager matched to the
#     installed R version.
#   - PoissonSeq is installed from the archived CRAN mirror on GitHub because
#     it is no longer distributed through current CRAN/Bioconductor releases.
#   - GitHub packages are installed using remotes.
#   - Pipeline-specific minimum package versions are enforced only when an
#     explicit minimum has been defined.
#
# Required bootstrap packages:
#   BiocManager
#   remotes
#
# Reproducibility:
#   This script is the single source of truth for package requirements used by
#   the pipeline. Downstream scripts should not install packages independently.
#
# Usage:
#   Run once from a fresh R session before executing the analysis pipeline:
#
#     source("R/00_packages_setup.R")
#
# ==============================================================================


# ==============================================================================
# 1. DEPENDENCIES
# ==============================================================================

# This script is itself responsible for installing the dependencies required
# downstream. Only base R functionality is assumed at startup.


# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 CRAN repository
# ------------------------------------------------------------------------------

CRAN_REPOSITORY <- "https://cloud.r-project.org"


# ------------------------------------------------------------------------------
# 2.2 CRAN packages
# ------------------------------------------------------------------------------

# BiocManager and remotes are included explicitly because they are required to
# install packages from Bioconductor and GitHub, respectively.

cran_packages <- c(
  "BiocManager",
  "remotes",
  "here",
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "purrr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "circlize",
  "pROC",
  "rstatix",
  "car",
  "RColorBrewer",
  "gridExtra",
  "officer",
  "flextable",
  "rlang"
)


# ------------------------------------------------------------------------------
# 2.3 Bioconductor packages
# ------------------------------------------------------------------------------

bioc_packages <- c(
  "edgeR",
  "DESeq2",
  "EDASeq",
  "RUVSeq",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "reactome.db",
  "ReactomePA",
  "enrichplot"
)


# ------------------------------------------------------------------------------
# 2.4 GitHub packages
# ------------------------------------------------------------------------------

# PoissonSeq is no longer distributed through current CRAN/Bioconductor
# releases. DANA documentation recommends installation from the archived CRAN
# GitHub mirror.
#
# Order is intentional:
#   PoissonSeq is installed before DANA because DANA may require it.

github_packages <- c(
  PoissonSeq = "cran/PoissonSeq",
  DANA       = "LXQin/DANA",
  miRPM      = "sergio30po/miRPM"
)


# ------------------------------------------------------------------------------
# 2.5 Pipeline-specific minimum versions
# ------------------------------------------------------------------------------

# Only package versions that are explicitly required by the current pipeline
# are constrained here.
#
# Dependencies without a project-specific minimum are managed by their package
# managers.

minimum_versions <- c(
  rlang = "1.3.0",
  miRPM = "0.2.0"
)


# Runtime dependencies only.
#
# Avoid installing all Suggests/Enhances automatically because they can greatly
# increase environment size without being required by the analysis.

INSTALL_DEPENDENCIES <- c(
  "Depends",
  "Imports",
  "LinkingTo"
)


# ==============================================================================
# 3. INPUTS
# ==============================================================================

# No biological data are read by this script.
#
# Software inputs are obtained from:
#
#   CRAN         -> cran_packages
#   Bioconductor -> bioc_packages
#   GitHub       -> github_packages
#
# Repository location is resolved after package installation using here::here().


# ==============================================================================
# 4. VALIDATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Package availability
# ------------------------------------------------------------------------------

package_installed <- function(pkg) {
  
  pkg %in% rownames(
    utils::installed.packages()
  )
}


# ------------------------------------------------------------------------------
# 4.2 Package loadability
# ------------------------------------------------------------------------------

package_loadable <- function(pkg) {
  
  isTRUE(
    tryCatch(
      
      suppressWarnings(
        requireNamespace(
          pkg,
          quietly = TRUE
        )
      ),
      
      error = function(e) {
        FALSE
      }
    )
  )
}


# ------------------------------------------------------------------------------
# 4.3 Minimum required version
# ------------------------------------------------------------------------------

minimum_version_for <- function(pkg) {
  
  required <- unname(
    minimum_versions[pkg]
  )
  
  if (
    length(required) == 0L ||
    is.na(required) ||
    !nzchar(required)
  ) {
    return(NULL)
  }
  
  required
}


# ------------------------------------------------------------------------------
# 4.4 Version validation
# ------------------------------------------------------------------------------

package_version_ok <- function(pkg) {
  
  required <- minimum_version_for(pkg)
  
  if (is.null(required)) {
    return(TRUE)
  }
  
  
  installed <- tryCatch(
    utils::packageVersion(pkg),
    error = function(e) NULL
  )
  
  
  if (is.null(installed)) {
    return(FALSE)
  }
  
  
  utils::compareVersion(
    as.character(installed),
    required
  ) >= 0
}


# ------------------------------------------------------------------------------
# 4.5 Complete package validation
# ------------------------------------------------------------------------------

package_ok <- function(pkg) {
  
  package_loadable(pkg) &&
    package_version_ok(pkg)
}


# ------------------------------------------------------------------------------
# 4.6 GitHub repository validation
# ------------------------------------------------------------------------------

github_repository_for <- function(pkg) {
  
  if (!pkg %in% names(github_packages)) {
    return(NA_character_)
  }
  
  unname(
    github_packages[[pkg]]
  )
}


# ==============================================================================
# 5. ANALYSIS / PROCESSING
# ==============================================================================

# No biological analysis occurs in this script.
#
# "Processing" in this context refers exclusively to software installation and
# environment validation.


# ------------------------------------------------------------------------------
# 5.1 Generic installation function
# ------------------------------------------------------------------------------

install_package <- function(
    pkg,
    source = c(
      "cran",
      "bioc",
      "github"
    ),
    repo = NULL
) {
  
  source <- match.arg(source)
  
  
  # Package is already suitable for the pipeline.
  if (package_ok(pkg)) {
    message(
      "OK: ",
      pkg
    )
    
    return(
      invisible(TRUE)
    )
  }
  
  
  installed <- package_installed(pkg)
  
  loadable <- if (installed) {
    package_loadable(pkg)
  } else {
    FALSE
  }
  
  version_ok <- if (installed && loadable) {
    package_version_ok(pkg)
  } else {
    FALSE
  }
  
  
  # ---------------------------------------------------------------------------
  # Installation status message
  # ---------------------------------------------------------------------------
  
  if (!installed) {
    
    message(
      "Installing ",
      pkg,
      " [",
      source,
      "]..."
    )
    
  } else if (!loadable) {
    
    message(
      "Repairing ",
      pkg,
      " [",
      source,
      "]: installed but not loadable..."
    )
    
  } else if (!version_ok) {
    
    message(
      "Updating ",
      pkg,
      " [",
      source,
      "] to meet minimum version ",
      minimum_version_for(pkg),
      "..."
    )
    
  } else {
    
    message(
      "Reinstalling ",
      pkg,
      " [",
      source,
      "]..."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # CRAN
  # ---------------------------------------------------------------------------
  
  if (source == "cran") {
    
    utils::install.packages(
      pkg,
      repos = CRAN_REPOSITORY,
      dependencies = INSTALL_DEPENDENCIES
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Bioconductor
  # ---------------------------------------------------------------------------
  
  if (source == "bioc") {
    
    if (!requireNamespace(
      "BiocManager",
      quietly = TRUE
    )) {
      
      utils::install.packages(
        "BiocManager",
        repos = CRAN_REPOSITORY,
        dependencies = INSTALL_DEPENDENCIES
      )
    }
    
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE,
      force = installed && !package_ok(pkg)
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # GitHub
  # ---------------------------------------------------------------------------
  
  if (source == "github") {
    
    if (
      is.null(repo) ||
      !nzchar(repo)
    ) {
      
      stop(
        "A GitHub repository must be supplied for package '",
        pkg,
        "'.",
        call. = FALSE
      )
    }
    
    
    if (!requireNamespace(
      "remotes",
      quietly = TRUE
    )) {
      
      utils::install.packages(
        "remotes",
        repos = CRAN_REPOSITORY,
        dependencies = INSTALL_DEPENDENCIES
      )
    }
    
    
    remotes::install_github(
      repo,
      dependencies = INSTALL_DEPENDENCIES,
      upgrade = "never",
      force = TRUE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Post-installation validation
  # ---------------------------------------------------------------------------
  
  if (!package_loadable(pkg)) {
    
    stop(
      "Package '",
      pkg,
      "' is installed but cannot be loaded.\n",
      "Restart R and rerun R/00_packages_setup.R.\n",
      "If the error identifies another package, repair that dependency first.",
      call. = FALSE
    )
  }
  
  
  if (!package_version_ok(pkg)) {
    
    stop(
      "Package '",
      pkg,
      "' does not meet the minimum required version.\n",
      "Required: ",
      minimum_version_for(pkg),
      "\nInstalled: ",
      as.character(
        utils::packageVersion(pkg)
      ),
      call. = FALSE
    )
  }
  
  
  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 5.2 Install CRAN packages
# ------------------------------------------------------------------------------

message(
  "\n=== CRAN packages ==="
)

for (pkg in cran_packages) {
  
  install_package(
    pkg = pkg,
    source = "cran"
  )
}


# ------------------------------------------------------------------------------
# 5.3 Install Bioconductor packages
# ------------------------------------------------------------------------------

message(
  "\n=== Bioconductor packages ==="
)

for (pkg in bioc_packages) {
  
  install_package(
    pkg = pkg,
    source = "bioc"
  )
}


# ------------------------------------------------------------------------------
# 5.4 Install GitHub packages
# ------------------------------------------------------------------------------

message(
  "\n=== GitHub packages ==="
)

for (pkg in names(github_packages)) {
  
  install_package(
    pkg = pkg,
    source = "github",
    repo = github_packages[[pkg]]
  )
}


# ==============================================================================
# 6. OUTPUTS
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 Validate project root
# ------------------------------------------------------------------------------

if (!requireNamespace(
  "here",
  quietly = TRUE
)) {
  
  stop(
    "Package 'here' could not be loaded after installation.",
    call. = FALSE
  )
}


project_root <- normalizePath(
  here::here(),
  winslash = "/",
  mustWork = TRUE
)


required_project_dirs <- c(
  file.path(
    project_root,
    "R"
  ),
  file.path(
    project_root,
    "data"
  )
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


# ------------------------------------------------------------------------------
# 6.2 Results directory
# ------------------------------------------------------------------------------

results_root <- file.path(
  project_root,
  "results"
)


dir.create(
  results_root,
  recursive = TRUE,
  showWarnings = FALSE
)


if (!dir.exists(results_root)) {
  
  stop(
    "Could not create results directory: ",
    results_root,
    call. = FALSE
  )
}


# ==============================================================================
# 7. QC / AUDIT
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Final package validation
# ------------------------------------------------------------------------------

packages_to_report <- unique(
  c(
    cran_packages,
    bioc_packages,
    names(github_packages)
  )
)


final_package_status <- vapply(
  packages_to_report,
  package_ok,
  logical(1)
)


if (any(!final_package_status)) {
  
  failed_packages <- packages_to_report[
    !final_package_status
  ]
  
  stop(
    "Package setup validation failed for:\n",
    paste(
      "-",
      failed_packages,
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 7.2 Package source helper
# ------------------------------------------------------------------------------

declared_package_source <- function(pkg) {
  
  if (pkg %in% cran_packages) {
    return("CRAN")
  }
  
  if (pkg %in% bioc_packages) {
    return("Bioconductor")
  }
  
  if (pkg %in% names(github_packages)) {
    return("GitHub")
  }
  
  NA_character_
}


# ------------------------------------------------------------------------------
# 7.3 Package DESCRIPTION helper
# ------------------------------------------------------------------------------

package_description_field <- function(
    pkg,
    field
) {
  
  description <- tryCatch(
    utils::packageDescription(pkg),
    error = function(e) NULL
  )
  
  
  if (
    is.null(description) ||
    is.null(description[[field]]) ||
    length(description[[field]]) == 0L
  ) {
    return(NA_character_)
  }
  
  
  value <- as.character(
    description[[field]]
  )
  
  
  if (
    length(value) == 0L ||
    !nzchar(value[1])
  ) {
    return(NA_character_)
  }
  
  
  value[1]
}


# ------------------------------------------------------------------------------
# 7.4 Package audit table
# ------------------------------------------------------------------------------

package_audit <- data.frame(
  
  package = packages_to_report,
  
  version = vapply(
    packages_to_report,
    function(pkg) {
      
      as.character(
        utils::packageVersion(pkg)
      )
    },
    character(1)
  ),
  
  declared_source = vapply(
    packages_to_report,
    declared_package_source,
    character(1)
  ),
  
  repository = vapply(
    packages_to_report,
    function(pkg) {
      
      repo <- github_repository_for(pkg)
      
      if (is.na(repo)) {
        return("")
      }
      
      repo
    },
    character(1)
  ),
  
  minimum_version = vapply(
    packages_to_report,
    function(pkg) {
      
      minimum <- minimum_version_for(pkg)
      
      if (is.null(minimum)) {
        return("")
      }
      
      minimum
    },
    character(1)
  ),
  
  loadable = vapply(
    packages_to_report,
    package_loadable,
    logical(1)
  ),
  
  version_ok = vapply(
    packages_to_report,
    package_version_ok,
    logical(1)
  ),
  
  remote_sha = vapply(
    packages_to_report,
    function(pkg) {
      
      sha <- package_description_field(
        pkg,
        "RemoteSha"
      )
      
      if (is.na(sha)) {
        return("")
      }
      
      sha
    },
    character(1)
  ),
  
  remote_ref = vapply(
    packages_to_report,
    function(pkg) {
      
      ref <- package_description_field(
        pkg,
        "RemoteRef"
      )
      
      if (is.na(ref)) {
        return("")
      }
      
      ref
    },
    character(1)
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 7.5 Write package audit
# ------------------------------------------------------------------------------

package_versions_file <- file.path(
  results_root,
  "package_versions_setup.csv"
)


utils::write.csv(
  package_audit,
  package_versions_file,
  row.names = FALSE
)


# ==============================================================================
# 8. REPRODUCIBILITY
# ==============================================================================

# ------------------------------------------------------------------------------
# 8.1 Bioconductor version
# ------------------------------------------------------------------------------

bioconductor_version <- tryCatch(
  
  as.character(
    BiocManager::version()
  ),
  
  error = function(e) {
    NA_character_
  }
)


# ------------------------------------------------------------------------------
# 8.2 Environment summary
# ------------------------------------------------------------------------------

environment_summary <- data.frame(
  
  item = c(
    "timestamp_utc",
    "R_version",
    "R_platform",
    "Bioconductor_version",
    "project_root"
  ),
  
  value = c(
    format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    ),
    
    R.version.string,
    
    R.version$platform,
    
    bioconductor_version,
    
    project_root
  ),
  
  stringsAsFactors = FALSE
)


environment_file <- file.path(
  results_root,
  "package_environment_setup.csv"
)


utils::write.csv(
  environment_summary,
  environment_file,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8.3 sessionInfo()
# ------------------------------------------------------------------------------

session_info_file <- file.path(
  results_root,
  "sessionInfo_package_setup.txt"
)


writeLines(
  capture.output(
    sessionInfo()
  ),
  con = session_info_file
)


# ------------------------------------------------------------------------------
# 8.4 Final report
# ------------------------------------------------------------------------------

message(
  "\nPackage setup complete."
)

message(
  "Package audit: ",
  package_versions_file
)

message(
  "Environment audit: ",
  environment_file
)

message(
  "Session information: ",
  session_info_file
)

message(
  "Bioconductor version detected: ",
  bioconductor_version
)


print(
  package_audit,
  row.names = FALSE
)

print(
  environment_summary,
  row.names = FALSE
)
