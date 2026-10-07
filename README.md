# Pediatric appendicitis circulating miRNA analysis

Reproducible R workflow describing the analytical methodology used to study circulating microRNA (miRNA) profiles in pediatric patients undergoing surgery for suspected appendicitis.

## Analytical workflow

The repository implements the following methodological stages:

1. Input-data validation and harmonization of sequencing counts, clinical metadata, and sequencing metrics.
2. Reads-per-million (RPM) normalization using the complete miRNA library size before expression filtering.
3. Prespecified group-aware expression filtering based on abundance and prevalence.
4. Descriptive quality assessment using PCA, heatmaps, sample-level diagnostics, and a sequencing-based hemolysis proxy.
5. Normalization assessment using DANA-derived correlation metrics.
6. Rank-based RPM/miRPM differential-expression analysis.
7. Count-based differential-expression analysis with edgeR using TMM normalization and DESeq2 using median-of-ratios normalization.
8. Age- and sex-adjusted count-model sensitivity analyses.
9. Cross-method integration of differential-expression evidence.
10. Exploratory clinical association and ROC analyses.
11. Exploratory validated-target and pathway-enrichment analyses.

Benjamini-Hochberg false-discovery rate correction is used for multiple testing. Exploratory downstream analyses do not alter the prespecified differential-expression evidence hierarchy.

## Expression filtering

The primary expression-filtering rule retains a miRNA when it reaches:

- RPM >= 100;
- in at least 50% of samples;
- in at least one phenotype group.

Filtering is applied after RPM normalization from the complete library-size denominator.

## Differential-expression strategy

Five prespecified comparisons are supported by the workflow:

- appendicitis versus non-inflamed appendix;
- global phenotype comparison;
- mild appendicitis versus non-inflamed appendix;
- severe appendicitis versus non-inflamed appendix;
- severe versus mild appendicitis.

The workflow combines complementary rank-based and count-based approaches. Unadjusted count-based models are primary, with age- and sex-adjusted models used as sensitivity analyses.

## Repository structure

- `R/00_analysis_config.R` — shared configuration and analysis parameters.
- `R/00_packages_setup.R` — package checks and session setup.
- `R/00_utils.R` — common helper and validation functions.
- `R/01_data_loading_and_filtering.R` — input harmonization, RPM normalization, filtering, and cohort preparation.
- `R/02_descriptive_analysis.R` — descriptive QC, PCA, heatmaps, and hemolysis-related analyses.
- `R/03_DANA_analysis.R` — normalization assessment using DANA metrics.
- `R/04_DE_analysis_miRPM.R` — rank-based RPM/miRPM differential-expression analyses.
- `R/05_DE_validation_edgeR_DESeq2.R` — edgeR and DESeq2 analyses, including age/sex-adjusted models.
- `R/06_integrated_method_comparison.R` — cross-method evidence integration.
- `R/07_clinical_associations.R` — exploratory clinical, complication, and ROC analyses.
- `R/08_in_silico_analysis.R` — exploratory validated-target and functional-enrichment analyses.

## Public repository scope

This repository contains the reproducible analytical code and methodological configuration only. Study input data, generated result tables and figures, manuscript files, study-specific sample identifiers, and observed candidate/result values are intentionally not committed.

To reproduce the workflow with authorized study data, place the required input files under `data/` using the expected filenames and run the scripts in the order shown below. Generated `results/`, `figures/`, and manuscript/editorial files are excluded from version control.

## Configuration

The public configuration is methodological and does not hard-code study-specific sample identifiers or observed result counts. Input paths and analysis settings are resolved through the project configuration and environment variables where appropriate.

The main shared thresholds are defined centrally in `R/00_analysis_config.R` so that downstream scripts do not redefine them independently.

## Running the workflow

Scripts are run sequentially:

1. `R/00_packages_setup.R`
2. `R/00_analysis_config.R`
3. `R/00_utils.R`
4. `R/01_data_loading_and_filtering.R`
5. `R/02_descriptive_analysis.R`
6. `R/03_DANA_analysis.R`
7. `R/04_DE_analysis_miRPM.R`
8. `R/05_DE_validation_edgeR_DESeq2.R`
9. `R/06_integrated_method_comparison.R`
10. `R/07_clinical_associations.R`
11. `R/08_in_silico_analysis.R`

Each analytical stage performs its own validation and writes outputs to the configured local analysis directories.

## Data availability

Deidentified input data used for the study are available from the corresponding author upon reasonable request, subject to applicable ethical and institutional requirements.

## License

See `LICENSE`.

## Contact

Sergio Pérez Oliveira
