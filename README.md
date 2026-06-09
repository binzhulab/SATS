<div align="center">

<img src="docs/assets/sats-logo.png" alt="Signature Analyzer for Targeted Sequencing (SATS)" width="780">

![Version](https://img.shields.io/badge/version-1.0.10-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.1.0-276DC3)
![License](https://img.shields.io/badge/license-CC%20BY--NC%204.0-lightgrey)
![Tests](https://img.shields.io/badge/tests-testthat-green)

<p>
  <strong>National Cancer Institute (NCI) Web Tools</strong><br>
  <a href="https://analysistools.cancer.gov/mutational-signatures/#/catalog/STS">Interactive targeted-sequencing signature catalogue</a><br>
  <a href="https://analysistools.cancer.gov/mutational-signatures/#/refitting">Online signature refitting tool</a>
</p>

[User Guide](User_Guide_SATS_v1.0.10.md) | [User Guide PDF](User_Guide_SATS_v1.0.10.pdf) | [R Manual](SATS-manual.pdf) | [Project Webpage](docs/)

</div>

Signature Analyzer for Targeted Sequencing (SATS) is a panel-aware framework for mutational signature analysis in targeted sequencing data. Unlike tools developed primarily for whole-exome sequencing (WES) or whole-genome sequencing (WGS), SATS models panel-specific sequence context and mutation opportunity, enabling generation of SBS/DBS mutation count matrices from MAF-like mutation records or simple single-sample Variant Call Format (VCF) files, panel-context generation from assay-coordinate tables or Browser Extensible Data (BED) files, de novo signature extraction, mapping to tumor mutational burden (TMB)-normalized reference signatures, individual-tumor signature refitting and calculation of signature-attributed mutation burdens.

The accompanying manuscript applies SATS to 111,711 tumors from American Association for Cancer Research (AACR) Project GENIE (Genomics Evidence Neoplasia Information Exchange) to construct a real-world, panel-calibrated pan-cancer catalogue of targeted sequencing-derived mutational signatures. The package and repository support analysis of targeted-panel cohorts and use of the catalogue in settings where WES/WGS data are unavailable.

## Current Software Status

The current reviewer-response version is **SATS v1.0.10**. This update adds lightweight converters for preparing SATS inputs from simple single-sample VCF files and BED target-region files (`ReadVCFAsMutationRecord()` and `ReadBEDAsPanelInfo()`), while retaining preprocessing utilities for constructing matched mutation-count and panel-context matrices from MAF-like mutation records and panel annotations (`GenerateVMatrix()` and `GenerateLMatrix()`). The repository also includes an expanded executable user guide, regression tests for the main user-facing functions, a minimal Nextflow example, Dockerfile and GitHub Actions R-CMD-check workflow.

The GENIE version 13.0-public panels analyzed in the manuscript remain substantially smaller than exome-scale assays. Supplementary Table 3 lists 56 targeted panels with assay lengths from 0.05 Mb to 9.95 Mb, with a median of 1.47 Mb and no panels in the 10-50 Mb, 50-80 Mb or 80 Mb-WGS ranges. SATS is panel-size aware and can be adapted to larger targeted panels, but WES/WGS remains preferred when available for de novo discovery or low-burden rare signatures.

## Study and Catalogue Overview

SATS was developed using AACR Project GENIE version 13.0-public, a real-world targeted-sequencing cohort spanning clinical sequencing programs in North America and Europe.

<p align="center">
  <img width="780" alt="AACR Project GENIE participating center distribution across North America and Europe" src="docs/assets/genie-site-distribution.png"><br>
  <em><strong>Figure 1A.</strong> AACR Project GENIE participating centers and sample counts used for the targeted-sequencing mutational-signature catalogue. Center labels show the participating-center acronym, and numbers in parentheses indicate the number of tumors contributed by that center.</em>
</p>

<p align="center">
  <img width="900" alt="Sample size by cancer category in AACR Project GENIE version 13" src="docs/assets/genie-cancer-type-sample-size.png"><br>
  <em><strong>Figure 1B.</strong> Distribution of 111,711 tumors across the 23 cancer categories used for downstream targeted-sequencing mutational-signature analysis.</em>
</p>

The resulting SATS catalogue includes **26 single base substitution (SBS) signatures** and **12 double base substitution (DBS) signatures** detected from targeted sequencing data. Dot size indicates the proportion of tumors carrying each signature within a cancer category, and the stacked bars summarize signature-attributed mutation burden.

<p align="center">
  <img width="900" alt="SATS SBS signature catalogue across cancer categories" src="docs/assets/sats-sbs-catalogue.png"><br>
  <em><strong>Figure 2A.</strong> Pan-cancer catalogue of 26 single base substitution (SBS) signatures detected from AACR Project GENIE targeted sequencing data using SATS. Rows represent cancer categories and columns represent SBS signatures. Dot size indicates the proportion of tumors carrying each signature within a cancer category, and stacked bars summarize signature-attributed mutation burden on the targeted-sequencing scale.</em>
</p>

<p align="center">
  <img width="900" alt="SATS DBS signature catalogue across cancer categories" src="docs/assets/sats-dbs-catalogue.png"><br>
  <em><strong>Figure 2B.</strong> Pan-cancer catalogue of 12 double base substitution (DBS) signatures detected from AACR Project GENIE targeted sequencing data using SATS. Rows represent cancer categories and columns represent DBS signatures. Dot size indicates the proportion of tumors carrying each signature within a cancer category, and stacked bars summarize signature-attributed mutation burden on the targeted-sequencing scale.</em>
</p>

---

## Quick Install

The current source version in this repository is **SATS v1.0.10**. Before the v1.0.10 branch is merged into `main`, install the current branch explicitly:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(
    "binzhulab/SATS",
    ref = "software-reviewer-response-updates",
    subdir = "source",
    upgrade = "never"
)

library(SATS)

if (packageVersion("SATS") < "1.0.10" ||
    !"GenerateVMatrix" %in% getNamespaceExports("SATS")) {
    stop("This workflow requires SATS >= 1.0.10. Reinstall SATS and restart R.")
}
```

Alternatively, download `SATS_1.0.10.tar.gz` and install the source archive:

```bash
R CMD INSTALL SATS_1.0.10.tar.gz
```

SATS was formerly available from the [Comprehensive R Archive Network (CRAN)](https://CRAN.R-project.org/package=SATS). CRAN currently lists the package as archived, so the GitHub source installation above is the recommended route for the current version.

---

## National Cancer Institute (NCI) Web Tools

- [Interactive targeted-sequencing signature catalogue](https://analysistools.cancer.gov/mutational-signatures/#/catalog/STS): SBS and DBS signature frequencies, etiologies and cancer-type patterns from the targeted-sequencing catalogue.
- [Online signature refitting tool](https://analysistools.cancer.gov/mutational-signatures/#/refitting): web interface for targeted sequencing signature refitting.

---

## Workflow

SATS separates panel-context generation, de novo signature detection, signature mapping, individual-tumor refitting and burden calculation.

> Targeted sequencing data should not be analyzed as if all tumors shared the same mutation opportunity. SATS models the panel context directly, then estimates signatures and burdens on the targeted-sequencing scale.

<p align="center">
  <img width="900" alt="SATS workflow schematic" src="https://github.com/binzhulab/SATS/assets/51965629/64b226ef-58c1-4fc5-aca1-2be4c4a7cf6b">
</p>

1. **Prepare matched mutation-count and panel-context matrices** from MAF-like mutation records, simple single-sample VCF files, panel coordinates, BED target-region files and sample-panel annotations using `ReadVCFAsMutationRecord()`, `ReadBEDAsPanelInfo()`, `GenerateVMatrix()` and `GenerateLMatrix()`.
2. **Detect de novo signatures** using panel-adjusted opportunity counts.
3. **Map reference signatures** with `MappingSignature()` and Catalogue of Somatic Mutations in Cancer (COSMIC) TMB-normalized signatures.
4. **Refit and estimate burdens** with `EstimateSigActivity()` and `CalculateSignatureBurdens()`.

---

## Usage and Examples

The full executable workflow is maintained in the [User Guide](User_Guide_SATS_v1.0.10.md) and [User Guide PDF](User_Guide_SATS_v1.0.10.pdf). It includes:

- converting simple single-sample VCF and BED files into SATS-compatible mutation-record and panel-coordinate tables;
- generating matched `V` and `L` matrices from MAF-like mutation records and panel-coordinate tables;
- checking and aligning sample IDs and mutation-context rows between `V` and `L`;
- selecting the initial `signeR()` discovery strategy: use individual samples directly for small cohorts, for example fewer than 100 samples, and use 100-sample pooled profiles for very large cohorts, such as analyses with about 10,000 tumors;
- mapping de novo profiles to TMB-normalized COSMIC reference signatures with `MappingSignature()`;
- estimating signature activities and mutation burdens with `EstimateSigActivity()` and `CalculateSignatureBurdens()`.

Keeping the detailed code in one guide avoids duplicated examples and makes the workflow easier to validate end to end.

---

## Repository Layout

- [`source/`](source/): current R package source, including preprocessing functions for MAF-like mutation records and simple single-sample VCF/BED input preparation.
- [`SATS_1.0.10.tar.gz`](SATS_1.0.10.tar.gz): source archive for the current version.
- [`User_Guide_SATS_v1.0.10.md`](User_Guide_SATS_v1.0.10.md): current user guide.
- [`SATS-manual.pdf`](SATS-manual.pdf): function-level R manual.
- [`Generating_L/`](Generating_L/): panel-context generation helper scripts and example panel files.
- [`nextflow/example1/`](nextflow/example1/): minimal Nextflow example for `GeneratePanelSize()`.
- [`docs/`](docs/): static project webpage for GitHub Pages.
- [`old_versions/`](old_versions/): older package archives.

---

## Software Quality

Unit tests are provided in `source/tests/testthat/` for the main user-facing functions, including `ReadVCFAsMutationRecord()`, `ReadBEDAsPanelInfo()`, `GenerateVMatrix()`, `GenerateLMatrix()`, `GeneratePanelSize()`, `CalculateSignatureBurdens()` and `EstimateSigActivity()`. The tests use simulated package data and small single-sample VCF, BED and MAF-like mutation-record examples.

Run tests locally:

```bash
cd source
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Run a local package check:

```bash
R CMD check source --no-manual
```

The repository also includes a GitHub Actions workflow, `.github/workflows/R-CMD-check.yaml`, that runs `R CMD check` on Linux, macOS and Windows. A Dockerfile is provided for building an R environment with SATS and its core genomic dependencies installed.

---

## Input Expectations

SATS accepts either summarized mutation-count and panel-context matrices, MAF-like mutation-record tables with panel-coordinate data frames, or simple single-sample VCF/BED input files that are converted with `ReadVCFAsMutationRecord()` and `ReadBEDAsPanelInfo()`. The VCF/BED converters are intended for standard targeted-panel input preparation. Complex VCF normalization, multi-sample genotype parsing, tumor-normal genotype interpretation, phasing and representation of complex events should be handled upstream when needed.

The row order of the mutation catalogue matrix `V`, panel-context matrix `L` and reference signature matrix `W` must match. For SBS analyses, the `SBS_order` argument controls mutation-type ordering only; the COSMIC reference-signature version used for mapping is controlled separately by `MappingSignature(COSMICv=...)`, with `"v3.4"` as the current default.

---

## Citation

If you use SATS or the targeted-sequencing mutational-signature catalogue, please cite:

Lee et al., "A real-world pan-cancer catalogue of mutational signatures from 111,711 tumors" (submitted).
