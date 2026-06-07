# SATS News

## SATS 1.0.9

- Added `GenerateVMatrix()` to generate SBS96 or DBS78 mutation-count matrices
  from MAF-like mutation-record tables.
- Extended `GenerateLMatrix()` so users can prepare an `L` matrix directly from
  a panel-coordinate table and clinical sample table.
- Added small SBS and DBS mutation-record examples and unit tests for the new
  preprocessing workflow.
- Added named-matrix alignment checks so `V`, `L`, `W` and `H` are reordered
  when they contain the same sample/context/signature names and rejected when
  named IDs differ.
- Clarified that SATS supports MAF-like mutation records and panel-coordinate
  tables, but does not directly parse raw VCF or BED files.

## SATS 1.0.8

- Updated the current source package under `source/`.
- Added unit tests for `CalculateSignatureBurdens()`, `EstimateSigActivity()` and
  `GeneratePanelSize()` using simulated data and stored expected results.
- Added a GitHub Actions workflow for multi-platform `R CMD check`.
- Added a minimal Nextflow example for calling `GeneratePanelSize()` from a
  workflow runner.
- Updated installation, Docker, citation and web-resource documentation to align
  with the revised SATS manuscript and the current source package.
- Retained older package archives under `old_versions/` for reproducibility.

## Earlier versions

Older package archives are available under `old_versions/`.
