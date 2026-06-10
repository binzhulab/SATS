# SATS Containerized Nextflow Workflow Example

This example runs a compact SATS workflow from bundled package data. It
validates matched `V`, `L`, `W` and `H` matrices with `ValidateSATSInputs()`,
maps example TMB-normalized profiles with `MappingSignature()`, estimates
sample-level signature activities with `EstimateSigActivity()`, calculates
signature burdens with `CalculateSignatureBurdens()` and writes reproducible
output files.

The workflow intentionally does not run full de novo signature extraction. That
step can require additional method-specific dependencies and should be inserted
upstream when users have cohort-derived de novo profiles. Here, known reference
profiles are used as compact stand-ins for de novo profiles so the workflow
tests the SATS mapping, refitting and burden-estimation steps in a stable
containerized template.

## Build the SATS Docker Image

Run from the repository root:

```bash
docker build -t sats:1.0.10 .
```

## Run with Docker-Enabled Nextflow

Run from the repository root:

```bash
nextflow run nextflow/sats_workflow/main.nf -profile docker \
  --outdir nextflow_results/sats_workflow \
  --n_samples 25
```

If the image uses a different tag, override it through the root
`nextflow.config` parameter:

```bash
nextflow run nextflow/sats_workflow/main.nf -profile docker \
  --sats_container sats:1.0.10 \
  --outdir nextflow_results/sats_workflow
```

## Expected Outputs

The workflow writes the following files to `nextflow_results/sats_workflow/`:

- `sats_mapping_results.csv`
- `sats_activity_matrix.csv`
- `sats_signature_burdens.csv`
- `sats_workflow_outputs.rds`
- `sats_workflow_summary.txt`

This example is intended as a reproducible workflow template rather than a
full clinical production pipeline or an nf-core-compliant pipeline.
