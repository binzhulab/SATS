# SATS Nextflow Example

This directory provides a minimal Nextflow workflow showing how to call the SATS
`GeneratePanelSize()` function from a reproducible workflow runner.

The example expects an input `.rda` file containing an object named
`genomic_information`. The object should include the columns required by
`GeneratePanelSize()`, including chromosome, start position, end position and
assay identifier fields. The input file may also include the optional objects
`Class`, `SBS_order` and `ref.genome`; if these are not provided, the workflow
uses `Class = "SBS"`, `SBS_order = "COSMIC"` and `ref.genome = "hg19"`.

Run the example from the repository root:

```bash
nextflow run nextflow/example1/main.nf \
  --outdir nextflow_results/example1 \
  --outfile outfile.rda
```

To run with the repository Docker configuration, first build the local SATS
image, then enable the Docker profile:

```bash
docker build -t sats:1.0.10 .
nextflow run nextflow/example1/main.nf -profile docker \
  --outdir nextflow_results/example1 \
  --outfile outfile.rda
```

The workflow writes an `.rda` file containing the generated panel-size matrix.
If SATS is already installed in the R environment, the helper script uses the
installed package. Otherwise, it installs the current GitHub source package into
a local R library within the Nextflow work directory.

This example is intended as a template for workflow integration. Production
pipelines should pin the SATS version, use an explicit R/container environment,
and provide project-specific input validation.
