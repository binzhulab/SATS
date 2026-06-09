#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.n_samples = 25
params.outdir = "nextflow_results/sats_workflow"
params.script = "${projectDir}/bin/run_sats_workflow.R"

Channel.fromPath(params.script).set { script_ch }

process SATS_WORKFLOW_EXAMPLE {
  tag "n_samples=${params.n_samples}"
  publishDir params.outdir, mode: 'copy'

  input:
    path script_file

  output:
    path "sats_mapping_results.csv"
    path "sats_activity_matrix.csv"
    path "sats_signature_burdens.csv"
    path "sats_workflow_summary.txt"
    path "sats_workflow_outputs.rds"

  script:
    """
    Rscript ${script_file} --n-samples ${params.n_samples} --out-prefix sats
    """
}

workflow {
  if (params.n_samples < 1) {
    exit 1, "Parameter --n_samples must be at least 1"
  }

  SATS_WORKFLOW_EXAMPLE(script_ch)
}
