#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.infile  = "./data/infile.rda"
params.outfile = "outfile.rda"
params.script  = "./bin/run_GeneratePanelSize.R"

// Define a channel for the input file
Channel.fromPath(params.infile).set { infile_ch }
Channel.fromPath(params.script).set { script_ch }

process runAnalysis {
  input: 
    path infile
    path script_file

  output:
    path "${params.outfile}"

  script:
    """
    Rscript ${script_file} ${infile} ${params.outfile}
    """
}

workflow {
    if (!params.infile) {
      exit 1, "input file (--infile) not specified"
    }
    if (!params.outfile) {
      exit 1, "output file (--outfile) not specified"
    }

    runAnalysis(infile_ch, script_ch)
}
