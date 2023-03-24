#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WorkflowMain.initialise(workflow, params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTALKER } from './workflows/variantalker'

process help {
    label "help"

    script:
    """
    echo "Usage: nextflow run main.nf -c nextflow.config"
    echo ""
    echo "This pipeline performs whole exome sequencing clinical annotation (snp, indel, cnv)."
    echo ""
    echo "Options:"
    echo "  --somatic.input_snp_indel <input_snp_indel>  A tab-separated file with two columns: tumor_type and path to VCF file containing somatic SNPs/Indels." 
    echo "  --germline.input_snp_indel <input_snp_indel> The input germline SNP/INDEL VCF file."
    echo "  --somatic.input_cnv <input_cnv>              The input somatic CNV file."
    echo "  --tumoronly                                  Flag to indicate that the input is tumor-only (default: false)."
    echo "  --pipeline <pipeline>                        The pipeline used to generate the input data (options: 'Sarek', 'DRAGEN' or 'Iontorrent' (only SNP/INDEL), default: 'Sarek')."
    echo ""
    """
}

workflow {
    VARIANTALKER ()
}