#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowMain.initialise(workflow, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, params.fasta, params.funcotator_somatic_db, params.funcotator_germline_db, 
    params.target, params.annovar_db, params.annovar_software_folder, 
    params.cancervar_evidence_file, params.intervar_evidence_file
    ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTALKER } from './workflows/variantalker'
include { BIOMARKERS; CLONAL_TMB } from './workflows/biomarkers'

process help {
    label "help"

    script:
    """
    echo "Usage: nextflow run main.nf -c nextflow.config"
    echo ""
    echo "This pipeline performs whole exome sequencing clinical annotation (snp, indel, cnv)."
    echo ""
    echo "Options:"
    echo "  --input samplesheet.csv                      CSV file with 4 columns: patient, tumor_tissue, sample_file, sample_type." 
    echo "                                               Available sample_type: somatic, germline and cnv. For somatic and germline" 
    echo "                                               vcf.gz file is accepted for sample_file. For CNV, either cnr from CNVkit or" 
    echo "                                               vcf.gz file from dragen is accepted." 
    echo "                                               For available tumor type, checkout https://github.com/zhanyinx/variantalker/blob/main/README.md" 
    echo "  --analysis <analysis>                        Type of analysis: annotation (default), biomarkers, clonal_tmb." 
    echo "  --pipeline <pipeline>                        The pipeline used to generate the input data (options: 'Sarek', 'DRAGEN'"
    echo "                                               or 'Iontorrent' (only SNP/INDEL), default: 'Sarek')." 
    echo "  --tumoronly                                  Flag to indicate that the input is tumor-only (default: false)."
    echo ""
    """
}

workflow {
    if (!params.analysis || params.analysis.isEmpty() || params.analysis == "annotation") {
        VARIANTALKER ()
    }else if(params.analysis == "biomarkers"){
        println("Biomarkers BETA version!")
        BIOMARKERS()
    }else if(params.analysis == "clonal_tmb"){
        println("clonal tmb BETA version, works only with conda env for now")
        CLONAL_TMB()
    }else{
        println("Provide one of the following analysis: annotation, biomarkers, clonal_tmb")
        println("e.g. nextflow run main.nf --analysis annotation")
    }
    
}