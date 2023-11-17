/*
All processes to extract biomarkers from maf files
*/

/******************************** RNA sequencing data ********************************/

// process to extract transcript per milion (tpm) from rnaseq count data from dragen
process extract_tpm{
    cpus 1
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}", mode: "copy"

    tag "extract_signatures_tpm"

    input: 
        tuple val(patient), path(counts)
    output:
        file("${patient}.rna.tpm.csv")
    script:
    """
        gene_expression.py -i ${counts} -o ${patient}.rna.tpm.csv
    """
}


/******************************** DNA sequencing data ********************************/

// process to extract mutational signature, total and nmd-escapees tumor mutational burden from maf file 
process calculate_tmb_signature{
    cpus 1
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}/", mode: "copy"
    tag "tmb calculation"

    input:
        tuple val(patient), path(maf)
    output:
        file(maf)
        file("tmb_signatures.${patient}.txt")
    script:

    """
    extract_signatures.R -i ${maf} \
        -g ${params.build} \
        -o signatures.txt

    calculate_tmb.py -m ${maf} \
        -t ${params.target} \
        -o tmb.txt

    cat signatures.txt tmb.txt > tmb_signatures.${patient}.txt
    """
}

