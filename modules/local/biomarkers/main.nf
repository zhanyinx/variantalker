/*
All processes to extract biomarkers from maf files
*/

/******************************** RNA sequencing data ********************************/

// process to extract transcript per milion (tpm) from rnaseq count data from dragen
process extract_tpm{
    fair true
    cpus 1
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}", mode: "copy"

    tag "extract_signatures_tpm"

    input: 
        tuple val(patient), path(counts)
    output:
        tuple val(patient), file("${patient}.rna.tpm.csv")
    script:
    """
        gene_expression.py -i ${counts} -o ${patient}.rna.tpm.csv
    """
}


/******************************** DNA sequencing data ********************************/

// process to extract mutational signature, total and nmd-escapees tumor mutational burden from maf file 
process calculate_tmb_signature{
    fair true
    cpus 1
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}/", mode: "copy"
    tag "tmb calculation"

    input:
        tuple val(patient), path(maf)
    output:
        tuple val(patient), file(maf)
        tuple val(patient), file("tmb_signatures.${patient}.txt")
    script:

    """
    mkdir input
    cp ${maf} input/

    extract_signatures.py -i ./input \
    -g ${params.build_alt_name} -c ${params.cosmic_version}\
    --cosmic_group ${params.cosmic_group} \
    -o signatures.txt

    calculate_tmb.py -m ${maf} \
        -t ${params.target} \
        --nmd_scores ${params.nmd_db} \
        -o tmb.txt

    cat signatures.txt tmb.txt > tmb_signatures.${patient}.txt
    """
}

