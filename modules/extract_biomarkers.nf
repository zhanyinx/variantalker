/*
All processes to extract biomarkers from maf files
*/

/******************************** RNA sequencing data ********************************/

// process to extract transcript per milion (tpm) from rnaseq count data from dragen
process extract_tpm{
    cpus 5
    maxRetries = 2
    memory { 10.GB * task.attempt }
    publishDir "${params.output}/${params.date}/biomarkers/${counts.simpleName}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${counts.simpleName}/biomarkers/", mode: "copy"

    tag "extract_signatures_tpm"

    input: 
        file(counts)
    output:
        file("${counts.simpleName}.rna.tpm.csv")
    script:
    """
        gene_expression.py -i ${counts} -o ${counts.simpleName}.rna.tpm.csv
    """
}


/******************************** DNA sequencing data ********************************/

// process to extract mutational signature, total and nmd-escapees tumor mutational burden from maf file 
process calculate_tmb_signature{
    cpus 2
    memory "2 G"
    publishDir "${params.output}/${params.date}/biomarkers/${maf.simpleName}/", mode: "copy"
    // publishDir "${params.output}/${params.date}/${maf.simpleName}/biomarkers/", mode: "copy"
    tag "tmb calculation"

    input:
        file(maf)
    output:
        file(maf)
        file("tmb_signatures.${maf.simpleName}.txt")
    script:

    """
    extract_signatures.R -i ${maf} \
        -g ${params.build} \
        -o signatures.txt

    calculate_tmb.py -m ${maf} \
        -t ${params.target} \
        -o tmb.txt

    cat signatures.txt tmb.txt > tmb_signatures.${maf.simpleName}.txt
    """
}

process ascat_calling{
    tag "ascat_calling"
    input:
        file(samp)
        file(conf)
    output:
        path("*.clonalTMB.txt")
    script:
    """     
        header="\$(awk -F ',' 'BEGIN{getline; printf "%s", \$1; for(i=2;i<=NF-1;i++) printf ",%s", \$i; print "";}' ${samp})" # get header for sarek input file

        # loop over patients
        for patient in `awk -F ',' 'BEGIN{getline}{print \$1}' ${samp} | sort | uniq`; do 
            # create patient specific input file and run nf sarek with ascat 
            echo \${header} > tmp.csv
            awk -F ',' '{if(\$1 == "'"\$patient"'") {printf "%s", \$1; for(i=2;i<=NF-1;i++) printf ",%s", \$i; print "";}}' ${samp} >> tmp.csv 

            if [[ ${params.build} == "hg19" ]]; then
                cp ${params.target} intervals.bed
                sed -i 's/chr//g' intervals.bed
                nextflow run nf-core/sarek -r 3.1.2 \
                -profile singularity \
                --input tmp.csv \
                -resume \
                -c ${conf} \
                --genome "GATK.GRCh37" \
                --step variant_calling \
                --intervals intervals.bed
            else
                nextflow run nf-core/sarek -r 3.1.2 \
                -profile singularity \
                --input tmp.csv \
                -resume \
                -c ${conf} \
                --genome "GATK.GRCh38" \
                --step variant_calling \
                --intervals ${params.target}
            fi

            # extract cellularity, tumor/normal sample name and name from dragen
            cellularity="\$(awk 'BEGIN{getline; nrow=NF}{if(\$1 == "'"\$patient"'") cellularity = \$nrow}END{print cellularity}' ${samp})"
            tumor="\$(awk '{if(\$1 == "'"\$patient"'" && \$3 == 1) out = \$4}END{print out}' ${samp})"
            normal="\$(awk '{if(\$1 == "'"\$patient"'" && \$3 == 0) out = \$4}END{print out}' ${samp})"
            name=\$tumor
            #name="\$(awk '{if(\$1 == "'"\$patient"'" && \$3 == 1) out = \$4}END{print out}' ${samp})"
            #name="\$(basename \$name | cut -d'.' -f-1)"
            
            # create pyclone input file
            ascat_file="\$(ls results/variant_calling/ascat/\${tumor}_vs_\${normal}/*cnvs.txt)"

            # create output folder
            if ! [ -d ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/ ]; then
                mkdir -p ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/
            fi

            # if ascat files do not exist, can't calculate clonal TMB, set it to NA
            if [ -z \${ascat_file} ]; then
                echo "Clonal TMB: NA" > \${name}.clonalTMB.txt
                cp \${name}.clonalTMB.txt ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/
                continue
            fi

            maf_file="\$(ls ${launchDir}/${params.output}/${params.date}/annotation/somatic/\$name/\$name*maf)"

            if ! [ -z \$cellularity ]; then
                create_input4pyclone.py -as \${ascat_file} -c \${cellularity} -m \${maf_file} -o \${patient}.pyclone.tsv
            else
                ascat_tumor_estimate="\$(ls results/variant_calling/ascat/\${tumor}_vs_\${normal}/*.purityploidy.txt)"
                create_input4pyclone.py -as \${ascat_file} -m \${maf_file} -o \${patient}.pyclone.tsv -ac \${ascat_tumor_estimate}
            fi

            # run pyclone
            pyclone-vi fit -i \${patient}.pyclone.tsv -o \${patient}.pyclone.h5 -c 40 -d beta-binomial -r 20
            pyclone-vi write-results-file -i \${patient}.pyclone.h5 -o \${patient}.pyclone.output.tsv

            # extract clonal tmb
            awk 'BEGIN{
                max=-99; ncount=0; getline
            }{
                if(\$4>max){
                    max = \$4; ncount = 1;
                }else if(\$4 == max){
                    ncount++
                }
            }END{print "Clonal TMB: ", ncount}' \${patient}.pyclone.output.tsv > \${name}.clonalTMB.txt

            cp \${name}.clonalTMB.txt ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/

        done
    """

}
