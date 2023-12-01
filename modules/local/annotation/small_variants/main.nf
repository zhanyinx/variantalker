// List of functions and processes used to annotate snp and indel from whole exome sequences

// add civic data to vcf
process add_civic{
    cpus 10
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }

    tag "civic2vcf"

    input:
        tuple val(meta), file(vcf), file(index)
    output:
        tuple val(meta), file("${meta.patient}.vcf.gz"), file("${meta.patient}.vcf.gz.tbi")
    script:
    """
    zcat ${vcf} > appo.vcf

    # split on chromosomes
    awk '{if(\$1~/^#/) print \$0}' appo.vcf > header
    for chr in `awk '{if(!(\$1~/^#/)) print \$1}' appo.vcf | sort | uniq`; do
        cp header tmp_\$chr
    done
    awk '{if(!(\$1~/^#/)) print \$0 >> "tmp_"\$1 }'  appo.vcf

    cpu=0
    for file in `ls tmp_*`; do
        civicpy annotate-vcf --input-vcf  \$file --output-vcf out.\$file --reference ${params.build_alt_name}  --include-status accepted --include-status submitted > log.\$file 2>&1 &
        let cpu=cpu+1

        if [ \$cpu -eq $task.cpus ]; then
            cpu=0
            wait
        fi
    done
    wait

    if [ -f final.vcf ]; then
        rm final.vcf
    fi

    for file in `ls out*`; do
        awk '{if(!(\$1~/^#/)) print \$0; else print \$0 > "header"}' \$file >> final.vcf
    done

    sort -k1,1V -k2,2n final.vcf > final.sorted.vcf
    sed -i 's/ /_/g' final.sorted.vcf
    cat header final.sorted.vcf > ${meta.patient}.vcf
    rm ${vcf} ${index}
    bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
    tabix -p vcf ${meta.patient}.vcf.gz
    """

}


// alpha_missense
process add_alpha_missense{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 3.GB * task.attempt }
    tag "alphamissense"

    input:
        tuple val(meta), file(maf), file(vcf)
    output:
        tuple val(meta), file("${maf.baseName}.missense.maf"), file("${meta.patient}.vcf")
    script:
    """
        if ! [ -f ${params.genomes[params.build].alpha_missense} ]; then
            echo "${params.genomes[params.build].alpha_missense} alpha Missense file does not exist." > .command.err
            exit 125
        fi

        zcat ${params.genomes[params.build].alpha_missense} > tmp
        awk -F'\\t' 'BEGIN{fn=0; count=0}{
            if(FNR==1) fn++
            if(fn==1){
                if(\$1~/^#/){
                    print \$0
                }else if(\$1=="Hugo_Symbol"){
                    print "# Alpha Missense file ""'"${params.genomes[params.build].alpha_missense}"'"
                    printf "%s\\t%s\\t%s\\n", \$0, "am_pathogenicity", "am_class"
                }else{
                    converter[\$5,\$6] = count
                    string[\$5,\$6]=\$5\$6\$12\$13
                    line[count]=\$0
                    bool[\$5,\$6]++
                    count++
                }
            }

            if(fn==2){
                if(bool[\$1,\$2] > 0){
                    if(string[\$1,\$2] == \$1\$2\$3\$4){
                        am_pathogenicity[converter[\$1,\$2]] = \$9
                        am_class[converter[\$1,\$2]] = \$10
                    }
                }
            }
        }END{for(i=0;i<count;i++) printf "%s\\t%s\\t%s\\n", line[i], am_pathogenicity[i], am_class[i]}' ${maf} tmp > ${maf.baseName}.missense.maf
        rm tmp
    """
}

// filter maf file
process filter_maf{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"
    tag "filtermaf"

    input:
        tuple val(meta), file(maf), file(vcf)
    output:
        tuple val(meta), file("${meta.patient}.maf")
        tuple file("filtered*.pass.tsv"), file("filtered*.nopass.tsv")
        file("${meta.patient}.vcf")
    script:
    """
        touch filtered.${meta.patient}.maf.pass.tsv filtered.${meta.patient}.maf.nopass.tsv

        filter_variants.py -m ${maf} \
         -o ${meta.patient}.maf \
         --filter_intervar "${params.filter_intervar}" \
         --filter_cancervar "${params.filter_cancervar}" \
         --filter_renovo "${params.filter_renovo}" \
         --sample_type ${meta.sample_type} \
         --min_depth ${params.filter_min_depth} \
         --vaf_threshold_germline ${params.filter_vaf_threshold_germline} \
         --vaf_threshold ${params.filter_vaf_threshold} \
         --filter_genes_germline ${params.filter_genes_germline} \
         --filter_genes_somatic ${params.filter_genes_somatic} \
         --filter_variant_classification "${params.filter_var_classification}" \
         --filter_civic "${params.filter_civic_evidence_level}"

    """
    
}


