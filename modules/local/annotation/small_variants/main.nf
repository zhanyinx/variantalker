// List of functions and processes used to annotate snp and indel from whole exome sequences

process split_chunks{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 2.GB * task.attempt }
    container "docker://yinxiu/civicpy:v1.0"
    tag "split_chunk"

    input:
        tuple val(meta), file(vcf)
    output:
        tuple val(meta), file("${meta.patient}_chunk_*.vcf.gz")
    script:
    """
        zcat ${vcf} > appo.vcf
        nfile=`awk 'BEGIN{counts = 0}{if(\$1==chr && \$2==pos){counts++;}else{counts=0;}; if(\$0~/^#/){print \$0 > "header";} else if(\$7=="PASS"){print \$0 > "body_"counts".vcf"}; pos=\$2; chr=\$1}END{print counts}' appo.vcf`
        for i in `seq 0 \$nfile`; do
            split -l ${params.chunk_size} body_\$i.vcf ${meta.patient}_chunk_\$i
        done

        for i in `ls ${meta.patient}_chunk_*`; do
            cat header \$i > \$i.vcf
            bgzip -c \$i.vcf > \$i.vcf.gz
        done
    """
}

// add civic data to vcf
process add_civic{
    cpus 3
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    container "docker://yinxiu/civicpy:v1.0"

    tag "civic2vcf"

    input:
        tuple val(meta), val(chunk_index), file(vcf) 
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.vcf.gz"), file("${meta.patient}.vcf.gz.tbi")
    script:
    """
    zcat ${vcf} > appo.vcf
    nline=`wc -l appo.vcf | awk '{print \$1}'`
    check=0
    while [ \$nline -ne \$check ]; do

        # split header and body
        awk '{if(\$1~/^#/) print \$0 > "header"; else print \$0 > "body"}' appo.vcf
        nlines=`wc -l body | awk '{print int(\$1/"'${task.cpus}'" + 1)}'`
        split -l \$nlines body chunks

        for file in `ls chunks*`; do
            awk '{print \$0}' header \$file > tmp
            mv tmp \$file
        done

        export CIVICPY_CACHE_FILE="/home/.civicpy/cache.pkl"
        export CIVICPY_CACHE_TIMEOUT_DAYS=${params.civic_cache_timeout_days}        
        

        for file in `ls chunks*`; do
            civicpy annotate-vcf --input-vcf  \$file --output-vcf out.\$file --reference ${params.build_alt_name}  --include-status accepted --include-status submitted > log.\$file 2>&1 &
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
        rm ${vcf}
        bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
        tabix -p vcf ${meta.patient}.vcf.gz
        check=`wc -l ${meta.patient}.vcf | awk '{print \$1-1}'`
    done
    """

}


// alpha_missense
process add_alpha_missense{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 3.GB * task.attempt }
    tag "alphamissense"
    container "docker://ubuntu:20.04"

    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf)
    output:
        tuple val(meta), val(chunk_index), file("${chunk_index}.missense.maf"), file(vcf)
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
        }END{for(i=0;i<count;i++) printf "%s\\t%s\\t%s\\n", line[i], am_pathogenicity[i], am_class[i]}' ${maf} tmp > ${chunk_index}.missense.maf
        rm tmp
    """
}

// filter maf file
process filter_maf{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    tag "filtermaf"
    container "docker://yinxiu/gatk:latest"
    
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf)
    output:
        tuple val(meta), file("${chunk_index}.maf"), file("filtered*.pass.tsv"), file("filtered*.nopass.tsv"), file(vcf)
    script:
    """
        touch filtered.${chunk_index}.maf.pass.tsv filtered.${chunk_index}.maf.nopass.tsv

        filter_variants.py -m ${maf} \
         -o ${chunk_index}.maf \
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

process merge_chunks{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 2.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"
    container "docker://ubuntu:20.04"
    tag "merge_chunk"

    input:
        tuple val(meta), path(mafs), path(pass), path(nopass), path(vcfs)
    output:
        tuple val(meta), file("${meta.patient}.maf"), file("${meta.patient}.pass.tsv"), file("${meta.patient}.nopass.tsv"), file("${meta.patient}.vcf"), emit: all
        tuple val(meta), path("${meta.patient}.vcf"), emit: vcf
    script:
    """
        sorted_mafs=`echo ${mafs} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_vcfs=`echo ${vcfs} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_pass=`echo ${pass} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_nopass=`echo ${nopass} | tr ' ' '\\n' | sort | tr '\\n' ' '`

        awk 'BEGIN{fn=0}{if(FNR==1) fn++; if(fn==1){print \$0} else if(!(\$0~/^#/) && (\$1!="Hugo_Symbol")){print \$0}}' \${sorted_mafs} > ${meta.patient}.maf
        awk 'BEGIN{fn=0}{if(FNR==1) {fn++; if(fn>1) getline}; print \$0}' \${sorted_nopass} > ${meta.patient}.nopass.tsv
        awk 'BEGIN{fn=0}{if(FNR==1) {fn++; if(fn>1) getline}; print \$0}' \${sorted_pass} > ${meta.patient}.pass.tsv
        awk 'BEGIN{fn=0}{if(FNR==1) fn++; if(fn==1){print \$0} else if(!(\$0~/^#/)){print \$0}}' \${sorted_vcfs} > ${meta.patient}.vcf
    """
}

