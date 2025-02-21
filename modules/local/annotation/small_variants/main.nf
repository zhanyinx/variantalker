// List of functions and processes used to annotate snp and indel from whole exome sequences

process split_chunks{
    input:
        tuple val(meta), file(vcf)
    output:
        tuple val(meta), file("${meta.patient}_chunk_*.vcf.gz")
    script:
    """
        zcat ${vcf} > appo.vcf
        nfile=`awk 'BEGIN{
                    counts = 0; 
                    nfile=0
                }{
                    if(\$0~/^#/){
                        print \$0 > "header";
                    } else {
                        if(\$1==chr && \$2==pos && \$7 == "PASS"){
                            counts++; 
                            if(counts>nfile) nfile=counts
                        }else{
                            counts=0
                        }
                        if(\$7=="PASS"){  
                            print \$0 > "body_"counts".vcf"; 
                            pos=\$2; 
                            chr=\$1
                        }
                    }
                }END{print nfile}' appo.vcf`
                
        for i in `seq 0 \$nfile`; do
            split -l ${params.chunk_size} body_\$i.vcf ${meta.patient}_chunk_\$i
        done

        for i in `ls ${meta.patient}_chunk_*`; do
            cat header \$i > \$i.vcf
            bgzip -c \$i.vcf > \$i.vcf.gz
        done
    """
}

process run_funcotator{
    input:
        tuple val(meta), val(chunk_index), file(vcf)
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.maf"), file("${meta.patient}.vcf")
    script:
    """
    tabix -p vcf ${vcf}

    trascript_params=""
    if [ -f ${params.transcript_list} ]; then
        trascript_params=" --transcript-list ${params.transcript_list}"
    fi

    if [ "${meta.sample_type}" == "germline" ]; then
        extra_params=" --annotation-default Matched_Norm_Sample_Barcode:${meta.patient} --data-sources-path ${params.funcotator_germline_db}"
    else
        normal="\$(zcat ${vcf} | grep 'normal_sample'  | cut -d'=' -f2)"
        tumor="\$(zcat ${vcf} | grep 'tumor_sample'  | cut -d'=' -f2)"
        extra_params=" --annotation-default Matched_Norm_Sample_Barcode:\${normal} --annotation-default Tumor_Sample_Barcode:\${tumor} --annotation-default Tumor_type:${meta.tumor_tissue} --data-sources-path ${params.funcotator_somatic_db}"
    fi

    gatk Funcotator \
        -L ${params.funcotator_target} \
        -R ${params.fasta} \
        -V ${vcf} \
        -O ${meta.patient}.maf \
        --remove-filtered-variants true \
        --output-file-format MAF \
        --ref-version ${params.build} \
        --transcript-selection-mode ${params.transcript_selection} \
        --splice-site-window-size ${params.splice_site_window_size} \
        --interval-padding ${params.target_padding} \
        \$trascript_params \$extra_params

    zcat ${vcf} > ${meta.patient}.vcf

    # check completeness of maf file
    nrow_maf=\$(wc -l ${meta.patient}.maf | cut -f1 -d' ')
    nrow_vcf=\$(wc -l ${meta.patient}.vcf | cut -f1 -d' ')
    diff=\$((\$nrow_maf-\$nrow_vcf))
    if [ \$diff -lt 5 ]; then
        echo "Funcotator maf file is incomplete!"
        exit
    fi
    """
}

process add_guidelines_escat{
    cpus 1
    memory "1 G"
    container "docker://yinxiu/gatk:latest"
    errorStrategy 'retry'
    maxRetries = 3
    tag "add_guidelines_escat"
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(guidelines), file(grl_p), file(config)
    output:
        tuple val(meta), val(chunk_index), file("${chunk_index}.guidelines.maf"), file("${chunk_index}.vcf")
    script:
    """
    if [ "${meta.sample_type}" == "somatic" ]; then
        extraopts=" -t ${meta.tumor_tissue}"
    else
        extraopts=" --germline"
    fi
    add_guidelines_and_escat_to_maf.py -m ${maf} \
        -c ${guidelines} \
        -cc ${config} \
        -o ${chunk_index}.guidelines.maf \
        -p ${params.projectid} \
        --escat ${params.escat_db} \
        -d ${params.date} \$extraopts

    cp ${vcf} ${chunk_index}.vcf
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
    memory { 1.GB * task.attempt }
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

