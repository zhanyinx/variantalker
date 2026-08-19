process standardize_vcf {
    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("${meta.patient}.vcf")

    script:
    """
    if [[ ${vcf} == *.gz ]]; then
        zcat ${vcf} > input.vcf
    else
        cp ${vcf} input.vcf
    fi

    awk '{if(\$1~/^#/) print \$0; else if(\$7=="PASS") print \$0}' input.vcf > pass.vcf

    bcftools norm -m-any --check-ref -w -f ${params.fasta} pass.vcf -Ou \
        | bcftools view -e 'ALT="*"' -Ov \
        | awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' > ${meta.patient}.vcf

    input_variants=\$(grep -vc '^#' input.vcf)
    output_variants=\$(grep -vc '^#' ${meta.patient}.vcf)

    echo "Sample: ${meta.patient}"
    echo "Input variants:  \${input_variants}"
    echo "Output variants: \${output_variants}"

    if [[ \$output_variants -eq 0 ]]; then
        echo "WARNING: standardization produced an empty VCF for sample ${meta.patient}; it will be reported with header-only outputs." >&2
    fi

    if ! grep -q '^#CHROM' ${meta.patient}.vcf; then
        echo "ERROR: missing #CHROM header in ${meta.patient}.vcf" >&2
        exit 1
    fi

    bcftools view -h ${meta.patient}.vcf > /dev/null
    """
}

process handle_empty_sample {

    publishDir {
        "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.custom_id ?: meta.patient}"
    }, mode: "copy"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), file("${meta.patient}.maf"), file("${meta.patient}.pass.tsv"), file("${meta.patient}.nopass.tsv"), file("${meta.patient}.vcf")

    script:
    """
        echo "WARNING: sample ${meta.patient} has no variants after standardization; emitting header-only outputs." >&2

        cp ${vcf} ${meta.patient}.vcf

        skip_civic_flag=""
        if [ "${params.skip_civic}" == "true" ]; then
            skip_civic_flag="--skip_civic"
        fi

        write_empty_outputs.py \
            --sample_type ${meta.sample_type} \
            --patient ${meta.patient} \
            \$skip_civic_flag
    """
}

process split_chunks {

    input:
        tuple val(meta), file(vcf)

    output:
        tuple val(meta), file("${meta.patient}_chunk_*.vcf")

    script:
    """
    balanced_chunker.py \
        -i ${vcf} \
        -p ${meta.patient} \
        -c ${params.chunk_size}

    for f in ${meta.patient}_chunk_*.vcf; do
        bcftools sort \$f -Ov -o tmp.vcf
        mv tmp.vcf \$f
    done
    """
}

process split_by_database{
    container "docker://yinxiu/gatk:latest"
    input: 
        tuple val(meta), val(chunk_index), file(vcf)
        val safe_funcotator_germline_db
        val safe_funcotator_somatic_db
    output:
        tuple val(meta), val(chunk_index),  file("*_noskip.vcf"), file("*_full.vcf"), file("*_skip.maf")

    script:
    """
    if [ "${meta.sample_type}" == "germline" ]; then
        funcotator_db="${safe_funcotator_germline_db}"
        annovar_config="${params.intervar_init}"
    else
        funcotator_db="${safe_funcotator_somatic_db}"
        annovar_config="${params.cancervar_init}"
    fi

    # Build PostgreSQL connection parameters if provided
    pg_params=""

    if [ -n "${params.pg_dsn}" ] && [ "${params.pg_dsn}" != "null" ]; then
        pg_params="--dsn '${params.pg_dsn}' --pg-database '${params.pg_database}'"
    else
        pg_params="--pg-host ${params.pg_host} \
                --pg-port ${params.pg_port} \
                --pg-user ${params.pg_user} \
                --pg-password ${params.pg_password} \
                --pg-database ${params.pg_database}"
    fi

    python ${projectDir}/bin/database/split_annotated_variants.py \
        --vcf ${vcf} \
        --config \${annovar_config} \
        --funcotator \${funcotator_db} \
        --technology ${params.technology} \
        --annovar_version ${params.annovar_version} \
        --sample_type ${meta.sample_type} \
        --genome_build ${params.build_alt_name} \
        --tissue_type ${meta.tumor_tissue} \
        \$pg_params

    newname=`basename ${vcf} | sed 's/.vcf/_full.vcf/g'`
    cp ${vcf} \${newname}
    """
}

process add_guidelines_escat_to_funcotator{
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(vcf_full), file(maf_skip), file(guidelines), file(grl_p), file(config)
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}_${chunk_index}.guidelines.maf"), file("${meta.patient}_${chunk_index}.vcf"), file(vcf_full), file(maf_skip)
    script:
    """
    if [ "${meta.sample_type}" == "somatic" ]; then
        extraopts=" -t ${meta.tumor_tissue}"
    else
        extraopts=" --germline"
    fi

    if [ "${params.skip_civic}" == "true" ]; then
            skip_civic="--skip_civic"
        else
            skip_civic=""
    fi

    add_guidelines_escat_to_funcotator.py -m ${maf} \
        -c ${guidelines} \
        -cc ${config} \
        -o ${meta.patient}_${chunk_index}.guidelines.maf \
        --escat ${params.escat_db} \$extraopts \$skip_civic

    cp ${vcf} ${meta.patient}_${chunk_index}.vcf
    """
}

process add_alpha_missense{
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(vcf_full), file(maf_skip)
        val safe_alpha_missense
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}_${chunk_index}.missense.maf"), file(vcf), file(vcf_full), file(maf_skip)
    script:
    """
        if ! [ -f ${safe_alpha_missense} ]; then
            echo "${safe_alpha_missense} alpha Missense file does not exist." > .command.err
            exit 125
        fi

        # Handle both compressed and uncompressed files
        if [[ ${safe_alpha_missense} == *.gz ]]; then
            zcat ${safe_alpha_missense} > tmp
        else
            cat ${safe_alpha_missense} > tmp
        fi
        awk -F'\\t' 'BEGIN{fn=0; count=0}{
            if(FNR==1) fn++
            if(fn==1){
                if(\$1~/^#/){
                    print \$0
                }else if(\$1=="Hugo_Symbol"){
                    print "# Alpha Missense file ""'"${safe_alpha_missense}"'"
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
        }END{for(i=0;i<count;i++) printf "%s\\t%s\\t%s\\n", line[i], am_pathogenicity[i], am_class[i]}' ${maf} tmp > ${meta.patient}_${chunk_index}.missense.maf
        rm tmp
    """
}

process insert_annotations {

    input:
        tuple val(sample_type), val(tumor_tissue), path(maf)
        val safe_funcotator_germline_db
        val safe_funcotator_somatic_db

    when:
        params.pg_database &&
        (
            params.pg_dsn ||
            (
                params.pg_host &&
                params.pg_port &&
                params.pg_user &&
                params.pg_password
            )
        )

    script:
    """
    if [ "${sample_type}" == "germline" ]; then
        funcotator_db="${safe_funcotator_germline_db}"
        annovar_config="${params.intervar_init}"
    else
        funcotator_db="${safe_funcotator_somatic_db}"
        annovar_config="${params.cancervar_init}"
    fi

    awk 'BEGIN{fn=0}{if(FNR==1) fn++; if(fn==1){if(\$1=="Hugo_Symbol") print \$0}if(!(\$0~/^#/) && \$1!="Hugo_Symbol"){print \$0}}' ${maf} > tmp

    pg_params=""

    if [ -n "${params.pg_dsn}" ] && [ "${params.pg_dsn}" != "null" ]; then
        pg_params="--dsn '${params.pg_dsn}' --pg-database '${params.pg_database}'"
    else
        pg_params="--pg-host ${params.pg_host} \
                   --pg-port ${params.pg_port} \
                   --pg-user ${params.pg_user} \
                   --pg-password ${params.pg_password} \
                   --pg-database ${params.pg_database}"
    fi

    python ${projectDir}/bin/database/insert_annotations.py -m tmp \
        --config \${annovar_config} \
        --funcotator \${funcotator_db} \
        --technology ${params.technology} \
        --annovar_version ${params.annovar_version} \
        --sample_type ${sample_type} \
        --genome_build ${params.build_alt_name} \
        --tissue_type ${tumor_tissue} \
        \$pg_params
    """
}

process generate_clinical_summary{
    input:
        tuple val(meta), val(chunk_index), path(maf), file(vcf), file(vcf_full), file(maf_skip)
    output:
        tuple val(meta), val(chunk_index), file("${maf}_clinsumm.maf"), file(vcf), file(vcf_full), file(maf_skip)
    script:
    """
        generate_clinical_summary.py -m ${maf} -o ${maf}_clinsumm.maf
    """
}

process filter_maf{
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(vcf_full), file(maf_skip)
        val safe_funcotator_germline_db
        val safe_funcotator_somatic_db
    output:
        tuple val(meta), file("${chunk_index}.maf"), file("filtered*.pass.tsv"), file("filtered*.nopass.tsv"), file(vcf_full)
    script:
    """
        touch filtered.${chunk_index}.maf.pass.tsv filtered.${chunk_index}.maf.nopass.tsv
        if [ "${meta.sample_type}" == "germline" ]; then
            funcotator_db="${safe_funcotator_germline_db}"
            annovar_config="${params.intervar_init}"
        else
            funcotator_db="${safe_funcotator_somatic_db}"
            annovar_config="${params.cancervar_init}"
        fi

        if [ "${params.skip_civic}" == "true" ]; then
            skip_civic="--skip_civic"
        else
            skip_civic=""
        fi

        if [ "${params.skip_pathogenic_retention}" == "true" ]; then
            skip_pathogenic_retention="--skip_pathogenic"
        else
            skip_pathogenic_retention=""
        fi

        filter_variants.py -m ${maf} \
         -o ${chunk_index}.maf \
         --filter_intervar "${params.filter_intervar}" \
         --filter_cancervar "${params.filter_cancervar}" \
         --filter_renovo "${params.filter_renovo}" \
         --filter_escat "${params.filter_escat}" \
         --filter_clinvar "${params.filter_clinvar}" \
         --sample_type ${meta.sample_type} \
         --min_depth ${params.filter_min_depth} \
         --vaf_threshold_germline ${params.filter_vaf_threshold_germline} \
         --vaf_threshold ${params.filter_vaf_threshold} \
         --filter_genes_germline ${params.filter_genes_germline} \
         --filter_genes_somatic ${params.filter_genes_somatic} \
         --filter_variant_classification "${params.filter_var_classification}" \
         --filter_civic "${params.filter_civic_evidence_level}" \
         --funcotator \${funcotator_db} \
         --technology ${params.technology} \
         --genome_build ${params.build_alt_name} \
         --config \$annovar_config \
         --annovar_version ${params.annovar_version} \
         --projectid ${params.projectid} \
         --maf_skip ${maf_skip} \$skip_civic \$skip_pathogenic_retention

    """
}

process merge_chunks{
    
    publishDir {
        "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.custom_id ?: meta.patient}"
    }, mode: "copy"
    
    input:
        tuple val(meta), path(mafs), path(pass), path(nopass), file(vcf_full)
    output:
        tuple val(meta), file("${meta.patient}.maf"), file("${meta.patient}.pass.tsv"), file("${meta.patient}.nopass.tsv"), file("${meta.patient}.vcf"), emit: all
        tuple val(meta), path("${meta.patient}.vcf"), emit: vcf
    script:
    """
        sorted_mafs=`echo ${mafs} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_vcfs=`echo ${vcf_full} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_pass=`echo ${pass} | tr ' ' '\\n' | sort | tr '\\n' ' '`
        sorted_nopass=`echo ${nopass} | tr ' ' '\\n' | sort | tr '\\n' ' '`

        clean_final_maf.py -m \${sorted_mafs} -o ${meta.patient}.maf
        awk 'BEGIN{fn=0}{if(FNR==1) {fn++; if(fn>1) getline}; print \$0}' \${sorted_nopass} > tmp.nopass.tsv
        awk 'BEGIN{fn=0}{if(FNR==1) {fn++; if(fn>1) getline}; print \$0}' \${sorted_pass} > tmp.pass.tsv
        awk 'BEGIN{fn=0}{if(FNR==1) fn++; if(fn==1){print \$0} else if(!(\$0~/^#/)){print \$0}}' \${sorted_vcfs} > tmp.vcf
        order_vcf_pass_nopass.py -v tmp.vcf \
        -p tmp.pass.tsv \
        -n tmp.nopass.tsv \
        -vo ${meta.patient}.vcf \
        -po ${meta.patient}.pass.tsv \
        -no ${meta.patient}.nopass.tsv
    """
}