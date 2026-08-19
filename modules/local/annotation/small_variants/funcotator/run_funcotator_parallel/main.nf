process run_funcotator_parallel {

    input:
        tuple val(group_key),
            val(meta),
            val(chunk_index),
            val(mini_index),
            file(vcf),
            file(vcf_full),
            file(maf_skip)        
        val safe_funcotator_germline_db
        val safe_funcotator_somatic_db
        val safe_transcript_list

    output:
    tuple val(group_key),
        val(meta),
        val(chunk_index),
        val(mini_index),
        file("${meta.patient}_${mini_index}.maf"),
        file("${meta.patient}_${mini_index}.vcf"),
        file(vcf_full),
        file(maf_skip)
        
    script:
    """
    bgzip -f ${vcf}
    tabix -f -p vcf ${vcf}.gz

    trascript_params=""
    if [ -f "${safe_transcript_list}" ]; then
        trascript_params=" --transcript-list ${safe_transcript_list}"
    fi

    if [ "${meta.sample_type}" == "germline" ]; then

        extra_params=" --annotation-default Matched_Norm_Sample_Barcode:${meta.patient} --data-sources-path ${safe_funcotator_germline_db}"

    else

        normal="\$(zcat -f ${vcf}.gz | grep -m1 'normal_sample' | cut -d'=' -f2)"
        tumor="\$(zcat -f ${vcf}.gz | grep -m1 'tumor_sample' | cut -d'=' -f2)"

        extra_params=" \
            --annotation-default Matched_Norm_Sample_Barcode:\${normal} \
            --annotation-default Tumor_Sample_Barcode:\${tumor} \
            --annotation-default Tumor_type:${meta.tumor_tissue} \
            --data-sources-path ${safe_funcotator_somatic_db}"
    fi

    gatk Funcotator \
        -L ${params.funcotator_target} \
        -R ${params.fasta} \
        -V ${vcf}.gz \
        -O ${meta.patient}_${mini_index}.maf \
        --remove-filtered-variants true \
        --output-file-format MAF \
        --ref-version ${params.build} \
        --transcript-selection-mode ${params.transcript_selection} \
        --splice-site-window-size ${params.splice_site_window_size} \
        --interval-padding ${params.target_padding} \
        \$trascript_params \
        \$extra_params

    zcat -f ${vcf}.gz > ${meta.patient}_${mini_index}.vcf

    nrow_maf=\$(wc -l ${meta.patient}_${mini_index}.maf | cut -f1 -d' ')
    nrow_vcf=\$(wc -l ${meta.patient}_${mini_index}.vcf | cut -f1 -d' ')

    diff=\$((\$nrow_maf-\$nrow_vcf))

    if [ \$diff -lt 5 ]; then
        echo "Funcotator maf file is incomplete!"
        exit 1
    fi
    """
}