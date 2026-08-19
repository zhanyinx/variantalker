process merge_funcotator_chunks {

    input:
        tuple val(meta), val(chunk_index), path(mafs), path(vcfs), file(vcf_full), file(maf_skip)

    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.maf"), file("${meta.patient}.vcf"), file(vcf_full), file(maf_skip)

    script:
    """
    first=1
    for vcf in \$(echo ${vcfs} | tr ' ' '\\n' | sort -V); do
        if [ \$first -eq 1 ]; then
            cat \$vcf > ${meta.patient}.vcf
            first=0
        else
            grep -v '^#' \$vcf >> ${meta.patient}.vcf
        fi
    done

    first=1
    for maf in \$(echo ${mafs} | tr ' ' '\\n' | sort -V); do
        if [ \$first -eq 1 ]; then
            cat \$maf > ${meta.patient}.maf
            first=0
        else
            grep -v '^#' \$maf | tail -n +2 >> ${meta.patient}.maf
        fi
    done
    """
}