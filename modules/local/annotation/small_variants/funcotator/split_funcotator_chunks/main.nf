process split_funcotator_chunks {

    input:
        tuple val(meta), val(chunk_index), file(vcf), file(vcf_full), file(maf_skip)
        val funcotator_n_chunks

    output:
        tuple val(meta), val(chunk_index), path("${meta.patient}_n_parts.txt"), path("${meta.patient}_part_*.vcf"), file(vcf_full), file(maf_skip)

    script:
    """
    split_threshold=50

    zcat -f ${vcf} | grep '^#' > header.vcf
    zcat -f ${vcf} | grep -v '^#' > body.vcf

    nvars=\$(wc -l body.vcf | cut -d' ' -f1)

    if [ \$nvars -le \$split_threshold ]; then

        n_parts=1

        cat header.vcf body.vcf > ${meta.patient}_part_0000.vcf

    else

        lines_per_chunk=\$(( (\$nvars + ${funcotator_n_chunks} - 1) / ${funcotator_n_chunks} ))

        split \
            -l \$lines_per_chunk \
            -d -a 4 \
            body.vcf \
            part_

        for f in part_*; do

            idx=\$(basename \$f | sed 's/part_//')

            cat header.vcf \$f > ${meta.patient}_part_\${idx}.vcf

        done

        n_parts=\$(ls ${meta.patient}_part_*.vcf | wc -l)

    fi

    echo \$n_parts > ${meta.patient}_n_parts.txt
    """
}