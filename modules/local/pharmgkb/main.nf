process liftover19to38 {
    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("*.hg38.norm.vcf"), emit: vcf
    script:
        """
            gatk LiftoverVcf \
            -I ${vcf} \
            -O ${vcf.baseName}.hg38.vcf \
            -CHAIN ${params.pharmgkb_hg38_chain} \
            -REJECT rejected.vcf.gz \
            -R ${params.pharmgkb_hg38_fasta} \
            --RECOVER_SWAPPED_REF_ALT true

            bcftools norm -f ${params.pharmgkb_hg38_fasta} -m -both ${vcf.baseName}.hg38.vcf -o ${vcf.baseName}.hg38.norm.vcf
        """
}

process pharmgkb {
    fair true
    maxRetries 3
    publishDir {
        "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.custom_id ?: meta.patient}"
    }, mode: "copy"

    when:
        !params.skip_pharmgkb

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*.html")

    script:
        def xmx = Math.floor(task.memory.toGiga() * 0.75) as int

        """
        set -euo pipefail

        if [[ "${vcf}" == *.gz ]]; then
            gunzip -c "${vcf}" > ${meta.patient}.vcf
        elif [[ "${vcf}" == "${meta.patient}.vcf" ]]; then
            :
        else
            cp "${vcf}" ${meta.patient}.vcf
        fi

        pharmcat_vcf_preprocessor -vcf ${meta.patient}.vcf

        java -Xmx${xmx}g -jar /pharmcat/pharmcat.jar \\
            -vcf ${meta.patient}.preprocessed.vcf.bgz
        """
}