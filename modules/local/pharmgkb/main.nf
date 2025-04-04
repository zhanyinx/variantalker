process pharmgkb {
    fair true
    cpus 1
    maxRetries = 3
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"
    container "docker://pgkb/pharmcat:2.15.5"

    tag "pharmgkb"    

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("*.html")
    script:
        """
            if [[ $vcf == *.gz ]]; then
                gunzip -c $vcf > ${vcf.baseName}.vcf
                java -jar /pharmcat/pharmcat.jar -vcf ${vcf.baseName}.vcf
            else
                java -jar /pharmcat/pharmcat.jar -vcf ${vcf}
            fi
        """
}