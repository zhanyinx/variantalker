/*
List of functions and processes used to generate the report
*/

process report_raw{
    cpus 2
    memory "5 G"
    publishDir "${params.output}/${params.date}/reports/${maf.simpleName}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${maf.simpleName}/reports/", mode: "copy"
    tag "report raw"

    input:
        file(maf)
        file(tmb)
    output:
        file("report.${maf.simpleName}.raw.csv")
    script:
    """
    report_raw.py -m ${maf} \
        -t ${tmb} \
         -o report.${maf.simpleName}.raw.csv \
    """
}