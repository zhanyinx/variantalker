/*
List of functions and processes used to generate the report
*/

process reporter{
    cpus 1
    memory "5 G"
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}", mode: "copy"
    tag "report"

    input:
        tuple val(patient), val(maf), val(maf_germline), val(variant_signatures), val(msi),  val(clonal_tmb), val(tmb), val(rna), val(cnv)
    output:
        file("${patient}.report.html")
    script:
    """
    # copy reporter
    cp ${projectDir}/bin/report.Rmd report.Rmd
    today=`date  '+%Y%m%d'`
    sed -i 's/DATE/'\$today'/g' report.Rmd
    
    Rscript -e "rmarkdown::render('report.Rmd', 'html_document' , output_file= '${patient}.report.html', params=list( patientid = '${patient}', variants_somatic_file = '${maf}', variants_germline_file = '${maf_germline}', clonal_tmb_file = '${clonal_tmb}', variant_signatures_file = '${variant_signatures}',dragen_tmb_file = '${tmb}', dragen_msi_file = '${msi}', rna_biomarker = '${rna}',  cnv_file = '${cnv}'))"
    """
}