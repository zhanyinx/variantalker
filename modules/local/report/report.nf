/*
List of functions and processes used to generate the report
*/

process reporter{
    cpus 1
    memory "5 G"
    publishDir "${params.outdir}/${params.date}/biomarkers/${patient}", mode: "copy"
    tag "report"
    container "docker://yinxiu/reporter:latest"

    input:
        tuple val(patient), val(maf), val(maf_germline), val(variant_signatures), val(msi),  val(clonal_tmb), val(tmb), val(rna), val(cnv), val(metrics), val(hrd)
    output:
        file("${patient}.report.html")
    script:
    """
    # copy reporter
    cp ${projectDir}/bin/report.Rmd report.Rmd
    today=`date  '+%Y%m%d'`
    sed -i 's/DATE/'\$today'/g' report.Rmd

    # get annotation filters
    annotation_filters=`echo "${params}" | sed 's/, /\\n/g' | grep filter | awk '{printf "%s|", \$0}' | sed "s/'//g"`

    # Generate the biomarkers report
    Rscript -e "rmarkdown::render('report.Rmd', 'html_document' , output_file= '${patient}.report.html', params=list( patientid = '${patient}', annotation_filters='\${annotation_filters}', metrics_file = '${metrics}', variants_somatic_file = '${maf}', variants_germline_file = '${maf_germline}', clonal_tmb_file = '${clonal_tmb}', variant_signatures_file = '${variant_signatures}',dragen_tmb_file = '${tmb}', dragen_msi_file = '${msi}', rna_biomarker = '${rna}',  cnv_file = '${cnv}', cnv_gene_file='${params.cnv_genes_keep}', hrd_file='${hrd}'))"
    """
}