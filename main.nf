/*
Nextflow pipeline for whole exome sequencing clinical annotation (snp, indel, cnv) and QC report generation
Author(s): Yinxiu Zhan
Contact: yinxiu.zhan@ieo.it
This software is distributed without any guarantee under the terms of the GNU General
Public License, either Version 2, June 1991 or Version 3, June 2007.

*/

nextflow.enable.dsl=2

params.date = new java.util.Date().format('yyMMdd')
// Basic configuration 
if (params.tumoronly){
    params.cancervar.input_type = "VCF"
}

// Modified version of cancervar and intervar; if defined in your config, it will be overwritten
params.cancervar_folder = "$projectDir/resources/CancerVar"
params.intervar_folder = "$projectDir/resources/InterVar"
params.cancervar_init = "$projectDir/resources/configs/config.init.CancerVar"
params.intervar_init = "$projectDir/resources/configs/config.init.intervar"


include {extractInfo; fixvcf; rename_somatic_vcf; somatic_annotate_snp_indel; filter_variants; normalise_rename_germline_vcf; germline_annotate_snp_indel; germline_renovo_annotation} from './modules/snp_indel_annotation.nf'
include {extractBam; extractInfos; access; cnvkit; cnvkit_call; annotate_cnv} from './modules/cnv_annotation.nf'


process help {
    label "help"

    script:
    """
    echo "Usage: nextflow run main.nf -c nextflow.config"
    echo ""
    echo "This pipeline performs whole exome sequencing clinical annotation (snp, indel, cnv)."
    echo ""
    echo "Options:"
    echo "  --somatic.input_snp_indel <input_snp_indel>  A tab-separated file with two columns: tumor_type and path to VCF file containing somatic SNPs/Indels." 
    echo "  --germline.input_snp_indel <input_snp_indel> The input germline SNP/INDEL VCF file."
    echo "  --somatic.input_cnv <input_cnv>              The input somatic CNV file."
    echo "  --tumoronly                                  Flag to indicate that the input is tumor-only (default: false)."
    echo "  --pipeline <pipeline>                        The pipeline used to generate the input data (options: 'Sarek', 'DRAGEN' or 'Iontorrent' (only SNP/INDEL), default: 'Sarek')."
    echo ""
    """
}

workflow{

    // Workflow for snp and indel variant annotation
    if (params.somatic.input_snp_indel){
        somatic_snp_indel_input_ch = extractInfo(file(params.somatic.input_snp_indel))
        fixvcf(somatic_snp_indel_input_ch)
        rename_somatic_vcf(fixvcf.out)
        somatic_annotate_snp_indel(rename_somatic_vcf.out)
    }

    if (params.germline.input_snp_indel){
        germline_snp_indel_input_ch = channel.fromPath(params.germline.input_snp_indel)
        if (params.pipeline.toUpperCase() == "SAREK") {
            filter_variants(germline_snp_indel_input_ch)
            normalise_rename_germline_vcf(filter_variants.out)
            germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
            germline_renovo_annotation(germline_annotate_snp_indel.out)
        }
        else{
            normalise_rename_germline_vcf(germline_snp_indel_input_ch)
            germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
            germline_renovo_annotation(germline_annotate_snp_indel.out)
        }
    }


    // workflow for somatic cnv annotation
    if (params.somatic.input_cnv){
        if (params.pipeline.toUpperCase() == "SAREK") {
            input_sarek_ch = extractBam(file(params.somatic.input_cnv))
            (genderMap, statusMap, input_sarek_ch) = extractInfos(input_sarek_ch)
            input_sarek_ch.branch{
                tumor: statusMap[it[0], it[1]] == 1
                normal: statusMap[it[0], it[1]] == 0
            }.set{result}
            access()
            cnvkit(result.tumor, result.normal, access.out)
            cnvkit_call(cnvkit.out)
            annotate_cnv(cnvkit_call.out)
        }
        if (params.pipeline.toUpperCase() == "DRAGEN"){
            somatic_cnv_input_ch = channel.fromPath(params.somatic.input_cnv)
            annotate_cnv(somatic_cnv_input_ch)
        }
    }

}

