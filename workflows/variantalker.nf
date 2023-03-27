/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.date = new java.util.Date().format('yyMMdd')

// Modified version of cancervar and intervar; if defined in your config, it will be overwritten
params.cancervar_folder = "$projectDir/resources/CancerVar"
params.intervar_folder = "$projectDir/resources/InterVar"
params.cancervar_init = "$projectDir/resources/configs/config.init.CancerVar"
params.intervar_init = "$projectDir/resources/configs/config.init.intervar"
params.cancervar_db = "$projectDir/resources/CancerVar/cancervardb"
params.intervar_db = "$projectDir/resources/InterVar/intervardb"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {extractInfo; fixvcf; rename_somatic_vcf; somatic_annotate_snp_indel; filter_variants; normalise_rename_germline_vcf; germline_annotate_snp_indel; germline_renovo_annotation} from '../modules/snp_indel_annotation.nf'
include {extractBam; extractInfos; access; cnvkit; cnvkit_call; annotate_cnv} from '../modules/cnv_annotation.nf'



workflow VARIANTALKER{

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

