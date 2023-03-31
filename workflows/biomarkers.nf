

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.biomarkers.dna.input_mafs || params.biomarkers.dna.input_mafs.isEmpty()) {
    params.annotation_output = "${params.output}/*/annotation/somatic/*/*maf"
}else{
    params.annotation_output = params.biomarkers.dna.input_mafs
}

params.date = new java.util.Date().format('yyMMdd')

include {extract_tpm; calculate_tmb_signature; ascat_calling} from '../modules/extract_biomarkers.nf'
include {report_raw} from '../modules/report.nf'

workflow BIOMARKERS {

    // tmb and mutational signatures
    maf_files_ch = channel.fromPath(params.annotation_output, checkIfExists: true)         
    calculate_tmb_signature(maf_files_ch)
    // report_raw(calculate_tmb_signature.out[0], calculate_tmb_signature.out[1])

    // trascript per million from RNAseq
    if (params.biomarkers.rna.input_counts){
        if (params.pipeline.toUpperCase() == "DRAGEN"){
            rna_input_ch = channel.fromPath(params.biomarkers.rna.input_counts)
            extract_tpm(rna_input_ch)
        }
    }
}

workflow CLONAL_TMB{
    // clonal tmb
    if (params.pipeline.toUpperCase() == "DRAGEN"){
        if (params.biomarkers.ascat.input){
            sample=Channel.from(file(params.biomarkers.ascat.input))
            if (params.biomarkers.ascat.config && !params.biomarkers.ascat.config.isEmpty()){
                config=Channel.from(file(params.biomarkers.ascat.config))
            }else{
                config=Channel.from(file("$projectDir/resources/configs/ascat.sarek.config"))
            }
            ascat_calling(sample, config)
        }
    }   
}