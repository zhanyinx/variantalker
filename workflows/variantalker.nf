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

if(params.build == "hg38"){
    params.build_alt_name = "GRCh38"
}else if(params.build == "hg19"){
    params.build_alt_name = "GRCh37"
}

if (!params.cancervar_evidence_file || params.cancervar_evidence_file.isEmpty()) {
    params.cancervar_evidence_file = "None"
}

if (!params.intervar_evidence_file || params.intervar_evidence_file.isEmpty()) {
    params.intervar_evidence_file = "None"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {filter_maf; add_civic; add_alpha_missense; fixvcf} from '../modules/local/annotation/small_variants/main.nf'
include {filter_maf as filter_maf_germline} from '../modules/local/annotation/small_variants/main.nf'
include {add_alpha_missense as add_alpha_missense_germline} from '../modules/local/annotation/small_variants/main.nf'

include {somatic_annotate_snp_indel} from '../modules/local/annotation/small_variants/somatic/main.nf'
include {filter_variants; normalise_rename_germline_vcf; germline_annotate_snp_indel; germline_renovo_annotation;} from '../modules/local/annotation/small_variants/germline/main.nf'
include {cnvkit_call; annotate_cnv} from '../modules/local/annotation/cnv/main.nf'

// extract channels from input annotation sample sheet 
def extract_csv(csv_file, sample_type) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
            if (numberOfLinesInSampleSheet == 1){
                def requiredColumns = ["patient", 'tumor_tissue', 'sample_file', 'sample_type']
                def headerColumns = line
                if (!requiredColumns.every { headerColumns.contains(it) }) {
                    log.error "Header missing or CSV file does not contain all of the required columns in the header: ${requiredColumns}"
                    System.exit(1)
                }
            }
        }
        
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Provided SampleSheet has less than two lines. Provide a samplesheet with header and at least a sample."
            System.exit(1)
        }
    }

    Channel.from(csv_file)
        .splitCsv(header: true)
        .map{ row ->
            if (row.sample_type == sample_type) {
                row
            } else {
                null
            }
        }
        .filter { row -> row != null }
        .map { row ->
            def meta = [:]

            meta.tumor_tissue = row.tumor_tissue
            meta.patient = row.patient
            meta.sample_type = row.sample_type

            if(sample_type == "somatic"){
                return [ meta, row.sample_file ]
            }else{
                return [ meta, row.sample_file ]
            }
        }
}

workflow VARIANTALKER{

    ch_somatic = extract_csv(file(params.input), "somatic")
    ch_germline = extract_csv(file(params.input), "germline")
    ch_cnv = extract_csv(file(params.input), "cnv")

    // ********** Workflow for snp and indel variant annotation **********
    // somatic 
    fixvcf(ch_somatic)
    add_civic(fixvcf.out)
    somatic_annotate_snp_indel(add_civic.out)
    add_alpha_missense(somatic_annotate_snp_indel.out)
    filter_maf(add_alpha_missense.out)

    // germline
    if (params.pipeline.toUpperCase() == "SAREK") {
        filter_variants(ch_germline)
        normalise_rename_germline_vcf(filter_variants.out)
        germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
        germline_renovo_annotation(germline_annotate_snp_indel.out)
        add_alpha_missense_germline(germline_renovo_annotation.out)
        filter_maf_germline(add_alpha_missense_germline.out)
    }
    else{
        normalise_rename_germline_vcf(ch_germline)
        germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
        germline_renovo_annotation(germline_annotate_snp_indel.out)
        add_alpha_missense_germline(germline_renovo_annotation.out)
        filter_maf_germline(add_alpha_missense_germline.out)
    }

    // workflow for somatic cnv annotation
    if (params.pipeline.toUpperCase() == "SAREK") {
        cnvkit_call(ch_cnv)
        annotate_cnv(cnvkit_call.out)
    }
    if (params.pipeline.toUpperCase() == "DRAGEN"){
        annotate_cnv(ch_cnv)
    }

    // biomarkers

}
