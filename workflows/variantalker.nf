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

include {fixvcf; rename_somatic_vcf; somatic_annotate_snp_indel; filter_variants; normalise_rename_germline_vcf; germline_annotate_snp_indel; germline_renovo_annotation} from '../modules/local/annotation/main.nf'
include {cnvkit_call; annotate_cnv} from '../modules/local/cnv/main.nf'

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
            if(sample_type == "somatic"){
                [row.patient, row.tumor_tissue, row.sample_file]
            }else{
                [row.patient, row.sample_file]
            }
        }
}

workflow VARIANTALKER{

    ch_somatic = extract_csv(file(params.input), "somatic")
    ch_germline = extract_csv(file(params.input), "germline")
    ch_cnv = extract_csv(file(params.input), "cnv")

    // Workflow for snp and indel variant annotation
    fixvcf(ch_somatic)
    rename_somatic_vcf(fixvcf.out)
    somatic_annotate_snp_indel(rename_somatic_vcf.out)

    if (params.pipeline.toUpperCase() == "SAREK") {
        filter_variants(ch_germline)
        normalise_rename_germline_vcf(filter_variants.out)
        germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
        germline_renovo_annotation(germline_annotate_snp_indel.out)
    }
    else{
        normalise_rename_germline_vcf(ch_germline)
        germline_annotate_snp_indel(normalise_rename_germline_vcf.out)
        germline_renovo_annotation(germline_annotate_snp_indel.out)
    }

    // workflow for somatic cnv annotation
    if (params.pipeline.toUpperCase() == "SAREK") {
        cnvkit_call(ch_cnv)
        annotate_cnv(cnvkit_call.out)
    }
    if (params.pipeline.toUpperCase() == "DRAGEN"){
        annotate_cnv(ch_cnv)
    }
}
