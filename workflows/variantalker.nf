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

// annotation
include {split_chunks; run_funcotator; add_guidelines_escat; add_alpha_missense; filter_maf; merge_chunks} from '../modules/local/annotation/small_variants/main.nf'
include {standardize_somatic_vcf; add_somatic_civic; run_somatic_cancervar; } from '../modules/local/annotation/small_variants/somatic/main.nf'
include {filter_germline_variants; standardize_germline_vcf; run_germline_intervar; run_germline_renovo} from '../modules/local/annotation/small_variants/germline/main.nf'
include {pharmgkb} from '../modules/local/pharmgkb/main.nf'
include {cnvkit_call; annotate_cnv} from '../modules/local/annotation/cnv/main.nf'

//biomarker
include {extract_tpm; calculate_tmb_signature} from '../modules/local/biomarkers/main.nf'

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

    /* Standardize vcf */
    standardize_somatic_vcf(ch_somatic) 
    if (params.pipeline.toUpperCase() == "SAREK") {
        filter_germline_variants(ch_germline)
        standardize_germline_vcf(filter_germline_variants.out)
    }
    else{
        standardize_germline_vcf(ch_germline)
    }
    vcf_ch = standardize_somatic_vcf.out.mix(standardize_germline_vcf.out)

    // split into chunks
    split_chunks(vcf_ch)
    chunks = split_chunks.out.map{ it->
        def meta = it[0]
        def vcf_files = it[1] instanceof List ? it[1] : [it[1]]
        return vcf_files.collect { vcf_file ->                    
            return [ meta, vcf_file.simpleName, vcf_file]
        }
    }.flatMap { it }

    // Annotations
    add_somatic_civic(chunks.filter{it -> it[0].sample_type == 'somatic'}) // CIViC
    run_somatic_cancervar(add_somatic_civic.out) // cancervar annotation on somatic

    merged_ch = add_somatic_civic.out.mix(chunks.filter{it -> it[0].sample_type == 'germline'})
    run_funcotator(merged_ch) // funcotator annotation
    run_germline_intervar(chunks.filter{it -> it[0].sample_type == 'germline'}) // intervar annotation on germline
    combined = run_funcotator.out.combine(run_somatic_cancervar.out.mix(run_germline_intervar.out), by: [0,1])
    add_guidelines_escat(combined)
    run_germline_renovo(add_guidelines_escat.out.filter(it -> it[0].sample_type == 'germline'))
    merged_ch = add_guidelines_escat.out.filter(it -> it[0].sample_type == 'somatic').mix(run_germline_renovo.out)
    add_alpha_missense(merged_ch)
    filter_maf(add_alpha_missense.out)
    merge_chunks(filter_maf.out.groupTuple(by: 0))

    // pharmGKB only for hg38
    if (params.build == "hg38"){
        pharmgkb(merge_chunks.out.vcf)
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
