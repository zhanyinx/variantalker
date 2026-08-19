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

def safe_funcotator_somatic_db =
    params.funcotator_somatic_db == null ?
    "${params.database_dir}/${params.build}/funcotator_dataSources.v1.8.2024s" :
    params.funcotator_somatic_db
def safe_funcotator_germline_db =
    params.funcotator_germline_db == null ?
    "${params.database_dir}/${params.build}/funcotator_dataSources.v1.8.2024g" :
    params.funcotator_germline_db
def safe_annovar_db =
    params.annovar_db == null ?
    "${params.database_dir}/humandb" :
    params.annovar_db
def safe_annovar_software_folder =
    params.annovar_software_folder == null ?
    "${params.database_dir}/annovar_software" :
    params.annovar_software_folder
def safe_transcript_list =
    params.transcript_list == null ?
    "${projectDir}/resources/mane.transcript.gencode.v43.${params.build}.txt" :
    params.transcript_list
def safe_alpha_missense =
    params.alpha_missense == null ?
    "${params.database_dir}/${params.build}/alpha_missense/${params.build}/AlphaMissense_${params.build}.231106.tsv.gz" :
    params.alpha_missense

def safe_cancervar_evidence_file =
    (!params.cancervar_evidence_file || params.cancervar_evidence_file.isEmpty()) ?
    "None" :
    params.cancervar_evidence_file
def safe_intervar_evidence_file =
    (!params.intervar_evidence_file || params.intervar_evidence_file.isEmpty()) ?
    "None" :
    params.intervar_evidence_file

log.info "Funcotator somatic database chosen: ${safe_funcotator_somatic_db}"
log.info "Funcotator germline database chosen: ${safe_funcotator_germline_db}"
log.info "Annovar database chosen: ${safe_annovar_db}"
log.info "Annovar software folder chosen: ${safe_annovar_software_folder}"
log.info "Transcript list chosen: ${safe_transcript_list}"
log.info "CancerVar evidence file chosen: ${safe_cancervar_evidence_file}"
log.info "InterVar evidence file chosen: ${safe_intervar_evidence_file}"
log.info "Alpha missense file chosen: ${safe_alpha_missense}"
log.info "Funcotator ncores: ${params.funcotator_ncores}"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// annotation
include {
    standardize_vcf
    handle_empty_sample
    split_chunks
    split_by_database
    add_guidelines_escat_to_funcotator
    add_alpha_missense
    insert_annotations
    filter_maf
    generate_clinical_summary
    merge_chunks
} from '../modules/local/annotation/small_variants/main.nf'

include {
    fix_somatic_header
    add_somatic_civic
    run_somatic_cancervar
} from '../modules/local/annotation/small_variants/somatic/main.nf'

include {
    filter_germline_variants
    run_germline_intervar
    run_germline_renovo
} from '../modules/local/annotation/small_variants/germline/main.nf'

include {
    pharmgkb
    liftover19to38
} from '../modules/local/pharmgkb/main.nf'

include {
    cnvkit_call
    annotate_cnv
} from '../modules/local/annotation/cnv/main.nf'

include { FUNCOTATOR } from '../subworkflows/funcotator_parallel.nf'


// Return true if the VCF has at least one variant record (a non-'#' line).
// Reads lazily and stops at the first record so large VCFs aren't fully read.
def vcfHasVariants(vcf) {
    def has = false
    vcf.withReader { reader ->
        def line = reader.readLine()
        while (line != null) {
            if (!line.startsWith('#')) {
                has = true
                break
            }
            line = reader.readLine()
        }
    }
    return has
}


def extract_csv(csv_file, sample_types = []) {

    file(csv_file).withReader('UTF-8') { reader ->
        def line
        def numberOfLinesInSampleSheet = 0
        def headerColumns = []
        def customIds = [] as Set

        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++

            if (numberOfLinesInSampleSheet == 1) {
                def requiredColumns = ["patient", "tumor_tissue", "sample_file", "sample_type"]
                headerColumns = line.split(',') as List

                if (!requiredColumns.every { headerColumns.contains(it) }) {
                    log.error "Header missing or CSV file does not contain all required columns: ${requiredColumns}"
                    System.exit(1)
                }
            } else {
                def values = line.split(',', -1)
                def row = [headerColumns, values].transpose().collectEntries()

                if (row.custom_id) {
                    if (customIds.contains(row.custom_id)) {
                        log.error "Duplicated custom_id in samplesheet: ${row.custom_id}"
                        System.exit(1)
                    }
                    customIds.add(row.custom_id)
                }
            }
        }

        if (numberOfLinesInSampleSheet < 2) {
            log.error "Provided SampleSheet has less than two lines."
            System.exit(1)
        }
    }

    Channel.from(csv_file)
        .splitCsv(header: true)
        .filter { row -> sample_types.isEmpty() || row.sample_type in sample_types }
        .map { row ->
            def meta = [:]

            meta.patient = row.patient
            meta.sample_type = row.sample_type

            def supportedTumorTypes = [
                'Adrenal_Gland', 'Bile_Duct', 'Bladder', 'Blood', 'Bone',
                'Bone_Marrow', 'Brain', 'Breast', 'Cancer_all', 'Cervix',
                'Colorectal', 'Esophagus', 'Eye', 'Head_and_Neck',
                'Inflammatory', 'Intrahepatic', 'Kidney', 'Liver',
                'Lung', 'Lymph_Nodes', 'Nervous_System', 'Other',
                'Ovary', 'Pancreas', 'Pleura', 'Prostate', 'Skin',
                'Soft_Tissue', 'Stomach', 'Testis', 'Thymus',
                'Thyroid', 'Uterus'
            ] as Set

            if (row.sample_type == 'germline') {

                // Artificial tissue used by downstream logic and DB caching
                meta.tumor_tissue = 'Germline'

            } else if (row.sample_type == 'somatic') {

                def tissue = row.tumor_tissue?.trim()

                if (!tissue) {
                    log.error "Missing tumor_tissue for somatic sample '${row.patient}'"
                    System.exit(1)
                }

                if (!supportedTumorTypes.contains(tissue)) {
                    log.error "Unsupported tumor_tissue '${tissue}' for somatic sample '${row.patient}'"
                    System.exit(1)
                }

                meta.tumor_tissue = tissue

            }

            if (row.containsKey('custom_id') && row.custom_id) {
                meta.custom_id = row.custom_id
            }

            return [ meta, row.sample_file ]
        }
}

workflow VARIANTALKER {

    ch_variants = extract_csv(file(params.input), ["somatic", "germline"])
    ch_cnv = extract_csv(file(params.input), ["cnv"])

    def split = ch_variants.branch {
        somatic  : it[0].sample_type == "somatic"
        germline : it[0].sample_type == "germline"
    }

    if (params.pipeline.toUpperCase() == "SAREK") {

        filter_germline_variants(split.germline)
        ch_variants = split.somatic.mix(filter_germline_variants.out)

    } else {

        fix_somatic_header(split.somatic)
        ch_variants = fix_somatic_header.out.mix(split.germline)

    }

    standardize_vcf(ch_variants)

    // Samples with no variants after standardization (empty input, or nothing
    // PASS / surviving normalization) are reported with header-only outputs
    // instead of crashing the run. Everything else proceeds as usual.
    def std = standardize_vcf.out.branch { meta, vcf ->
        nonempty: vcfHasVariants(vcf)
        empty:    !vcfHasVariants(vcf)
    }

    handle_empty_sample(std.empty)

    vcf_ch = std.nonempty

    // split into chunks
    split_chunks(vcf_ch)

    chunks = split_chunks.out
        .map { it ->
            def meta = it[0]
            def vcf_files = it[1] instanceof List ? it[1] : [it[1]]

            vcf_files.collect { vcf_file ->
                def chunk_index = vcf_file.simpleName.replaceFirst("^${meta.patient}_", "")
                [ meta, chunk_index, vcf_file ]
            }
        }
        .flatMap()

    split_by_database(
        chunks,
        safe_funcotator_germline_db,
        safe_funcotator_somatic_db
    )

    // Split chunks into those that still need annotation (noskip.vcf has variants)
    // and those fully served by the cache (noskip.vcf is header-only). The latter
    // bypass the whole annotation chain and feed straight into clinical summary.
    split_db = split_by_database.out.branch { meta, chunk_index, noskip, full, skip ->
        annotate: vcfHasVariants(noskip)
        cached:   !vcfHasVariants(noskip)
    }

    // Cached chunks: shape matches add_alpha_missense.out -> (meta, chunk_index, maf, vcf, vcf_full, maf_skip)
    // maf = skip.maf (already fully annotated, canonical order); maf_skip = empty placeholder.
    cached_ch = split_db.cached.map { meta, chunk_index, noskip, full, skip ->
        [ meta, chunk_index, skip, file("${projectDir}/assets/EMPTY_FILE_vcf"), full, file("${projectDir}/assets/EMPTY_FILE_maf") ]
    }

    // Annotations
    add_somatic_civic(
        split_db.annotate.filter { it -> it[0].sample_type == 'somatic' }
    )

    run_somatic_cancervar(
        add_somatic_civic.out,
        safe_annovar_db,
        safe_annovar_software_folder,
        safe_cancervar_evidence_file
    )

    run_germline_intervar(
        split_db.annotate.filter { it -> it[0].sample_type == 'germline' },
        safe_annovar_db,
        safe_annovar_software_folder,
        safe_intervar_evidence_file
    )

    merged_ch = add_somatic_civic.out.mix(
        split_db.annotate.filter { it -> it[0].sample_type == 'germline' }
    )

    FUNCOTATOR(
        merged_ch,
        safe_funcotator_germline_db,
        safe_funcotator_somatic_db,
        safe_transcript_list,
        params.funcotator_ncores,
    )

    combined = FUNCOTATOR.out.combine(
        run_somatic_cancervar.out.mix(run_germline_intervar.out),
        by: [0,1]
    )

    add_guidelines_escat_to_funcotator(combined)

    run_germline_renovo(
        add_guidelines_escat_to_funcotator.out.filter { it -> it[0].sample_type == 'germline' },
        safe_annovar_db,
        safe_annovar_software_folder
    )

    merged_ch = add_guidelines_escat_to_funcotator.out
        .filter { it -> it[0].sample_type == 'somatic' }
        .mix(run_germline_renovo.out)

    add_alpha_missense(
        merged_ch,
        safe_alpha_missense
    )

    annotation_ch = add_alpha_missense.out
        .map {
            def meta = it[0]
            def maf = it[2]
            return [ [meta.sample_type, meta.tumor_tissue], maf ]
        }
        .groupTuple(by: 0)
        .map { key, maf ->
            def sample_type = key[0]
            def tumor_tissue = key[1]
            return [ sample_type, tumor_tissue, maf ]
        }

    insert_annotations(
        annotation_ch,
        safe_funcotator_germline_db,
        safe_funcotator_somatic_db
    )

    // Freshly-annotated chunks plus cache-only chunks converge here.
    generate_clinical_summary(add_alpha_missense.out.mix(cached_ch))

    filter_maf(
        generate_clinical_summary.out,
        safe_funcotator_germline_db,
        safe_funcotator_somatic_db
    )

    merge_chunks(
        filter_maf.out.groupTuple(by: 0).map { it ->
            def meta = it[0]
            def maf = it[1]
            def pass = it[2]
            def nopass = it[3]
            def vcf_full = it[4]
            return [ meta, maf, pass, nopass, vcf_full ]
        }
    )

    germline_vcf_ch = vcf_ch.filter { meta, vcf ->
        meta.sample_type == "germline"
    }

    // pharmGKB
    if (params.build == "hg38") {
        pharmgkb(germline_vcf_ch)
    } else if (params.build == "hg19") {
        if (params.pharmgkb_hg38_fasta && params.pharmgkb_hg38_chain) {
            liftover19to38(germline_vcf_ch)
            pharmgkb(liftover19to38.out.vcf)
        } else {
            println "Skipping liftover + pharmgkb because 'pharmgkb_hg38_fasta' or 'pharmgkb_hg38_chain' are missing."
        }
    } else {
        error "Invalid 'build' parameter! Must be 'hg19' or 'hg38'."
    }

    // workflow for somatic cnv annotation
    if (params.pipeline.toUpperCase() == "SAREK") {
        cnvkit_call(ch_cnv)
        annotate_cnv(cnvkit_call.out)
    } else if (params.pipeline.toUpperCase() == "DRAGEN") {
        annotate_cnv(ch_cnv)
    }
}