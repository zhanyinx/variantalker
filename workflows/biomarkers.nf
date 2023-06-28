

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.date = new java.util.Date().format('yyMMdd')
params.ascat_genome = params.build

include {extract_tpm; calculate_tmb_signature; clonal_tmb} from '../modules/local/biomarkers/main.nf'
include {ASCAT; generate_ascat_loci; generate_ascat_alleles; generate_ascat_rt; generate_ascat_gc} from '../modules/local/ascat/main.nf'
include {SAMTOOLS_CONVERT as BAM_TO_CRAM} from '../modules/nf-core/samtools/convert/main.nf'
include {generate_pyclone; pyclone} from '../modules/local/pyclone/main.nf'
//include {report_raw} from '../modules/report.nf'



workflow BIOMARKERS {

    ch_rna = extract_csv(file(params.input), "rna")
    ch_dna = extract_csv(file(params.input), "dna")
    // tmb and mutational signatures
    calculate_tmb_signature(ch_dna)

    // trascript per million from RNAseq
    if (params.pipeline.toUpperCase() == "DRAGEN"){
        extract_tpm(ch_rna)
    }
}

workflow CLONAL_TMB{
    if (params.pipeline.toUpperCase() == "DRAGEN"){

        // create channel from input genome files
        ch_loci = Channel.from(file(params.genomes[params.build].ascat_loci))
        ch_alleles = Channel.from(file(params.genomes[params.build].ascat_alleles))
        ch_rt = Channel.from(file(params.genomes[params.build].ascat_loci_rt))
        ch_gc = Channel.from(file(params.genomes[params.build].ascat_loci_gc))

        // extract channels from input file
        input_sample = extract_csv_clonal_tmb(file(params.input))
        input_variant_calling_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
        }

        // convert to cram
        BAM_TO_CRAM(input_variant_calling_convert.bam, params.fasta, "${params.fasta}.fai")
        cram_variant_calling = Channel.empty().mix(BAM_TO_CRAM.out.alignment_index, input_variant_calling_convert.cram)
        cram_variant_calling_status = cram_variant_calling.branch{
                normal: it[0].status == 0
                tumor:  it[0].status == 1
            }
        
        // generate ascat db
        generate_ascat_gc(ch_gc)
        generate_ascat_rt(ch_rt)
        generate_ascat_loci(ch_loci)
        generate_ascat_alleles(ch_alleles)

        // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross.cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample
                meta.maf        = tumor[1].maf
                meta.cellularity = tumor[1].cellularity

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
            }

        ASCAT(cram_variant_calling_pair, generate_ascat_alleles.out.first(), generate_ascat_loci.out.first(), Channel.fromPath(params.target).first(), Channel.fromPath(params.fasta).first(), generate_ascat_gc.out.first(), generate_ascat_rt.out.first())
        generate_pyclone(ASCAT.out.cram.groupTuple())
        pyclone(generate_pyclone.out)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// extract channels from input biomarkers sample sheet 
def extract_csv(csv_file, sample_type) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
            if (numberOfLinesInSampleSheet == 1){
                def requiredColumns = ["patient", 'sample_file', 'sample_type']
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
            if(sample_type == "rna"){
                [row.patient, row.sample_file]
            }else if(sample_type == "dna"){
                [row.patient, row.sample_file]
            }
        }
}

def extract_csv_clonal_tmb(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, samplesheet_line_count = 0;
        while ((line = reader.readLine()) != null) {samplesheet_line_count++}
        if (samplesheet_line_count < 2) {
            error("Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines.")
        }
    }

    // Additional check of sample sheet:
    // 1. If params.step == "mapping", then each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            if (!(row.patient && row.sample)) {
                error("Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'.")
            }
            else if (row.patient.contains(" ") || row.sample.contains(" ")) {
                error("Invalid value in csv file. Values for 'patient' and 'sample' can not contain space.")
            }
            [ [ row.patient.toString(), row.sample.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 1) sample_count_tumor++
        else sample_count_normal++

        // Two checks for ensuring that the pipeline stops with a meaningful error message if
        // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
        // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
        if ((sample_count_normal == sample_count_all) && !params.build_only_index) { // In this case, the sample-sheet contains no tumor-samples
            def tools_tumor_asked = ['ascat']
            tools_tumor.each{ tool ->
                if (params.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
            }
        } else if ((sample_count_tumor == sample_count_all)) {  // In this case, the sample-sheet contains no normal/germline-samples
            def tools_requiring_normal_samples = ['ascat']
            def requested_tools_requiring_normal_samples = []
            tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
                if (params.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
            }
            if (!requested_tools_requiring_normal_samples.isEmpty()) {
                error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
            }
        }


        // start from BAM
        if (row.lane && row.bam) {
            if (!row.bai) {
                error("BAM index (bai) should be provided.")
            }
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def bai         = row.bai ? file(row.bai,   checkIfExists: true) : []
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.sample}_${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            meta.num_lanes  = num_lanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'bam'
            meta.maf = row.maf
            meta.cellularity = row.cellularity

            meta.size       = 1 // default number of splitted fastq

            if(bam.getExtension() != 'bam'){
                error("A column with name 'bam' was specified, but it contains a file not ending on '.bam': " + bam.getName())
            }

            return [ meta, bam, bai ]

        // recalibration
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)

            meta.data_type  = 'cram'
            meta.maf = row.maf
            meta.cellularity = row.cellularity

            if(cram.getExtension() != 'cram'){
                error("A column with name 'cram' was specified, but it contains a file not ending on '.cram': " + cram.getName())
            }
            return [ meta, cram, crai ]

        }  else {
            error("Missing or unknown field in csv file header. Please check your samplesheet")
        }
    }
}

