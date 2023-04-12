

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.date = new java.util.Date().format('yyMMdd')

include {extract_tpm; calculate_tmb_signature; ascat_calling} from '../modules/extract_biomarkers.nf'
include {report_raw} from '../modules/report.nf'

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
    // clonal tmb
    if (params.pipeline.toUpperCase() == "DRAGEN"){
        sample=Channel.from(file(params.input))
        if (params.biomarkers_ascat_config && !params.biomarkers_ascat_config.isEmpty()){
            config=Channel.from(file(params.biomarkers_ascat_config))
        }else{
            config=Channel.from(file("$projectDir/resources/configs/ascat.sarek.config"))
        }
        ascat_calling(sample, config)
    }   
}