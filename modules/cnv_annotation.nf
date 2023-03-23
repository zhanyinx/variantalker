/*
List of functions and processes used to [call] and annotate cnv from whole exome sequences
[] = optional, performed if sarek is set to true

[1.] access: create access file for cnvkit
[2.] cnvkit + cnvkit_call: call cnv from whole exome sequences
3. annotate_cnv: annotate cnv using classifyCNV
*/

/* 
Extract bam file from tsvFile 
tsvFile = file containing list of bam files
format of tsvFile: *_no_duplicate

*/
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension."
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension."

            [idPatient, gender, status, idSample, bamFile, baiFile]
    }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}"
    return file(it)
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}"
    return it
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

process access{

    output:
        file("access${params.build}.bed")
    script:
        """
        cnvkit.py access ${params.fasta} -o "access${params.build}.bed"
        """
}

process cnvkit{
    cpus 5
    memory "5 G"

    input:
        tuple val(patientid), val(sampleid), file(tumorbam), file(tumorbai)
        tuple val(patientid1), val(sampleid1), file(normbam), file(normbai)
        file("access${params.build}.bed")

    output:
        file("results/${tumorbam.baseName}.cnr")
    script:
        """
        cnvkit.py batch ${tumorbam} \
        --normal ${normbam} \
        --targets ${params.target} \
        --annotate ${params.cnvkit.refFlat} \
        --fasta ${params.fasta} \
        --access access${params.build}.bed \
        --output-reference ${patientid}.cnn \
        --output-dir results/ \
        --diagram --scatter \
        -p 5
        """

}

process cnvkit_call{
    publishDir "${params.output}/${params.date}/intermediate_files/cnvkit/cnv", mode: "copy"
    cpus 5
    memory "5 G"

    input:
        file(cnr)

    output:
        file("${cnr.baseName}.call.cnv")
    script:
        if(!params.cnvkit.cellularity)
        """
        cnvkit.py call ${cnr} --drop-low-coverage -m threshold --t=${params.cnvkit.threshold} -o ${cnr.baseName}.call.cnv
        """
        else
        """
        cnvkit.py call ${cnr} -y -m clonal --purity ${params.cnvkit.cellularity} -o ${cnr.baseName}.call.cnv
        """
}

process annotate_cnv {
    cpus 5
    memory "5 G"
    publishDir "${params.output}/${params.date}/annotation/somatic/${input.simpleName}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${input.simpleName}/annotation/somatic/", mode: "copy"

    input:
        file(input)
    output:
        file("${input.simpleName}.cnv.annotated.tsv")
    script:
    if (params.pipeline.toUpperCase() == "DRAGEN")
        """
        echo ${params}
        name="\$(basename ${input} | sed 's,vcf\\.gz,bed,g')"
        zcat ${input} | awk '{if(\$7=="PASS") print \$3,\$5}' > \$name
        sed -i 's/DRAGEN:GAIN://g' \$name
        sed -i 's/DRAGEN:LOSS://g' \$name
        sed -i 's/:/\\t/g' \$name
        sed -i 's/-/\\t/g' \$name
        sed -i 's/<//g' \$name
        sed -i 's/>//g' \$name
        awk '{print \$1, \$2, \$3 +1, \$4}' \$name > appo
        mv appo \$name
        sed -i 's/ /\\t/g' \$name
        python ${params.classifyCNV.folder}/ClassifyCNV.py --infile \$name --GenomeBuild ${params.build} --cores 5 --outdir tmp
        mv tmp/Scoresheet.txt ${input.simpleName}.cnv.annotated.tsv

        """
    else if (params.pipeline.toUpperCase() == "SAREK")
        """
        awk 'BEGIN{getline;}{if(\$4!="Antitarget"){if(\$6>3) print  \$1, \$2,\$3, "DUP"; if(\$6<1) print \$1, \$2,\$3, "DEL"}}' ${input} > appo
        sed -i 's/ /\\t/g' appo
        python ${params.classifyCNV.folder}/ClassifyCNV.py --infile appo --GenomeBuild ${params.build} --cores 5 --outdir tmp
        mv tmp/Scoresheet.txt ${input.simpleName}.cnv.annotated.tsv
        """
}