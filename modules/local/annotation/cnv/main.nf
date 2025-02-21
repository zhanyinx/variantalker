process cnvkit_call{
    publishDir "${params.outdir}/${params.date}/intermediate_files/cnvkit/cnv", mode: "copy"
    cpus 5
    memory "5 G"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.10--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.10--pyhdfd78af_0' }"
    
    input:
        tuple val(meta), path(cnr)

    output:
        tuple val(meta), path("${meta.patient}.call.cnv")
    script:
        if(!params.cnvkit_cellularity || params.cnvkit_cellularity.isEmpty())
        """
        cnvkit.py call ${cnr} --drop-low-coverage -m threshold --t=${params.cnvkit_threshold} -o ${meta.patient}.call.cnv
        """
        else
        """
        cnvkit.py call ${cnr} -y -m clonal --purity ${params.cnvkit_cellularity} -o ${meta.patient}.call.cnv
        """
}

process annotate_cnv {
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"

    input:
        tuple val(meta), path(input)
    output:
        file("${meta.patient}.cnv.annotated.tsv")
    script:
    if (params.pipeline.toUpperCase() == "DRAGEN")
        """
        name="\$(basename ${input} | sed 's,vcf\\.gz,bed,g')"
        zcat ${input} | awk '{if(\$7=="PASS") print \$3,\$5, \$10}' > \$name
        sed -i 's/DRAGEN:GAIN://g' \$name
        sed -i 's/DRAGEN:LOSS://g' \$name
        sed -i 's/:/\\t/g' \$name
        sed -i 's/-/\\t/g' \$name
        sed -i 's/<//g' \$name
        sed -i 's/>//g' \$name
        awk '{print \$1, \$2, \$3 +1, \$4, \$6}' \$name > appo
        mv appo \$name
        sed -i 's/ /\\t/g' \$name
        python ${params.classifyCNV_folder}/ClassifyCNV.py --infile \$name --GenomeBuild ${params.build} --cores 5 --outdir tmp
        awk 'BEGIN{fn=0}{if(FNR==1) fn++; if(fn==1){score[\$1,\$2,\$3] = \$5}; if(fn==2){if(header==0) {print \$0, "logRatio"; header++}else{print \$0, score[\$2,\$3,\$4]}}}' \$name tmp/Scoresheet.txt > ${meta.patient}.cnv.annotated.tsv
        """
    else if (params.pipeline.toUpperCase() == "SAREK")
        """
        awk 'BEGIN{getline;}{if(\$4!="Antitarget"){if(\$6>3) print  \$1, \$2,\$3, "DUP"; if(\$6<1) print \$1, \$2,\$3, "DEL"}}' ${input} > appo
        sed -i 's/ /\\t/g' appo
        python ${params.classifyCNV_folder}/ClassifyCNV.py --infile appo --GenomeBuild ${params.build} --cores 5 --outdir tmp
        mv tmp/Scoresheet.txt ${meta.patient}.cnv.annotated.tsv
        """
}
