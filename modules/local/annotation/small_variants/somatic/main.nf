process standardize_somatic_vcf{
    cpus 1
    memory "1 G"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("${meta.patient}.vcf.gz")
    script:
    if (params.pipeline.toUpperCase() != "SAREK")
        if (params.tumoronly)
            """
            if [[ $vcf != *.gz ]]; then 
                bgzip $vcf 
                mv ${vcf}.gz $vcf
            fi
            
            fix_vcf_header4funcotator.sh -i $vcf -o "tmp.vcf" -t
            
            # normalise and split multiallelic
            bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf

            # remove 0 allele frequency and wt genotype
            awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${meta.patient}.vcf
            bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
            """
        else
            """
            if [[ $vcf != *.gz ]]; then 
                bgzip $vcf 
                mv ${vcf}.gz $vcf
            fi
            fix_vcf_header4funcotator.sh -i $vcf -o "tmp.vcf" 
            bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf
            awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${meta.patient}.vcf
            bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
            """


    else
        """
        if [[ $vcf != *.gz ]]; then 
            bgzip $vcf 
            mv ${vcf}.gz $vcf
        fi
        zcat $vcf > tmp.vcf
        bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf
        awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${meta.patient}.vcf
        bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
        """
}

process run_somatic_funcotator{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    tag "somatic_funcotator"
    container "docker://yinxiu/gatk:latest"
    input:
        tuple val(meta), val(chunk_index), file(vcf), file(index)
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.maf"), file("${meta.patient}.vcf")
    script:
    """
    normal="\$(zcat ${vcf} | grep 'normal_sample'  | cut -d'=' -f2)"
    tumor="\$(zcat ${vcf} | grep 'tumor_sample'  | cut -d'=' -f2)"
    trascript_params=""
    if [ -f ${params.transcript_list} ]; then
        trascript_params=" --transcript-list ${params.transcript_list}"
    fi

    gatk Funcotator \
        -L ${params.funcotator_target} \
        -R ${params.fasta} \
        -V ${vcf} \
        -O ${meta.patient}.maf \
        --annotation-default Matched_Norm_Sample_Barcode:\${normal} \
        --annotation-default Tumor_Sample_Barcode:\${tumor} \
        --annotation-default Tumor_type:${meta.tumor_tissue} \
        --remove-filtered-variants true \
        --output-file-format MAF \
        --data-sources-path ${params.funcotator_somatic_db}\
        --ref-version ${params.build} \
        --transcript-selection-mode ${params.transcript_selection} \
        --splice-site-window-size ${params.splice_site_window_size} \
        --interval-padding ${params.target_padding} \
        \$trascript_params

    zcat ${vcf} > ${meta.patient}.vcf

    # check completeness of maf file
    nrow_maf=\$(wc -l ${meta.patient}.maf | cut -f1 -d' ')
    nrow_vcf=\$(wc -l ${meta.patient}.vcf | cut -f1 -d' ')
    diff=\$((\$nrow_maf-\$nrow_vcf))
    if [ \$diff -lt 5 ]; then
        echo "Funcotator maf file is incomplete!"
        exit
    fi
    """
}

process run_cancervar{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    tag "cancervar"
    container "docker://yinxiu/gatk:latest"
    input:
        tuple val(meta), val(chunk_index), file(vcf), file(index)
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.cancervar"), file("${meta.patient}.grl_p"), file("config.init")
    script:
    """
    # cancervar call
    zcat ${vcf} > file.vcf
    cp ${params.cancervar_init} config.init
    sed -i "s,INPUTYPE,${params.cancervar_input_type},g" config.init
    sed -i "s,BUILD,${params.build},g" config.init
    sed -i "s,INPUTFILE,file.vcf,g" config.init
    sed -i "s,OUTFILE,cancervar,g" config.init
    sed -i "s,ANNOVARDB,${params.annovar_db},g" config.init
    sed -i "s,ANNOVAR,${params.annovar_software_folder},g" config.init
    sed -i "s,CANCERVARDB,${params.cancervar_db},g" config.init
    sed -i "s,SPLICE_WINDOW,${params.splice_site_window_size},g" config.init

    if ! [ ${params.cancervar_evidence_file} == "None" ]; then
        if ! [ -f ${params.cancervar_evidence_file} ]; then
            echo "${params.cancervar_evidence_file} cancervar evidence file does not exist!"
            exit
        fi
        sed -i "s,None,${params.cancervar_evidence_file},g" config.init
    fi

    python ${params.cancervar_folder}/CancerVar.py -c config.init --cancer_type=${meta.tumor_tissue}

    # merge cancervar and funcotator
    if [ ${params.tumoronly} == "true" ]; then
        cancervar_file="cancervar.${params.build}_multianno.txt.cancervar"
    else
        cancervar_file="\$(ls cancervar.${vcf.simpleName}*.${params.build}_multianno.txt.cancervar)"
    fi
    cp \${cancervar_file} ${meta.patient}.cancervar
    grl_p=`echo \$cancervar_file | sed 's/\\.cancervar/\\.grl_p/g'`
    cp \$grl_p ${meta.patient}.grl_p
    """
}

process add_guidelines_escat{
    cpus 1
    memory "1 G"
    container "docker://yinxiu/gatk:latest"
    errorStrategy 'retry'
    maxRetries = 3
    tag "add_guidelines_escat"
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(cancervar), file(grl_p), file(config)
    output:
        tuple val(meta), val(chunk_index), file("${chunk_index}.guidelines.maf"), file("${chunk_index}.vcf")
    script:
    """
    add_guidelines_and_escat_to_maf.py -m ${maf} \
        -c ${cancervar} \
        -cc ${config} \
        -o ${chunk_index}.guidelines.maf \
        -t ${meta.tumor_tissue} \
        -p ${params.projectid} \
        --escat ${params.escat_db} \
        -d ${params.date} 

    cp ${vcf} ${chunk_index}.vcf
    """
}