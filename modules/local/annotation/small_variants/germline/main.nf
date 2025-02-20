process filter_germline_variants {
    cpus 1
    memory "4 G"
    tag "filter_variants"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), path(vcf) 
    output:
        tuple val(meta), path("*.vcf.gz") 
    script:
        """
        tabix -p vcf ${vcf}
        gatk VariantFiltration \
            -V ${vcf} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O filtered.vcf.gz
        """
}

process standardize_germline_vcf {
    cpus 1
    memory "1 G"
    tag "rename_and_index"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("*.vcf1.gz")
    script:
        """
        if [[ $vcf != *.gz ]]; then 
            bgzip $vcf 
            mv ${vcf}.gz $vcf
        fi
        # name="\$(zcat $vcf | awk '{if(\$1 =="#CHROM"){print \$NF} }')"
        zcat $vcf > tmp.vcf
        rm $vcf
        bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf
        awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > tmp2.vcf
        bgzip -c tmp2.vcf > ${meta.patient}.vcf1.gz
        tabix -p vcf ${meta.patient}.vcf1.gz
        """
}


process run_germline_intervar{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    tag "cancervar"
    container "docker://yinxiu/gatk:latest"
    input:
        tuple val(meta), val(chunk_index), file(vcf)
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.intervar"), file("${meta.patient}.grl_p"), file("config.init")
    script:
    """
    # Intervar call
    zcat ${vcf} > file.vcf
    cp ${params.intervar_init} config.init
    sed -i "s,INPUTYPE,${params.intervar_input_type},g" config.init
    sed -i "s,BUILD,${params.build},g" config.init
    sed -i "s,INPUTFILE,file.vcf,g" config.init
    sed -i "s,OUTFILE,intervar,g" config.init
    sed -i "s,ANNOVARDB,${params.annovar_db},g" config.init
    sed -i "s,ANNOVAR,${params.annovar_software_folder},g" config.init
    sed -i "s,INTERVARDB,${params.intervar_db},g" config.init
    sed -i "s,SPLICE_WINDOW,${params.splice_site_window_size},g" config.init

    if ! [ ${params.intervar_evidence_file} == "None" ]; then
        if ! [ -f ${params.intervar_evidence_file} ]; then
            echo "${params.intervar_evidence_file} intervar evidence file does not exist!"
            exit
        fi
        sed -i "s,None,${params.intervar_evidence_file},g" config.init
    fi

    python ${params.intervar_folder}/Intervar.py -c config.init

    intervarfile="intervar.${params.build}_multianno.txt.intervar"
    cp \${intervarfile} ${meta.patient}.intervar
    grl_p=`echo \$intervarfile | sed 's/\\.intervar/\\.grl_p/g'`
    cp \$grl_p ${meta.patient}.grl_p
    """
}


process run_germline_renovo{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    tag "vcf2maf"
    container "docker://yinxiu/renovo:latest"

    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf)
    output:
        tuple val(meta), val(chunk_index), file("${maf.baseName}.renovo.maf"), file(vcf)
    script:
    """
        python ${params.renovo_path}/ReNOVo.py \
        -p . -a ${params.annovar_software_folder} \
        -d ${params.annovar_db} \
        -b ${params.build} 
        
        add_renovo_to_maf.py -m ${maf} \
         -r ReNOVo_output/${vcf.baseName}_ReNOVo_and_ANNOVAR_implemented.txt \
         -o ${maf.baseName}.renovo.maf
    """
}
