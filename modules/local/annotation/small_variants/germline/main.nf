process filter_variants {
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

process normalise_rename_germline_vcf {
    cpus 1
    memory "1 G"
    tag "rename_and_index"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("*.vcf1.gz"), file("*.vcf1.gz.tbi")
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


process germline_annotate_snp_indel{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 8.GB * task.attempt }
    tag "vcf2maf"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), file(vcf), file(index) 
    output:
        tuple val(meta), file("${meta.patient}.small_mutations.intervar.escat.maf"), file("${meta.patient}.vcf")
    script:
    """
    zcat ${vcf} | awk '{if(\$7 == "PASS") print \$0; if( (\$0 ~/^#/) ) print \$0}' > ${meta.patient}.vcf

    cat ${meta.patient}.vcf | awk 'BEGIN{counts = 0}{ if(\$1==chr && \$2==pos){counts++;}else{counts=0;}; if(\$0~/^#/){print \$0 > "header";} else {print \$0 > "tmp_"counts".vcf"}; pos=\$2; chr=\$1}'
    
    for file in `ls tmp_*vcf`; do
        cat header \$file > file.vcf
        bgzip -c file.vcf > file.vcf.gz
        tabix -p vcf file.vcf.gz

        # # GATK funcotator
        # gatk Funcotator \
        #     -L ${params.funcotator_target} \
        #     -R ${params.fasta} \
        #     -V file.vcf.gz \
        #     -O ${meta.patient}.maf \
        #     --annotation-default Matched_Norm_Sample_Barcode:${meta.patient} \
        #     --remove-filtered-variants true \
        #     --output-file-format MAF \
        #     --data-sources-path ${params.funcotator_germline_db} \
        #     --ref-version ${params.build} \
        #     --transcript-selection-mode ${params.transcript_selection} \
        #     --splice-site-window-size ${params.splice_site_window_size} \
        #     --interval-padding ${params.target_padding}

        check=1
        while [[ \$check -ne 0 ]]
        do
            # GATK funcotator
            gatk Funcotator \
                -L ${params.funcotator_target} \
                -R ${params.fasta} \
                -V file.vcf.gz \
                -O ${meta.patient}.maf \
                --annotation-default Matched_Norm_Sample_Barcode:${meta.patient} \
                --remove-filtered-variants true \
                --output-file-format MAF \
                --data-sources-path ${params.funcotator_germline_db}\
                --ref-version ${params.build} \
                --transcript-selection-mode ${params.transcript_selection} \
                --splice-site-window-size ${params.splice_site_window_size} \
                --interval-padding ${params.target_padding}

            check=\$?
            if [[ \$check -ne 0 ]]; then
                grep -nr "ERROR GencodeFuncotationFactory" .command.err | awk '{split(\$14,res,"-"); print res[1]}'  | sort | uniq | awk '{split(\$0,res,":"); print res[1],res[2]}'  > remove
                awk 'BEGIN{fn=0; count=1; a=0}{
                    if(FNR==1) fn++;

                    if(fn==1){
                        chr[count] = \$1
                        pos[count] = \$2
                        count++
                    }    

                    if(fn==2){
                        bool=1
                        for(i=1;i<count;i++){
                            if(\$1 == chr[i] && \$2 == pos[i]){
                                bool=0
                                break
                            }
                        }
                        if(bool==1) print \$0
                    }
                }' remove file.vcf > appo
                mv appo file.vcf
                bgzip -c file.vcf > file.vcf.gz
                tabix -p vcf file.vcf.gz
            fi
        done

        # intervar call
        cp ${params.intervar_init} config.init
        sed -i "s,INPUTYPE,${params.intervar_input_type},g" config.init
        sed -i "s,BUILD,${params.build},g" config.init
        sed -i "s,INPUTFILE,file.vcf,g" config.init
        sed -i "s,OUTFILE,intervar,g" config.init
        sed -i "s,ANNOVARDB,${params.annovar_db},g" config.init
        sed -i "s,ANNOVAR,${params.annovar_software_folder},g" config.init
        sed -i "s,INTERVARDB,${params.intervar_db},g" config.init
        if ! [ ${params.intervar_evidence_file} == "None" ]; then
            if ! [ -f ${params.intervar_evidence_file} ]; then
                echo "${params.intervar_evidence_file} intervar evidence file does not exist!"
                exit
            fi
            sed -i "s,None,${params.intervar_evidence_file},g" config.init
        fi


        python ${params.intervar_folder}/Intervar.py -c config.init

        # merge intervar and funcotator
        # intervarfile="\$(ls intervar.${vcf.simpleName}*.${params.build}_multianno.txt.intervar)"
        intervarfile="intervar.${params.build}_multianno.txt.intervar"
        add_guidelines_and_escat_to_maf.py -m ${meta.patient}.maf \
            -c \${intervarfile} \
            -cc config.init \
            -o tmp \
            -p ${params.projectid} \
            -d ${params.date} \
            --escat ${params.escat_db} \
            --germline
        
        if ! [ -f ${meta.patient}.small_mutations.intervar.escat.maf ]; then
            cp tmp ${meta.patient}.small_mutations.intervar.escat.maf
        else
            awk '{if(!(\$0 ~/^#/) && \$1!="Hugo_Symbol") print \$0}' tmp >> ${meta.patient}.small_mutations.intervar.escat.maf
        fi

    done
    """
}


process germline_renovo_annotation{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    tag "vcf2maf"
    container "docker://yinxiu/renovo:latest"

    input:
        tuple val(meta), file(maf), file(vcf)
    output:
        tuple val(meta), file("${maf.baseName}.renovo.maf"), file("${meta.patient}.vcf")
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
