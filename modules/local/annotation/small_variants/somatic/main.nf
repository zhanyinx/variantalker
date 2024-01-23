process fixvcf{
    cpus 1
    memory "1 G"
    container "docker://yinxiu/gatk:latest"

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), file("${meta.patient}.vcf.gz"), file("${meta.patient}.vcf.gz.tbi")
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
            tabix -p vcf ${meta.patient}.vcf.gz
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
            tabix -p vcf ${meta.patient}.vcf.gz
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
        tabix -p vcf ${meta.patient}.vcf.gz
        """
}


process somatic_annotate_snp_indel{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 8.GB * task.attempt }
    tag "vcf2maf"
    container "docker://yinxiu/gatk:latest"
    input:
        tuple val(meta), file(vcf), file(index)
    output:
        tuple val(meta), file("${meta.patient}.small_mutations.cancervar.escat.maf"), file("${meta.patient}.vcf")
    script:
    """
    normal="\$(zcat ${vcf} | grep 'normal_sample'  | cut -d'=' -f2)"
    tumor="\$(zcat ${vcf} | grep 'tumor_sample'  | cut -d'=' -f2)"
    zcat ${vcf} | awk '{if(\$7 == "PASS") print \$0; if( (\$0 ~/^#/) ) print \$0}' > ${meta.patient}.vcf

    awk 'BEGIN{counts = 0}{ if(\$1==chr && \$2==pos){counts++;}else{counts=0;}; if(\$0~/^#/){print \$0 > "header";} else {print \$0 > "tmp_"counts".vcf"}; pos=\$2; chr=\$1}' ${meta.patient}.vcf

    for file in `ls tmp_*vcf`; do
        cat header \$file > file.vcf
        bgzip -c file.vcf > file.vcf.gz
        tabix -p vcf file.vcf.gz

        # GATK funcotator
        # gatk Funcotator \
        #    -L ${params.funcotator_target} \
        #    -R ${params.fasta} \
        #    -V file.vcf.gz \
        #    -O ${meta.patient}.maf \
        #    --annotation-default Matched_Norm_Sample_Barcode:\${normal} \
        #    --annotation-default Tumor_Sample_Barcode:\${tumor} \
        #    --annotation-default Tumor_type:${meta.tumor_tissue} \
        #    --remove-filtered-variants true \
        #    --output-file-format MAF \
        #    --data-sources-path ${params.funcotator_somatic_db}\
        #    --ref-version ${params.build} \
        #    --transcript-selection-mode ${params.transcript_selection} \
        #    --splice-site-window-size ${params.splice_site_window_size}  \
        #    --interval-padding ${params.target_padding}

        
        check=1
        while [[ \$check -ne 0 ]]
        do
            # GATK funcotator
            gatk Funcotator \
                -L ${params.funcotator_target} \
                -R ${params.fasta} \
                -V file.vcf.gz \
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


        # cancervar call
        cp ${params.cancervar_init} config.init
        sed -i "s,INPUTYPE,${params.cancervar_input_type},g" config.init
        sed -i "s,BUILD,${params.build},g" config.init
        sed -i "s,INPUTFILE,file.vcf,g" config.init
        sed -i "s,OUTFILE,cancervar,g" config.init
        sed -i "s,ANNOVARDB,${params.annovar_db},g" config.init
        sed -i "s,ANNOVAR,${params.annovar_software_folder},g" config.init
        sed -i "s,CANCERVARDB,${params.cancervar_db},g" config.init

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
        add_guidelines_and_escat_to_maf.py -m ${meta.patient}.maf \
            -c \${cancervar_file} \
            -cc config.init \
            -o tmp \
            -t ${meta.tumor_tissue} \
            -p ${params.projectid} \
            --escat ${params.escat_db} \
            -d ${params.date} 
        
        if ! [ -f ${meta.patient}.small_mutations.cancervar.escat.maf ]; then
            cp tmp ${meta.patient}.small_mutations.cancervar.escat.maf
        else
            awk '{if(!(\$0 ~/^#/) && \$1!="Hugo_Symbol") print \$0}' tmp >> ${meta.patient}.small_mutations.cancervar.escat.maf
        fi

    done
    """
}