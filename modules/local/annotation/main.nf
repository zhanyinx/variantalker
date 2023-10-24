// List of functions and processes used to annotate snp and indel from whole exome sequences

// fix vcf header in case of dragen vcf, remove split multi allelic and remove 0 af multi allelic from ION 
process fixvcf{
    cpus 1
    memory "1 G"

    input:
        tuple val(patient), val(tumor), path(vcf)
    output:
        tuple val(patient), val(tumor), file("${patient}.vcf.fix.gz"), file("${patient}.vcf.fix.gz.tbi")
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
            awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${patient}.vcf
            bgzip -c ${patient}.vcf > ${patient}.vcf.fix.gz
            tabix -p vcf ${patient}.vcf.fix.gz
            """
        else
            """
            if [[ $vcf != *.gz ]]; then 
                bgzip $vcf 
                mv ${vcf}.gz $vcf
            fi
            fix_vcf_header4funcotator.sh -i $vcf -o "tmp.vcf" 
            bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf
            awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${patient}.vcf
            bgzip -c ${patient}.vcf > ${patient}.vcf.fix.gz
            tabix -p vcf ${patient}.vcf.fix.gz
            """


    else
        """
        if [[ $vcf != *.gz ]]; then 
            bgzip $vcf 
            mv ${vcf}.gz $vcf
        fi
        zcat $vcf > tmp.vcf
        bcftools norm -m-any --check-ref -w -f ${params.fasta} tmp.vcf -o tmp1.vcf
        awk '{if(!(\$8~/AF=0;/) && !(\$NF~/0\\/0/)) print \$0}' tmp1.vcf > ${patient}.vcf
        bgzip -c ${patient}.vcf > ${patient}.vcf.fix.gz
        tabix -p vcf ${patient}.vcf.fix.gz
        """
}

// rename vcf file to have standard naming system.
process rename_somatic_vcf {
    cpus 1
    memory "1 G"
    tag "rename"

    input:
        tuple val(patient), val(tumor), file(vcf), file(index)
    output:
        tuple val(patient), val(tumor), file("*.vcf.gz"), file("*.vcf.gz.tbi")
    script:
        """
        tumor="\$(zcat ${vcf} | grep 'tumor_sample'  | cut -d'=' -f2)"
        cp $vcf \${tumor}.vcf.gz
        cp $index \${tumor}.vcf.gz.tbi
        """
}


// annotate vcf with funcotator and cancervar and output to maf format
process somatic_annotate_snp_indel{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 8.GB * task.attempt }
    publishDir "${params.output}/${params.date}/annotation/somatic/${patient}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${vcf.simpleName}/annotation/somatic/", mode: "copy"
    tag "vcf2maf"

    input:
        tuple val(patient), val(tumor_type), file(vcf), file(index)
    output:
        file("${patient}.small_mutations.cancervar.escat.maf")
        file("filtered.${patient}.small_mutations.cancervar.escat.maf.nopass.tsv")
        file("filtered.${patient}.small_mutations.cancervar.escat.maf.pass.tsv")
    script:
    """
    normal="\$(zcat ${vcf} | grep 'normal_sample'  | cut -d'=' -f2)"
    tumor="\$(zcat ${vcf} | grep 'tumor_sample'  | cut -d'=' -f2)"
    zcat ${vcf} | awk '{if(\$7 == "PASS") print \$0; if( (\$0 ~/^#/) ) print \$0}' > ${patient}.vcf

    awk 'BEGIN{counts = 0}{ if(\$1==chr && \$2==pos){counts++;}else{counts=0;}; if(\$0~/^#/){print \$0 > "header";} else {print \$0 > "tmp_"counts".vcf"}; pos=\$2; chr=\$1}' ${patient}.vcf

    for file in `ls tmp_*vcf`; do
        cat header \$file > file.vcf
        bgzip -c file.vcf > file.vcf.gz
        tabix -p vcf file.vcf.gz

        # GATK funcotator
        # gatk Funcotator \
        #    -L ${params.funcotator_target} \
        #    -R ${params.fasta} \
        #    -V file.vcf.gz \
        #    -O ${patient}.maf \
        #    --annotation-default Matched_Norm_Sample_Barcode:\${normal} \
        #    --annotation-default Tumor_Sample_Barcode:\${tumor} \
        #    --annotation-default Tumor_type:${tumor_type} \
        #    --remove-filtered-variants true \
        #    --output-file-format MAF \
        #    --data-sources-path ${params.funcotator_somatic_db}\
        #    --ref-version ${params.build} \
        #    --transcript-selection-mode ${params.transcript_selection} \
        #    --interval-padding ${params.target_padding}

        
        check=1
        while [[ \$check -ne 0 ]]
        do
            # GATK funcotator
            gatk Funcotator \
                -L ${params.funcotator_target} \
                -R ${params.fasta} \
                -V file.vcf.gz \
                -O ${patient}.maf \
                --annotation-default Matched_Norm_Sample_Barcode:\${normal} \
                --annotation-default Tumor_Sample_Barcode:\${tumor} \
                --annotation-default Tumor_type:${tumor_type} \
                --remove-filtered-variants true \
                --output-file-format MAF \
                --data-sources-path ${params.funcotator_somatic_db}\
                --ref-version ${params.build} \
                --transcript-selection-mode ${params.transcript_selection} \
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

        python ${params.cancervar_folder}/CancerVar.py -c config.init --cancer_type=${tumor_type}

        # merge cancervar and funcotator
        if [ ${params.tumoronly} == "true" ]; then
            cancervar_file="cancervar.${params.build}_multianno.txt.cancervar"
        else
            cancervar_file="\$(ls cancervar.${vcf.simpleName}*.${params.build}_multianno.txt.cancervar)"
        fi
        add_guidelines_and_escat_to_maf.py -m ${patient}.maf \
            -c \${cancervar_file} \
            -cc config.init \
            -o tmp \
            -t ${tumor_type} \
            -p ${params.projectid} \
            -d ${params.date} 
        
        if ! [ -f ${patient}.small_mutations.cancervar.escat.maf ]; then
            cp tmp ${patient}.small_mutations.cancervar.escat.maf
        else
            awk '{if(!(\$0 ~/^#/) && \$1!="Hugo_Symbol") print \$0}' tmp >> ${patient}.small_mutations.cancervar.escat.maf
        fi

    done
    mv ${patient}.small_mutations.cancervar.escat.maf tmp.maf 
    filter_variants.py -m tmp.maf \
         -o ${patient}.small_mutations.cancervar.escat.maf \
         --somatic \
         -md ${params.filter_min_depth} \
         -vt ${params.filter_vaf_threshold} \
         --filter_cancervar "${params.filter_cancervar}" \
         --filter_genes_somatic ${params.filter_genes_somatic} \
         --filter_variant_classification "${params.filter_var_classification}"
    """
}


// For germline mutation annotation 

process filter_variants {
    cpus 1
    memory "4 G"
    tag "filter_variants"

    input:
        tuple val(patient), path(vcf) 
    output:
        tuple val(patient), path("*.vcf.gz") 
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

    input:
        tuple val(patient), path(vcf)
    output:
        tuple val(patient), file("*.vcf1.gz"), file("*.vcf1.gz.tbi")
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
        bgzip -c tmp2.vcf > ${patient}.vcf1.gz
        tabix -p vcf ${patient}.vcf1.gz
        """
}


process germline_annotate_snp_indel{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 8.GB * task.attempt }
    // publishDir "${params.output}/${params.date}/annotation/germline/${vcf.simpleName}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${vcf.simpleName}/annotation/germline/", mode: "copy"
    tag "vcf2maf"

    input:
        tuple val(patient), file(vcf), file(index) 
    output:
        tuple file("${patient}.small_mutations.intervar.escat.maf"), file("${patient}.vcf")
    script:
    """
    zcat ${vcf} | awk '{if(\$7 == "PASS") print \$0; if( (\$0 ~/^#/) ) print \$0}' > ${patient}.vcf

    cat ${patient}.vcf | awk 'BEGIN{counts = 0}{ if(\$1==chr && \$2==pos){counts++;}else{counts=0;}; if(\$0~/^#/){print \$0 > "header";} else {print \$0 > "tmp_"counts".vcf"}; pos=\$2; chr=\$1}'
    
    for file in `ls tmp_*vcf`; do
        cat header \$file > file.vcf
        bgzip -c file.vcf > file.vcf.gz
        tabix -p vcf file.vcf.gz

        # # GATK funcotator
        # gatk Funcotator \
        #     -L ${params.funcotator_target} \
        #     -R ${params.fasta} \
        #     -V file.vcf.gz \
        #     -O ${patient}.maf \
        #     --annotation-default Matched_Norm_Sample_Barcode:${patient} \
        #     --remove-filtered-variants true \
        #     --output-file-format MAF \
        #     --data-sources-path ${params.funcotator_germline_db} \
        #     --ref-version ${params.build} \
        #     --transcript-selection-mode ${params.transcript_selection} \
        #     --interval-padding ${params.target_padding}

        check=1
        while [[ \$check -ne 0 ]]
        do
            # GATK funcotator
            gatk Funcotator \
                -L ${params.funcotator_target} \
                -R ${params.fasta} \
                -V file.vcf.gz \
                -O ${patient}.maf \
                --annotation-default Matched_Norm_Sample_Barcode:${patient} \
                --remove-filtered-variants true \
                --output-file-format MAF \
                --data-sources-path ${params.funcotator_germline_db}\
                --ref-version ${params.build} \
                --transcript-selection-mode ${params.transcript_selection} \
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
        add_guidelines_and_escat_to_maf.py -m ${patient}.maf \
            -c \${intervarfile} \
            -cc config.init \
            -o tmp \
            -p ${params.projectid} \
            -d ${params.date} \
            --germline
        
        if ! [ -f ${patient}.small_mutations.intervar.escat.maf ]; then
            cp tmp ${patient}.small_mutations.intervar.escat.maf
        else
            awk '{if(!(\$0 ~/^#/) && \$1!="Hugo_Symbol") print \$0}' tmp >> ${patient}.small_mutations.intervar.escat.maf
        fi

    done
    """
}


process germline_renovo_annotation{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.output}/${params.date}/annotation/germline/${maf.simpleName}", mode: "copy"
    // publishDir "${params.output}/${params.date}/${maf.simpleName}/annotation/germline/", mode: "copy"
    tag "vcf2maf"

    input:
        tuple file(maf), file(vcf) 
    output:
        file("${maf.baseName}.renovo.maf")
        file("filtered.${maf.baseName}.renovo.maf.pass.tsv")
        file("filtered.${maf.baseName}.renovo.maf.nopass.tsv")
    script:
    """
        python ${params.renovo_path}/ReNOVo.py \
        -p . -a ${params.annovar_software_folder} \
        -d ${params.annovar_db} \
        -b ${params.build} 
        
        add_renovo_to_maf.py -m ${maf} \
         -r ReNOVo_output/${vcf.baseName}_ReNOVo_and_ANNOVAR_implemented.txt \
         -o ${maf.baseName}

        filter_variants.py -m ${maf.baseName} \
         -o ${maf.baseName}.renovo.maf \
         --filter_intervar "${params.filter_intervar}" \
         --filter_renovo "${params.filter_renovo}" \
         --germline \
         -md ${params.filter_min_depth} \
         -vtg ${params.filter_vaf_threshold_germline} \
         --filter_genes_germline ${params.filter_genes_germline} \
         --filter_variant_classification "${params.filter_var_classification}"

    """
}