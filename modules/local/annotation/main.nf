// List of functions and processes used to annotate snp and indel from whole exome sequences

// add civic data to vcf
process add_civic{
    cpus 10
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }

    tag "civic2vcf"

    input:
        tuple val(meta), file(vcf), file(index)
    output:
        tuple val(meta), file("${meta.patient}.vcf.gz"), file("${meta.patient}.vcf.gz.tbi")
    script:
    """
    zcat ${vcf} > appo.vcf

    # split on chromosomes
    awk '{if(\$1~/^#/) print \$0}' appo.vcf > header
    for chr in `awk '{if(!(\$1~/^#/)) print \$1}' appo.vcf | sort | uniq`; do
        cp header tmp_\$chr
    done
    awk '{if(!(\$1~/^#/)) print \$0 >> "tmp_"\$1 }'  appo.vcf

    cpu=0
    for file in `ls tmp_*`; do
        civicpy annotate-vcf --input-vcf  \$file --output-vcf out.\$file --reference ${params.build_alt_name}  --include-status accepted --include-status submitted > log.\$file 2>&1 &
        let cpu=cpu+1

        if [ \$cpu -eq $task.cpus ]; then
            cpu=0
            wait
        fi
    done
    wait

    if [ -f final.vcf ]; then
        rm final.vcf
    fi

    for file in `ls out*`; do
        awk '{if(!(\$1~/^#/)) print \$0; else print \$0 > "header"}' \$file >> final.vcf
    done

    sort -k1,1V -k2,2n final.vcf > final.sorted.vcf
    sed -i 's/ /_/g' final.sorted.vcf
    cat header final.sorted.vcf > ${meta.patient}.vcf
    rm ${vcf} ${index}
    bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
    tabix -p vcf ${meta.patient}.vcf.gz
    """

}


// alpha_missense
process add_alpha_missense{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 3.GB * task.attempt }
    tag "alphamissense"

    input:
        tuple val(meta), file(maf), file(vcf)
    output:
        tuple val(meta), file("${maf.baseName}.missense.maf"), file("${meta.patient}.vcf")
    script:
    """
        if ! [ -f ${params.genomes[params.build].alpha_missense} ]; then
            echo "${params.genomes[params.build].alpha_missense} alpha Missense file does not exist." > .command.err
            exit 125
        fi

        zcat ${params.genomes[params.build].alpha_missense} > tmp
        awk -F'\\t' 'BEGIN{fn=0; count=0}{
            if(FNR==1) fn++
            if(fn==1){
                if(\$1~/^#/){
                    print \$0
                }else if(\$1=="Hugo_Symbol"){
                    print "# Alpha Missense file ""'"${params.genomes[params.build].alpha_missense}"'"
                    printf "%s\\t%s\\t%s\\n", \$0, "am_pathogenicity", "am_class"
                }else{
                    converter[\$5,\$6] = count
                    string[\$5,\$6]=\$5\$6\$12\$13
                    line[count]=\$0
                    bool[\$5,\$6]++
                    count++
                }
            }

            if(fn==2){
                if(bool[\$1,\$2] > 0){
                    if(string[\$1,\$2] == \$1\$2\$3\$4){
                        am_pathogenicity[converter[\$1,\$2]] = \$9
                        am_class[converter[\$1,\$2]] = \$10
                    }
                }
            }
        }END{for(i=0;i<count;i++) printf "%s\\t%s\\t%s\\n", line[i], am_pathogenicity[i], am_class[i]}' ${maf} tmp > ${maf.baseName}.missense.maf
        rm tmp
    """
}

// filter maf file
process filter_maf{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 2
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"
    tag "filtermaf"

    input:
        tuple val(meta), file(maf), file(vcf)
    output:
        tuple val(meta), file("${meta.patient}.maf"), file("${meta.patient}.vcf")
        tuple file("filtered*.pass.tsv"), file("filtered*.nopass.tsv")
    script:
    """
        filter_variants.py -m ${maf} \
         -o ${meta.patient}.maf \
         --filter_intervar "${params.filter_intervar}" \
         --filter_cancervar "${params.filter_cancervar}" \
         --filter_renovo "${params.filter_renovo}" \
         --sample_type ${meta.sample_type} \
         --min_depth ${params.filter_min_depth} \
         --vaf_threshold_germline ${params.filter_vaf_threshold_germline} \
         --vaf_threshold ${params.filter_vaf_threshold} \
         --filter_genes_germline ${params.filter_genes_germline} \
         --filter_genes_somatic ${params.filter_genes_somatic} \
         --filter_variant_classification "${params.filter_var_classification}" \
         --filter_civic "${params.filter_civic_evidence_level}"

    """
    
}

// *************************************************************** For germline mutation annotation ***************************************************************
// fix vcf header in case of dragen vcf, remove split multi allelic and remove 0 af multi allelic from ION 
process fixvcf{
    cpus 1
    memory "1 G"

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


// annotate vcf with funcotator and cancervar and output to maf format
process somatic_annotate_snp_indel{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 8.GB * task.attempt }
    tag "vcf2maf"

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
            -d ${params.date} 
        
        if ! [ -f ${meta.patient}.small_mutations.cancervar.escat.maf ]; then
            cp tmp ${meta.patient}.small_mutations.cancervar.escat.maf
        else
            awk '{if(!(\$0 ~/^#/) && \$1!="Hugo_Symbol") print \$0}' tmp >> ${meta.patient}.small_mutations.cancervar.escat.maf
        fi

    done
    """
}


// *************************************************************** For germline mutation annotation ***************************************************************

process filter_variants {
    cpus 1
    memory "4 G"
    tag "filter_variants"

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

