
process generate_pyclone{
    publishDir "${params.output}/${params.date}/biomarkers/${patient}/pyclone", mode: "copy"
    input:
        tuple val(patient), val(sex), path(mafs), val(cellularity), path(crams), path(crais), path(pluriploidy), path(cnvs)
    output:
       tuple val(patient), path("*.pyclone.tsv")
    script:
    """
        if [ -f appo ]; then
            rm appo
        fi

        IFS=' ' read -r -a mafs <<< "$mafs"
        IFS=' ' read -r -a cellularity <<< "$cellularity"
        IFS=' ' read -r -a crams <<< "$crams"
        IFS=' ' read -r -a crais <<< "$crais"
        IFS=' ' read -r -a pluriploidy <<< "$pluriploidy"
        IFS=' ' read -r -a cnvs <<< "$cnvs"
        IFS=' ' read -r -a sexes <<< "$sex"

        for index in \${!mafs[@]}; do
            patientmaf=\${mafs[\$index]}
            awk -F '\\t' '{if(\$10=="SNP") print \$5"-"\$6"-"\$7"-"\$11"-"\$13"-"\$1"-"\$9"-"\$10}' \$patientmaf >> appo
        done

        cat appo | sort | uniq > list_unique_mutations
        sed -i 's/-/ /g' list_unique_mutations
        line=`wc -l list_unique_mutations | awk '{print \$1}'`

        echo "mutation_id sample_id ref_counts alt_counts normal_cn major_cn minor_cn tumour_content" > ${patient}.pyclone.tsv
        sed -i 's/ /\\t/g' ${patient}.pyclone.tsv

        re='^[0-9]+\$'

        for index in \${!mafs[@]}; do
            echo "Chromosome Start_Position End_Position Tumor_Sample_Barcode t_ref_count t_alt_count Genome_Change" > tmp_maf
            tumor_sample=`basename \${mafs[\$index]}`
            for idx in `seq 1 \$line`; do
                tumor_bam=\${crams[\$index]}
                strings=`awk '{if(NR=="'"\$idx"'"+0.) print \$1":"\$2"-"\$3}' list_unique_mutations`
                start_end_chr=`awk '{if(NR=="'"\$idx"'"+0.) print \$1,\$2,\$3}' list_unique_mutations`
                alt=`awk '{if(NR=="'"\$idx"'"+0.) print \$5}' list_unique_mutations`
                altcounts=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} |  cut -f 5 | tr '[a-z]' '[A-Z]' | fold -w 1 | sort | uniq -c | awk 'BEGIN{ee=0}{if(\$2=="'"\$alt"'") {print \$1; ee=99}}END{if(ee==0) print "0"}'`

                ref=`awk '{if(NR=="'"\$idx"'"+0.) print \$4}' list_unique_mutations`
                ref_detected=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} | awk '{print \$3}'`
                refcounts=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} | awk '{print \$4 - "'"\$altcounts"'"+0.}'`
                echo "\$start_end_chr \${tumor_sample} \$refcounts \$altcounts \$strings.\$ref.\$alt" >> tmp_maf

            done

            sed -i 's, ,\\t,g' tmp_maf

            sex=`echo \${sexes[\$index]}  | sed 's,\\[,,g' | sed 's,\\],,g' | sed 's/,//g'`
            if [[ \${cellularity[\$index]} =~ \$re ]] ; then
                cell=`echo \${cellularity[\$index]}   | sed 's,\\[,,g' | sed 's,\\],,g' | sed 's/,//g'`
                create_input4pyclone.py -as \${cnvs[\$index]} -c \${cell} -m tmp_maf -o ${patient}.pyclone.tsv -s \$sex
            else
                create_input4pyclone.py -as \${cnvs[\$index]} -m tmp_maf -o ${patient}.pyclone.tsv -ac \${pluriploidy[\$index]} -s \$sex
            fi
        done

    """
}


process pyclone{
    publishDir "${params.output}/${params.date}/biomarkers/${patient}/", mode: "copy"
    input:
        tuple val(patient), path(pyclone)
    output:
       tuple val(patient), path("${patient}.clonalTMB.txt")
       file("${patient}.pyclone.output.tsv")
    script:
    """
        # run pyclone
        pyclone-vi fit -i ${pyclone} -o ${patient}.pyclone.h5 -c ${params.num_clusters} -d beta-binomial -r ${params.num_restarts}
        pyclone-vi write-results-file -i ${patient}.pyclone.h5 -o ${patient}.pyclone.output.tsv
        
        # extract clonal tmb
        awk 'BEGIN{
            max=-99; ncount=0; getline
        }{
            if(\$4>max){
                max = \$4; ncount = 1;
            }else if(\$4 == max){
                ncount++
            }
        }END{print "Clonal TMB: ", ncount}' ${patient}.pyclone.output.tsv > ${patient}.clonalTMB.txt
    """
}