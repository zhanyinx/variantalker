process standardize_somatic_vcf{
    cpus 1
    memory "1 G"
    container "docker://yinxiu/gatk:latest"
    tag "standardize_vcf"

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

// add civic data to vcf
process add_somatic_civic{
    cpus 3
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    container "docker://yinxiu/civicpy:v1.0"

    tag "civic2vcf"

    input:
        tuple val(meta), val(chunk_index), file(vcf) 
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.vcf.gz")
    script:
    """
    zcat ${vcf} > appo.vcf
    nline=`wc -l appo.vcf | awk '{print \$1}'`
    check=0
    while [ \$nline -ne \$check ]; do

        # split header and body
        awk '{if(\$1~/^#/) print \$0 > "header"; else print \$0 > "body"}' appo.vcf
        nlines=`wc -l body | awk '{print int(\$1/"'${task.cpus}'" + 1)}'`
        split -l \$nlines body chunks

        for file in `ls chunks*`; do
            awk '{print \$0}' header \$file > tmp
            mv tmp \$file
        done

        export CIVICPY_CACHE_FILE="/home/.civicpy/cache.pkl"
        export CIVICPY_CACHE_TIMEOUT_DAYS=${params.civic_cache_timeout_days}        
        

        for file in `ls chunks*`; do
            civicpy annotate-vcf --input-vcf  \$file --output-vcf out.\$file --reference ${params.build_alt_name}  --include-status accepted --include-status submitted > log.\$file 2>&1 &
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
        rm ${vcf}
        bgzip -c ${meta.patient}.vcf > ${meta.patient}.vcf.gz
        check=`wc -l ${meta.patient}.vcf | awk '{print \$1-1}'`
    done
    """

}

process run_somatic_cancervar{
    cpus 1
    errorStrategy 'retry'
    maxRetries = 3
    memory { 4.GB * task.attempt }
    tag "cancervar"
    container "docker://yinxiu/gatk:latest"
    input:
        tuple val(meta), val(chunk_index), file(vcf)
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

