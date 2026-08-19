process fix_somatic_header {

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("fixed.vcf")

    script:
    """
    file="${vcf}"

    if [[ ${vcf} != *.gz ]]; then
        bgzip ${vcf}
        file="${vcf}.gz"
    fi

    ncol=\$(zcat \$file | awk -F '\t' '{if(!(\$1~/^#/)) {print NF; exit}}')

    if [ \$ncol -eq 10 ]; then
        fix_vcf_header4funcotator.sh -i \$file -o fixed.vcf -t
    else
        fix_vcf_header4funcotator.sh -i \$file -o fixed.vcf
    fi
    """
}

process add_somatic_civic{
    cpus 3
    errorStrategy 'retry'
    memory { 4.GB * task.attempt }
    container "docker://yinxiu/civicpy:v1.0"

    tag "civic2vcf"

    input:
        tuple val(meta), val(chunk_index), file(vcf), file(vcf_full), file(maf_skip) 
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.civic.vcf"), file(vcf_full), file(maf_skip) 
    script:
    """
    if [ ${params.skip_civic} == "true" ]; then
        cp ${vcf} ${meta.patient}.civic.vcf
    else
        cp ${vcf} appo.vcf
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
            cat header final.sorted.vcf > ${meta.patient}.civic.vcf
            check=`wc -l ${meta.patient}.civic.vcf | awk '{print \$1-1}'`
        done
    fi
    """
}

process run_somatic_cancervar{
    input:
        tuple val(meta), val(chunk_index), file(vcf), file(vcf_full), file(maf_skip)
        val safe_annovar_db
        val safe_annovar_software_folder
        val safe_cancervar_evidence_file
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.cancervar"), file("${meta.patient}.grl_p"), file("config.init")
    script:
    """
    # cancervar call
    cat ${vcf} > file.vcf
    cp ${params.cancervar_init} config.init
    sed -i "s,BUILD,${params.build},g" config.init
    sed -i "s,INPUTFILE,file.vcf,g" config.init
    sed -i "s,OUTFILE,cancervar,g" config.init
    sed -i "s,ANNOVARDB,${safe_annovar_db},g" config.init
    sed -i "s,ANNOVAR,${safe_annovar_software_folder},g" config.init
    sed -i "s,CANCERVARDB,${params.cancervar_db},g" config.init
    sed -i "s,SPLICE_WINDOW,${params.splice_site_window_size},g" config.init

    if ! [ "${safe_cancervar_evidence_file}" == "None" ]; then
        if ! [ -f "${safe_cancervar_evidence_file}" ]; then
            echo "${safe_cancervar_evidence_file} cancervar evidence file does not exist!"
            exit 1
        fi
        sed -i "s,None,${safe_cancervar_evidence_file},g" config.init
    fi

    ncol=`awk -F '\t' '{if(!(\$1~/^#/)) {print NF; exit}}' file.vcf`
    
    if [ \$ncol -eq 10 ]; then
       cancervar_input_type="VCF"
    else
       cancervar_input_type="VCF_m"
    fi

    sed -i "s,INPUTYPE,\$cancervar_input_type,g" config.init

    python ${params.cancervar_folder}/CancerVar.py -c config.init --cancer_type=${meta.tumor_tissue}

    # merge cancervar and funcotator
    if [ \$cancervar_input_type == "VCF" ]; then
        cancervar_file="cancervar.${params.build}_multianno.txt.cancervar"
    else
        cancervar_file="\$(ls cancervar.${vcf.simpleName}*.${params.build}_multianno.txt.cancervar)"
    fi
    cp \${cancervar_file} ${meta.patient}.cancervar
    grl_p=`echo \$cancervar_file | sed 's/\\.cancervar/\\.grl_p/g'`
    cp \$grl_p ${meta.patient}.grl_p
    """
}



