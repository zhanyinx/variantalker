process filter_germline_variants {
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

process run_germline_intervar {
    input:
        tuple val(meta), val(chunk_index), file(vcf), file(vcf_full), file(maf_skip)
        val safe_annovar_db
        val safe_annovar_software_folder
        val safe_intervar_evidence_file
    output:
        tuple val(meta), val(chunk_index), file("${meta.patient}.intervar"), file("${meta.patient}.grl_p"), file("config.init")
    script:
    """
    # Intervar call
    cat ${vcf} > file.vcf
    cp ${params.intervar_init} config.init
    sed -i "s,INPUTYPE,VCF,g" config.init
    sed -i "s,BUILD,${params.build},g" config.init
    sed -i "s,INPUTFILE,file.vcf,g" config.init
    sed -i "s,OUTFILE,intervar,g" config.init
    sed -i "s,ANNOVARDB,${safe_annovar_db},g" config.init
    sed -i "s,ANNOVAR,${safe_annovar_software_folder},g" config.init
    sed -i "s,INTERVARDB,${params.intervar_db},g" config.init
    sed -i "s,SPLICE_WINDOW,${params.splice_site_window_size},g" config.init

    if ! [ "${safe_intervar_evidence_file}" == "None" ]; then
        if ! [ -f "${safe_intervar_evidence_file}" ]; then
            echo "${safe_intervar_evidence_file} intervar evidence file does not exist!"
            exit 1
        fi
        sed -i "s,None,${safe_intervar_evidence_file},g" config.init
    fi

    python ${params.intervar_folder}/Intervar.py -c config.init

    intervarfile="intervar.${params.build}_multianno.txt.intervar"
    cp \${intervarfile} ${meta.patient}.intervar
    grl_p=`echo \$intervarfile | sed 's/\\.intervar/\\.grl_p/g'`
    cp \$grl_p ${meta.patient}.grl_p
    """
}

process run_germline_renovo {
    input:
        tuple val(meta), val(chunk_index), file(maf), file(vcf), file(vcf_full), file(maf_skip)
        val safe_annovar_db
        val safe_annovar_software_folder
    output:
        tuple val(meta), val(chunk_index), file("${maf.baseName}.renovo.maf"), file(vcf), file(vcf_full), file(maf_skip)
    script:
    """
        unset R_LIBS
        unset R_LIBS_USER
        unset R_LIBS_SITE

        export HOME=/tmp
        export R_PROFILE_USER=/dev/null
        export R_ENVIRON_USER=/dev/null
        
        mkdir tmp
        mv ${vcf_full} tmp/

        python ${params.renovo_path}/ReNOVo.py \
        -p . -a ${safe_annovar_software_folder} \
        -d ${safe_annovar_db} \
        -b ${params.build}
        
        add_renovo_to_maf.py -m ${maf} \
         -r ReNOVo_output/${vcf.baseName}_ReNOVo_and_ANNOVAR_implemented.txt \
         -o ${maf.baseName}.renovo.maf
         
        mv tmp/${vcf_full} ./
        rm -r tmp
    """
}