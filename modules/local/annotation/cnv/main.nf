process cnvkit_call{
    publishDir "${params.outdir}/${params.date}/intermediate_files/cnvkit/cnv", mode: "copy"
    cpus 5
    memory "5 G"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.10--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.10--pyhdfd78af_0' }"
    
    input:
        tuple val(meta), path(cnr)

    output:
        tuple val(meta), path("${meta.patient}.call.cnv")
    script:
        if(!params.cnvkit_cellularity || params.cnvkit_cellularity.isEmpty())
        """
        cnvkit.py call ${cnr} --drop-low-coverage -m threshold --t=${params.cnvkit_threshold} -o ${meta.patient}.call.cnv
        """
        else
        """
        cnvkit.py call ${cnr} -y -m clonal --purity ${params.cnvkit_cellularity} -o ${meta.patient}.call.cnv
        """
}

process annotate_cnv {
    publishDir "${params.outdir}/${params.date}/annotation/${meta.sample_type}/${meta.patient}", mode: "copy"

    input:
        tuple val(meta), path(input)
    output:
        file("${meta.patient}.cnv.annotated.tsv")
    script:
    if (params.pipeline.toUpperCase() == "DRAGEN")
        """
        name="\$(basename "${input}" | sed 's,vcf\\.gz,bed,g')"

        zcat "${input}" | awk -v OFS="\\t" '
        BEGIN {
            header_found = 0
        }

        # Check if row "#CHROM" contains required columns and save their indices
        (/^#CHROM/) {
            for (i = 1; i <= NF; i++) {
                if (\$i == "ID") col_id = i
                else if (\$i == "ALT") col_alt = i
                else if (\$i == "FORMAT") {
                    col_format = i
                    col_sample = i + 1  # SAMPLE column right after FORMAT
                }
            }

            # Check if all required columns were found
            if (!col_id || !col_alt || !col_format || !col_sample) {
                printf("ERROR: incomplete header or missing columns\\n") > "/dev/stderr"
                exit 1
            }

            # Save header with only relevant columns
            print "ID", "ALT", "FORMAT", "SAMPLE"
            header_found = 1
            next
        }

        # Ignore rows with metadata "##"
        /^##/ { next }

        # Filter rows with FILTER=="PASS"
        (header_found && \$7 == "PASS") {
            print \$col_id, \$col_alt, \$col_format, \$col_sample
        }' > "\$name"

        awk -v fname="\$name" '
        BEGIN {
            OFS = "\\t"
        }

        NR == 1 {
            # Find FORMAT and SAMPLE columns
            for (i = 1; i <= NF; i++) {
                if (\$i == "FORMAT") {
                    col_format = i
                    col_sample = i + 1
                    break
                }
            }

            if (!col_format) {
                printf("ERROR [%s]: no FORMAT column found in header\\n", fname) > "/dev/stderr"
                exit 1
            }

            next  # skip header
        }

        {
            # Split FORMAT and SAMPLE
            n_key = split(\$col_format, keys, ":")
            n_val = split(\$col_sample, values, ":")

            if (n_key != n_val) {
                printf("WARNING [%s]: row %d -> %d keys but %d values\\n", fname, NR, n_key, n_val) > "/dev/stderr"
            }

            out = ""
            m = (n_key < n_val) ? n_key : n_val
            for (i = 1; i <= m; i++) {
                out = out keys[i] "=" values[i] "\\t"
            }
            sub(/\\t\$/, "", out)  # remove trailing tab

            # Rebuild row without FORMAT and SAMPLE
            for (i = 1; i <= NF; i++) {
                if (i == col_format) continue
                else if (i == col_sample) printf "%s%s", out, (i < NF ? OFS : "")
                else printf "%s%s", \$i, (i < NF ? OFS : "")
            }
            printf "\\n"
        }' "\$name" > tmp && mv tmp "\$name"

        sed -i 's/DRAGEN:[^:]*://g' "\$name"
        sed -i 's/:/\\t/g; s/-/\\t/g; s/[<>]//g' "\$name"

        awk -F"\\t" '{
            cnf = ""

            # Look for CNF first
            for (i = 5; i <= NF; i++) {
                if (\$i ~ /^CNF=/) {
                    split(\$i, arr, "=")
                    cnf = arr[2]
                    break
                }
            }

            # If CNF not found, try CN
            if (cnf == "") {
                for (i = 5; i <= NF; i++) {
                    if (\$i ~ /^CN=/) {
                        split(\$i, arr, "=")
                        cnf = arr[2]
                        break
                    }
                }
            }

            # If CNF & CN not found, try SM
            if (cnf == "") {
                for (i = 5; i <= NF; i++) {
                    if (\$i ~ /^SM=/) {
                        split(\$i, arr, "=")
                        cnf = arr[2]
                        break
                    }
                }
            }

            # Check CNF/CN/SM presence
            if (cnf == "" || cnf == ".") {
                printf("WARNING: CNF/CN/SM not found on line %d\\n", NR) > "/dev/stderr"
                cnf = NA
            }

            # Compute logR only if cnf is numeric
            if (cnf ~ /^[0-9.]+\$/) {
                logR = (cnf == 0) ? -17.600000 : log(cnf)/log(2)
            } else {
                logR = "NA"
            }

            # Print output (same format as before)
            printf "%s\\t%s\\t%s\\t%s\\t%.6f\\t%.6f\\n", \$1, \$2, \$3+1, \$4, logR, cnf
        }' "\$name" > aaa

        mv aaa "\$name"

        sed -i 's/ /\\t/g' "\$name"

        : > bbb
        awk '{if (\$4=="DEL" || \$4=="DUP") print > "bbb"}' "\$name"
        mv bbb "\$name"

        python ${params.classifyCNV_folder}/ClassifyCNV.py --infile \$name --GenomeBuild ${params.build} --cores 5 --outdir tmp
        awk 'BEGIN{OFS="\\t"; fn=0}
        FNR==NR {
            key = \$1 FS \$2 FS \$3;
            logr[key] = \$5;
            cnf[key]  = \$6;
            next
        }
        FNR==1 && NR!=FNR {
            print \$0, "logRatio", "CNF";
            next
        }
        {
            key = \$2 FS \$3 FS \$4;
            print \$0, logr[key], cnf[key]
        }' \$name tmp/Scoresheet.txt > ${meta.patient}.cnv.annotated.tsv
        """
    else if (params.pipeline.toUpperCase() == "SAREK")
        """
        awk 'BEGIN{getline;}{if(\$4!="Antitarget"){if(\$6>3) print  \$1, \$2,\$3, "DUP"; if(\$6<1) print \$1, \$2,\$3, "DEL"}}' ${input} > ccc
        sed -i 's/ /\\t/g' ccc
        python ${params.classifyCNV_folder}/ClassifyCNV.py --infile ccc --GenomeBuild ${params.build} --cores 5 --outdir tmp
        mv tmp/Scoresheet.txt ${meta.patient}.cnv.annotated.tsv
        """
}