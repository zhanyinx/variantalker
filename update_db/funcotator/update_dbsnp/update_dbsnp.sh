#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : update_dbsnp.sh -d db_dir [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Get the latest dbsnp database for funcotator"
    echo "---------------"
    echo "OPTIONS"
    echo "   -d|--db db_dir : directory of the current dbsnp vcf. Funcotator folder with hg38 and hg19 within it."
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--db") set -- "$@" "-d" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

db_dir=""

while getopts ":d:h" OPT
do
    case $OPT in
        d) db_dir=$OPTARG;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done


if [ $# -lt 2 ]
then
    usage
    exit
fi

if ! [ -d $db_dir ]; then
    echo "$db_dir directory does not exist!"
    exit
fi

if ! [ -d "$db_dir/hg38" ]; then
    echo "WARNING!! $db_dir does not contain hg38 directory!"
fi

if ! [ -d "$db_dir/hg19" ]; then
    echo "WARNING!! $db_dir does not contain hg19 directory!"
fi

# get info
curl https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/ 2>/dev/null > appo
awk '{if($0 ~/^<a href/) print $0}' appo | cut -d '"' -f2 > list_files


######################################## hg19 ########################################
latest_version=`awk '{print $1}' list_files | grep -v "tbi" | grep -v "md5" | grep GCF | sort | awk 'BEGIN{getline; print $0}'`
curl https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}.md5 > ${latest_version}.md5
latest_md5sum=`awk '{print $1}' ${latest_version}.md5`

# check if current file is different from latest
check=0
if ! [ -f "$db_dir/hg19/${latest_version}.md5" ]; then
    let check=check+1
else
    current_md5sum=`awk '{print $1}' $db_dir/hg19/${latest_version}.md5`

    if ! [[ $latest_md5sum == $current_md5sum ]]; then
        let check=check+1
    fi
fi

# if current file is different from former
if ! [ $check -eq 0 ]; then

    # get data
    wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}
    name=`echo $latest_version | sed 's,\.gz,,g'`

    # rename chromosomes
    zcat  $latest_version |  awk 'BEGIN{getline; prev = $1}{if($1!=prev && !($1~/^#/)) print $1; prev=$1}' | grep "NC_0" > primary_chromosomes.txt
    awk 'BEGIN{for(i=1;i<=22;i++) print "chr"i; print "chrX"; print "chrY"; print "chrM"}' > ucsc_primary_chromosomes.txt
    command=`paste primary_chromosomes.txt ucsc_primary_chromosomes.txt | awk '{printf "%s; ", "gsub(/"$1"/,\""$2"\")"}'`
    zcat  $latest_version | awk '{'"$command"'; print}' > $name

    # get version
    version=`head -n 10 $name | grep dbSNP_BUILD_ID | awk -F '=' '{print $2}'`

    bgzip -c $name > $name.gz
    tabix -p vcf $name.gz
    rm $name

    mkdir -p dbsnp/hg19/
    mv $name.gz $name.gz.tbi dbsnp/hg19/
    mv ${latest_version}.md5 dbsnp/hg19/

    echo "name = dbSNP
version = b$version 
src_file = $name.gz
origin_location = https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}
preprocessing_script = update_dnsnp.sh

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
# cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = vcf

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key =

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column =

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter =

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols =

# Required field for locatableXSV files.
# The 0-based index of the column containing the contig for each row
contig_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the start position for each row
start_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the end position for each row
end_column =" > dbSNP.config

    mv dbSNP.config dbsnp/hg19/

    # cleaning
    rm appo primary_chromosomes.txt ucsc_primary_chromosomes.txt
fi

######################################## hg38 ########################################
latest_version=`awk '{print $1}' list_files | grep -v "tbi" | grep -v "md5" | grep GCF | sort | awk 'BEGIN{getline;getline; print $0}'`
curl https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}.md5 > ${latest_version}.md5
latest_md5sum=`awk '{print $1}' ${latest_version}.md5`

# check if current file is different from latest
check=0
if ! [ -f "$db_dir/hg38/${latest_version}.md5" ]; then
    let check=check+1
else
    current_md5sum=`awk '{print $1}' $db_dir/hg38/${latest_version}.md5`

    if ! [[ $latest_md5sum == $current_md5sum ]]; then
        let check=check+1
    fi
fi

# if current file is different from former
if ! [ $check -eq 0 ]; then

    # get data
    wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}
    name=`echo $latest_version | sed 's,\.gz,,g'`

    # rename chromosomes
    zcat  $latest_version |  awk 'BEGIN{getline; prev = $1}{if($1!=prev && !($1~/^#/)) print $1; prev=$1}' | grep "NC_0" > primary_chromosomes.txt
    awk 'BEGIN{for(i=1;i<=22;i++) print "chr"i; print "chrX"; print "chrY"; print "chrM"}' > ucsc_primary_chromosomes.txt
    command=`paste primary_chromosomes.txt ucsc_primary_chromosomes.txt | awk '{printf "%s; ", "gsub(/"$1"/,\""$2"\")"}'`
    zcat  $latest_version | awk '{'"$command"'; print}' > $name

    # get version
    version=`head -n 10 $name | grep dbSNP_BUILD_ID | awk -F '=' '{print $2}'`

    bgzip -c $name > $name.gz
    tabix -p vcf $name.gz
    rm $name

    mkdir -p dbsnp/hg38/
    mv $name.gz $name.gz.tbi dbsnp/hg38/
    mv ${latest_version}.md5 dbsnp/hg38/

    echo "name = dbSNP
version = b$version 
src_file = $name.gz
origin_location = https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/${latest_version}
preprocessing_script = update_dnsnp.sh

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
# cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = vcf

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key =

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column =

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter =

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols =

# Required field for locatableXSV files.
# The 0-based index of the column containing the contig for each row
contig_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the start position for each row
start_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the end position for each row
end_column =" > dbSNP.config

    mv dbSNP.config dbsnp/hg38/

    rm list_files appo primary_chromosomes.txt ucsc_primary_chromosomes.txt
fi



