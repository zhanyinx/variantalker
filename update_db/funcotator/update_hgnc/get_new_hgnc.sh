#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : get_new_hgnc.sh [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Get the latest HGNC database for funcotator"
    echo "---------------"
    echo "OPTIONS"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done


while getopts ":h" OPT
do
    case $OPT in
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


today=$(date +%b%d%Y)

# get data
curl 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_locus_type&col=gd_locus_group&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_date_mod&col=gd_date_sym_change&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=family.id&col=family.name&col=gd_ccds_ids&col=gd_vega_ids&col=md_eg_id&col=md_mim_id&col=md_refseq_id&col=md_prot_id&col=md_ensembl_id&col=md_ucsc_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit%27'  > appo

# rename header
awk 'BEGIN{getline; print "HGNC ID\tApproved Symbol\tApproved Name\tStatus\tLocus Type\tLocus Group\tPrevious Symbols\tPrevious Name\tSynonyms\tName Synonyms\tChromosome\tDate Modified\tDate Symbol Changed\tDate Name Changed\tAccession Numbers\tEnzyme IDs\tEntrez Gene ID\tEnsembl Gene ID\tPubmed IDs\tRefSeq IDs\tGene Family ID\tGene Family Name\tCCDS IDs\tVega ID\tEntrez Gene ID(supplied by NCBI)\tOMIM ID(supplied by OMIM)\tRefSeq(supplied by NCBI)\tUniProt ID(supplied by UniProt)\tEnsembl ID(supplied by Ensembl)\tUCSC ID(supplied by UCSC)"}{print $0}' appo > hgnc_$today.tsv

# create config
echo "name = HGNC
version = $today
src_file = hgnc_$today.tsv
origin_location = https://www.genenames.org/cgi-bin/download
preprocessing_script = downloadHgncDataSource.sh

# Whether this data source is for the b37 reference.
# Required and defaults to false.
isB37DataSource = false

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
#       cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = simpleXSV

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for GENCODE files.
# NCBI build version (either hg19 or hg38):
ncbi_build_version =

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key = GENE_NAME

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column = 1

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter =\t

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols = true

# Required field for locatableXSV files.
# The 0-based index of the column containing the contig for each row
contig_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the start position for each row
start_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the end position for each row
end_column =" > hgnc.config

# create new database
mkdir -p hgnc/hg38
mv hgnc.config hgnc_$today.tsv hgnc/hg38
cp -r hgnc/hg38 hgnc/hg19