#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

# Defence in depth, not the authority: the Python caller validates the tables this script
# produces before it lets either the annotation script or the config start naming them, and that
# is what decides whether an update is adopted. What this adds is that the script stops at the
# first real failure instead of grinding on and leaving something that looks plausible.
#
# `pipefail` is the half that matters, and `set -e` alone would not have caught the failure this
# was added for. Both `date=` assignments below take their value from a pipeline, and a
# pipeline's status is its LAST command's -- so when `zcat` fails on a truncated or absent
# download, `tr` still succeeds, the assignment succeeds, `$date` is empty, and the script
# carries on to build `clinvar_hg38_.vcf`. `pipefail` is what surfaces `zcat`'s failure.
set -eo pipefail

function usage {
    echo -e "usage : update_clinvar_annovar.sh -v VT -o OUTPUT -n name [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Generate clinvar latest version for annovar."
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -v|--vt VT : path to vt software"
    echo "   -n|--name NAME : database name (e.g. clinvar_date)"
    echo "   -o|--output OUTPUT : folder where to save databases (typically the humandb from annovar)"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--vt") set -- "$@" "-v" ;;
      "--name") set -- "$@" "-n" ;;
      "--output") set -- "$@" "-o" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

vt=""
output=""
name=""

while getopts ":v:n:o:h" OPT
do
    case $OPT in
        v) vt=$OPTARG;;
        n) name=$OPTARG;;
        o) output=$OPTARG;;
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

# These three refusals used to be a bare `exit`, which exits with the status of the last command
# run -- an `echo` -- so the script reported SUCCESS for having declined to do anything at all.
if [ $# -lt 6 ]
then
    usage
    exit 1
fi

if ! [ -d "$vt" ]; then
    echo "$vt folder does not exist!"
    exit 1
fi

if ! [ -d "$output" ]; then
    echo "$output folder does not exist! if this is expected, make it first: mkdir -p $output"
    exit 1
fi


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "=== ClinVar Annovar Database Update ==="
echo "VT path: $vt"
echo "Database name: $name"
echo "Output directory: $output"
echo ""

# Process hg38 (GRCh38)
echo "Processing GRCh38..."
wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O clinvar_hg38.vcf.gz
if [ ! -f "clinvar_hg38.vcf.gz" ]; then
    echo "ERROR: Failed to download GRCh38 ClinVar VCF"
    exit 1
fi

date=$(zcat clinvar_hg38.vcf.gz | grep "fileDate" | cut -d "=" -f2- | tr -d '-')
filename="clinvar_hg38_${date}.vcf"

# decompose multiple allelic data
"$vt/vt" decompose clinvar_hg38.vcf.gz -o "$filename"

if [ ! -f "$filename" ]; then
    echo "ERROR: Failed to decompose GRCh38 VCF"
    exit 1
fi

# Remove problematic YT lines before GATK processing
grep -v $'\tYT\t' "$filename" > "${filename}.tmp" && mv "${filename}.tmp" "$filename"

# extract information using gatk tools
gatk VariantsToTable -F CHROM -F POS -F POS -F REF -F ALT -F ALLELEID -F CLNDN -F CLNDISDB -F CLNREVSTAT -F CLNSIG -V "$filename" -O "hg38_${name}.txt"

if [ ! -f "hg38_${name}.txt" ]; then
    echo "ERROR: Failed to create hg38_${name}.txt with GATK"
    exit 1
fi

# Fix escaped characters - use perl for more reliable replacement
perl -i -pe 's/provided,_/provided\\x2c_/g' "hg38_${name}.txt"

# remodel header
awk 'BEGIN{getline; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#Chr","Start","End","Ref","Alt","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG"}{
	print $0
}' "hg38_${name}.txt" > "hg38_${name}.txt.tmp"
mv "hg38_${name}.txt.tmp" "hg38_${name}.txt"

# remove chr prefix
perl -i -pe 's/^chr//' "hg38_${name}.txt"

# index 
perl "${SCRIPT_DIR}/utils/index_annovar.pl" "hg38_${name}.txt" 1000 > "hg38_${name}.txt.idx"

# clean intermediate files
rm -f clinvar_hg38.vcf.gz "$filename"

echo "  SUCCESS: Created hg38_${name}.txt"

# Process hg19 (GRCh37)
echo "Processing GRCh37..."
wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz -O clinvar_hg19.vcf.gz
if [ ! -f "clinvar_hg19.vcf.gz" ]; then
    echo "ERROR: Failed to download GRCh37 ClinVar VCF"
    exit 1
fi

date=$(zcat clinvar_hg19.vcf.gz | grep "fileDate" | cut -d "=" -f2- | tr -d '-')
filename="clinvar_hg19_${date}.vcf"

# decompose multiple allelic data
"$vt/vt" decompose clinvar_hg19.vcf.gz -o "$filename"

if [ ! -f "$filename" ]; then
    echo "ERROR: Failed to decompose GRCh37 VCF"
    exit 1
fi

# Remove problematic YT lines before GATK processing
grep -v $'\tYT\t' "$filename" > "${filename}.tmp" && mv "${filename}.tmp" "$filename"

# extract information using gatk tools
gatk VariantsToTable -F CHROM -F POS -F POS -F REF -F ALT -F ALLELEID -F CLNDN -F CLNDISDB -F CLNREVSTAT -F CLNSIG -V "$filename" -O "hg19_${name}.txt"

if [ ! -f "hg19_${name}.txt" ]; then
    echo "ERROR: Failed to create hg19_${name}.txt with GATK"
    exit 1
fi

# Fix escaped characters - use perl for more reliable replacement
perl -i -pe 's/provided,_/provided\\x2c_/g' "hg19_${name}.txt"

# remodel header
awk 'BEGIN{getline; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#Chr","Start","End","Ref","Alt","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG"}{
        print $0
}' "hg19_${name}.txt" > "hg19_${name}.txt.tmp"
mv "hg19_${name}.txt.tmp" "hg19_${name}.txt"

# remove chr prefix
perl -i -pe 's/^chr//' "hg19_${name}.txt"

# index 
perl "${SCRIPT_DIR}/utils/index_annovar.pl" "hg19_${name}.txt" 1000 > "hg19_${name}.txt.idx"

# clean intermediate files
rm -f clinvar_hg19.vcf.gz "$filename"

echo "  SUCCESS: Created hg19_${name}.txt"

# Move final files to output directory
mv "hg19_${name}.txt" "hg19_${name}.txt.idx" "hg38_${name}.txt" "hg38_${name}.txt.idx" "$output/"

echo ""
echo "=== ClinVar Annovar Update Complete ==="
echo "Files created in $output:"
echo "  - hg19_${name}.txt"
echo "  - hg19_${name}.txt.idx"
echo "  - hg38_${name}.txt"
echo "  - hg38_${name}.txt.idx"
exit 0