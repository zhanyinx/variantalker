#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


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

if [ $# -lt 6 ]
then
    usage
    exit
fi

if ! [ -d $vt ]; then
    echo "$vt folder does not exist!"
    exit
fi

if ! [ -d $output ]; then
    echo "$output folder does not exist! if this is expected, make it first: mkdir -p $output"
    exit
fi


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# download data
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
date=`zcat clinvar.vcf.gz | grep "fileDate" | cut -d "=" -f2- | sed 's/-//g'`
filename=clinvar_$date.vcf

# decompose multiple allelic data
$vt/vt decompose clinvar.vcf.gz -o $filename


# extract information using gatk tools
# name=`basename $filename | sed 's,\.vcf,,g'`
sed -i '/\tYT\t/d' $filename
gatk VariantsToTable -F CHROM -F POS -F POS -F REF -F ALT -F ALLELEID -F CLNDN -F CLNDISDB -F CLNREVSTAT -F CLNSIG -V $filename -O hg38_$name.txt
sed -i 's/provided,_/provided\\x2c_/g' hg38_$name.txt

# remodel header
awk 'BEGIN{getline; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#Chr","Start","End","Ref","Alt","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG"}{

	print $0
}' hg38_$name.txt > a
mv a hg38_$name.txt

# remove chr
sed -i 's/chr//' hg38_$name.txt

# index 
perl ${SCRIPT_DIR}/utils/index_annovar.pl hg38_$name.txt 1000 > hg38_$name.txt.idx

# clean
rm clinvar*




# download data
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
date=`zcat clinvar.vcf.gz | grep "fileDate" | cut -d "=" -f2- | sed 's/-//g'`
filename=clinvar_$date.vcf

# decompose multiple allelic data
$vt/vt decompose clinvar.vcf.gz -o $filename


# extract information using gatk tools
# name=`basename $filename | sed 's,\.vcf,,g'`
sed -i '/\tYT\t/d' $filename
gatk VariantsToTable -F CHROM -F POS -F POS -F REF -F ALT -F ALLELEID -F CLNDN -F CLNDISDB -F CLNREVSTAT -F CLNSIG -V $filename -O hg19_$name.txt
sed -i 's/provided,_/provided\\x2c_/g' hg19_$name.txt

# remodel header
awk 'BEGIN{getline; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#Chr","Start","End","Ref","Alt","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG"}{

        print $0
}' hg19_$name.txt > a
mv a hg19_$name.txt

# remove chr
sed -i 's/chr//' hg19_$name.txt

# index 
perl ${SCRIPT_DIR}/utils/index_annovar.pl hg19_$name.txt 1000 > hg19_$name.txt.idx

# clean
rm clinvar*

mv hg19_$name.txt hg19_$name.txt.idx hg38_$name.txt hg38_$name.txt.idx $output