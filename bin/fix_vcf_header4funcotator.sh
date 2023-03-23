#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : fix_vcf_header4funcotator.sh -i INPUT -o OUTPUT [-t] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Dragen vcf file does not have normal and tumor sample label in the header. 
	  This creates problem in reading allele frequency in funcotator. 
	  This script adds normal and tumor sample label to the vcf file, and reindex it."
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT : VCF.gz file as input"
    echo "   -o|--output OUTPUT : output VCF.gz file"
    echo "   [-t|--tumoronly]: if defined, dragen vcf file is tumor only file"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--output") set -- "$@" "-o" ;;
      "--tumoronly")   set -- "$@" "-t" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

input=""
output=""
tumoronly=0

while getopts ":i:o:th" OPT
do
    case $OPT in
        i) input=$OPTARG;;
        o) output=$OPTARG;;
        t) tumoronly=1;;
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

if [ $# -lt 4 ]
then
    usage
    exit
fi

if ! [ -f $input ]; then
    echo "$input does not exist!"
    exit
fi

if [ $tumoronly -eq 0 ]; then

	zcat $input | awk '{
		if($1 =="#CHROM"){
			print "##tumor_sample="$NF; 
			print "##normal_sample="$(NF-1); 
		} 
		print $0;
	}' > $output

else
	zcat $input | awk '{
		if($1 =="#CHROM"){
			print "##tumor_sample="$NF; 
		} 
		print $0;
	}' > $output

fi
