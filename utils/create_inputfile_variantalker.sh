#!/bin/bash
    
## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : $(basename $0) -i INPUT -t TISSUE_SAMPLE [-o OUTPUT] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Create input file for variantalker annotation. See https://github.com/zhanyinx/variantalker#input"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT: folder that contains folders with *hard-filtered*vcf.gz and *cnv.vcf.gz files"
    echo "   -t|--tissue SAMPLE_TISSUE: Sample tissue."
    echo "   [-o|--output OUTPUT]: output name of the csv file, default out.csv"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--tissue") set -- "$@" "-t" ;;
      "--output") set -- "$@" "-o" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

path=""
tissue=""
output="out.csv"

while getopts ":i:t:o:h" OPT
do
    case $OPT in
        i) path=$OPTARG;;
        t) tissue=$OPTARG;;
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

if [ $# -lt 4 ]
then
    usage
    exit
fi


path=`realpath $path`

if ! [ -d $path ]; then
    echo "$path directory does not exist!"
    exit
fi



echo "patient,tumor_tissue,sample_file,sample_type" > $output

for file in `ls -d $path/*/*hard*vcf.gz | grep -v germline`; do

	dir=`dirname $file`
	name=`basename $dir`
	echo "$name,$tissue,$file,somatic" >> $output
	

done

for file in `ls -d $path/*/*hard*vcf.gz | grep germline`; do

        dir=`dirname $file`
        name=`basename $dir`
        echo "$name,$tissue,$file,germline" >> $output
        

done

for file in `ls -d $path/*/*cnv*vcf.gz`; do
        dir=`dirname $file`
        name=`basename $dir`
        echo "$name,,$file,cnv" >> $output
done
