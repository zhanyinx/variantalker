#!/bin/bash
    
## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : $(basename $0) -i INPUT -a ANNOTATION_FOLDER [-o OUTPUT] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Create input file for variantalker biomarker pipeline. See https://github.com/zhanyinx/variantalker#input"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT: path to variantalker annotation folder. It should contain germline and somatic folders."
    echo "   [-d|--dragen_folder DRAGEN: path to folder that contains dragen folders: should contain *.tmb.metrics.csv and *.microsat_output.json]"
    echo "   [-r|--dragen_rna_folder DRAGEN_RNA: path to folder that contains dragen RNA folders: should contain *.sf]"
    echo "   [-o|--output OUTPUT]: output name of the csv file, default out.csv"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--dragen_folder") set -- "$@" "-d" ;;
      "--dragen_rna_folder") set -- "$@" "-r" ;;
      "--output") set -- "$@" "-o" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

input=""
dragen_folder="random_folder"
dragen_rna_folder="random_folder"
output="out.csv"


while getopts ":i:d:r:o:h" OPT
do
    case $OPT in
        i) input=$OPTARG;;
        d) dragen_folder=$OPTARG;;
        r) dragen_rna_folder=$OPTARG;;
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

if [ $# -lt 2 ]
then
    usage
    exit
fi
input=`realpath $input`
if ! [ -d $input ]; then
    echo "$input directory does not exist!"
    exit
fi


echo "patient,sample_file,sample_type" > $output

# append dragen msi and tmb if they exist
if [ -d $dragen_folder ]; then
    echo "Appending dragen msi, tmb, coverage, hrd"
    abspath=`realpath $dragen_folder`
    for file in `ls $abspath/*/*.tmb.metrics.csv`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,tmb" >> $output
    done

    for file in `ls $abspath/*/*.microsat_output.json`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,msi" >> $output
    done

    for file in `ls $abspath/*/*.target_bed_coverage* | grep -v normal`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,coverage" >> $output
    done

    for file in `ls $abspath/*/*hrd*`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,hrd" >> $output
    done
fi

# append rna if they exist
if [ -d $dragen_rna_folder ]; then
    echo "Appending dragen rna files"
    abspath=`realpath $dragen_rna_folder`
    for file in `ls $abspath/*/*.sf`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,rna" >> $output
    done
fi

# Append germline maf samples
if [ -d $input/germline ]; then
    echo "Appending germline samples...."
    for file in `ls $input/germline/*/*maf`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        echo "$patient,$absfile,variant_germline" >> $output
    done
fi

# append somatic maf samples and rename corresponding germline
if [ -d $input/somatic ]; then
    echo "Appending somatic samples...."
    for file in `ls $input/somatic/*/*maf`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        patient_normal=`awk -F '\t' '{if(!($0~/^#/) && ($1!="Hugo_Symbol")) print $17}' $absfile | sort | uniq`
        echo "$patient,$absfile,variant_somatic" >> $output
        if ! [ -z $patient_normal ]; then
            sed -i 's;'${patient_normal},';'${patient},';g' $output
        fi
    done
fi

if [ -d $input/cnv ]; then
    for file in `ls $input/cnv/*/*cnv*tsv`; do
        absfile=`realpath $file`
        foldername=`dirname $absfile`
        patient=`basename $foldername`
        patient=`basename $patient | sed 's,_CNV,,g'`
        if ls $input/germline/${patient}* 1> /dev/null 2>&1; then
            echo "skilling $patient CNV"
        else
            echo "$patient,$absfile,cnv" >> $output
        fi
    done
fi