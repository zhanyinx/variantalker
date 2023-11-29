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
    echo "Create input file for clonal tmb analysis. See https://github.com/zhanyinx/clonal_evolution#input"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT: folder that contains folders with 2 cram files, one of them must be *_tumor.cram."
    echo "   -a|--annotation_folder ANNOTATION: folder with output from variantalker annotation. e.g. path2/somatic"
    echo "   [-o|--output OUTPUT]: output name of the csv file, default out.csv"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--annotation_folder") set -- "$@" "-a" ;;
      "--output") set -- "$@" "-o" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

input=""
annotation_folder=""
output="out.csv"
while getopts ":i:o:a:h" OPT
do
    case $OPT in
        i) input=$OPTARG;;
        a) annotation_folder=$OPTARG;;
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

if ! [ -d $input ]; then
    echo "$input directory does not exist!"
    exit
fi


if ! [ -d $annotation_folder ]; then
    echo "$annotation_folder directory does not exist!"
    exit
fi


selected_folders=""
for dir in `ls -d $input/*`; do 
	nline=`ls $dir/*cram | grep -v "evidence" | wc -l`; 
	if [ $nline -eq 2 ]; then
		selected_folders+=" $dir"
	fi
done

if ! [ -f $output ];then
	echo "patient,sex,status,sample,lane,cram,crai,cellularity,maf" > $output
fi

for dir in $selected_folders; do

	patient=`basename $dir`
	file=`ls ${annotation_folder}/$patient/*maf`
	normal=`ls $dir/*cram | grep -v tumor | grep -v "evidence"`
	normalbai=`ls $normal*i`
	tumor=`ls $dir/*cram | grep tumor | grep -v evidence`
	tumorbai=`ls $tumor*i`
	sex="XX"
	echo "$patient,$sex,1,tumor,1,$tumor,$tumorbai,,$file" >> $output
	echo "$patient,$sex,0,normal,1,$normal,$normalbai,," >> $output
done


echo "IMPORTANT: remember to change the sex in the clonal tmb input file!!!"
