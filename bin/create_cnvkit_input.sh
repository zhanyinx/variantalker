#!/bin/bash
    
## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : create_cnvkit_input.sh -i INPUT  [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Create tsv file for calling and annotating cnv from sarek"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT : folder where sarek was run; for multiple runs, put the common basefolder path"
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

input=""

while getopts ":i:h" OPT
do
    case $OPT in
        i) input=$OPTARG;;
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

if ! [ -d $input ]; then
    echo "$input directory does not exist!"
    exit
fi

if [ -f cnvkit_input.tsv ]; then
    rm cnvkit_input.tsv
fi

for file in `find $input -name duplicates_marked_no_table.tsv`; do 
	basepath=`realpath $file | sed 's/results.*//g'`; 
	
	cat $file >> cnvkit_input.tsv;
	sed -i 's,\.\/,'$basepath',g' cnvkit_input.tsv 
done
