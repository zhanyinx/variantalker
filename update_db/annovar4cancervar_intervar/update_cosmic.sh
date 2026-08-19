#!/bin/bash

if ! [ $# -eq 3 ]; then
        echo "Usage: update_cosmic.sh version passcode output_folder"
        exit
fi

if ! [ -d $3 ]; then
        echo "Output folder does not exist, creating it.."
        mkdir -p $3
fi

version=$1
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v$version/CosmicMutantExport.tsv.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicMutantExport.tsv.gz

link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v$version/VCF/CosmicCodingMuts.vcf.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicCodingMuts.vcf.gz

link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v$version/VCF/CosmicNonCodingVariants.vcf.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicNonCodingVariants.vcf.gz

gunzip CosmicNonCodingVariants.vcf.gz CosmicCodingMuts.vcf.gz CosmicMutantExport.tsv.gz

mutantfile="CosmicMutantExport.tsv"
codingfile="CosmicCodingMuts.vcf"
noncodingfile="CosmicNonCodingVariants.vcf"


if ! [ -f $mutantfile ]; then
	echo "CosmicMutantExport.tsv does not exist. ($mutantfile)"
	exit
fi

if ! [ -f $codingfile ]; then
	echo "CosmicMutantExport.tsv does not exist. ($codingfile)"
	exit
fi

if ! [ -f $noncodingfile ]; then
	echo "CosmicMutantExport.tsv does not exist. ($noncodingfile)"
	exit
fi


perl ${SCRIPT_DIR}/utils/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicCodingMuts.vcf > hg38_cosmic$version.txt
perl ${SCRIPT_DIR}/utils/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicNonCodingVariants.vcf >> hg38_cosmic$version.txt

bedtools sort -i hg38_cosmic$version.txt > appo
mv appo hg38_cosmic$version.txt
perl ${SCRIPT_DIR}/utils/index_annovar.pl hg38_cosmic$version.txt 1000 > hg38_cosmic$version.txt.idx


rm $mutantfile $codingfile $noncodingfile


link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v$version/CosmicMutantExport.tsv.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicMutantExport.tsv.gz

link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v$version/VCF/CosmicCodingMuts.vcf.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicCodingMuts.vcf.gz

link=`curl -H "Authorization: Basic $2" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v$version/VCF/CosmicNonCodingVariants.vcf.gz`
url=`echo ${link::-2} | cut -d"\"" -f4-`
curl ${url} --output CosmicNonCodingVariants.vcf.gz

gunzip CosmicNonCodingVariants.vcf.gz CosmicCodingMuts.vcf.gz CosmicMutantExport.tsv.gz




if ! [ -f $mutantfile ]; then
        echo "CosmicMutantExport.tsv does not exist. ($mutantfile)"
        exit
fi

if ! [ -f $codingfile ]; then
        echo "CosmicMutantExport.tsv does not exist. ($codingfile)"
        exit
fi

if ! [ -f $noncodingfile ]; then
        echo "CosmicMutantExport.tsv does not exist. ($noncodingfile)"
        exit
fi


perl ${SCRIPT_DIR}/utils/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicCodingMuts.vcf > hg19_cosmic$version.txt
perl ${SCRIPT_DIR}/utils/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicNonCodingVariants.vcf >> hg19_cosmic$version.txt

bedtools sort -i hg19_cosmic$version.txt > appo
mv appo hg19_cosmic$version.txt
perl ${SCRIPT_DIR}/utils/index_annovar.pl hg19_cosmic$version.txt 1000 > hg19_cosmic$version.txt.idx


rm $mutantfile $codingfile $noncodingfile

mv hg19_cosmic$version.txt hg19_cosmic$version.txt.idx hg38_cosmic$version.txt hg38_cosmic$version.txt.idx $3/
