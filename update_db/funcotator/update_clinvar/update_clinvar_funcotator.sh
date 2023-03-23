#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir clinvar
mkdir clinvar/hg19
mkdir clinvar/hg38

cd clinvar/hg38
	# download data
	wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

	# get date of the db
	date=`zcat clinvar.vcf.gz | grep "fileDate" | cut -d "=" -f2- | sed 's/-//g'`
	filename=clinvar_$date.vcf

	# convert to UCSC annotation of chromosomes
	gunzip clinvar.vcf.gz
	sed -i '/\tYT\t/d' clinvar.vcf 
	awk '{if($0~/^#/) print $0; else print "chr"$0}' clinvar.vcf > appo
	mv appo clinvar.vcf
	mv clinvar.vcf clinvar_$date.vcf

	# update config file
	cp $SCRIPT_DIR/config/clinvar_vcf.config .
	sed -i 's/DATE/'$date'/g' clinvar_vcf.config

	#index vcf
	gatk IndexFeatureFile -I clinvar_$date.vcf

cd ../../clinvar/hg19

	# same as above but for hg19
	wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
	date=`zcat clinvar.vcf.gz | grep "fileDate" | cut -d "=" -f2- | sed 's/-//g'`
	filename=clinvar_$date.vcf

	gunzip clinvar.vcf.gz
	sed -i '/\tYT\t/d' clinvar.vcf 
	awk '{if($0~/^#/) print $0; else print "chr"$0}' clinvar.vcf > appo
	mv appo clinvar.vcf
	mv clinvar.vcf clinvar_$date.vcf
	cp $SCRIPT_DIR/config/clinvar_vcf.config clinvar_vcf.config

	sed -i 's/DATE/'$date'/g' clinvar_vcf.config
	gatk IndexFeatureFile -I clinvar_$date.vcf

cd ../../
