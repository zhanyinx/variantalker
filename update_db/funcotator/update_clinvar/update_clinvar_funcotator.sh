#!/bin/bash

# Stop at the first real failure instead of grinding on. Defence in depth only: the authority on
# whether the live database is replaced is the Python-side validation in `update_clinvar`, because
# no shell-side fix can reach that decision -- it is made after this script has exited.
#
# `pipefail` is the load-bearing half, and `set -e` alone would not have caught the bug this
# script is famous for. The `date=` assignment below takes its status from the LAST command in its
# pipeline (`sed`, which succeeds), so under `set -e` alone a failed `zcat` sails past and leaves
# `$date` empty -- which is what produces a 0-byte `clinvar_.vcf` and a blank-version config that
# Funcotator loads as a valid data source with no records. `pipefail` makes the pipeline report
# `zcat`'s failure so `set -e` can act on it.
set -eo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Default build
BUILD="both"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --build)
            BUILD="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# `-p`, not a bare `mkdir`: with `set -e` above, a leftover `clinvar/` from a crashed run would
# otherwise abort every subsequent run before it downloaded anything. A leftover is not a reason to
# refuse -- the caller clears the stage before invoking this script, and a stale data file that
# survived both would be caught by the staged-database validation, which refuses two files matching
# the same pattern as an ambiguous database.
mkdir -p clinvar


if [[ "$BUILD" == "hg38" || "$BUILD" == "both" ]]; then
mkdir -p clinvar/hg38
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
cd ../..
fi

if [[ "$BUILD" == "hg19" || "$BUILD" == "both" ]]; then
mkdir -p clinvar/hg19
cd clinvar/hg19

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
cd ../..
fi
