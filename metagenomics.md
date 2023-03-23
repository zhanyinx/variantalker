# metagenomics

**1. Installation - directly through bioconda causes dependency cashes**

```bash
conda create --name mpa python=3.7
conda activate mpa
pip install metaphlan
conda install -c bioconda bowtie2
metaphlan --install
```

Additional information [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0).

**2. Get read files**

* Example SRA files from [here](https://www.hmpdacc.org/HMASM/)
* Raw fastq read files from clinical sample

**3. Getting clades**

```bash
metaphlan <INPUT_FASTA> --input_type fasta --nproc <NCORES> > <OUPUT_TXT>
```

**4. Raw clade analysis**

* `bowtie2out.txt` file - information on unique sequence markers. Shouldn't be of importance
* `.txt` file - organism abundances listed as one clade per line, tab-separated from the clade's percent abundance
  * `#` Header with metadata
  * Prefixes indicate taxonomic level: `Kingdom: k__, Phylum: p__, Class: c__, Order: o__, Family: f__, Genus: g__, Species: s__`
* Extraction with `grep s__ FILE | cut -f1 | sed 's/|/\n/g'

**5. Visual clade analysis**

* In py3.x environment
  * Clone `hclust2` [here](https://github.com/SegataLab/hclust2)
* Create py2.7 environment
  * Clone `export2graphlan` [here](https://github.com/SegataLab/export2graphlan)
  * Clone `graphlan` [here](https://github.com/biobakery/graphlan)
  * Copy `hclust2.py` file from cloned repo
* Heatmap
  * Create a species only abundance table and generate heatmap with `python hclust2/hclust2.py -i INPUT_SPECIES.txt -o FILE_PLOT.png` or just create a custom plot
* Cladogram
  * `tail -n +2 ABUNDANCE_TABLE.txt | cut -f1,3- > FILE.reformatted.txt`
  * Activate py2.7 environment
  * `python export2graphlan/export2graphlan.py --skip_rows 1 -i FILE.reformatted.txt --tree FILE.tree.txt --annotation FILE.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1`
  * Activate py3.x environment
  * `python graphlan/graphlan_annotate.py --annot FILE.annot.txt FILE.tree.txt merged_abundance.xml`
  * `python graphlan/graphlan.py --dpi 300 merged_abundance.xml merged_abundance.png --external_legends`

Additional information [here](https://github.com/biobakery/biobakery/wiki/metaphlan3).
