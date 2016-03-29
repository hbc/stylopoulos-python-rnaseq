This is a re-analysis of some old data looking at changes in the gut of the python before and after feeding. The original authors of the paper
didn't make some of the processed data available, but the raw data is up so we are regenerate it. We can also improve on their analysis,
the metadata I saw in the SRA project description described some suboptimal steps. The project on SRA is SRP051827. I think we also get
a side benefit of having a much nicer annotated genome than they had.

1. First downloaded the data from SRA:

```r
source('http://bioconductor.org/biocLite.R')
biocLite('SRAdb')
library(SRAdb)
biocLite('DBI')
library(DBI)
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(), srafile)
getSRAfile('SRP051827', con)
```

2. Dump out FASTQ files from the stupid SRA files.

```bash
for file in *.sra
do
    fastq-dump $file
done
```

3. Convert the Excel file about which file is what to something we can more easily use. I did this by making the Excel
file a simpler format and exporting to CSV. The original file, the simplified file and the CSV are in the metadata directory.
I ended up joining this with the metadata file from SRA to match up the filenames to the samples. The result of this 
is in python-fixed.csv in the metadata directory.

4. Combine the technical replicates into one file for each sample:

```bash
bcbio_prepare_samples.py --csv python-fixed.csv --out merged
mv python-fixed-merged.csv ../metadata/python-paper-redo.csv
```

5. Set up the bcbio run:

```bash
bcbio_nextgen.py -w template ./fastrnaseq-template.yaml metadata/python-paper-redo.csv data/merged/
```

6. Run bcbio:

```bash
bcbio_nextgen.py --timeout 6000 --tag python -t ipython -s slurm -q general -n 96 /n/regal/hsph_bioinfo/bcbio_nextgen/galaxy/bcbio_system.yaml ../config/python-paper-redo.yaml
```

