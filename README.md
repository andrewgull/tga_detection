# Description

This is the code for detecting gene amplifications in ultra-deep Nanopore long read sequencing data at frequencies as low as $10^{-5}$.
Full study and method description can be found [in this publication](link).

The pipeline is written using [Snakemake](https://snakemake.github.io/) v8.23.1.

# Dependencies

Given that you have Snakemake installed, all other dependencies will be installed automatically when you run the pipeline.

List of dependencies:

 - [Snakemake](https://snakemake.readthedocs.io/en/stable/) v8.23.1
 - [R](https://www.r-project.org/) v4.4.0
 - [dplyr](https://dplyr.tidyverse.org/) v1.1.4
 - [readr](https://readr.tidyverse.org/) v2.1.5
 - [purrr](https://purrr.tidyverse.org/) v1.0.2
 - [tidyr](https://tidyr.tidyverse.org/) v1.3.1
 - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) v2.70.1
 - [Python](https://www.python.org/) v3.12.10
 - [pandas](https://pandas.pydata.org/) v1.5.3 
 - [OpenPyXL](https://openpyxl.readthedocs.io/en/stable/) v3.1.5
 - [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/) v2.12.0
 - [SeqKit](https://bioinf.shenwei.me/seqkit/) v2.0.0
 - [Filtlong](https://github.com/rrwick/Filtlong) v0.2.1
 - [bedtools](https://bedtools.readthedocs.io/en/latest/) v2.30.0
 - [gzip](https://www.gzip.org/) v1.10
 - [pigz](https://zlib.net/pigz/) v2.6


# Input & output

All the input files are described in `config.yaml`:

- a TSV file with sample names and their corresponding FASTQ files;
- FASTA file with the plasmid sequence;
- FASTA file with the blaSHV gene sequence;
- the output table name (without extension);
- a TSV file with parameters for the analysis (see the section below)

# Usage

## Test run

If you want to test the pipeline with a small dataset, you can use the `test.yaml` configuration file.

```bash
snakemake --use-conda --cores <number_of_cores> --configfile test.yaml
```

## Full run

To run the pipeline on all the samples, just replace the config file name with `config.yaml`:

```bash
snakemake --use-conda --cores <number_of_cores> --configfile config.yaml
```

You might need to edit `config/samples.tsv` to include actual paths to the FASTQ files (default is `resources/reads/`)

# Parameters

The parameters for the analysis are specified in the `params.tsv` file. The parameters include:

- minimum read length: reads shorter than this are discarded from the analysis.
- fr_red_start and fr_red_end: the start and end positions of the flanking region 1.
- fr_green_start and fr_green_end: the start and end positions of the flanking region 2.
- ru_start and ru_end: the start and end positions of the repeat unit (IS element).
- bla_start and bla_end: the start and end positions of the blaSHV gene.
- format: blast output format.
- n_fr_aligns: the number of flanking region alignments to consider.
- n_bla_aligns: the number of blaSHV gene alignments to consider.
- min_fr_len: the minimum length of the flanking region to consider.
- min_identity: the minimum identity BLAST hits to consider.
- max_e_value: the maximum E-value of BLAST hits to consider.
- min_ru_len: the minimum length of the repeat unit to consider.
- max_dist: the maximum distance between the BLAST hits to consider.
- dist_to_end: the distance from the end of the read to the end of the repeat unit to consider.
- max_cn: the maximum copy number of the repeat unit to consider.
- increment: increase in length of the DNA segment with each new blaSHV copy.
- base_len: length of a blaSHV gene for expected  copy number calculation.
- dist: the distance between BLAST hits to use in `bedtools merge`.

# Rule graph

![DAG](images/rulegraph.png)
