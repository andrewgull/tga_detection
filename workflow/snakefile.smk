from snakemake.io import expand, directory

config = "config.yaml"
flanking_region = [1, 2]

rule all:
    input:
        expand("results/CNV/{replicate}_cnv_histogram.pdf", replicate=config['replicates'])

rule filter_reads:
    input:
        "resources/reads/{replicate}/CNV_reads_all.fastq.gz"
    output:
        "results/reads/{replicate}/CNV_reads_all.fastq.gz"
    params: min_len=3500
    threads: 10
    log: "results/logs/{replicate}_filtering.log"
    conda: "envs/filtlong.yaml"
    shell:
        "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule fq2fasta:
    input:
        "results/reads/{replicate}/CNV_reads_all.fastq.gz"
    output:
        "results/reads/{replicate}/CNV_reads_all.fasta"
    threads: 10
    log: "results/logs/{replicate}_converting.log"
    conda: "envs/seqkit.yaml"
    shell:
        "seqkit fq2fa -j {threads} {input} 1> {output} 2> {log}"

rule blast_db:
    input:
        "results/reads/{replicate}/CNV_reads_all.fasta"
    output:
        directory("results/blastdb_{replicate}")
    log: "results/logs/{replicate}_makeblastdb.log"
    conda: "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl -out {output}/FAU50052 -logfile {log}"

rule blast_fr1:
    input: 
        fr="resources/flanking_regions/region_1.fasta",
        db="results/blastdb_{replicate}"
    output: 
        "results/tables/{replicate}_FR1_blast.tsv"
    params: fmt=6, nalns=1000000
    threads: 10
    log: "results/logs/{replicate}_blast_fr1.log"
    conda: "envs/blast.yaml"
    shell:
        "blastn -query {input.fr} -db {input.db}/FAU50052 -outfmt {params.fmt} -num_threads {threads} "
        "-num_alignments {params.nalns} 1> {output} 2> {log}"

rule blast_fr2:
    input:
        fr="resources/flanking_regions/region_2.fasta",
        db="results/blastdb_{replicate}"
    output:
        "results/tables/{replicate}_FR2_blast.tsv"
    params: fmt=6, nalns=1000000
    threads: 10
    log: "results/logs/{replicate}_blast_fr2.log"
    conda: "envs/blast.yaml"
    shell:
        "blastn -query {input.fr} -db {input.db}/FAU50052 -outfmt {params.fmt} -num_threads {threads} "
        "-num_alignments {params.nalns} 1> {output} 2> {log}"

rule make_histogram:
    input:
        script="workflow/scripts/parse_blast.R",
        tables=expand("results/tables/{replicate}_FR{fr}_blast.tsv", fr=flanking_region)
    output: 
        hist="results/CNV/{replicate}_cnv_histogram.pdf",
        table="results/CNV/{replicate}_cnv_counts.csv"
    params: hit_len=250, unit_len=3500
    log: "results/logs/{replicate}_parse_tables.log"
    conda: "envs/rscripts.yaml"
    shell:
        "Rscript {input.script} -i {input.tables[0]} -e {input.tables[1]} -s {output.hist} -t {output.table} -l {params.hit_len} -u {params.unit_len} &> {log}"
