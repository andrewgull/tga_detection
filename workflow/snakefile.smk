from snakemake.io import expand, multiext

#SAMPLES = ["FAU50052"]
#config = "config.yaml"
flanking_region = [1, 2]

rule all:
    input:
        "results/CNV/FAU50052_cnv_histogram.pdf"

rule filter:
    input:
        "resources/reads/FAU50052_pass_all.fastq.gz"
    output:
        "results/reads/FAU50052_pass_all_filt.fastq.gz"
    params: min_len=3500
    threads: 10
    log: "results/logs/filtering.log"
    conda: "envs/filtlong.yaml"
    shell:
        "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule convert:
    input:
        "results/reads/FAU50052_pass_all_filt.fastq.gz"
    output:
        "results/reads/FAU50052_pass_all_filt.fasta"
    threads: 10
    log: "results/logs/converting.log"
    conda: "envs/seqkit.yaml"
    shell:
        "seqkit fq2fa -j {threads} {input} 1> {output} 2> {log}"

rule blast_db:
    input:
        "results/reads/FAU50052_pass_all_filt.fasta"
    output:
        directory("results/blastdb")
    log: "results/logs/makeblastdb.log"
    conda: "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl -out {output}/FAU50052 -logfile {log}"

rule blast_fr1:
    input: 
        fr="resources/reads/flanking_region_1.fasta",
        db="results/blastdb"
    output: 
        "results/tables/FAU50052_FR1_blast.tsv"
    params: fmt=6, nalns=1000000
    threads: 10
    log: "results/logs/blast_fr1.log"
    conda: "envs/blast.yaml"
    shell:
        "blastn -query {input.fr} -db {input.db}/FAU50052 -outfmt {params.fmt} -num_threads {threads} "
        "-num_alignments {params.nalns} 1> {output} 2> {log}"

rule blast_fr2:
    input:
        fr="resources/reads/flanking_region_2.fasta",
        db="results/blastdb"
    output:
        "results/tables/FAU50052_FR2_blast.tsv"
    params: fmt=6, nalns=1000000
    threads: 10
    log: "results/logs/blast_fr2.log"
    conda: "envs/blast.yaml"
    shell:
        "blastn -query {input.fr} -db {input.db}/FAU50052 -outfmt {params.fmt} -num_threads {threads} "
        "-num_alignments {params.nalns} 1> {output} 2> {log}"

rule parse_tables:
    input:
        script="workflow/scripts/parse_blast.R",
        tables=expand("results/tables/FAU50052_FR{fr}_blast.tsv", fr=flanking_region)
    output: 
        hist="results/CNV/FAU50052_cnv_histogram.pdf",
        table="results/CNV/FAU50052_cnv_counts.csv"
    params: hit_len=250, unit_len=3500
    log: "results/logs/parse_tables.log"
    conda: "envs/rscripts.yaml"
    shell:
        "Rscript {input.script} -i {input.tables[0]} -e {input.tables[1]} -h {output.hist} -t {output.table} -l {params.hit_len} -u {params.unit_len} &> {log}"
