from snakemake.io import expand, directory, touch

rule all:
    input:
        expand("results/final/{sample}_all.done", sample=config['samples'])

rule filter_reads:
    input: "resources/reads/{sample}/reads_all.fastq.gz"
    output: "results/reads/{sample}/reads_all.fastq.gz"
    threads: 2
    log: "results/logs/{sample}_filtlong.log"
    conda: "filtlong-env"
    params: config['min_read_len']
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule read_length_histogram:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/plots/{sample}_read_length_histogram.png"
    threads: 2
    log: "results/logs/{sample}_seqkit_length.log"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads}"

rule read_quality_histogram:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/reads_stat/{sample}_read_quality_histogram.png"
    threads: 2
    log: "results/logs/{sample}_seqkit_quality.log"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads} -f MeanQual"

rule fq2fasta:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/reads/{sample}/reads_all.fasta"
    threads: 2
    log: "results/logs/{sample}_seqkit_fq2fa.log"
    conda: "seqkit-env"
    shell: "seqkit fq2fa -j {threads} {input} 1> {output} 2> {log}"

rule create_fr_red:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/{sample}_fr_red.fa"
    threads: 2
    log: "results/logs/{sample}_seqkit_subseq_fr_red.log"
    conda: "seqkit-env"
    params: start=['fr_red_start'], end=['fr_red_end']
    shell: "seqkit subseq -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule create_fr_green:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/{sample}_fr_green.fa"
    threads: 2
    log: "results/logs/{sample}_seqkit_subseq_fr_green.log"
    conda: "seqkit-env"
    params: start=['fr_green_start'], end=['fr_green_end']
    shell: "seqkit subseq -j {threads} -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule blast_red:
    input: fr="results/flanking_regions/{sample}_fr_red.fa", rd="results/reads/{sample}/reads_all.fasta"
    output: "results/tables/{sample}_blast_red.tsv"
    params: fmt=config['format'], nalns=config['n_aligns']
    log: "results/logs/{sample}_blast_red.log"
    conda: "blast-env"
    shell: "blastn -query {input.fr} -subject {input.rd} -outfmt {params.fmt} -num_alignments {params.nalns} 1> {output} 2> {log}"

rule blast_green:
    input: fr="results/flanking_regions/{sample}_fr_green.fa", rd="results/reads/{sample}/reads_all.fasta"
    output: "results/tables/{sample}_blast_green.tsv"
    params: fmt=config['format'], nalns=config['n_aligns']
    log: "results/logs/{sample}_blast_green.log"
    conda: "blast-env"
    shell: "blastn -query {input.fr} -subject {input.rd} -outfmt {params.fmt} -num_alignments {params.nalns} 1> {output} 2> {log}"

rule final:
    input: blast_red="results/tables/{sample}_blast_red.tsv", 
            blast_green="results/tables/{sample}_blast_green.tsv",
            len_hist="results/plots/{sample}_read_length_histogram.png",
            qual_hist="results/reads_stat/{sample}_read_quality_histogram.png"

    output: touch("results/final/{sample}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow fininshed, no errors")
