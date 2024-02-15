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
    params: start = config['fr_red_start'], end = config['fr_red_end']
    shell: "seqkit subseq -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule create_fr_green:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/{sample}_fr_green.fa"
    threads: 2
    log: "results/logs/{sample}_seqkit_subseq_fr_green.log"
    conda: "seqkit-env"
    params: start = config['fr_green_start'], end = config['fr_green_end']
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
    params: fmt = config['format'], nalns = config['n_fr_aligns']
    log: "results/logs/{sample}_blast_green.log"
    conda: "blast-env"
    shell: "blastn -query {input.fr} -subject {input.rd} -outfmt {params.fmt} -num_alignments {params.nalns} 1> {output} 2> {log}"

rule filter_flanking_regions_blast:
    input: script="workflow/scripts/filter_fr_hits.R", red = "results/tables/{sample}_blast_red.tsv", green = "results/tables/{sample}_blast_green.tsv"
    output: "results/tables/{sample}_blast_joined.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    conda: "rscripts-env"
    params: identity = config['min_identity'], e_val = config['max_e_value'], length = config['min_hit_len']
    shell: "Rscript {input.script} -r {input.red} -g {input.green} -i {params.identity} -e {params.e_val} -i {params.length} -o {output} &> {log}"

rule plot_distance_between_FRs:
    input: script="workflow/scripts/plot_fr_distances.R", blast="results/tables/{sample}_blast_joined.tsv"
    output: "results/plots/{sample}_FR_distances.png"
    log: "results/logs/{sample}_FR_distances.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.blast} -o {output} &> {log}"

rule plot_read_length_FR_distance:
    input: script="workflow/scripts/plot_read_length_and_FR_distances.R", blast="results/tables/{sample}_blast_joined.tsv", reads="results/reads/{sample}/reads_all.fasta"
    output: "results/plots/{sample}_reads_FR_distances.png"
    log: "results/logs/{sample}_reads_FR_distances.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -b {input.blast} -r {input.reads} -o {output} &> {log}"

rule blast_blaSHV:
    input: query="resources/genes/blaSHV.fa", subject="results/reads/{sample}/reads_all.fasta"
    output: "results/tables/{sample}_blast_blaSHV.tsv"
    log: "results/logs/{sample}_blast_blaSHV.log"
    conda: "blast-env"
    params: alns=config["n_bla_aligns"], fmt=config["format"]
    shell: "blast -query {input.query} -subject {input.subject} -outfmt {params.fmt} -num_alignments {params.alns} 1> {output} 2> {log}"

rule filter_blaSHV_hits:
    input: script="workflow/scripts/filter_blaSHV_blast.R", bla="results/tables/{sample}_blast_blaSHV.tsv", fr="results/tables/{sample}_blast_joined.tsv"
    output: "results/tables/{sample}_blast_blaSHV_filtered.tsv"
    log: "results/logs/{sample}_blast_blaSHV_filter.log"
    conda: "rscripts-env"
    params: e_val=config["max_e_value"]
    shell: "Rscript {input.script} -b {input.bla} -f {input.fr} -e {params.e_val} -o {output} $> {log}"

rule plot_blaSHV_counts:
    input: script="workflow/scripts/plot_blaSHV_counts.R", bla="results/tables/{sample}_blast_blaSHV_filtered.tsv"
    output: "results/plots/{sample}_blaSHV_counts.png"
    log: "results/logs/{sample}_blaSHV_counts.log"
    params: length=config["bla_len"]
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.bla} -l {params.length} -o {output} $> {log}"

rule final:
    input: blast_red="results/tables/{sample}_blast_red.tsv", 
            blast_green="results/tables/{sample}_blast_green.tsv",
            len_hist="results/plots/{sample}_read_length_histogram.png",
            qual_hist="results/reads_stat/{sample}_read_quality_histogram.png"

    output: touch("results/final/{sample}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow fininshed, no errors")
