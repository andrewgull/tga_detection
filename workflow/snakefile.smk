from snakemake.io import expand, directory, touch

rule all:
    input:
        expand("results/final/{sample}_all.done", sample=config['samples'])

rule filter_reads:
    input: "resources/reads/{sample}/reads_all.fastq.gz"
    output: "results/reads/{sample}/reads_all.fastq.gz"
    threads: 10
    log: "results/logs/{sample}_filtlong.log"
    benchmark: "results/benchmarks/filtlong_filter/{sample}.tsv"
    conda: "filtlong-env"
    params: min_len=config['min_read_len']
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule read_length_histogram:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/plots/{sample}/reads_length_histogram.png"
    threads: 10
    log: "results/logs/{sample}_seqkit_length.log"
    benchmark: "results/benchmarks/seqkit_length_filter/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads} &> {log}"

rule read_quality_histogram:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/plots/{sample}/reads_quality_histogram.png"
    threads: 10
    log: "results/logs/{sample}_seqkit_quality.log"
    benchmark: "results/benchmarks/seqkit_quality/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads} -f MeanQual &> {log}"

rule fq2fasta:
    input: "results/reads/{sample}/reads_all.fastq.gz"
    output: "results/reads/{sample}/reads_all.fasta"
    threads: 10
    log: "results/logs/{sample}_seqkit_fq2fa.log"
    benchmark: "results/benchmarks/seqkit_convert/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit fq2fa -j {threads} {input} 1> {output} 2> {log}"

rule split_fasta:
    input: "results/reads/{sample}/reads_all.fasta"
    # here you may need a smarter splitting, e.g. depending on the file size
    output: directory("results/reads_split/{sample}")
    threads: 10
    log: "results/logs/{sample}_seqkit_split.log"
    benchmark: "results/benchmarks/seqkit_split/{sample}.tsv"
    conda: "seqkit-env"
    params: parts=config['parts']
    shell: "seqkit split {input} -j {threads} -O {output} -p {params.parts} --two-pass -U &> {log}"

rule create_fr_red:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/fr_red.fa" # the same for every sample
    threads: 10
    log: "results/logs/seqkit_subseq_fr_red.log"
    conda: "seqkit-env"
    params: start = config['fr_red_start'], end = config['fr_red_end']
    shell: "seqkit subseq -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule create_fr_green:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/fr_green.fa" # the same for every sample
    threads: 10
    log: "results/logs/seqkit_subseq_fr_green.log"
    conda: "seqkit-env"
    params: start = config['fr_green_start'], end = config['fr_green_end']
    shell: "seqkit subseq -j {threads} -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule blast_red:
    input: query="results/flanking_regions/fr_red.fa",
           subject_dir="results/reads_split/{sample}"
    output: "results/tables/{sample}/blast_red.tsv"
    params: fmt=config['format'], n_alns=config['n_fr_aligns']
    log: "results/logs/{sample}_blast_red.log"
    benchmark: "results/benchmarks/blast_red/{sample}.tsv"
    conda: "blast-env"
    shell: "for F in {input.subject_dir}/*.fasta; do blastn -query {input.query} -subject $F -outfmt {params.fmt} -num_alignments {params.n_alns} 1> {output} 2> {log}; done"

rule blast_green:
    input: query="results/flanking_regions/fr_green.fa",
           subject_dir="results/reads/{sample}"
    output: "results/tables/{sample}/blast_green.tsv"
    params: fmt = config['format'], n_alns = config['n_fr_aligns']
    log: "results/logs/{sample}_blast_green.log"
    benchmark: "results/benchmarks/blast_green/{sample}.tsv"
    conda: "blast-env"
    shell: "for F in {input.subject_dir}/*.fasta; do blastn -query {input.query} -subject $F -outfmt {params.fmt} -num_alignments {params.n_alns} 1> {output} 2> {log}; done"

rule filter_flanking_regions_blast:
    input: script="workflow/scripts/filter_fr_hits.R",
           red = "results/tables/{sample}/blast_red.tsv",
           green = "results/tables/{sample}/blast_green.tsv"
    output: "results/tables/{sample}/blast_joined.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    benchmark: "results/benchmarks/filter_fr/{sample}.tsv"
    conda: "rscripts-env"
    params: identity = config['min_identity'], e_val = config['max_e_value'], length = config['min_hit_len']
    shell: "Rscript {input.script} -r {input.red} -g {input.green} -i {params.identity} -e {params.e_val} -l {params.length} -o {output} &> {log}"

rule plot_distance_between_FRs:
    input: script="workflow/scripts/plot_fr_distances.R",
           blast="results/tables/{sample}/blast_joined.tsv"
    output: "results/plots/{sample}/FR_distances.png"
    log: "results/logs/{sample}_FR_distances.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.blast} -o {output} &> {log}"

rule plot_read_length_FR_distance:
    input: script="workflow/scripts/plot_read_length_and_FR_distances.R",
           blast="results/tables/{sample}/blast_joined.tsv",
           reads="results/reads/{sample}/reads_all.fasta"
    output: "results/plots/{sample}/reads_FR_distances.png"
    log: "results/logs/{sample}_reads_FR_distances.log"
    conda: "biostrings-env"
    shell: "Rscript {input.script} -b {input.blast} -r {input.reads} -o {output} &> {log}"

rule blast_blaSHV:
    input: query="resources/genes/blaSHV.fa",
           subject_dir="results/reads_split/{sample}"
    output: "results/tables/{sample}/blast_blaSHV.tsv"
    log: "results/logs/{sample}_blast_blaSHV.log"
    benchmark: "results/benchmarks/blast_blaSHV/{sample}.tsv"
    conda: "blast-env"
    params: alns=config["n_bla_aligns"], fmt=config["format"]
    shell: "for F in {input.subject_dir}/*.fasta; do blastn -query {input.query} -subject $F -outfmt {params.fmt} -num_alignments {params.alns} 1> {output} 2> {log}; done"

rule filter_blaSHV_hits:
    input: script="workflow/scripts/filter_blaSHV_blast.R",
           bla="results/tables/{sample}/blast_blaSHV.tsv",
           fr="results/tables/{sample}/blast_joined.tsv"
    output: "results/tables/{sample}/blast_blaSHV_filtered.tsv"
    log: "results/logs/{sample}_blast_blaSHV_filter.log"
    conda: "rscripts-env"
    params: e_val=config["max_e_value"]
    shell: "Rscript {input.script} -b {input.bla} -f {input.fr} -e {params.e_val} -o {output} &> {log}"

rule plot_blaSHV_counts:
    input: script="workflow/scripts/plot_blaSHV_counts.R",
           bla="results/tables/{sample}/blast_blaSHV_filtered.tsv"
    output: "results/plots/{sample}/blaSHV_counts.png"
    log: "results/logs/{sample}_blaSHV_counts.log"
    params: length=config["bla_len"]
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.bla} -l {params.length} -o {output} &> {log}"

rule make_bed:
    input: script="workflow/scripts/make_bed.R",
           bla="results/tables/{sample}/blast_blaSHV_filtered.tsv"
    output: "results/bedfiles/{sample}/blaSHV_hits.bed"
    log: "results/logs/{sample}_make_bed.log"
    benchmark: "results/benchmarks/make_bed/{sample}.tsv"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.bla} -o {output} &> {log}"

rule cluster_blaSHV_hits:
    input: "results/bedfiles/{sample}/blaSHV_hits.bed"
    output: sorted="results/bedfiles/{sample}/blaSHV_hits_sorted.bed",
            merged="results/bedfiles/{sample}/blaSHV_hits_merged.bed"
    log: "results/logs/{sample}_bedtools_merge.log"
    benchmark: "results/benchmarks/bedtools_merge/{sample}.tsv"
    conda: "varcalling-env"
    params: dist=config["dist"]
    shell: "sort -k1,1 -k2,2n {input} > {output.sorted} && bedtools merge -i {output.sorted} -s -d {params.dist} > {output.merged} 2> {log}"

rule plot_bla_clusters:
    input: script="workflow/scripts/plot_blaSHV_counts_merged.R",
           bed="results/bedfiles/{sample}/blaSHV_hits_merged.bed",
           blast="results/tables/{sample}/blast_joined.tsv"
    output: "results/plots/{sample}/blaSHV_merged_counts.png"
    log: "results/logs/{sample}_blaSHV_counts.log"
    conda: "rscripts-env"
    params: length=config["bla_len"]
    shell: "Rscript {input.script} -i {input.bed} -b {input.blast} -l {params.length} -o {output} &> {log}"

rule final:
    input:  len_hist="results/plots/{sample}/reads_length_histogram.png",
            qual_hist="results/plots/{sample}/reads_quality_histogram.png",
            blast_join="results/tables/{sample}/blast_joined.tsv",
            plot_dist="results/plots/{sample}/FR_distances.png",
            plot_len_dist="results/plots/{sample}/reads_FR_distances.png",
            bla_counts="results/plots/{sample}/blaSHV_counts.png",
            clusters="results/bedfiles/{sample}/blaSHV_hits_merged.bed",
            plot_clust="results/plots/{sample}/blaSHV_merged_counts.png"
    output: touch("results/final/{sample}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow fininshed, no errors")
