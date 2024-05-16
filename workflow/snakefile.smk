from snakemake.io import expand, directory, touch, temp

rule all:
    input:
        expand("results/final/{sample}_all.done", sample=config['samples'])

rule merge_reads:
    input: "resources/reads_separate/{sample}"
    output: temp("resources/reads/{sample}/reads_all.fastq.gz")
    threads: 10
    log: "results/logs/{sample}_zcat.log"
    benchmark: "results/benchmarks/zcat_reads/{sample}.tsv"
    shell: "zcat {input}/*.gz | pigz -c -p {threads} 1> {output} 2> {log}"

# rule chop_reads:
#     input: "resources/reads/{sample}/reads_all.fastq.gz"
#     output: temp("results/reads/{sample}/reads_chopped.fastq.gz")
#     threads: 10
#     log: "results/logs/{sample}_porechop.log"
#     benchmark: "results/benchmarks/porechop_reads/{sample}.tsv"
#     conda: "porechop-env"
#     shell: "porechop -i {input} -o {output} --threads {threads} &> {log}"

rule filter_reads:
    input: "resources/reads/{sample}/reads_all.fastq.gz"
    output: temp("results/reads/{sample}/reads_filtered.fastq.gz")
    threads: 10
    log: "results/logs/{sample}_filtlong.log"
    benchmark: "results/benchmarks/filtlong_filter/{sample}.tsv"
    conda: "filtlong-env"
    params: min_len=config['min_read_len']
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule read_length_histogram:
    input: "results/reads/{sample}/reads_filtered.fastq.gz"
    output: "results/plots/{sample}/reads_length_histogram.png"
    threads: 10
    log: "results/logs/{sample}_seqkit_length.log"
    benchmark: "results/benchmarks/seqkit_length_filter/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads} &> {log}"

rule read_quality_histogram:
    input: "results/reads/{sample}/reads_filtered.fastq.gz"
    output: "results/plots/{sample}/reads_quality_histogram.png"
    threads: 10
    log: "results/logs/{sample}_seqkit_quality.log"
    benchmark: "results/benchmarks/seqkit_quality/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit watch {input} -O {output} -j {threads} -f MeanQual &> {log}"

rule fq2fasta:
    input: "results/reads/{sample}/reads_filtered.fastq.gz"
    output: "results/reads/{sample}/reads_filtered.fasta.gz"
    threads: 10
    log: "results/logs/{sample}_seqkit_fq2fa.log"
    benchmark: "results/benchmarks/seqkit_convert/{sample}.tsv"
    conda: "seqkit-env"
    shell: "seqkit fq2fa -j {threads} {input} | pigz -c -p {threads} 1> {output} 2> {log}"

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

rule make_blast_db:
    input: "results/reads/{sample}/reads_filtered.fasta.gz"
    output: directory("results/blast_databases/{sample}")
    threads: 10
    log: "results/logs/{sample}_blastdb.log"
    conda: "blast-env"
    shell: "pigz -p {threads} {input} | makeblastdb -in stdin -dbtype nucl -out {output}/blastdb &> {log}"

rule blast_red:
    input: query="results/flanking_regions/fr_red.fa", database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_red.tsv"
    threads: 10
    params: fmt=config['format'], n_alns=config['n_fr_aligns']
    log: "results/logs/{sample}_blast_red.log"
    benchmark: "results/benchmarks/blast_red/{sample}.tsv"
    conda: "blast-env"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule blast_green:
    input: query="results/flanking_regions/fr_green.fa", database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_green.tsv"
    threads: 10
    params: fmt = config['format'], n_alns = config['n_fr_aligns']
    log: "results/logs/{sample}_blast_green.log"
    benchmark: "results/benchmarks/blast_green/{sample}.tsv"
    conda: "blast-env"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule blast_repeat_unit:
    input: query="resources/flanking_regions/repeat_unit.fa",
           database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_repeat_unit.tsv"
    threads: 10
    log: "results/logs/{sample}_blast_repeat_unit.log"
    benchmark: "results/benchmarks/blast_repeat_unit/{sample}.tsv"
    conda: "blast-env"
    params: fmt=config["format"], n_alns=config["n_bla_aligns"]
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule filter_red_and_repeat_unit_blast:
    input: script="workflow/scripts/filter_red_repunit.R",
           red = "results/tables/{sample}/blast_red.tsv",
           repunit = "results/tables/{sample}/blast_repeat_unit.tsv"
    output: "results/tables/{sample}/blast_joined_red_repunit.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    benchmark: "results/benchmarks/filter_fr/{sample}.tsv"
    conda: "rscripts-env"
    params: identity = config['min_identity'], e_val = config['max_e_value'],
            length_fr = config['min_fr_len'], length_ru = config['min_ru_len'],
            distance = config['max_dist']
    shell: "Rscript {input.script} -r {input.red} -u {input.repunit} -i {params.identity} -e {params.e_val} "
           "-l {params.length_fr} -k {params.length_ru} -d {params.distance} -o {output} &> {log}"

# filter RED+RU+Oriented length
rule filter_min_orient_length:
    input: script = "workflow/scripts/filter_1820.R",
           table = "results/tables/{sample}/blast_joined_red_repunit.tsv"
    output: "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv"
    log: "results/logs/{sample}_filt_min_orient_len.log"
    conda: "rscripts-env"
    params: orient_dist = config['orient_dist']
    shell: "Rscript {input.script} -i {input.table} -d {params.orient_dist} -o {output} &> {log}"

# get read lengths which passed the previous filtering step
rule cn_reads_bins:
    input: script = "workflow/scripts/get_cn_read_counts.R",
           table = "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv",
           reads = "results/reads/{sample}/reads_filtered.fasta.gz"
    output: "results/tables/{sample}/number_reads_containing_CN.tsv"
    log: "results/logs/{sample}_filt_read_len.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -d {input.table} -r {input.reads} -o {output} &> {log}"

# filter GREEN
rule filter_flanking_regions:
    input: script="workflow/scripts/filter_fr_hits.R",
           red_ru = "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv",
           green = "results/tables/{sample}/blast_green.tsv"
    output: "results/tables/{sample}/blast_joined.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    conda: "rscripts-env"
    params: identity = config['min_identity'], e_val = config['max_e_value'], length = config['min_fr_len']
    shell: "Rscript {input.script} -r {input.red_ru} -g {input.green} -i {params.identity} -e {params.e_val} "
           "-l {params.length} -o {output} &> {log}"

rule blast_blaSHV:
    input: query="resources/genes/blaSHV.fa",
           database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_blaSHV.tsv"
    threads: 10
    log: "results/logs/{sample}_blast_blaSHV.log"
    benchmark: "results/benchmarks/blast_blaSHV/{sample}.tsv"
    conda: "blast-env"
    params: fmt=config["format"], n_alns=config["n_bla_aligns"]
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

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

rule make_bed_blaSHV_filtered:
    input: script="workflow/scripts/make_bed.R",
           bla="results/tables/{sample}/blast_blaSHV_filtered.tsv"
    output: "results/bedfiles/{sample}/blaSHV_hits.bed"
    log: "results/logs/{sample}_make_bed.log"
    benchmark: "results/benchmarks/make_bed/{sample}.tsv"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.bla} -o {output} &> {log}"

rule merge_blaSHV_filtered:
    input: "results/bedfiles/{sample}/blaSHV_hits.bed"
    output: sorted="results/bedfiles/{sample}/blaSHV_hits_sorted.bed",
            merged="results/bedfiles/{sample}/blaSHV_hits_merged.bed"
    log: "results/logs/{sample}_bedtools_merge.log"
    benchmark: "results/benchmarks/bedtools_merge/{sample}.tsv"
    conda: "varcalling-env"
    params: dist=config["dist"]
    shell: "sort -k1,1 -k2,2n {input} > {output.sorted} && bedtools merge -i {output.sorted} -s -d {params.dist} > {output.merged} 2> {log}"

rule blaSHV_counts:
    input: script="workflow/scripts/plot_blaSHV_counts_merged.R",
           bed="results/bedfiles/{sample}/blaSHV_hits_merged.bed",
           blast="results/tables/{sample}/blast_joined.tsv"
    output: plot="results/plots/{sample}/blaSHV_merged_counts.png",
            table="results/tables/{sample}/blaSHV_counts.tsv"
    log: "results/logs/{sample}_blaSHV_counts.log"
    conda: "rscripts-env"
    params: length=config["bla_len"]
    shell: "Rscript {input.script} -i {input.bed} -b {input.blast} -l {params.length} -p {output.plot} -a {output.table} &> {log}"

rule make_bed_blaSHV_all:
    input: script="workflow/scripts/make_bed.R",
           bla="results/tables/{sample}/blast_blaSHV.tsv"
    output: "results/bedfiles/{sample}/blaSHV_hits_all.bed"
    log: "results/logs/{sample}_make_bed_all.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} -i {input.bla} -o {output} &> {log}"

rule merge_blaSHV_all:
    input: "results/bedfiles/{sample}/blaSHV_hits_all.bed"
    output: sorted="results/bedfiles/{sample}/blaSHV_hits_all_sorted.bed",
            merged="results/bedfiles/{sample}/blaSHV_hits_all_merged.bed"
    log: "results/logs/{sample}_bedtools_merge_all.log"
    conda: "varcalling-env"
    params: dist=config["dist"]
    shell: "sort -k1,1 -k2,2n {input} > {output.sorted} && bedtools merge -i {output.sorted} -s -d {params.dist} > {output.merged} 2> {log}"

rule blaSHV_counts_all:
    input: script="workflow/scripts/plot_save_blaSHV_counts_all.R",
           bed="results/bedfiles/{sample}/blaSHV_hits_all_merged.bed"
    output: table="results/tables/{sample}/blaSHV_counts_all.tsv",
            plot="results/plots/{sample}/blaSHV_counts_all.png"
    log: "results/logs/{sample}_plot_save_blaSHV_counts_all.log"
    conda: "rscripts-env"
    shell: "Rscript {input.script} --input {input.bed} --plot {output.plot} --table {output.table} &> {log}"

rule final:
    input:  len_hist="results/plots/{sample}/reads_length_histogram.png",
            qual_hist="results/plots/{sample}/reads_quality_histogram.png",
            blast_join="results/tables/{sample}/blast_joined.tsv",
            plot_dist="results/plots/{sample}/FR_distances.png",
            plot_len_dist="results/plots/{sample}/reads_FR_distances.png",
            read_cn="results/tables/{sample}/number_reads_containing_CN.tsv",
            bla_counts="results/plots/{sample}/blaSHV_counts.png",
            save_counts="results/tables/{sample}/blaSHV_counts.tsv",
            table_counts_all="results/tables/{sample}/blaSHV_counts_all.tsv"

    output: touch("results/final/{sample}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
