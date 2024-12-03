from snakemake.io import expand, directory, touch, temp

rule all:
    input:
        expand("results/final/{sample}_all.done", sample=config['samples']), "results/tables/aggregate/frequencies_full_table.tsv"

rule merge_reads:
    input: "resources/reads_separate/{sample}"
    output: temp("resources/reads/{sample}/reads_all.fastq.gz")
    threads: 18
    log: "results/logs/{sample}_zcat.log"
    benchmark: "results/benchmarks/zcat_reads/{sample}.tsv"
    conda: "compress-env"
    container: "containers/compress.sif"
    shell: "zcat {input}/*.gz | pigz -c -p {threads} 1> {output} 2> {log}"

rule filter_reads:
    input: "resources/reads/{sample}/reads_all.fastq.gz"
    output: temp("results/reads/{sample}/reads_filtered.fastq.gz")
    threads: 18
    log: "results/logs/{sample}_filtlong.log"
    benchmark: "results/benchmarks/filtlong_filter/{sample}.tsv"
    conda: "filtlong-env"
    container: "containers/filtlong.sif"
    params: min_len=config['min_read_len']
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule fq2fasta:
    input: "results/reads/{sample}/reads_filtered.fastq.gz"
    output: "results/reads/{sample}/reads_filtered.fasta.gz"
    threads: 18
    log: "results/logs/{sample}_seqkit_fq2fa.log"
    benchmark: "results/benchmarks/seqkit_convert/{sample}.tsv"
    conda: "seqkit-env"
    container: "containers/seqkit.sif"
    shell: "seqkit fq2fa -j {threads} {input} | pigz -c -p {threads} 1> {output} 2> {log}"

rule create_fr_red:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/fr_red.fa" # the same for every sample
    threads: 18
    log: "results/logs/seqkit_subseq_fr_red.log"
    conda: "seqkit-env"
    container: "containers/seqkit.sif"
    params: start = config['fr_red_start'], end = config['fr_red_end']
    shell: "seqkit subseq -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule create_fr_green:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/fr_green.fa" # the same for every sample
    threads: 18
    log: "results/logs/seqkit_subseq_fr_green.log"
    conda: "seqkit-env"
    container: "containers/seqkit.sif"
    params: start = config['fr_green_start'], end = config['fr_green_end']
    shell: "seqkit subseq -j {threads} -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule make_blast_db:
    input: "results/reads/{sample}/reads_filtered.fasta.gz"
    output: directory("results/blast_databases/{sample}")
    threads: 18
    log: "results/logs/{sample}_blastdb.log"
    conda: "blast-env"
    container: "containers/blast.sif"
    shell: "pigz -c -d -p {threads} {input} | makeblastdb -in - -dbtype nucl -title blastdb -out {output}/blastdb &> {log}"

rule blast_red:
    input: query="results/flanking_regions/fr_red.fa", database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_red.tsv"
    threads: 18
    params: fmt=config['format'], n_alns=config['n_fr_aligns']
    log: "results/logs/{sample}_blast_red.log"
    benchmark: "results/benchmarks/blast_red/{sample}.tsv"
    conda: "blast-env"
    container: "containers/blast.sif"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule blast_green:
    input: query="results/flanking_regions/fr_green.fa", database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_green.tsv"
    threads: 18
    params: fmt = config['format'], n_alns = config['n_fr_aligns']
    log: "results/logs/{sample}_blast_green.log"
    benchmark: "results/benchmarks/blast_green/{sample}.tsv"
    conda: "blast-env"
    container: "containers/blast.sif"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule create_repeat_unit:
    input: "resources/plasmid/DA61218_plasmid.fa"
    output: "results/flanking_regions/repeat_unit.fa"
    log: "results/logs/seqkit_repunit.log"
    params: start=config["ru_start"], end=config["ru_end"]
    conda: "seqkit-env"
    container: "containers/seqkit.sif"
    shell:  "seqkit subseq -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule blast_repeat_unit:
    input: query="results/flanking_regions/repeat_unit.fa",
           database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_repeat_unit.tsv"
    threads: 18
    log: "results/logs/{sample}_blast_repeat_unit.log"
    conda: "blast-env"
    container: "containers/blast.sif"
    params: fmt=config["format"], n_alns=config["n_bla_aligns"]
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule filter_rr_ru_blast:
    input: red = "results/tables/{sample}/blast_red.tsv",
           repunit = "results/tables/{sample}/blast_repeat_unit.tsv"
    output: "results/tables/{sample}/blast_joined_red_repunit.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    params: identity = config['min_identity'], e_val = config['max_e_value'],
            length_fr = config['min_fr_len'], length_ru = config['min_ru_len'],
            distance = config['max_dist']
    script: "scripts/filter_red_repunit.R"

# filter RED+RU+Oriented length
rule filter_min_orient_length:
    input: "results/tables/{sample}/blast_joined_red_repunit.tsv"
    output: "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv"
    log: "results/logs/{sample}_filt_min_orient_len.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    params: base_len = config['base_len']
    script: "scripts/filter_1820.R"

# get counts of read possibly containing various CNs
rule cn_reads_bins:
    input: table = "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv",
           reads = "results/reads/{sample}/reads_filtered.fasta.gz"
    output: "results/tables/{sample}/number_reads_containing_CN.tsv"
    log: "results/logs/{sample}_cn_reads_bins.log"
    conda: "biostrings-env"
    container: "containers/biostrings.sif"
    params: max_cn = config['max_cn'], incr = config['increment'], base_len = config['base_len']
    shell: "scripts/get_cn_read_counts.R"

# filter GREEN
rule filter_flanking_regions:
    input: red_ru = "results/tables/{sample}/blast_joined_red_repunit_orient_len.tsv",
           green = "results/tables/{sample}/blast_green.tsv"
    output: "results/tables/{sample}/blast_joined.tsv"
    log: "results/logs/{sample}_blast_joined.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    params: identity = config['min_identity'], e_val = config['max_e_value'], length = config['min_fr_len']
    script: "scripts/filter_fr_hits.R"

rule blast_blaSHV:
    input: query="resources/genes/blaSHV.fa",
           database="results/blast_databases/{sample}"
    output: "results/tables/{sample}/blast_blaSHV.tsv"
    threads: 18
    log: "results/logs/{sample}_blast_blaSHV.log"
    benchmark: "results/benchmarks/blast_blaSHV/{sample}.tsv"
    conda: "blast-env"
    container: "containers/blast.sif"
    params: fmt=config["format"], n_alns=config["n_bla_aligns"]
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule filter_blaSHV_hits:
    input: bla="results/tables/{sample}/blast_blaSHV.tsv",
           fr="results/tables/{sample}/blast_joined.tsv"
    output: "results/tables/{sample}/blast_blaSHV_filtered.tsv"
    log: "results/logs/{sample}_blast_blaSHV_filter.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    params: e_val=config["max_e_value"]
    script: "scripts/filter_blaSHV_blast.R"

rule make_bed_blaSHV_filtered:
    input: "results/tables/{sample}/blast_blaSHV_filtered.tsv"
    output: "results/bedfiles/{sample}/blaSHV_hits.bed"
    log: "results/logs/{sample}_make_bed.log"
    benchmark: "results/benchmarks/make_bed/{sample}.tsv"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    script: "scripts/make_bed.R"

rule merge_blaSHV_filtered:
    input: "results/bedfiles/{sample}/blaSHV_hits.bed"
    output: sorted="results/bedfiles/{sample}/blaSHV_hits_sorted.bed",
            merged="results/bedfiles/{sample}/blaSHV_hits_merged.bed"
    log: "results/logs/{sample}_bedtools_merge.log"
    benchmark: "results/benchmarks/bedtools_merge/{sample}.tsv"
    conda: "bedtools-env"
    container: "containers/bedtools.sif"
    params: dist=config["dist"]
    shell: "sort -k1,1 -k2,2n {input} > {output.sorted} && "
           "bedtools merge -i {output.sorted} -s -d {params.dist} > {output.merged} 2> {log}"

rule blaSHV_counts:
    input: script="workflow/scripts/plot_blaSHV_counts_merged.R",
           bed="results/bedfiles/{sample}/blaSHV_hits_merged.bed",
           blast="results/tables/{sample}/blast_joined.tsv"
    output: plot="results/plots/{sample}/blaSHV_merged_counts.png",
            table="results/tables/{sample}/blaSHV_counts.tsv"
    log: "results/logs/{sample}_blaSHV_counts.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    params: length=config["bla_len"]
    shell: "Rscript {input.script} -i {input.bed} -b {input.blast} -l {params.length} -p {output.plot} -a {output.table} &> {log}"

rule frequency_calculation:
    input: bins="results/tables/{sample}/number_reads_containing_CN.tsv",
           bla="results/tables/{sample}/blaSHV_counts.tsv"
    output: "results/tables/{sample}/frequencies.tsv"
    log: "results/logs/{sample}_frequencies.log"
    conda: "rscripts-env"
    container: "containers/rscripts.sif"
    shell: "scripts/final_calculations.R"

rule aggregate_freq_tables:
    input: expand("results/tables/{sample}/frequencies.tsv", sample=config['samples'])
    output: tsv = "results/tables/aggregate/frequencies_full_table.tsv",
            xlsx = "results/tables/aggregate/frequencies_full_table.xlsx"
    log: "results/logs/aggregate_freq_tables.log"
    conda: "pandas-env"
    container: "containers/pandas.sif"
    script: "scripts/aggregate_frequency_tables.py"

rule final:
    input: freqs="results/tables/{sample}/frequencies.tsv"
    output: touch("results/final/{sample}_all.done")
    log: "results/logs/{sample}_final.log"
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
