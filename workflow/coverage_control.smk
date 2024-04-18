# this smk file is for blasting all control regions (1-4) on all samples
# basically it repeats the main snakefile.smk but to a much smaller extent

from snakemake.io import expand, touch

configfile: "workflow/config.yaml"

rule all_regions:
    input:
        expand("results/coverage_control/final/{region}_{sample}_all.done", sample=config['samples'],
            region=config['control_regions'])

module main_workflow:
    snakefile:
        "snakefile.smk"
    #config: config["main-workflow"]

# import rules that you're going to use from the main file
use rule merge_reads, filter_reads, fq2fasta, make_blast_db from main_workflow

rule blast_region:
    input: query="results/control_regions/{region}.fa", database="results/blast_databases/{sample}"
    output: "results/tables/{region}/blast_{region}_{sample}.tsv"
    threads: 10
    params: fmt=config['format'], n_alns=config['n_fr_aligns']
    log: "results/logs/blast_{region}_{sample}.log"
    benchmark: "results/benchmarks/blast_regions/blast_{region}_{sample}.tsv"
    conda: "blast-env"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule final_regions:
    input: "results/tables/{region}/blast_{region}_{sample}.tsv"
    output: touch("results/coverage_control/final/{region}_{sample}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")