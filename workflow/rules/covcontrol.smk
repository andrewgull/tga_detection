from snakemake.io import expand, directory, touch
import pandas as pd

rule all:
    input:
        expand("results/final/cc/{sample}_{region}_all.done", sample=config['samples'], region=config['regions']), "results/cc/coverage_control_table.tsv"

rule create_cr:
    input: "resources/plasmid/plasmid.fasta"
    output: "results/control_regions/{region}.fasta"
    threads: 18
    log: "results/logs/seqkit_subseq_{region}.log"
    conda: "seqkit-env"
    container: "containers/seqkit.sif"
    params: 
        start=lambda wildcards: config['regions'][wildcards.region]["start"],
        end=lambda wildcards: config['regions'][wildcards.region]["end"]
    shell: "seqkit subseq -j {threads} -r {params.start}:{params.end} {input} 1> {output} 2> {log}"

rule make_blast_db:
    input: "results/reads/{sample}/reads_filtered.fasta.gz"
    output: directory("results/blast_databases/{sample}")
    threads: 18
    log: "results/logs/{sample}_blastdb.log"
    conda: "blast-env"
    shell: "pigz -c -d -p {threads} {input} | makeblastdb -in - -dbtype nucl -title blastdb -out {output}/blastdb 2> {log}"

rule blast_region:
    input: 
        query="results/control_regions/{region}.fasta", 
        database="results/blast_databases/{sample}"
    output: "results/cc/{sample}/blast_{region}.tsv"
    threads: 18
    params: fmt=config['format'], n_alns=config['n_cc_alns']
    log: "results/logs/{sample}_{region}_blast.log"
    conda: "blast-env"
    shell: "blastn -query {input.query} -db {input.database}/blastdb -outfmt {params.fmt} "
           "-num_threads {threads} -num_alignments {params.n_alns} 1> {output} 2> {log}"

rule aggregate_сс_tables:
    input: expand("results/cc/{sample}/blast_{region}.tsv", sample=config['samples'], region=config['regions'])
    output: tsv = "results/cc/coverage_control_table.tsv"
            # xlsx = "results/cc/coverage_control_table.xlsx"
    run:
        dfs = []
        samples_count = len(input) / len(config['samples']) # there should be the same number of samples as input files
        samples = [key for key in config['samples'] for _ in range(int(samples_count))]
        regions = [key for key in config['regions']] * len(config['samples'])
        for sample, file, region in zip(samples, input, regions):
            df = pd.read_csv(file, sep='\t', header=None) # ensure no headers to avoid wrong concatenating
            df['sample'] = sample
            df[0] = region
            dfs.append(df)
        merged_df = pd.concat(dfs, axis=0, ignore_index=True)
        merged_df.to_csv(output.tsv, index=False, sep='\t')
        # requires openpyxl
        # merged_df.to_excel(output.xlsx, index=False, sheet_name='coverage')

rule final:
    input: cc="results/cc/{sample}/blast_{region}.tsv"
    output: touch("results/final/cc/{sample}_{region}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
