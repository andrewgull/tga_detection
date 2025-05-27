# script to concatenate tables with frequencies
import pandas as pd
import openpyxl

table_names = snakemake.input
output_tsv = snakemake.output['tsv']
output_xlsx = snakemake.output['xlsx']

# parse the sample names table config
samples = pd.read_csv(snakemake.config["sample_table"], sep="\t", dtype={"sample": str, "path": str})

dfs = []
for sample, file in zip(samples["sample"], table_names):
    df = pd.read_csv(file, sep='\t')
    df['sample'] = sample
    dfs.append(df)
merged_df = pd.concat(dfs, axis=0)
merged_df.to_csv(output_tsv, index=False, sep='\t')
merged_df.to_excel(output_xlsx, index=False, sheet_name='Frequencies')
