import pandas as pd


def get_sample_path(wildcards: snakemake.io.wildcards, data_frame: pd.DataFrame) -> str:
    """
    Get the path to the fastq file using its sample name.
    arguments:
        wildcards: wildcards object.
        data_frame: A DataFrame object containing sample names and paths.
    returns:
        str: The path to the fastq file.
    """
    row = data_frame[data_frame["sample"] == wildcards.sample]
    if not row.empty:
        return row["path"].values[0]
    else:
        raise ValueError(f"Sample '{wildcards.sample}' not found in the provided DataFrame.")