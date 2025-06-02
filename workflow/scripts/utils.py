import pandas as pd
from typing import Any


def get_sample_path(wildcards: Any, data_frame: pd.DataFrame) -> str:
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


def tsv2dict(file_path: str) -> dict:
    """
    Read a TSV file and convert it to a dictionary.
    arguments:
        file_path: Path to the TSV file.
    returns:
        dict: A dictionary with parameters as keys and their values.
    """
    df = pd.read_csv(
        file_path,
        sep="\t",
        dtype={"param": str, "value": str}
    )
    data_dict = df.set_index("param")["value"].to_dict()
    return data_dict