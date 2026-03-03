#!/usr/bin/env python3
"""Aggregate per-sample frequency tables into a single TSV and XLSX."""
import argparse
import pathlib
import sys

import openpyxl
import pandas as pd


# --- Functions ---


def main(
    sample_table: str,
    input_files: list[str],
    output_tsv: str,
    output_xlsx: str,
) -> None:
    samples = pd.read_csv(sample_table, sep="\t", dtype={"sample": str, "path": str})

    # Match each input file to its sample by the parent-directory component of the
    # path (pipeline convention: results/tables/{sample}/frequencies.tsv).
    # This avoids silent misalignment from relying on CLI argument order.
    sample_to_file: dict[str, str] = {}
    for file in input_files:
        name = pathlib.Path(file).parent.name
        if name not in samples["sample"].values:
            raise ValueError(
                f"Input file {file!r}: parent directory {name!r} does not match "
                "any sample in the sample table."
            )
        if name in sample_to_file:
            raise ValueError(f"Duplicate input file for sample {name!r}.")
        sample_to_file[name] = file

    missing = [s for s in samples["sample"] if s not in sample_to_file]
    if missing:
        raise ValueError(f"No input file found for samples: {missing}")

    dfs = []
    for sample in samples["sample"]:
        df = pd.read_csv(sample_to_file[sample], sep="\t")
        df["sample"] = sample
        dfs.append(df)

    merged_df = pd.concat(dfs, axis=0)
    merged_df.to_csv(output_tsv, index=False, sep="\t")
    merged_df.to_excel(output_xlsx, index=False, sheet_name="Frequencies")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate per-sample frequency tables into a single TSV and XLSX."
    )
    parser.add_argument(
        "--sample-table",
        required=True,
        metavar="TSV",
        help="Sample table TSV with 'sample' and 'path' columns.",
    )
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        dest="input_files",
        metavar="TSV",
        help="Per-sample frequency TSVs in the same order as rows in --sample-table.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        metavar="TSV",
        help="Path for the merged output TSV.",
    )
    parser.add_argument(
        "--output-xlsx",
        required=True,
        metavar="XLSX",
        help="Path for the merged output XLSX.",
    )
    return parser.parse_args(argv)


# --- Execution ---

if "snakemake" in globals():
    _stdout, _stderr = sys.stdout, sys.stderr
    try:
        with open(snakemake.log[0], "w") as _log:  # noqa: F821
            sys.stdout = _log
            sys.stderr = _log
            main(
                sample_table=snakemake.config["sample_table"],  # noqa: F821
                input_files=list(snakemake.input),  # noqa: F821
                output_tsv=snakemake.output["tsv"],  # noqa: F821
                output_xlsx=snakemake.output["xlsx"],  # noqa: F821
            )
    finally:
        sys.stdout = _stdout
        sys.stderr = _stderr
elif __name__ == "__main__":
    args = parse_args()
    main(
        sample_table=args.sample_table,
        input_files=args.input_files,
        output_tsv=args.output_tsv,
        output_xlsx=args.output_xlsx,
    )
