#!/usr/bin/env python3
"""Aggregate per-sample frequency tables into a single TSV and XLSX."""
import argparse
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

    dfs = []
    for sample, file in zip(samples["sample"], input_files):
        df = pd.read_csv(file, sep="\t")
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

try:
    snakemake  # noqa: F821 — injected by Snakemake at runtime
    with open(snakemake.log[0], "w") as _log:  # noqa: F821
        sys.stdout = _log
        sys.stderr = _log
        main(
            sample_table=snakemake.config["sample_table"],  # noqa: F821
            input_files=list(snakemake.input),  # noqa: F821
            output_tsv=snakemake.output["tsv"],  # noqa: F821
            output_xlsx=snakemake.output["xlsx"],  # noqa: F821
        )
except NameError:
    args = parse_args()
    main(
        sample_table=args.sample_table,
        input_files=args.input_files,
        output_tsv=args.output_tsv,
        output_xlsx=args.output_xlsx,
    )
