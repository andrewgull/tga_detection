#!/usr/bin/env python3
"""Copy Nanopore reads from """
import os
import pandas as pd
import shutil
import argparse
from tqdm import tqdm


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for copying multiplexed (i.e. with barcodes) read files from Argos to the local computer using cp.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("barcodes", metavar="<barcodes and paths file>",
                        help="a CSV file with the following columns: Barcode, Sample, Path\n"
                             "Barcode - integer representing sample's barcode\n"
                             "Sample - sample name, text, e.g. 76595_D_20g\n"
                             "Path - full path to the Sample's reads on Argos, without trailing slash\n"
                             "Commented lines will be skipped")
    parser.add_argument("destination", metavar="<destination prefix>",
                        help="a common prefix directory for the reads - "
                             "directory within the project where all the reads will be stored (with trailing slash)")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    csv_file_path = args.barcodes
    destination_prefix = args.destination

    # Read the CSV file with pandas, skipping commented lines
    df = pd.read_csv(
        csv_file_path,
        comment="#",  # Skip lines starting with "#" in case you want to slip some files use this
        usecols=["Barcode", "Sample", "Path"]  # Select specific columns
    )

    for index, row in df.iterrows():
        source_dir = row["Path"]
        row_name = row["Sample"]
        # Use destination directory from CSV
        destination_dir = os.path.join(destination_prefix, row["Sample"])
        print(f"Copying reads for {row_name}")

        # Loop through all files in the source directory
        for filename in tqdm(os.listdir(source_dir)):
            if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
                source_file = os.path.join(source_dir, filename)
                destination_file = os.path.join(destination_dir, filename)

                # Check if target directory exists, create if needed
                os.makedirs(destination_dir, exist_ok=True)

                # Copy the file
                shutil.copy2(source_file, destination_file)

print("Done")
