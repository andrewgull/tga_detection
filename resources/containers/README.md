# Docker Environments

This directory contains the Dockerfiles used to build the container environments for the Snakemake pipeline. These containers ensure reproducibility by locking down specific versions of software.

## Available Containers

| Image Name | Purpose | Key Software Included |
| :--- | :--- | :--- |
| **default** | General purpose utilities | `pigz`, `seqkit`, `filtlong`, `bedtools`, `blast`, `pandas`, `openpyxl` |
| **rscripts** | Statistical analysis steps | `R v4.4.3`, `tidyr v1.3.1`, `dplyr v1.1.4`, `readr v2.1.5`, `purrr v1.0.2` |
| **biostrings** | Sequence manipulation | `R v4.3`, `biostrings v2.70.1`, `dplyr v1.1.4`, `readr v2.1.5`, `purrr v1.0.2`|
§
## Prerequisites

To rebuild these images, you need:
1.  **Docker** (requires sudo/root privileges).
2.  **Apptainer** (formerly Singularity) to convert the images to `.sif` format.

## How to Build the Containers

We provide a helper script to automate the build process.

**1. Run from the Project Root:**
```bash
bash scripts/build_containers.sh
```

This script will use Dockerfiles located under this directory to build the containers and place them in `resources/apptainer/` directory where `snakemake` will look for them.