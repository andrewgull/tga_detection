#!/usr/bin/env bash
set -euo pipefail

CONTAINERS=( "default" "rscripts" "biostrings" )
OUTPUT_DIR="resources/apptainer"

# Verify singularity/apptainer is available
if ! command -v singularity &> /dev/null; then
    echo "Error: 'singularity' command not found. Please install Apptainer/Singularity." >&2
    exit 1
fi

# Check if the output directory exists, if not, create it
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi


for C in "${CONTAINERS[@]}"; do
    echo "Building Apptainer image for container: $C"
    # Note: The 'singularity build' invocation below for "${OUTPUT_DIR}/${C}.sif" from
    # "resources/containers/${C}/Apptainer.def" requires elevated privileges (either 
    # sudo or --fakeroot when the system supports unprivileged builds). 
    # --fakeroot is not universally available due to kernel/subordinate UID/GID mapping requirements.
    singularity build --force "${OUTPUT_DIR}/${C}.sif" "resources/containers/${C}/Apptainer.def"
    echo "Successfully built ${OUTPUT_DIR}/${C}.sif"
done

echo "All containers built successfully."