#!/usr/bin/bash
set -euo pipefail

CONTAINERS=( "default" "rscripts" "biostrings" )
OUTPUT_DIR="resources/apptainer"
# Check if the output directory exists, if not, create it
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi


for C in "${CONTAINERS[@]}"; do
    echo "Building Apptainer image for container: $C"
    # Build Apptainer/Singularity image from definition file
    singularity build "${OUTPUT_DIR}/${C}.sif" "resources/containers/${C}/Apptainer.def"
    echo "Successfully built ${OUTPUT_DIR}/${C}.sif"
done

echo "All containers built successfully."