#!/usr/bin/bash
set -euo pipefail

cleanup_tar() {
    for C in "${CONTAINERS[@]}"; do
        rm -f resources/containers/${C}.tar
    done
}

trap cleanup_tar EXIT

CONTAINERS=( "default" "rscripts" "biostrings" )
OUTPUT_DIR="resources/apptainer"

# Check if the output directory exists, if not, create it
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

for C in "${CONTAINERS[@]}"; do
    echo "Building SIF for container: $C"
    
    # Build Docker image
    docker build -t ${C} resources/containers/${C}
    
    # Save Docker image to tar archive
    docker save -o resources/containers/${C}.tar ${C}
    
    # Build Apptainer/Singularity image from Docker archive
    singularity build ${OUTPUT_DIR}/${C}.sif docker-archive://resources/containers/${C}.tar
    
    # Clean up tar file after successful build
    rm -f resources/containers/${C}.tar
    
    echo "Successfully built ${OUTPUT_DIR}/${C}.sif"
done

echo "All containers built successfully."