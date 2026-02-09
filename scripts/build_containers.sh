#!/usr/bin/bash
set -euo pipefail

# Mode: docker or apptainer
MODE=${1:-""}

if [[ "$MODE" != "docker" && "$MODE" != "apptainer" ]]; then
    echo "Usage: $0 {docker|apptainer}"
    exit 1
fi

CONTAINERS=( "default" "rscripts" "biostrings" )


if [[ "$MODE" == "apptainer" ]]; then
    OUTPUT_DIR="resources/apptainer"
    # Check if the output directory exists, if not, create it
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
        mkdir -p "$OUTPUT_DIR"
    fi
elif [[ "$MODE" == "docker" ]]; then
    OUTPUT_DIR="resources/docker"
    # Check if the output directory exists, if not, create it
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
        mkdir -p "$OUTPUT_DIR"
    fi
fi

for C in "${CONTAINERS[@]}"; do
    if [[ "$MODE" == "docker" ]]; then
        echo "Building Docker image for container: $C"
        docker build -t "tga_detection_${C}:latest" "resources/containers/${C}/Dockerfile"
        echo "Successfully built Docker image: tga_detection_${C}:latest"
    elif [[ "$MODE" == "apptainer" ]]; then
        echo "Building Apptainer image for container: $C"
        # Build Apptainer/Singularity image from definition file
        singularity build "${OUTPUT_DIR}/${C}.sif" "resources/containers/${C}/Apptainer.def"
        echo "Successfully built ${OUTPUT_DIR}/${C}.sif"
    fi
done

echo "All $MODE containers built successfully."