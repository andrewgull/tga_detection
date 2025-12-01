#!/usr/bin/bash

CONTAINERS=( "default" "rscripts" "biostrings" )

export APPTAINER_TMPDIR=/home/andrei/Data/singularity/tmp
export APPTAINER_CACHEDIR=/home/andrei/Data/singularity/cache

for C in "${CONTAINERS[@]}"; do
    echo "Building SIF for container: $C"
    
    # Build Docker image
    docker build -t ${C} workflow/docker/${C}
    
    # Save Docker image to tar archive
    docker save -o workflow/docker/${C}.tar ${C}
    
    # Build Apptainer/Singularity image from Docker archive
    singularity build resources/apptainer/${C}.sif docker-archive://workflow/docker/${C}.tar
    
    # Clean up tar file after successful build
    rm -f workflow/docker/${C}.tar
    
    echo "Successfully built resources/apptainer/${C}.sif"
done

echo "All containers built successfully."