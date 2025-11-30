ARG VERSION=4.4.3
ARG ENV_NAME=rscripts

# The build-stage image:
FROM mambaorg/micromamba:1.5.8 AS build

ARG ENV_NAME=rscripts

# Install tools using micromamba
# R requires the environment to be at the same path as it was built.
# We use a consistent path /opt/conda/envs/${ENV_NAME} for both build and runtime.
# We ensure the parent directory is writable so micromamba can create the environment.
USER root
RUN mkdir -p /opt/conda/envs && \
    chown mambauser:mambauser /opt/conda/envs
USER mambauser

RUN micromamba create -y -p /opt/conda/envs/${ENV_NAME} -c bioconda -c conda-forge \
    r-base=4.4.3 \
    r-dplyr=1.1.4 \
    r-readr=2.1.5 \
    r-purrr=1.0.2 \
    r-tidyr=1.3.1 \
    && micromamba clean --all --yes

# The runtime-stage image
FROM debian:bookworm-slim AS runtime

ARG VERSION
ARG ENV_NAME
ENV VERSION=${VERSION}
ENV PATH="/opt/conda/envs/${ENV_NAME}/bin:${PATH}"

# Create a non-root user
RUN groupadd -r appuser && useradd -r -g appuser appuser

# Copy the environment from the build stage
# We copy from the SAME path to the SAME path to ensure R works
COPY --from=build --chown=appuser:appuser /opt/conda/envs/${ENV_NAME} /opt/conda/envs/${ENV_NAME}

################## METADATA ######################
LABEL maintainer="andrei.guliaev@scilifelab.uu.se"
LABEL version=${VERSION}
LABEL env=${ENV_NAME}

# Switch to non-root user
USER appuser
