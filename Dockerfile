# Use a slim Python base image and install Julia manually.
FROM python:3.10-slim-bullseye

ARG JULIA_RELEASE=1.10
ARG JULIA_VERSION=1.10.4

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV SKIP_UPLOAD=true
ENV SKIP_WANDB=true
ENV JULIA_HISTORY=/workspace/.julia_history

# Install basic dependencies
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    curl ca-certificates nano wget \
    build-essential software-properties-common \
    git cmake ninja-build zip unzip \
    libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Julia
RUN curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_RELEASE}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz | \
    tar -C /usr/local -x -z --strip-components=1 -f -

# Install additional runtime dependencies
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    vim net-tools ffmpeg pkg-config && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the full repository so tests and docs are available inside the image.
COPY . /workspace

# Setup thread management for Julia
COPY setup_threads.sh /etc/profile.d/

# Install Python dependencies and pin numpy for ADRT compatibility.
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install -r /workspace/requirements.txt && \
    python3 -m pip install numpy==1.23.2

# Instantiate the Julia environment using the checked-in lockfile.
RUN julia --project=/workspace -e 'using Pkg; Pkg.instantiate()'

# Make the container start in an interactive, reviewer-friendly shell.
RUN chmod +x /workspace/entrypoint.sh

WORKDIR /workspace
ENTRYPOINT ["/bin/bash", "/workspace/entrypoint.sh"]

# Use the shell started by entrypoint.sh unless a command is provided.
