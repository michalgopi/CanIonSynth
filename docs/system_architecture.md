# System Architecture Overview

The Synthetic Tomography Data Generator is organized as a small Julia-centered application with Python helper scripts around the edges.

## Runtime Layers

### 1. Julia Entry Points

The user-facing generation layer lives in:

- `in_docker_organized/main_create_phantom_can.jl`
- `in_docker_organized/main_create_phantom_ionic_chamber.jl`

These scripts parse arguments, create or load parameter sets, generate voxelized phantoms, save outputs, and optionally trigger downstream processing such as Radon reconstruction and DICOM export.

### 2. Geometry Construction Layer

The lower-level geometry logic is split across:

- `in_docker_organized/get_geometry_main.jl`
- `in_docker_organized/get_rounded_bottom_b.jl`
- `in_docker_organized/volume_integration.jl`

This layer is responsible for:

- defining the can and ionic chamber primitives
- composing complex shapes from cylinders, ellipsoids, half-spheres, and related masks
- computing analytical volume estimates for comparison against voxel counts

### 3. Image IO And Interop Layer

`in_docker_organized/geometry_utils.jl` contains the image manipulation and export utilities used by both generators.

Key responsibilities include:

- writing NIfTI outputs with spatial metadata
- rotating or shifting masks when the geometry requires it
- attempting DICOM slice export through `nii2dcm`
- launching the Python DICOM-SEG conversion helper

### 4. Python Helper Layer

The Python scripts are focused on tasks that are easier to express with Python imaging libraries.

- `in_docker_organized/get_approximate_radon_inverse.py`: approximate projection and reconstruction helper
- `in_docker_organized/nifti_to_dicom_seg.py`: DICOM-SEG conversion helper
- `scripts/visualize_nifti.py`: visualization helper for manual inspection
- `scripts/generate_all_manual_tests.py`: manual test batch generator

### 5. Validation And Documentation Layer

- `tests/setup_env.jl`: environment bootstrap and `PyCall` binding
- `tests/run_tests.jl`: end-to-end test harness for both main generators
- `docs/`: user documentation and reviewer-facing project notes

### 6. Environment Layer

Reproducibility depends on:

- `Project.toml` and `Manifest.toml` for Julia
- `requirements.txt` for Python
- `Dockerfile` and `.devcontainer/devcontainer.json` for containerized development
- GitHub Actions workflows for tests and documentation builds

## Data Flow

The normal local execution path is:

1. an operator runs one of the Julia entry points
2. the script resolves configuration from arguments and optional JSON
3. Julia geometry helpers build the phantom and masks
4. the generator saves NIfTI outputs and parameter logs
5. optional Python helpers generate Radon-derived outputs or DICOM-SEG artifacts
6. the script zips the result directory and either keeps it locally or uploads it

## External Integrations

Several integrations are optional rather than required for local use.

- Weights and Biases logging can be disabled with `SKIP_WANDB=true`.
- Google Cloud Storage upload can be disabled with `SKIP_UPLOAD=true`.
- DICOM slice export is skipped when `nii2dcm` is unavailable.

## Coordinator Script Status

`in_docker_organized/coordinate_phantom_create.py` is still present, but it is tightly coupled to specific Google Cloud Storage locations and temporary filesystem paths. It is better viewed as a project-specific automation helper than as the primary documented interface for local users.
