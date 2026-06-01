# System Architecture Overview

CanIonSynth is organized as a small Julia-centered application with Python helper scripts around the edges.

## Runtime Layers

### 1. Julia Entry Points

The user-facing generation layer lives in:

- `in_docker_organized/main_create_phantom_can.jl`
- `in_docker_organized/main_create_phantom_ionic_chamber.jl`

These scripts parse arguments, create or load parameter sets, generate voxelized phantoms, save outputs, and optionally trigger downstream processing such as Radon reconstruction and DICOM export. All console output uses a shared `tlog` helper that prefixes every message with a `[YYYY-MM-DD HH:MM:SS] [INFO]` timestamp.

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

`in_docker_organized/geometry_utils.jl` contains the image manipulation and export utilities used by both generators. It also defines the `tlog` helper function (guarded by `@isdefined`) that all Julia files use for timestamped output.

Key responsibilities include:

- writing NIfTI outputs with spatial metadata
- rotating or shifting masks when the geometry requires it
- attempting DICOM slice export through `nii2dcm`
- launching the Python DICOM-SEG conversion helper

### 4. Python Helper Layer

The Python scripts are focused on tasks that are easier to express with Python imaging libraries. All Python helpers use the standard `logging` module configured with the `[YYYY-MM-DD HH:MM:SS] [INFO]` format.

- `in_docker_organized/radon_iradon_3d.py`: two-stage parallel Radon and inverse-Radon reconstruction helper; accepts `--n-theta` (projection angle count, default 10) and `--noise-level` (sinogram noise scale 0–1, default 0); invoked by the Julia entry points when `add_radon=true`
- `in_docker_organized/nifti_to_dicom_seg.py`: DICOM-SEG conversion helper
- `in_docker_organized/coordinate_phantom_create.py`: batch orchestration helper; reads a local JSON config file (path set by `CONFIG_JSON_PATH` env var, defaulting to `control_json_can_128.json`) and launches the Julia generators
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
6. the script zips the result directory and prints the local path

## External Integrations

The only optional external tool is `nii2dcm` for DICOM slice export. DICOM export is skipped gracefully when `nii2dcm` is unavailable.

## Coordinator Script

`in_docker_organized/coordinate_phantom_create.py` reads its configuration from a local JSON file rather than from any cloud storage. The config file path defaults to `in_docker_organized/control_json_can_128.json` and can be overridden with the `CONFIG_JSON_PATH` environment variable. If `json_folders_path` in the config points to a local directory, the script lists JSON files from that directory and runs one generation job per file.
