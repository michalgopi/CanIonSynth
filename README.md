# Synthetic Tomography Data Generator

Synthetic Tomography Data Generator creates synthetic industrial CT phantoms, segmentation masks, and optional reconstruction artifacts for metrology and reconstruction studies. The main workflows are implemented in Julia, with Python helper scripts for visualization, DICOM-SEG conversion, and approximate Radon or inverse-Radon processing.

## What The Repository Does

- Generates can phantoms with configurable wall thickness, rounded or flat bottoms, one or two fluids, meniscus effects, and optional internal objects.
- Generates ionic chamber phantoms with layered multi-material geometry and multiple head shapes.
- Saves the main outputs as NIfTI volumes and component masks.
- Optionally creates DICOM and DICOM-SEG outputs when the required external tooling is installed.
- Optionally runs a projection and reconstruction pass for reconstruction-oriented experiments.
- Supports reproducible runs from JSON configuration files and a validated automated test suite.

## Recommended Local Setup

The most predictable local setup today is Python 3.10 plus Julia 1.11 with an activated virtual environment before running `tests/setup_env.jl`. The checked-in `Manifest.toml` was generated with Julia 1.11.4.

### Windows PowerShell

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
julia --project=. tests/setup_env.jl
```

### Linux Or macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
julia --project=. tests/setup_env.jl
```

`tests/setup_env.jl` instantiates the Julia environment and rebuilds `PyCall` against the active virtual environment when `VIRTUAL_ENV` is set.

## Validate The Installation

Disable optional cloud logging and uploads for local work, then run the test suite.

### PowerShell

```powershell
$env:SKIP_WANDB = "true"
$env:SKIP_UPLOAD = "true"
$env:FIXED_UUIDS = "1"
julia --project=. tests/run_tests.jl
```

### Bash

```bash
export SKIP_WANDB=true
export SKIP_UPLOAD=true
export FIXED_UUIDS=1
julia --project=. tests/run_tests.jl
```

Build the documentation site with:

```bash
julia --project=docs docs/make.jl
```

## Generate Data

Run the main generators from the repository root.

For local runs, keep `SKIP_WANDB=true` and `SKIP_UPLOAD=true` set unless you have explicitly configured Weights and Biases and Google Cloud Storage uploads in your environment.

### Can Phantom

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 false false run-can-001 false false 0.0
```

### Can Phantom From JSON

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 false false run-can-json false false 0.0 tests/configs/can_phantom.json
```

### Ionic Chamber

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-001 false false 0.0
```

### Ionic Chamber From JSON

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-json false false 0.0 tests/configs/ionic_chamber_square.json
```

When `SKIP_UPLOAD=true`, the scripts keep the generated files locally and print the output directory and zip archive path.

## Reproducibility

- Julia dependencies are pinned in `Project.toml` and `Manifest.toml`.
- Python dependencies are listed in `requirements.txt`.
- JSON fixtures in `tests/configs/` provide reproducible reference parameter sets.
- `FIXED_UUIDS=1` makes test output naming deterministic.
- `SAVE_OUTPUTS_TO=<path>` preserves generated test artifacts before cleanup.
- `SKIP_WANDB=true` and `SKIP_UPLOAD=true` disable external services during local runs and CI.

## Optional External Tools

- `nii2dcm`: required only if you want DICOM slice export from generated NIfTI files.
- `gcloud`: required only if you want to use the built-in Google Cloud Storage upload paths.
- `wandb`: installed by default through `requirements.txt`, but can be disabled entirely with `SKIP_WANDB=true`.

## Container And Dev Container

This repository ships with:

- `Dockerfile` for a self-contained environment.
- `.devcontainer/devcontainer.json` for VS Code Dev Containers.

Build the image with:

```bash
docker build -t synthetic-tomo .
```

Open a shell in the container with:

```bash
docker run -it --rm synthetic-tomo
```

## Repository Layout

- `in_docker_organized/`: Julia generators and Python helper scripts.
- `tests/`: automated tests, JSON fixtures, and comparison helpers.
- `docs/`: user documentation, workflow notes, verification guides, and architecture notes.
- `scripts/`: utility scripts for visualization and manual test generation.
- `examples/`: convenience scripts and orchestration examples.
- `packages/`: vendored dependencies, including the patched `ImagePhantoms.jl` copy used by this project.

## Documentation

- [User Guide](docs/user_guide.md)
- [Workflow Summary](docs/workflow_summary.md)
- [System Architecture](docs/system_architecture.md)
- [Phantom Types](docs/phantom_types.md)
- [Functions Reference](docs/functions_reference.md)
- [Volume Calculation](docs/volume_calculation.md)
- [Manual Testing and Visualization](docs/manual_testing_visualization.md)
- [Manual Verification Guide](docs/manual_verification.md)
- [SoftwareX Submission Notes](SOFTWAREX_SUBMISSION.md)

## Citation And License

- Citation metadata is provided in `CITATION.cff`.
- Archive metadata is provided in `.zenodo.json`.
- The repository is licensed under the MIT License. See `LICENSE`.

## CI

GitHub Actions workflows are defined in `.github/workflows/ci.yml` and `.github/workflows/docs.yml`.
