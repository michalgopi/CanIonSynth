# User Guide

This guide explains how to install, validate, run, and reproduce CanIonSynth: Synthetic CT Volumetric Data Generator for Aerosol Cans and Ionization Chambers.

## Overview

The repository has two primary generation entry points:

- `in_docker_organized/main_create_phantom_can.jl`
- `in_docker_organized/main_create_phantom_ionic_chamber.jl`

The supporting Python code is used for visualization, approximate Radon or inverse-Radon processing, and optional DICOM-SEG conversion. For everyday local use, the Julia entry points above are the recommended interface.

## Prerequisites

- Python 3.10 is recommended.
- Julia 1.11 is recommended for the checked-in `Manifest.toml`.
- `nii2dcm` is optional and only needed for DICOM slice export.
- `gcloud` is optional and only needed if you want to use the built-in Google Cloud upload paths.

## Local Installation

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

`tests/setup_env.jl` activates the repository project, runs `Pkg.instantiate()`, and rebuilds `PyCall` against the active virtual environment when `VIRTUAL_ENV` is present.

## Containerized Installation

This repository also ships with:

- `Dockerfile` for a containerized environment
- `.devcontainer/devcontainer.json` for VS Code Dev Containers

Build the container manually with:

```bash
docker build -t canionsynth .
```

Run it with:

```bash
docker run -it --rm canionsynth
```

If you use VS Code, open the repository and reopen it in the Dev Container.

## Validate The Environment

Before generating larger datasets, disable optional services and run the test suite.

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

Build the documentation to verify the docs environment too:

```bash
julia --project=docs docs/make.jl
```

## Main Commands

Both generators accept the same command-line structure:

```text
julia --project=. <script> <dims> <add_radon> <variable_spacing> <uuid> <randomize> <add_smooth> <additive_noise> [json_path]
```

### Argument Reference

| Position | Name | Description |
| --- | --- | --- |
| 1 | `dims` | Output dimensions in `NxMxD` form, for example `64x64x64`. |
| 2 | `add_radon` | `true` to run the projection and reconstruction helper workflow. |
| 3 | `variable_spacing` | `true` to allow variable spacing behavior where the workflow supports it. |
| 4 | `uuid` | Run identifier used for output folder naming and seed derivation. |
| 5 | `randomize` | Enables extra randomized parameter branches where supported. |
| 6 | `add_smooth` | Applies Gaussian smoothing to the final phantom volume. |
| 7 | `additive_noise` | Adds Gaussian noise, for example `0.0` or `0.1`. |
| 8 | `json_path` | Optional JSON config path for reproducible parameter sets. |

For local generation, keep `SKIP_WANDB=true` and `SKIP_UPLOAD=true` set unless you have configured those external services on purpose.

## Running A Can Phantom

Generate a can phantom with internal defaults:

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 false false run-can-001 false false 0.0
```

Generate a reproducible can phantom from the checked-in example config:

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 false false run-can-json false false 0.0 tests/configs/can_phantom.json
```

Generate a can phantom with the Radon or inverse-Radon stage enabled:

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 true false run-can-radon false false 0.0
```

## Running An Ionic Chamber Phantom

Generate an ionic chamber with internal defaults:

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-001 false false 0.0
```

Generate a reproducible square-top ionic chamber from JSON:

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-square false false 0.0 tests/configs/ionic_chamber_square.json
```

Generate a ball-like ionic chamber from JSON:

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-ball false false 0.0 tests/configs/ionic_chamber_ball.json
```

Generate a lollipop-like ionic chamber from JSON:

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false run-ionic-lollipop false false 0.0 tests/configs/ionic_chamber_lolipop.json
```

## JSON Configuration

Use the JSON fixtures in `tests/configs/` whenever you need a documented, reproducible setup.

- `tests/configs/can_phantom.json`: rounded can with multiple internal features
- `tests/configs/ionic_chamber.json`: fully specified square-top chamber
- `tests/configs/ionic_chamber_ball.json`: ball-like chamber shape selector
- `tests/configs/ionic_chamber_lolipop.json`: lollipop-like chamber shape selector
- `tests/configs/ionic_chamber_square.json`: square-top chamber shape selector

If you provide a JSON file, keep the command-line `dims`, `uuid`, and processing flags aligned with the run you want to reproduce.

## Output Artifacts

When `SKIP_UPLOAD=true`, the generators keep the output directory locally and print its location.

Typical outputs include:

- `example_can.nii.gz` or `ionic_chamber.nii.gz`: main phantom volume
- `fluid_mask.nii.gz` or component-specific masks such as `air_bool.nii.gz`
- `argss.json` or `ionic_chamber_params.json`: parameter log for the run
- `after_radon.nii.gz` and `after_radon_plus_before.nii.gz` when `add_radon=true`
- DICOM folders and DICOM-SEG folders when `nii2dcm` and the related helper path are available
- a zip archive that mirrors the local output directory

## Environment Variables

The following environment variables are useful during local development and CI:

- `SKIP_WANDB=true`: disables Weights and Biases logging
- `SKIP_UPLOAD=true`: disables Google Cloud upload and keeps outputs locally
- `FIXED_UUIDS=1`: makes the test harness use deterministic UUIDs
- `SAVE_OUTPUTS_TO=<path>`: copies test outputs to a chosen directory before cleanup

## Examples And Helper Scripts

- `scripts/visualize_nifti.py`: creates a PNG visualization for a generated NIfTI volume
- `scripts/generate_all_manual_tests.py`: generates a manual validation set and matching PNGs in `manual_test_outputs/`
- `examples/*.sh`: POSIX shell convenience wrappers for a few common generation cases
- `examples/batch_generation.py`: small orchestration example for repeated runs

The POSIX shell scripts in `examples/` are convenient on Linux or macOS. On Windows, running the Julia entry points directly is usually simpler.

## Troubleshooting

### PyCall Uses The Wrong Python

Activate your virtual environment first, then rerun:

```bash
julia --project=. tests/setup_env.jl
```

### DICOM Conversion Is Skipped

The repository will still generate NIfTI outputs if `nii2dcm` is missing. Install the DICOM toolchain only if you need DICOM slice output.

### WandB Or Cloud Uploads Fail

For normal local work, set `SKIP_WANDB=true` and `SKIP_UPLOAD=true`.

### The Coordinator Script Is Too Environment-Specific

`in_docker_organized/coordinate_phantom_create.py` still assumes specific Google Cloud Storage paths and is best treated as a project-specific automation helper. For local use, run the Julia entry points directly.

### Julia Version Mismatch

If you see manifest-related warnings on Julia 1.10, switch to Julia 1.11. The current lockfile was generated with Julia 1.11.4.
