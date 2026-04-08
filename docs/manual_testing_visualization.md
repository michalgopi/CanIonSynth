# Manual Testing and Visualization Guide

This guide explains how to generate reviewable outputs and turn them into quick-look PNG visualizations.

## Prerequisites

Set up the project as described in `README.md`, then disable optional external services.

### PowerShell

```powershell
$env:SKIP_UPLOAD = "true"
$env:SKIP_WANDB = "true"
```

### Bash

```bash
export SKIP_UPLOAD=true
export SKIP_WANDB=true
```

## Fastest Manual Test Path

The easiest way to generate a review bundle is:

```bash
python scripts/generate_all_manual_tests.py
```

This creates `manual_test_outputs/` in the repository root and writes both `.nii.gz` outputs and matching `*_vis.png` quick-look images.

## Visualizing A Single NIfTI File

Use the helper below for any generated NIfTI file.

```bash
python scripts/visualize_nifti.py <path_to_nifti_file>
```

The script writes a `*_vis.png` file next to the input volume and includes sagittal, coronal, and axial slices plus an intensity histogram.

## Targeted Manual Runs

Run these commands from the repository root.

### Rounded Can From The Checked-In Config

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 false false manual-can false false 0.0 tests/configs/can_phantom.json
```

Expected files:

- `example_can.nii.gz`
- `fluid_mask.nii.gz`
- `argss.json`

### Can With Radon Or Reconstruction

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 64x64x64 true false manual-can-radon false false 0.0
```

Expected additional files:

- `after_radon.nii.gz`
- `after_radon_plus_before.nii.gz`

### Square-Top Ionic Chamber

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false manual-square false false 0.0 tests/configs/ionic_chamber_square.json
```

### Ball-Like Ionic Chamber

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false manual-ball false false 0.0 tests/configs/ionic_chamber_ball.json
```

### Lollipop-Like Ionic Chamber

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 64x64x64 false false manual-lollipop false false 0.0 tests/configs/ionic_chamber_lolipop.json
```

## Finding The Output Directory

When `SKIP_UPLOAD=true`, each generator prints:

```text
Output stored in: <path>
Zip stored in: <path>
```

Use that output directory with the visualization script.

Example:

```bash
python scripts/visualize_nifti.py <output_dir>/example_can.nii.gz
python scripts/visualize_nifti.py <output_dir>/ionic_chamber.nii.gz
```

## POSIX Example Scripts

The shell scripts in `examples/` are convenience wrappers for POSIX shells.

- `examples/generate_steel_can.sh`
- `examples/generate_aluminum_can.sh`
- `examples/generate_ionic_chamber.sh`

They are useful on Linux or macOS, but on Windows the direct Julia commands above are usually easier.

## What To Look For In The PNGs

- clear outer boundaries without missing slices
- expected internal objects or layered regions
- fluid masks aligned with visible fluid regions
- expected head shape for the ionic chamber variant
- reconstructed outputs that remain structurally related to the original phantom
