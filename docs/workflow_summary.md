# Synthetic Tomographic Data Generation Workflow

## Overview

This repository generates synthetic industrial CT phantoms and their derived artifacts for validation, metrology, and reconstruction studies. Two Julia entry points create the main datasets, while Python helpers handle a small number of post-processing and interchange tasks.

## Canonical Entry Points

Use these scripts directly from the repository root:

- `in_docker_organized/main_create_phantom_can.jl`
- `in_docker_organized/main_create_phantom_ionic_chamber.jl`

The Python coordinator `in_docker_organized/coordinate_phantom_create.py` still exists, but it assumes specific Google Cloud Storage paths and is best treated as a project-specific automation helper rather than the primary local interface.

## Can Phantom Workflow

The can workflow produces cylindrical container phantoms with internal fluids and optional objects.

### Inputs

- command-line arguments for dimensions and processing flags
- optional JSON configuration for reproducible geometry and material settings
- UUID string used for output folder naming and seed derivation

### Main Steps

1. Parse command-line arguments and optional JSON overrides.
2. Build the geometric object set in Julia using the geometry helper files.
3. Convert the geometry to voxelized masks and density volumes.
4. Compute numerical and analytical fluid volumes for validation.
5. Save the main volume and masks as NIfTI files.
6. Optionally run the approximate Radon and inverse-Radon helper script.
7. Optionally attempt DICOM and DICOM-SEG export.
8. Zip the output directory and optionally upload it to Google Cloud Storage.

### Typical Outputs

- `example_can.nii.gz`
- `fluid_mask.nii.gz`
- component masks such as `ball1_mask.nii.gz` or `pipe_mask.nii.gz`
- `argss.json`
- `after_radon.nii.gz` and `after_radon_plus_before.nii.gz` when `add_radon=true`

## Ionic Chamber Workflow

The ionic chamber workflow produces multi-layer, multi-material chamber phantoms.

### Inputs

- command-line arguments for dimensions and processing flags
- optional JSON configuration for exact shape, dimensions, and materials
- UUID string used for output folder naming and seed derivation

### Main Steps

1. Parse command-line arguments and JSON overrides.
2. Build a layered cylinder specification for the requested chamber variant.
3. Voxelize the chamber, material layers, and air cavity.
4. Compute analytical and numerical air-volume estimates.
5. Save the phantom and component masks as NIfTI files.
6. Optionally run the approximate Radon and inverse-Radon helper script.
7. Optionally create DICOM and DICOM-SEG outputs when the DICOM toolchain is available.
8. Zip the output directory and optionally upload it to Google Cloud Storage.

### Typical Outputs

- `ionic_chamber.nii.gz`
- `ionic_chamber_params.json`
- masks such as `air_bool.nii.gz`, `copper_el_bool.nii.gz`, and `combined_bool_array_electrode.nii.gz`
- `after_radon.nii.gz` when `add_radon=true`

## Optional Post-Processing

### Radon And Reconstruction

When `add_radon=true`, the Julia scripts call `in_docker_organized/get_approximate_radon_inverse.py` through the same Python interpreter used by `PyCall`.

### DICOM Output

The repository attempts DICOM export only when the required tooling is available.

- If `nii2dcm` is not installed, the NIfTI outputs are still produced.
- If no reference DICOM series exists, DICOM-SEG conversion is skipped instead of breaking the run.

## Output Lifecycle

During a local run with `SKIP_UPLOAD=true`:

- a temporary output directory is created
- the directory path is printed as `Output stored in: ...`
- a zip archive is created next to the output
- the files remain available for inspection

When uploads are enabled, the scripts upload the archive to hardcoded Google Cloud Storage locations and then delete the local outputs.

## Reproducibility Controls

- `tests/configs/` contains reproducible JSON fixtures.
- `SKIP_WANDB=true` disables optional experiment tracking.
- `SKIP_UPLOAD=true` keeps local outputs.
- `FIXED_UUIDS=1` makes the test harness deterministic.
- `SAVE_OUTPUTS_TO=<path>` preserves test outputs before cleanup.

## Validation Workflow

The automated validation path is:

1. activate a Python virtual environment
2. run `julia --project=. tests/setup_env.jl`
3. run `julia --project=. tests/run_tests.jl`
4. run `julia --project=docs docs/make.jl` when validating the documentation too
