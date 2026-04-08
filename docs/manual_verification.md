# Manual Verification Guide

This guide provides a reviewer-friendly checklist for validating the generated outputs beyond the automated test suite.

## Before You Start

Use a local environment with the optional external services disabled.

### PowerShell

```powershell
$env:SKIP_WANDB = "true"
$env:SKIP_UPLOAD = "true"
```

### Bash

```bash
export SKIP_WANDB=true
export SKIP_UPLOAD=true
```

Generate outputs with the commands from `docs/manual_testing_visualization.md`, then inspect them with `scripts/visualize_nifti.py` or a medical imaging viewer such as 3D Slicer or ITK-SNAP.

## Can Phantom Checklist

### Geometry

- verify the can walls form a clean cylindrical shell
- verify the lower geometry matches the intended flat or rounded-bottom configuration
- verify internal objects such as balls, pipe, and dispenser appear in plausible locations when enabled

### Fluids

- verify the fluid mask aligns with the visible fluid region
- verify dual-phase cases show two intensity regions when configured
- verify meniscus curvature appears near the wall region when the configuration enables it
- verify nonzero cut angles visibly tilt the fluid surface

### Files

- confirm `example_can.nii.gz` exists
- confirm `fluid_mask.nii.gz` exists
- confirm `argss.json` exists and contains the run parameters

## Ionic Chamber Checklist

### Shape Variants

- square-top: verify a flat head profile
- ball-like: verify a rounded or spherical head profile
- lollipop-like: verify a narrow neck feeding a wider head
- rounded-top: verify a smooth cap rather than a flat top

### Internal Layering

- verify the central electrode is centered
- verify the air cavity forms a distinct low-density region
- verify concentric material layers remain aligned with the chamber axis
- verify component masks such as `air_bool.nii.gz` and `combined_bool_array_electrode.nii.gz` align with the main chamber volume

### Files

- confirm `ionic_chamber.nii.gz` exists
- confirm `ionic_chamber_params.json` exists
- confirm expected component masks exist for the selected variant

## Radon And Reconstruction Checklist

For runs with `add_radon=true`:

- confirm `after_radon.nii.gz` exists
- for can phantoms, confirm `after_radon_plus_before.nii.gz` exists
- verify the reconstructed output is recognizably related to the source phantom while showing the expected reconstruction softening or artifacts

## Volume Verification Checklist

### Can Workflow

- inspect the console output for the analytical and numerical fluid-volume comparison
- inspect `argss.json` for the saved parameter and spacing information
- treat small relative differences as expected discretization error
- investigate large mismatches by checking spacing, masks, and complex geometry interactions

### Ionic Chamber Workflow

- inspect the console output for the analytical and numerical air-volume comparison
- confirm the reported difference stays in a plausible range for the chosen resolution

## Optional DICOM Verification

If `nii2dcm` is installed and DICOM export is enabled:

- confirm the output includes DICOM slice directories such as `example_can/` or `ionic_chamber/`
- confirm DICOM-SEG output directories exist for masks that were eligible for conversion
- confirm the DICOM series spatially matches the NIfTI source when opened in a viewer

If `nii2dcm` is not installed, a warning and NIfTI-only output is expected behavior rather than a failure.

## Reproducibility Checklist

- rerun the same JSON-driven command with a different UUID and confirm the geometry remains the same
- rerun the automated tests with `FIXED_UUIDS=1` and confirm the same set of expected files is produced
- preserve outputs with `SAVE_OUTPUTS_TO=<path>` when you need to compare runs side by side
