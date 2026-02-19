# User Guide

This guide provides detailed instructions on how to use the Synthetic Tomographic Data Generation repository.

## Table of Contents
1.  [Prerequisites & Installation](#prerequisites--installation)
2.  [Repository Structure](#repository-structure)
3.  [Running the Simulations](#running-the-simulations)
    *   [Can Phantom](#can-phantom)
    *   [Ionic Chamber Phantom](#ionic-chamber-phantom)
4.  [Configuration](#configuration)
    *   [Command Line Arguments](#command-line-arguments)
    *   [JSON Configuration](#json-configuration)
5.  [Testing & Verification](#testing--verification)
6.  [Output Artifacts](#output-artifacts)
7.  [Troubleshooting](#troubleshooting)

---

## Prerequisites & Installation

### Core Requirements
*   **Julia (v1.10+)**: The core simulation logic is written in Julia.
*   **Python (v3.10+)**: Used for DICOM conversion, logging, and visualization.
*   **Google Cloud CLI (`gcloud`)**: Required for uploading artifacts (unless skipped).

### Julia Dependencies
This project relies on several Julia packages, including a custom version of `ImagePhantoms`.

To install dependencies, you can use the provided setup script:
```bash
julia tests/setup_env.jl
```
Or manually install them:
```julia
using Pkg
Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git")
Pkg.add(["Sinograms", "ImageGeoms", "Unitful", "MIRTjim", "FFTW", "LazyGrids", "Plots", "PyCall", "JSON", "ImageFiltering", "Accessors"])
```

### Python Dependencies
Install the required Python packages:
```bash
pip install SimpleITK numpy wandb scikit-image adrt matplotlib nibabel
```

### System Tools
*   **nii2dcm**: Part of `dcmtk` or similar suites, used to convert NIfTI files to DICOM. If missing, the scripts will skip DICOM conversion but still generate NIfTI files.

---

## Repository Structure

*   **`in_docker_organized/`**: Contains the primary source code for simulations.
    *   `main_create_phantom_can.jl`: Entry point for Can Phantom generation.
    *   `main_create_phantom_ionic_chamber.jl`: Entry point for Ionic Chamber generation.
    *   `geometry_utils.jl`: Helper functions for geometry and file handling.
*   **`tests/`**: Test suite and configurations.
    *   `run_tests.jl`: Main test runner.
    *   `configs/`: JSON configuration files for reproducible runs.
    *   `comparisons/`: Visual comparisons of outputs.
*   **`docs/`**: Documentation files.

---

## Running the Simulations

The simulations are run via command-line Julia scripts located in `in_docker_organized/`.

### Can Phantom
Generates a cylindrical container with fluid, meniscus, and optional objects.

**Basic Command:**
```bash
julia in_docker_organized/main_create_phantom_can.jl <dims> <radon> <var_spacing> <uuid> <random> <smooth> <noise> [json_path]
```

**Example:**
```bash
# Generate a 128x128x128 phantom, no radon transform, random parameters
julia in_docker_organized/main_create_phantom_can.jl 128x128x128 false false my-run-001 true false 0.0
```

### Ionic Chamber Phantom
Generates a complex multi-layered electrode structure.

**Basic Command:**
```bash
julia in_docker_organized/main_create_phantom_ionic_chamber.jl <dims> <radon> <var_spacing> <uuid> <random> <smooth> <noise> [json_path]
```

**Example:**
```bash
# Generate from a specific JSON configuration
julia in_docker_organized/main_create_phantom_ionic_chamber.jl 32x32x32 false false my-run-002 false false 0.0 tests/configs/ionic_chamber.json
```

---

## Configuration

### Command Line Arguments

| Position | Name | Type | Description |
| :--- | :--- | :--- | :--- |
| 1 | `dims` | String | Dimensions format `NxMxD` (e.g., `128x128x128`). |
| 2 | `add_radon` | Bool | `true` to compute Radon transform (projections). |
| 3 | `variable_spacing` | Bool | Adjusts voxel spacing dynamically (mostly for Ionic Chamber). |
| 4 | `uuid` | String | Unique run identifier used for output folder naming. |
| 5 | `randomize` | Bool | `true` to randomize parameters (ignored if JSON provided). |
| 6 | `add_smooth` | Bool | `true` to apply Gaussian smoothing to the final image. |
| 7 | `additive_noise` | Float | Magnitude of Gaussian noise to add (e.g., `0.1`). |
| 8 | `json_path` | String | (Optional) Path to a JSON config file. Overrides `randomize`. |

### JSON Configuration
For reproducible research and testing, define parameters in a JSON file.

**Sample Can Config (`tests/configs/can_phantom.json`):**
```json
{
    "center_cylinder": [0.0, 0.0, 0.0],
    "bigger_cyl_size": [3.0, 3.0, 8.0],
    "cylinder_wall_thickness": 0.1,
    "rounded_bottom": true,
    "density_inside": 0.5,
    "dual_phase_percentage": 0.5
}
```

**Sample Ionic Chamber Config (`tests/configs/ionic_chamber_lolipop.json`):**
```json
{
    "lolipop_like": true,
    "total_len": 3.5,
    "graphite_density": 1.0,
    "dims": [64, 64, 64]
}
```

---

## Testing & Verification

The repository includes a robust test suite to verify that code changes do not break the physics or geometry of the phantoms.

### Running Tests
To run the full test suite (using small 32x32x32 volumes):
```bash
julia tests/run_tests.jl
```

### Verification Visualization
To visually inspect changes against a known reference, use the visualization tool:
```bash
python3 tests/visualize_comparison.py
```
This script looks for output in `tests/vis_out_ref` and compares it against `tests/vis_out_orig` (if generated), saving difference images to `tests/comparisons/`.

### Environment Variables for Testing
*   `SKIP_WANDB=true`: Disables Weights & Biases logging (useful for local testing).
*   `SKIP_UPLOAD=true`: Disables `gcloud` upload.
*   `FIXED_UUIDS=1`: Forces deterministic UUIDs for reproducible filenames.
*   `SAVE_OUTPUTS_TO=<path>`: Copies generated artifacts to a specified directory instead of deleting them.

---

## Output Artifacts

The simulation generates a folder named with the `uuid` containing:

*   **`*.nii.gz`**: NIfTI format volumes.
    *   `example_can.nii.gz` / `ionic_chamber.nii.gz`: The final phantom.
    *   `fluid_mask.nii.gz`: Binary mask of the fluid volume.
    *   `*mask.nii.gz`: Masks for internal components (balls, pipes, etc.).
*   **`argss.json`**: A log of the full parameter set used for generation.
*   **`DICOM/`**: (If `nii2dcm` is installed) Directory containing DICOM slices and SEG objects.

---

## Troubleshooting

*   **`nii2dcm` not found:** The scripts check for this tool. If missing, a warning is printed, and DICOM conversion is skipped. The NIfTI files are still generated.
*   **WandB errors:** Ensure you are logged in (`wandb login`) or use `SKIP_WANDB=true`.
*   **Package loading errors:** Run `julia tests/setup_env.jl` to ensure the registry and `ImagePhantoms` are correctly configured.
