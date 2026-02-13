# Synthetic Tomographic Data Generation Workflow

## Overview

This repository is designed for the generation of synthetic tomographic data for metrology purposes in industrial computed tomography (CT). It simulates complex 3D phantoms, performs CT projection (Radon transform), and saves the data in various formats (HDF5, NIfTI, DICOM) for further analysis and reconstruction.

The core logic is implemented in Julia, leveraging packages like `ImagePhantoms`, `Sinograms`, and `ImageGeoms`. Python scripts are used for specific tasks like DICOM conversion and inverse Radon transform approximations.

## Main Workflows

There are two primary workflows for generating phantoms:

### 1. Can Phantom Generation
Simulates a cylindrical container ("can") filled with fluid(s) and containing various objects.

**Script:** `in_docker_organized/main_create_phantom_can.jl`

**Key Features:**
*   **Container:** A cylindrical can with configurable dimensions, wall thickness, and bottom curvature (flat or rounded).
*   **Contents:**
    *   **Fluids:** Single or dual-phase fluids with configurable densities and meniscus.
    *   **Objects:** Pipes, dispensers, and spherical "balls" floating in the fluid.
    *   **Defects:** Simulated by placing objects or modifying densities.
*   **Geometry:** Defined using constructive solid geometry (CSG) principles with cylinders, ellipsoids, and half-spheres.

**Workflow Steps:**
1.  **Configuration:** parameters are read from command-line arguments or a JSON file.
2.  **Phantom Creation:** `random_can` or `json_based_can` functions generate the 3D phantom based on the configuration.
3.  **Volume Calculation:** Analytical and numerical volumes of the fluid are calculated for metrology verification.
4.  **Projection:** (Optional) CT projections are generated using `get_CTFAN_proj`.
5.  **Output:**
    *   Phantom and masks are saved as NIfTI files.
    *   DICOM SEG (Segmentation) objects are created.
    *   Data is zipped and uploaded to Google Cloud Storage.

### 2. Ionic Chamber Phantom Generation
Simulates an ionic chamber, a device used for measuring ionizing radiation, which has a complex internal structure.

**Script:** `in_docker_organized/main_create_phantom_ionic_chamber.jl`

**Key Features:**
*   **Structure:** A multi-layered cylindrical device with electrodes, insulation, and housing.
*   **Components:**
    *   **Central Electrode:** Copper or graphite.
    *   **Insulation:** Polyethylene layers.
    *   **Shielding:** Aluminum and graphite layers.
    *   **Air Cavity:** Defined by the gap between electrodes and housing.
*   **Variations:** Supports different chamber tip shapes:
    *   Square top
    *   Ball-like (spherical)
    *   Lollipop-like (flat head)
    *   Rounded top (standard)

**Workflow Steps:**
1.  **Configuration:** Parameters are generated randomly or read from a JSON file.
2.  **Phantom Creation:** `create_ionic_chamber_phantom` constructs the 3D volume layer by layer.
3.  **Volume Calculation:** Air cavity volume is calculated analytically and numerically.
4.  **Projection:** (Optional) CT projections are generated.
5.  **Output:**
    *   Phantom and component masks are saved as NIfTI and DICOM.
    *   Parameters are saved to JSON.
    *   Data is zipped and uploaded to Google Cloud Storage.

## Configuration

### Command Line Arguments

Both main scripts accept command-line arguments to control the simulation.

**Example for Can Phantom:**
```bash
julia in_docker_organized/main_create_phantom_can.jl <dims> <add_radon> <variable_spacing> <uuid> <randomize> <add_smooth> <additive_noise> [json_path]
```

*   `dims`: Dimensions of the output volume (e.g., "256x256x256").
*   `add_radon`: Boolean, whether to perform Radon transform.
*   `variable_spacing`: Boolean, (not widely used in can phantom, but present in ionic).
*   `uuid`: Unique identifier for the run.
*   `randomize`: Boolean, whether to randomize parameters (if JSON not provided).
*   `add_smooth`: Boolean, whether to apply Gaussian smoothing.
*   `additive_noise`: Float, amount of noise to add.
*   `json_path`: (Optional) Path to a JSON file with specific parameters.

**Example for Ionic Chamber:**
```bash
julia in_docker_organized/main_create_phantom_ionic_chamber.jl <dims> <add_radon> <variable_spacing> <uuid> <randomize> <add_smooth> <additive_noise> [json_path]
```

### JSON Configuration

For reproducible simulations, parameters can be supplied via a JSON file.

**Key JSON Parameters (Can Phantom):**
*   `center_cylinder`: [x, y, z] coordinates.
*   `bigger_cyl_size`: [radius, radius, height].
*   `cylinder_wall_thickness`: Thickness of the can wall.
*   `density_inside`: Density of the fluid.
*   `pipe_len`, `pipe_cross_section`: Dimensions of the internal pipe.
*   `rounded_bottom`: Boolean, enables rounded bottom geometry.
*   `dual_phase_percentage`: Float (0-1), ratio of the two fluid phases.

**Key JSON Parameters (Ionic Chamber):**
*   `graphite_density`, `copper_density`, etc.: Material densities.
*   `total_len`, `base_len`, `main_radius`: Overall dimensions.
*   `square_top`, `ball_like`, `lolipop_like`: Boolean flags for shape variants.
*   `air_thickness`: Thickness of the air cavity (critical for metrology).

## Output Structure

The output is typically organized in a temporary folder structure labeled with the UUID, then zipped and uploaded.

*   `*.nii.gz`: 3D volumetric data (NIfTI format).
    *   `example_can.nii.gz` / `ionic_chamber.nii.gz`: The main phantom image.
    *   `*_mask.nii.gz`: Binary masks for specific components (fluid, pipe, air, etc.).
*   `DICOM/`: Directory containing DICOM slices and DICOM SEG objects.
*   `argss.json` / `ionic_chamber_params.json`: Configuration parameters used for generation.
*   `*.h5`: (In some versions) HDF5 files containing raw data and metadata.

## Dependencies

*   **Julia:** `ImagePhantoms`, `Sinograms`, `ImageGeoms`, `MIRTjim`, `Unitful`, `HDF5`, `PyCall`.
*   **Python:** `SimpleITK`, `numpy`, `wandb` (Weights & Biases for logging), `pydicom` (via `nifti_to_dicom_seg.py`).
*   **System:** `dcmtk` (for `nii2dcm` and `dcm2niix` utilities).
