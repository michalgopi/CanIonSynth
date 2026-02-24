# System Architecture Overview

The Synthetic Tomography Data Generation System is designed to create realistic phantom models for industrial computed tomography (CT).

## Core Components

1.  **Coordinate Script (`coordinate_phantom_create.py`)**
    *   **Role**: Orchestrator.
    *   **Function**: Reads JSON configuration, selects the appropriate phantom script (Can or Ionic Chamber), passes parameters, and manages the execution flow.
    *   **Modes**: Folder-based generation (processing multiple JSONs) or Batch generation (generating N random phantoms).

2.  **Generation Scripts (Julia)**
    *   **`main_create_phantom_can.jl`**:
        *   Builds can-based phantoms.
        *   Features: Fluid layers, meniscus, balls, pipes, dispensers, steel (flat) vs aluminum (rounded) bottoms.
    *   **`main_create_phantom_ionic_chamber.jl`**:
        *   Builds ionic chamber phantoms.
        *   Features: Multi-layer electrode structures (Copper, Insulator, Aluminum, Graphite), various shapes (Ball, Cylinder, Lollipop).

3.  **CT Simulation (`get_approximate_radon_inverse.py`)**
    *   **Role**: Simulation of acquisition and reconstruction.
    *   **Function**: Performs forward projection (Radon transform) to create sinograms and inverse reconstruction (Filtered Back-Projection) to simulate the CT imaging chain.

4.  **DICOM Conversion (`save_dicom.py`, `nifti_to_dicom_seg.py`)**
    *   **Role**: Data formatting.
    *   **Function**: Converts NIfTI volumetric data into standard DICOM slices and DICOM-SEG objects, preserving medical metadata (Window/Level, Rescale Slope/Intercept).

5.  **Docker Environment**
    *   **Role**: Reproducibility.
    *   **Function**: Encapsulates all dependencies (Julia 1.10, Python 3.10, ADRT, SimpleITK) to ensure consistent execution across platforms.

## Workflow

1.  **Configuration**: User provides a JSON file or parameters.
2.  **Orchestration**: `coordinate_phantom_create.py` parses config and calls Julia.
3.  **Generation**: Julia scripts create geometric primitives, apply boolean operations, and generate voxel grids.
4.  **Post-Processing**:
    *   Noise injection (Gaussian/Uniform).
    *   Smoothing (Gaussian blur).
    *   (Optional) Radon Transform cycle.
5.  **Output**: NIfTI files, DICOM series, and JSON metadata are saved to disk.
