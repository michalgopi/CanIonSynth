# Synthetic Tomography Data Generator

This repository creates synthetic phantoms and data for tomographic reconstruction. It generates NIfTI files and optionally DICOM and DICOM-SEG files.

## Setup

### Prerequisites

*   Julia 1.10+
*   Python 3.10+

### Installation

1.  **Install Python dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

2.  **Instantiate Julia environment:**

    ```bash
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    ```

## Running Tests

To run the full test suite (which also verifies the environment setup):

```bash
julia --project=. tests/run_tests.jl
```

## Running Generation Scripts

You can run the generation scripts directly. Ensure you activate the project environment.

**Example: Generate Can Phantom**

```bash
cd in_docker_organized
julia --project=.. main_create_phantom_can.jl 32x32x32 false false <uuid> false false 0.0
```

Arguments:
1.  Dimensions (e.g., `32x32x32`)
2.  Add Radon transform (`true`/`false`)
3.  Variable spacing (`true`/`false`)
4.  UUID (can be random string)
5.  Randomize (`true`/`false`)
6.  Add smoothing (`true`/`false`)
7.  Additive noise level (float)

## Repository Structure

*   `in_docker_organized/`: Main generation scripts (`main_create_phantom_can.jl`, etc.).
*   `tests/`: Test suite and configurations.
*   `packages/`: Vendored dependencies (e.g., patched `ImagePhantoms.jl`).
*   `Project.toml` / `Manifest.toml`: Julia dependencies lock files.
*   `requirements.txt`: Python dependencies.

## CI/CD

The repository uses GitHub Actions for CI. See `.github/workflows/ci.yml`.
