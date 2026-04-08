# Synthetic Tomography Data Generator

Synthetic Tomography Data Generator produces synthetic industrial CT datasets for phantom-based metrology and reconstruction workflows. The repository combines Julia geometry generation with Python-based format conversion and reconstruction helpers.

## Key Capabilities

- Generate parameterized can phantoms with fluids, meniscus effects, internal objects, and optional defects.
- Generate parameterized ionic chamber phantoms with layered materials and multiple tip geometries.
- Export NIfTI volumes and masks, with optional DICOM and DICOM-SEG conversion when external tools are available.
- Optionally run forward and inverse tomographic steps for reconstruction-oriented experiments.
- Reproduce runs from JSON configurations and validate the workflows with an automated test suite.

## Repository Layout

- `in_docker_organized/`: Main Julia and Python generation workflows.
- `tests/`: Automated tests, JSON fixtures, and comparison helpers.
- `docs/`: User guide, workflow notes, architecture, and verification material.
- `examples/`: Example entry points and helper scripts.
- `packages/`: Vendored dependencies, including the patched `ImagePhantoms.jl` copy used by the project.
- `Project.toml` and `Manifest.toml`: Julia environment definition and lockfile.
- `requirements.txt`: Python dependencies.
- `Dockerfile` and `.devcontainer/`: Reproducible containerized environment.

## Requirements

- Julia 1.10+
- Python 3.10+
- Optional system tools: `nii2dcm` for DICOM conversion and `gcloud` for uploads

## Installation

1. Install Python dependencies:

```bash
pip install -r requirements.txt
```

2. Instantiate the Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

3. If you want a fully reproducible environment, build the provided container instead of configuring the host manually:

```bash
docker build -t synthetic-tomo .
```

## Quick Start

Run the main test suite:

```bash
julia --project=. tests/run_tests.jl
```

Generate a can phantom from the source directory:

```bash
julia --project=. in_docker_organized/main_create_phantom_can.jl 32x32x32 false false run-001 false false 0.0
```

Generate an ionic chamber phantom:

```bash
julia --project=. in_docker_organized/main_create_phantom_ionic_chamber.jl 32x32x32 false false run-002 false false 0.0 tests/configs/ionic_chamber.json
```

## Reproducibility

- The Julia dependency graph is locked in `Manifest.toml`.
- Python dependencies are listed in `requirements.txt`.
- A containerized environment is provided through `Dockerfile` and `.devcontainer/`.
- Test runs can be made deterministic with `FIXED_UUIDS=1`.
- Upload and experiment tracking can be disabled with `SKIP_UPLOAD=true` and `SKIP_WANDB=true`.

## Output Artifacts

Typical runs produce:

- `*.nii.gz` phantom volumes and component masks
- `argss.json` or `ionic_chamber_params.json` parameter logs
- optional projection or reconstruction artifacts
- optional DICOM and DICOM-SEG outputs when external tooling is installed

## Documentation

- [User Guide](docs/user_guide.md)
- [Workflow Summary](docs/workflow_summary.md)
- [System Architecture](docs/system_architecture.md)
- [Phantom Types](docs/phantom_types.md)
- [Volume Calculation](docs/volume_calculation.md)
- [Manual Testing and Visualization](docs/manual_testing_visualization.md)
- [Manual Verification Guide](docs/manual_verification.md)
- [SoftwareX Submission Notes](SOFTWAREX_SUBMISSION.md)

## Citation

Citation metadata for software archiving and paper preparation is provided in `CITATION.cff`.

## License

This repository is licensed under the MIT License. See `LICENSE`.

## Publishing Notes

The repository includes SoftwareX-oriented citation and archive metadata and a repository-level open-source license. The remaining publication steps are release-oriented rather than repository-structure blockers.

## CI

GitHub Actions workflows are defined in `.github/workflows/ci.yml` and `.github/workflows/docs.yml`.
