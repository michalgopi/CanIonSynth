# Synthetic Tomography Documentation

This documentation covers installation, execution, validation, and reproducibility for the Synthetic Tomography Data Generator.

## Start Here

1. Read the repository `README.md` for the shortest setup path.
2. Use the [User Guide](user_guide.md) for detailed installation and usage instructions.
3. Use the [Workflow Summary](workflow_summary.md) and [System Architecture](system_architecture.md) when you need to understand how the repository is organized internally.

## Common Commands

Set up the environment:

```bash
pip install -r requirements.txt
julia --project=. tests/setup_env.jl
```

Run the automated tests:

```bash
julia --project=. tests/run_tests.jl
```

Build the documentation site:

```bash
julia --project=docs docs/make.jl
```

## Documentation Map

- [User Guide](user_guide.md)
- [Workflow Summary](workflow_summary.md)
- [System Architecture](system_architecture.md)
- [Phantom Types](phantom_types.md)
- [Functions Reference](functions_reference.md)
- [Volume Calculation](volume_calculation.md)
- [Manual Testing and Visualization](manual_testing_visualization.md)
- [Manual Verification Guide](manual_verification.md)
