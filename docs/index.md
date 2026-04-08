# Synthetic Tomography Documentation

This documentation covers the Synthetic Tomography Data Generator, a Julia and Python workflow for generating synthetic industrial CT phantoms, derived reconstruction artifacts, and validation outputs.

## Main Pages

- [User Guide](user_guide.md)
- [Workflow Summary](workflow_summary.md)
- [Functions Reference](functions_reference.md)
- [System Architecture](system_architecture.md)
- [Volume Calculation](volume_calculation.md)

## Validation

The repository includes automated tests for both can and ionic chamber generation workflows. Run them from the repository root with:

```bash
julia --project=. tests/run_tests.jl
```
