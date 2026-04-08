# Function Reference

This document is a code map for the project files that matter most when you need to extend, debug, or review the repository.

## Main Generation Entry Points

### `in_docker_organized/main_create_phantom_can.jl`

Role: top-level can phantom workflow.

Important functions:

- `json_based_can(ig, json_path, is_2d=true, is_debug=false)`: builds a can phantom from a JSON configuration.
- `random_can(ig, image_size_cm, is_2d, seed, spacing, is_debug=false)`: builds a randomized can parameter set.
- `create_boolean_density_dict(name_type_list, density_map)`: converts geometry objects into masks and density volumes.
- `get_random_can_uploaded(is_2d, seed, uuid=nothing)`: full generation, export, zip, and optional upload workflow.
- `calculate_and_save_fluid_volume(seed, is_2d=false, uuid=nothing)`: wrapper around the main can workflow.

### `in_docker_organized/main_create_phantom_ionic_chamber.jl`

Role: top-level ionic chamber workflow.

Important functions:

- `json_params(args_json_path, uuid=nothing, dims=(128, 128, 128))`: loads or completes a chamber parameter dictionary from JSON.
- `generate_random_params(uuid, dims, randomize=true, variable_spacing=false)`: builds a random parameter set.
- `create_ionic_chamber_phantom(params)`: assembles the chamber geometry and voxelized output.
- `build_image(base_params, specs, angle, spacing, dims)`: converts cylinder specifications into a final chamber volume.
- `get_random_chamber(dims, uuid, temp_fold, variable_spacing, randomize)`: full generation, export, zip, and optional upload workflow.

## Geometry Construction Files

### `in_docker_organized/get_geometry_main.jl`

Role: can geometry definition and placement helpers.

Important functions:

- `empty_cylinder_with_half_sphere_bottom_p(...)`: builds the main can geometry, fluid regions, and internal objects.
- `create_top_cut_objects(...)`: creates the cut geometry used for fluid surfaces and related features.
- `find_ball_positions(...)`: computes sphere placement inside the can.
- `volume_of_elliptical_cylinder(pipe_cross_section, overlay_length)`: helper for analytical volume calculations.

### `in_docker_organized/get_rounded_bottom_b.jl`

Role: rounded-bottom geometry support.

Important functions:

- `create_sphere(...)`: helper for curved bottom elements.
- `create_torus(...)`: helper for toroidal rounded-bottom geometry.
- `calculate_torus_volumes(...)`: analytical support for rounded-bottom volume estimates.
- `get_rounded_bottom(...)`: builds the lower can geometry for rounded-bottom variants.

## Utility And Interop Files

### `in_docker_organized/geometry_utils.jl`

Role: shared image manipulation and export helpers.

Important functions:

- `rotation3d(image, axis, theta)`: rotates SimpleITK images.
- `save_sitk_image_as_dicom(img, output_folder)`: attempts DICOM slice export through `nii2dcm`.
- `save_nifti_with_meta(arr, cast_to_uint8, spacing, output_path)`: writes NIfTI outputs with spatial metadata.
- `save_mask_as_nifti(mask, output_path, spacing)`: writes boolean masks as NIfTI.
- `convert_nifti_to_dicom_seg(nifti_path, reference_dicom_path, output_folder, reference_nifti_path="")`: wraps the Python DICOM-SEG helper.
- `move_image(...)`: translates masks in voxel space while respecting physical spacing.

### `in_docker_organized/volume_integration.jl`

Role: analytical and numerical volume comparison for the can workflow.

Important functions:

- `compute_fluid_volume_in_can(...)`: baseline analytical and numerical volume comparison.
- `compute_fluid_volume_in_can_v2(...)` and `compute_fluid_volume_in_can_v3(...)`: intermediate calculation variants.
- `compute_accurate_fluid_volume(...)`: higher-fidelity analytical comparison logic.
- `compute_accurate_fluid_volume_fixed(...)`: current detailed implementation used by the main workflow.

## Python Helpers

### `in_docker_organized/get_approximate_radon_inverse.py`

Role: approximate projection and reconstruction helper used when `add_radon=true`.

Important functions:

- `iadrt_cg(...)`: iterative inverse helper
- `next_power_of_two(n)`: array sizing helper
- `pad_to_power_of_two(arr)`: preprocessing helper

### `in_docker_organized/nifti_to_dicom_seg.py`

Role: converts a NIfTI segmentation mask into a DICOM-SEG output using a reference DICOM series.

Important functions:

- `parse_args()`
- `validate_inputs(...)`
- `read_nifti_mask(...)`
- `load_reference_dicoms(...)`
- `ensure_spacing_consistency(...)`
- `convert_nifti_to_dicom_seg(...)`

### `scripts/visualize_nifti.py`

Role: saves a PNG with orthogonal slices and an intensity histogram for a NIfTI file.

Important function:

- `visualize_nifti(nifti_path, output_path=None)`

## Testing And Documentation Entry Points

### `tests/setup_env.jl`

Role: instantiates the Julia project and binds `PyCall` to the active virtual environment.

### `tests/run_tests.jl`

Role: end-to-end integration test harness for both main generators. The harness also checks that expected output files are produced and optionally preserves them with `SAVE_OUTPUTS_TO`.

### `docs/make.jl`

Role: builds the documentation site from the Markdown files in `docs/`.
