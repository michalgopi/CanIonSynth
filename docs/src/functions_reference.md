# Function Reference

This document provides detailed documentation for the functions in the repository, organized by source file.

## 1. `in_docker_organized/main_create_phantom_can.jl`

This script is the entry point for generating "can" phantoms. It defines functions to create, configure, and save the phantom data.

### `json_based_can(ig, json_path, is_2d=true, is_debug=false)`
Generates a "can" phantom based on parameters loaded from a JSON file.
*   **Arguments:**
    *   `ig`: `ImageGeom` object defining the image geometry (dimensions, spacing).
    *   `json_path`: Path to the JSON configuration file.
    *   `is_2d`: Boolean, if true, forces certain angles to 0 (for 2D slice generation).
    *   `is_debug`: Boolean, enables debug print statements.
*   **Returns:** A tuple containing:
    *   `vol`: Dictionary with volume information (e.g., fluid volume).
    *   `args`: List of tuples representing the phantom parameters.
    *   `density_map`: Dictionary mapping object names to their densities.
    *   `name_type_list`: List of object names and types.
    *   Boolean flags for geometry features (double bottom, rounded bottom, etc.).

### `random_can(ig, image_size_cm, is_2d, seed, spacing, is_debug=false)`
Generates a "can" phantom with randomized parameters within defined constraints.
*   **Arguments:**
    *   `ig`: `ImageGeom` object.
    *   `image_size_cm`: Tuple of image dimensions in cm.
    *   `is_2d`: Boolean for 2D mode.
    *   `seed`: Integer seed for random number generation.
    *   `spacing`: Voxel spacing (dx, dy, dz).
    *   `is_debug`: Boolean for debug output.
*   **Returns:** Same structure as `json_based_can`.

### `get_cylinder_bool_mask(ob_el)`
Converts a cylinder object (or list of objects) into a boolean mask.
*   **Arguments:** `ob_el`: Phantom object(s).
*   **Returns:** 3D BitArray (mask).

### `get_half_s_bool(ob_el)`
Converts a half-sphere object into a boolean mask.
*   **Arguments:** `ob_el`: Phantom object(s).
*   **Returns:** 3D BitArray (mask).

### `get_temp_folder()`
Creates a temporary directory for storing intermediate files.
*   **Returns:** Path to the temporary directory.

### `save_args(args, vol, main_folder)`
Saves the generation arguments and volume data to a `argss.json` file in the specified folder.
*   **Arguments:**
    *   `args`: List of parameter tuples.
    *   `vol`: Volume data.
    *   `main_folder`: Output directory path.

### `create_boolean_density_dict(name_type_list, density_map)`
Constructs a dictionary of phantom objects with their boolean masks and density-scaled float masks.
*   **Arguments:**
    *   `name_type_list`: List of (name, object, type) tuples.
    *   `density_map`: Dictionary of densities.
*   **Returns:** Dictionary where keys are object names and values are tuples `(bool_mask, float_mask, original_object)`.

### `add_noise_at(phantoms_dict, key, density_map)`
Adds Gaussian noise to a specific object's region in the phantom.
*   **Arguments:**
    *   `phantoms_dict`: Dictionary of phantom objects.
    *   `key`: Key of the object to add noise to.
    *   `density_map`: Density map.
*   **Returns:** Float array with noise added to the object's region.

### `get_random_can_uploaded(is_2d, seed, uuid=nothing)`
Main workflow function to generate, save, and upload a random can phantom.
*   **Arguments:**
    *   `is_2d`: Boolean flag.
    *   `seed`: Random seed.
    *   `uuid`: Unique identifier (optional).
*   **Returns:** Tuple `(numerical_vol, analytical_vol, fluid_volume_numerical_my)` representing different volume calculations.

### `calculate_and_save_fluid_volume(seed, is_2d=false, uuid=nothing)`
Wrapper function that calls `get_random_can_uploaded`.

## 2. `in_docker_organized/main_create_phantom_ionic_chamber.jl`

This script generates "ionic chamber" phantoms.

### `cylinder_interval(c_z, l_z)`
Helper to calculate the z-range of a cylinder.
*   **Returns:** `(z_start, z_end)`.

### `overlaps(c1, l1, c2, l2)`
Checks if two z-intervals overlap.
*   **Returns:** Boolean.

### `compute_adjusted_density(new_desired, new_z, new_lz, existing)`
Calculates the effective density for a new cylinder to account for overlaps with existing cylinders (density superposition).
*   **Returns:** Adjusted density value.

### `create_cylinder(base_x, base_y, base_bottom_z, spec, angle, name)`
Creates a cylinder parameter NamedTuple from a specification dictionary.
*   **Returns:** NamedTuple with cylinder parameters.

### `adjust_densities!(cylinders)`
Iteratively adjusts the densities of a list of cylinders to ensure the final superposition matches the desired densities.
*   **Returns:** Vector of adjusted cylinder NamedTuples.

### `build_image(base_params, specs, angle, spacing, dims)`
Constructs the full phantom image from base parameters and component specifications.
*   **Returns:** Tuple containing `(cylinders_prim, named_res, cylinder_densities, vol_res)`.

### `create_ionic_chamber_phantom(params)`
Main function to generate the ionic chamber phantom. It determines the specific geometry (square top, ball-like, etc.) and assembles the components.
*   **Arguments:** `params`: Dictionary of parameters.
*   **Returns:** Tuple containing the final image array, raw image, masks, and volume data.

### `save_ionic_chamber_params(params, filename)`
Saves the parameter dictionary to a JSON file.

### `json_params(args_json_path, uuid=nothing, dims=(128, 128, 128))`
Loads parameters from a JSON file, filling in defaults where necessary.

### `generate_random_params(uuid, dims, randomize=true, variable_spacing=false)`
Generates a random parameter set for an ionic chamber, selecting a random type (square top, ball-like, etc.) and randomizing dimensions.
*   **Returns:** Parameter dictionary.

### `get_random_chamber(dims, uuid, temp_fold, variable_spacing, randomize)`
Main workflow function to generate, save, and upload an ionic chamber phantom.

## 3. `in_docker_organized/get_geometry_main.jl`

Contains the geometry definitions for the "can" phantom.

### `ionic_chamber_p(...)`
Defines the geometry of a parameterized ionic chamber using basic shapes. (Note: This seems to be an alternative or older implementation compared to the main script).

### `empty_cylinder_with_half_sphere_bottom_p(...)`
Constructs the components of a "can" phantom. This is a complex function that defines:
*   The main cylinder body.
*   The bottom (flat or rounded with half-spheres).
*   The fluid (single or dual phase) with meniscus.
*   The cut plane (tilt) at the top of the fluid.
*   Internal objects (pipe, dispenser).
*   Defects or "cuts".
*   **Returns:** A tuple containing the list of objects, volume dictionary, and individual component masks.

### `volume_of_elliptical_cylinder(pipe_cross_section, overlay_length)`
Calculates the volume of an elliptical cylinder.

## 4. `in_docker_organized/geometry_utils.jl`

Utility functions for image manipulation and saving.

### `matrix_from_axis_angle(a)`
Computes a rotation matrix from an axis-angle representation.

### `resample(image, transform)`
Resamples a SimpleITK image using a given transform.

### `rotation3d(image, axis, theta)`
Rotates a 3D image around a specific axis by `theta` degrees.

### `save_sitk_image_as_dicom(img, output_folder)`
Saves a SimpleITK image as a series of DICOM files. Uses `nii2dcm` utility.

### `get_per_slice_reconstruction(arr)`
Performs a slice-by-slice reconstruction (inverse Radon transform) using `skimage.transform.iradon`.

### `save_2d_sitk_image_as_dicom(immm, output_path)`
Saves a 2D SimpleITK image as a DICOM file, handling intensity scaling.

### `save_nifti_with_meta(arr, cast_to_uint8, spacing, output_path)`
Saves an array as a NIfTI file with specified metadata (spacing, direction, origin).

### `save_mask_as_nifti(mask, output_path, spacing)`
Wrapper to save a boolean mask as a NIfTI file.

### `convert_nifti_to_dicom_seg(nifti_path, reference_dicom_path, output_folder, reference_nifti_path="")`
Converts a NIfTI mask to DICOM-SEG format by calling a Python script (`nifti_to_dicom_seg.py`).

### `move_image(array, axis, distance, direction, spacing)`
Shifts a binary image along a specified axis by a given distance (in cm).

## 5. `in_docker_organized/volume_integration.jl`

Functions for analytical volume calculation of the fluid in the can phantom.

### `compute_fluid_volume_in_can(phantoms_dict, spacing, first_ball, second_ball, rounded_bottom, params)`
Calculates the fluid volume using both numerical (sum of voxels) and analytical (geometric formulas) methods.
It accounts for:
*   Base cylinder volume.
*   Top cut plane (tilt).
*   Bottom geometry (rounded or flat).
*   Subtracted objects (balls, pipe overlap).
*   **Returns:** Tuple `(numerical_vol, analytical_vol)`.

### `compute_accurate_fluid_volume(phantoms_dict, spacing, first_ball, second_ball, rounded_bottom, params, ig)`
A more precise version of volume calculation, handling complex interactions between components.

### `compute_accurate_fluid_volume_fixed(...)`
Fixed version of the accurate volume calculation logic.

## 6. `CtFanArc_params.jl`

Defines CT geometry and projection functions.

### `get_CTFAN_proj(ob, is_high_res=false)`
Configures a Fan Beam CT geometry and generates projections (sinogram) for a given object.
*   **Arguments:**
    *   `ob`: Phantom object(s) to project.
    *   `is_high_res`: Boolean, selects between standard and high-resolution geometry settings.
*   **Returns:** Tuple `(proj_arc, rg)`, where `proj_arc` is the projection data and `rg` is the geometry object.
