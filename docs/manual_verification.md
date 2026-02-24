# Manual Verification Guide

This guide explains how to manually verify the functionalities of the Synthetic Tomography Data Generation System, based on the system design and implementation details.

## 1. Verify Can Phantoms

### Steel Can (Flat Bottom)
*   **Action**: Generate a can phantom with parameters conducive to a flat bottom.
    *   Set `rounded_bottom` (or equivalent parameter in config) to `false`.
    *   Set `cylinder_wall_thickness` to approx `0.02` cm.
*   **Verification**:
    *   **Visual**: Open the generated NIfTI file (e.g., `example_can.nii.gz`) in a viewer (e.g., ITK-SNAP, 3D Slicer).
    *   **Check**: In the coronal or sagittal view, verify the bottom of the can is flat (sharp corners) relative to the walls.
    *   **Check**: Verify the wall thickness is thin relative to the can size.

### Aluminum Can (Rounded Bottom)
*   **Action**: Generate a can phantom with parameters for a rounded bottom.
    *   Set `rounded_bottom` to `true`.
    *   Set `cylinder_wall_thickness` to approx `0.03` cm.
*   **Verification**:
    *   **Visual**: Open the generated NIfTI file.
    *   **Check**: In the coronal/sagittal view, verify the bottom has a curved transition (torus geometry) between the base and the walls.
    *   **Check**: Look for the "dome" or "outer sphere" curvature at the very bottom center.

### Fluid and Meniscus
*   **Action**: Ensure `menisc_cut_height` is > 0 in the configuration.
*   **Verification**:
    *   **Visual**: Inspect the fluid surface.
    *   **Check**: Verify the fluid surface is curved (concave) near the walls, representing the meniscus.
    *   **Check**: If `x_cut_angle` or `y_cut_angle` are non-zero, verify the fluid surface is tilted.

### Dual-Phase Fluid
*   **Action**: Enable `two_fluids=true` (or ensure default is used).
*   **Verification**:
    *   **Visual**: Inspect the coronal view.
    *   **Check**: Identify two distinct intensity levels within the fluid volume. The lower layer (secondary fluid) should have lower density (intensity) than the upper layer.

## 2. Verify Ionic Chambers

### Ball-shaped Chamber
*   **Action**: Run generation with `ball_like=true` (or `rand_ver` mapping to ball shape).
*   **Verification**:
    *   **Visual**: Open `ionic_chamber.nii.gz`.
    *   **Check**: Verify the top section is spherical or ellipsoidal.
    *   **Check**: Identify concentric layers: Copper (center) -> Insulator -> Aluminum -> Insulator -> Aluminum -> Graphite (outer).

### Rounded-top Cylinder
*   **Action**: Run generation for the default cylindrical type with rounded top.
*   **Verification**:
    *   **Visual**: Check for a cylindrical base capped with a half-ellipsoid.

### Lollipop-shaped Chamber
*   **Action**: Run generation with `lollipop_like=true`.
*   **Verification**:
    *   **Visual**: Check for a thin stem connected to a wider, flat cylindrical head.

### Standardized Sizes
*   **Action**: Enable `new_flat_sizes=true` and vary `rand_ver` (1, 2, or 3).
*   **Verification**:
    *   **Quantitative**: Measure the dimensions of the generated chamber.
        *   Ver 1: Radius ~20mm, Height ~48mm.
        *   Ver 2: Radius ~15mm, Height ~42.9mm.
        *   Ver 3: Radius ~10mm, Height ~32.5mm.

## 3. Verify Radon Transform

*   **Action**: Set `add_radon=true` in the configuration.
*   **Verification**:
    *   **Files**: Check for the existence of `after_radon.nii.gz` and `after_radon_plus_before.nii.gz`.
    *   **Visual**: Open `after_radon.nii.gz`. It should look like a reconstructed CT image (possibly with streak artifacts or noise characteristic of filtered back-projection).
    *   **Visual**: Open `after_radon_plus_before.nii.gz`. It should show the average of the original phantom and the reconstructed one.

## 4. Verify Volume Calculations

*   **Action**: Run any generation script and capture the standard output (stdout) or check `argss.json`.
*   **Verification**:
    *   **Check Logs**: Look for "CYLINDER VOLUME DEBUG" or "Main execution completed with calculated fluid volumes".
    *   **Compare**: Check the "Numerical" volume (voxel count) vs "Analytical" volume.
    *   **Criterion**: The relative difference should be small (< 5% typically, though complex geometries might be higher). Large discrepancies indicate potential issues in mask generation or analytical formulas.

## 5. Verify Output Formats

### DICOM Generation
*   **Action**: Ensure `save_dicom.py` is executed (usually automatic if dependencies are met).
*   **Verification**:
    *   **Files**: Check for a directory (e.g., `example_can/`) containing `.dcm` files.
    *   **Metadata**: Open a `.dcm` file and check tags like `Modality` (should be "CT"), `RescaleSlope`, `RescaleIntercept`.

### Segmentation Masks
*   **Action**: Check the output directory.
*   **Verification**:
    *   **Files**: Verify existence of `fluid_mask.nii.gz`, `pipe_mask.nii.gz`, etc.
    *   **Visual**: Overlay these masks on the main phantom image to verify alignment.
