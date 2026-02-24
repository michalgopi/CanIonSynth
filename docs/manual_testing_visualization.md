# Manual Testing and Visualization Guide

This guide describes how to run a comprehensive set of manual tests to generate various phantom types and visualize the resulting NIfTI images using the provided helper tools.

## Prerequisites

Ensure you have set up the environment as described in the README:

```bash
pip install -r requirements.txt
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

**Note:** Generation scripts rely on environment variables to skip cloud uploads during testing:
```bash
export SKIP_UPLOAD=true
export SKIP_WANDB=true
```

## Helper Script: `visualize_nifti.py`

Since this environment is headless (no GUI), we provide a script to generate visualization images (PNG) from the output NIfTI files.

**Usage:**
```bash
python3 scripts/visualize_nifti.py <path_to_nifti_file>
```
This will create a `*_vis.png` image in the same directory as the NIfTI file, showing Axial, Coronal, and Sagittal views along with an intensity histogram.

---

## Test Sets

The following sections provide commands to generate specific phantom configurations. Run these commands from the root of the repository.

**Important:** The generation scripts store output in a temporary directory and print the path at the end (e.g., `Output stored in: /tmp/...`). You will need to use this path for visualization.

### 1. Steel Can (Flat Bottom)

Generates a can with a flat bottom, characteristic of steel cans.

**Command:**
```bash
cd in_docker_organized
# 64x64x64 resolution for speed. Use 128 or 256 for high detail.
# Arguments: dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]
# We use a temporary config to enforce flat bottom
echo '{"rounded_bottom": false, "cylinder_wall_thickness": 0.02}' > temp_steel.json
julia --project=.. main_create_phantom_can.jl 64x64x64 false false "test_steel" false false 0.0 temp_steel.json
rm temp_steel.json
cd ..
```

**Verify:**
1.  **Output:** Look for the line `Output stored in: ...` in the console output. Let's assume it is `$OUTPUT_DIR`.
2.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/example_can.nii.gz
    ```
3.  **Check Image:** Open `$OUTPUT_DIR/example_can_vis.png`. Verify the bottom corners are sharp (90 degrees) in the Coronal/Sagittal views.

### 2. Aluminum Can (Rounded Bottom)

Generates a can with a rounded/domed bottom and two fluid layers.

**Command:**
```bash
cd in_docker_organized
echo '{"rounded_bottom": true, "cylinder_wall_thickness": 0.035, "two_fluids": true}' > temp_alum.json
julia --project=.. main_create_phantom_can.jl 64x64x64 false false "test_alum" false false 0.0 temp_alum.json
rm temp_alum.json
cd ..
```

**Verify:**
1.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/example_can.nii.gz
    ```
3.  **Check Image:** Verify the bottom has a curved transition (torus shape) and a central dome. Look for two distinct gray levels in the fluid.

### 3. Can with Radon Transform

Simulates CT acquisition (Sinogram generation) and reconstruction.

**Command:**
```bash
cd in_docker_organized
# add_radon = true
julia --project=.. main_create_phantom_can.jl 64x64x64 true false "test_radon" false false 0.0
cd ..
```

**Verify:**
1.  **Output:** Check for `after_radon.nii.gz` in the output directory.
2.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/after_radon.nii.gz
    ```
3.  **Check Image:** The image should look like a reconstructed CT scan, potentially with characteristic streak artifacts or noise compared to the "perfect" phantom.

### 4. Ionic Chamber - Ball Shaped

Generates a spherical ionic chamber.

**Command:**
```bash
cd in_docker_organized
echo '{"ball_like": true, "lollipop_like": false}' > temp_ball.json
julia --project=.. main_create_phantom_ionic_chamber.jl 64x64x64 false false "test_ball" false false 0.0 temp_ball.json
rm temp_ball.json
cd ..
```

**Verify:**
1.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/ionic_chamber.nii.gz
    ```
3.  **Check Image:** Confirm the top part of the object is spherical.

### 5. Ionic Chamber - Lollipop Shaped

Generates a chamber with a thin stem and wide head.

**Command:**
```bash
cd in_docker_organized
echo '{"lollipop_like": true}' > temp_lollipop.json
julia --project=.. main_create_phantom_ionic_chamber.jl 64x64x64 false false "test_lollipop" false false 0.0 temp_lollipop.json
rm temp_lollipop.json
cd ..
```

**Verify:**
1.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/ionic_chamber.nii.gz
    ```
3.  **Check Image:** Identify the thin neck/stem connecting to the head.

### 6. Ionic Chamber - Standardized Sizes

Generates a standardized chamber (Version 1: approx 60cm³).

**Command:**
```bash
cd in_docker_organized
echo '{"new_flat_sizes": true, "rand_ver": 1}' > temp_std.json
julia --project=.. main_create_phantom_ionic_chamber.jl 64x64x64 false false "test_std" false false 0.0 temp_std.json
rm temp_std.json
cd ..
```

**Verify:**
1.  **Visualize:**
    ```bash
    python3 scripts/visualize_nifti.py $OUTPUT_DIR/ionic_chamber.nii.gz
    ```
3.  **Check Image:** Verify it has a cylindrical shape with flat ends.

## Using the Example Scripts

For convenience, you can use the scripts in `examples/` which handle the temporary directory finding for you:

```bash
export SKIP_UPLOAD=true
export SKIP_WANDB=true
./examples/generate_steel_can.sh
```

This will run the generation and print the exact command to visualize the output.
