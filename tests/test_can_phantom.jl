
using Pkg
using Pkg;
Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git");

using Sinograms: SinoPar, rays, plan_fbp, fbp, fbp_sino_filter, CtFanArc, CtFanFlat, Window, Hamming, fdk, ct_geom_plot3, project_bdd, backproject_bdd
using ImageGeoms: ImageGeom, fovs, MaskCircle, axesf
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom, Object, spectrum, Cylinder, cylinder, ellipsoid_parameters, ellipsoid, Ellipsoid, Cuboid
using Unitful: mm, unit, °, cm
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Plots: plot, plot!, scatter!, default, gui
using PyCall
using Revise, Statistics
using UUIDs
using JSON
using Dates, Logging
using ImagePhantoms
import ImagePhantoms
using Accessors
using ImageFiltering
using Random

# Mock WandB
struct MockWandB
    init::Function
    log::Function
end
const wandb = MockWandB((args...; kwargs...) -> nothing, (args...; kwargs...) -> nothing)

# Import external packages with pyimport_conda
sitk = pyimport_conda("SimpleITK", "simpleitk")
np = pyimport_conda("numpy", "numpy")
# wandb = pyimport_conda("wandb", "wandb") # Disabled for testing
skimage = pyimport_conda("skimage", "skimage")
adrt = pyimport_conda("adrt", "adrt")

# wandb.init(project="synth") # Disabled
isinteractive() ? jim(:prompt, true) : prompt(:draw)

# Adjust includes to point to the correct location relative to tests/
include.([
    "../in_docker_organized/get_geometry_main.jl",
    "../in_docker_organized/geometry_utils.jl"
])

include("../in_docker_organized/get_rounded_bottom_b.jl")
const global ImagePhantoms.IRREGULARITY = 0.3

# Import volume calculation function
include("../in_docker_organized/volume_integration.jl")

# Test Configuration
max_len = 10.0
max_radius = 5.0
is_2d = false # 3D for test
dims = (32, 32, 32) # Small dimensions for fast execution
max_z = (max_len) * 1.6
max_y_x = max_radius
spacing = (max_y_x / dims[1], max_y_x / dims[2], max_z / dims[3])
spacing = (maximum(spacing), maximum(spacing), maximum(spacing))

ig = ImageGeom(dims=dims, deltas=(spacing[1]cm, spacing[2]cm, spacing[3]cm))

# Modified function to skip JSON and use defaults or random
function test_random_can(ig, image_size_cm, is_2d, seed, spacing, is_debug=false)
    # ... (Reuse logic from random_can, but simplified if needed, or just call random_can if we can import it)
    # Since random_can is in main_create_phantom_can.jl which is a script, we can't easily import it without running it.
    # So we will copy the necessary parts or include the file if it was a module.
    # Given the structure, we'll redefine a minimal version or assume we can include the main file if it was cleaner.
    # BUT, main_create_phantom_can.jl has top-level execution code.

    # Best approach: Copy the critical functions `random_can`, `get_random_can_uploaded` (modified) here.
    # To avoid huge duplication, I will rely on the `includet` above which brings in geometry functions.
    # I need `random_can` logic.

    # Let's define a minimal test workflow that calls `empty_cylinder_with_half_sphere_bottom_p` directly.

    center_cylinder = (0.0, 0.0, 0.0)
    bigger_cyl_size = [2.0, 2.0, 5.0]
    cylinder_wall_thickness = 0.1
    cylinder_bottom_curvature = 0.5
    cylinder_top_curvature = 0.5
    angle = 0.0
    density_inside = 0.5
    pipe_len = 3.0
    pipe_cross_section = (0.5, 0.5)
    pipe_density = 0.8
    dispenser_len = 1.0
    dispenser_cross_section = (0.6, 0.6)
    dispenser_density = 0.9
    len_cut = 1.0
    menisc_radius = 2.0
    dual_phase_percentage = 0.5
    density_inside_b = 0.4
    x_cut_angle = 0.0
    y_cut_angle = 0.0
    cylinder_bottom_curvature_b = 0.6
    rel_size_bottom_curvature = 0.5
    menisc_radius_mult = 1.0
    menisc_cut_height = 0.5
    menisc_cut_total_height = 1.0
    rounded_bottom = true
    curvature = 0.1
    tube_ball_fraction = 0.8
    rel_dist = 0.5
    torus_height_factor = 0.1
    r_xy_small_factor = 0.2
    r_z_small_factor = 0.02

    vol, density_map, name_type_list = empty_cylinder_with_half_sphere_bottom_p(
        center_cylinder,
        bigger_cyl_size,
        cylinder_wall_thickness,
        cylinder_bottom_curvature,
        cylinder_top_curvature,
        angle,
        density_inside,
        pipe_len,
        pipe_cross_section,
        pipe_density,
        dispenser_len,
        dispenser_cross_section,
        dispenser_density,
        len_cut,
        menisc_radius,
        dual_phase_percentage,
        density_inside_b,
        x_cut_angle,
        y_cut_angle,
        cylinder_bottom_curvature_b,
        rel_size_bottom_curvature,
        menisc_radius_mult,
        menisc_cut_height,
        menisc_cut_total_height,
        rounded_bottom,
        curvature,
        tube_ball_fraction,
        rel_dist,
        torus_height_factor, r_xy_small_factor, r_z_small_factor
    )

    return vol, density_map, name_type_list
end

function run_test()
    println("Starting Can Phantom Test...")
    seed = 1234

    vol, density_map, name_type_list = test_random_can(ig, (dims[1] * spacing[1], dims[2] * spacing[2], dims[3] * spacing[3]), is_2d, seed, spacing)

    println("Phantom generated.")

    # Generate dictionary of phantoms
    phantoms_dict = create_boolean_density_dict(name_type_list, density_map)

    # Create result tensor
    result_tensor = zeros(dims)

    # Simple summation for testing (ignoring complex masking logic for speed)
    for (key, val) in phantoms_dict
        result_tensor += val[2]
    end

    println("Tensor created with size: ", size(result_tensor))

    # Save output
    output_dir = "tests/output_can"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    output_path = joinpath(output_dir, "test_can.nii.gz")
    save_nifti_with_meta(result_tensor, false, spacing, output_path)

    println("Saved NIfTI to ", output_path)

    if isfile(output_path)
        println("Test PASSED: Output file exists.")
    else
        println("Test FAILED: Output file not found.")
    end
end

run_test()
