#!/usr/bin/env julia
__precompile__(false)

# Parse arguments for dims and add_radon
dims_str = split(ARGS[1], "x")
dims = (parse(Int, dims_str[1]),
        parse(Int, dims_str[2]),
        parse(Int, dims_str[3]))


add_radon = parse(Bool, ARGS[2])

variable_spacing = parse(Bool, ARGS[3])

uuid = ARGS[4]

randomize = parse(Bool, ARGS[5])

add_smooth = parse(Bool, ARGS[6])
additive_noise = parse(Float64, ARGS[7])



args_json_path = length(ARGS) >= 8 ? ARGS[8] : " "


print("\n dims $dims add_radon $add_radon add_smooth $add_smooth additive_noise $additive_noise \n")

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
using Random  # Add this line to import the Random module

# Import external packages with pyimport_conda
sitk = pyimport_conda("SimpleITK", "simpleitk")
np = pyimport_conda("numpy", "numpy")

if !haskey(ENV, "SKIP_WANDB")
    wandb = pyimport_conda("wandb", "wandb")
    wandb.init(project="synth")
else
    struct MockWandB
        init::Function
        log::Function
    end
    wandb = MockWandB((args...; kwargs...) -> nothing, (args...; kwargs...) -> nothing)
end

try
    skimage = pyimport_conda("skimage", "skimage")
    adrt = pyimport_conda("adrt", "adrt")
catch e
    println("Warning: Could not import skimage or adrt. Proceeding without them.")
end

isinteractive() ? jim(:prompt, true) : prompt(:draw)
includet.([
    "get_geometry_main.jl"
    "geometry_utils.jl"])

includet("get_rounded_bottom_b.jl")
const global ImagePhantoms.IRREGULARITY = 0.3

# Import volume calculation function
includet("volume_integration.jl")





#get size of the generated image

max_len = 10.0
max_radius = 5.0
is_2d = true


# add_radon = false
# dims = (128, 128, 128)

# dims = (1024,1024 , 500)
#maximal size in z dim
max_z = (max_len) * 1.6
#maximal size in x or y  dim
max_y_x = max_radius
#now based on dims - so number of voxels and max dimensions we can calculate the spacing
spacing = (max_y_x / dims[1], max_y_x / dims[2], max_z / dims[3])
spacing = (maximum(spacing), maximum(spacing), maximum(spacing))

ig = ImageGeom(dims=dims, deltas=(spacing[1]cm, spacing[2]cm, spacing[3]cm))
# global_logger(lg)


function json_based_can(ig, json_path, is_2d=true, is_debug=false)
    # Load parameters from the JSON file
    params = JSON.parsefile(json_path)

    # Define the constants that were in the original function
    MAX_TILT_ANGLE = 0.05

    # Extract parameters from JSON with defaults if not present
    center_cylinder = get(params, "center_cylinder", (0.0, 0.0, 0.0))
    bigger_cyl_size = get(params, "bigger_cyl_size", [1.0, 1.0, 1.0])
    cylinder_wall_thickness = get(params, "cylinder_wall_thickness", 0.025)
    cylinder_bottom_curvature = get(params, "cylinder_bottom_curvature", 0.05 * bigger_cyl_size[3])
    cylinder_bottom_curvature_b = get(params, "cylinder_bottom_curvature_b", cylinder_bottom_curvature * 2.0)
    rel_size_bottom_curvature = get(params, "rel_size_bottom_curvature", 0.5)
    cylinder_top_curvature = get(params, "cylinder_top_curvature", 0.1 * bigger_cyl_size[3])

    # Extract angles
    angle = get(params, "angle", 0.0)
    x_cut_angle = get(params, "x_cut_angle", 0.0)
    y_cut_angle = get(params, "y_cut_angle", 0.0)

    # Force angles to 0 if 2D mode is enabled
    if is_2d
        x_cut_angle = 0.0
        y_cut_angle = 0.0
    end

    # Extract pipe parameters
    pipe_len = get(params, "pipe_len", 0.8 * (bigger_cyl_size[3] - cylinder_bottom_curvature + cylinder_top_curvature))
    pipe_cross_section = get(params, "pipe_cross_section", (0.15 * bigger_cyl_size[1], 0.15 * bigger_cyl_size[1]))
    pipe_density = get(params, "pipe_density", 0.3)

    # Extract dispenser parameters
    dispenser_len = get(params, "dispenser_len", 0.2 * bigger_cyl_size[3])
    dispenser_cross_section = get(params, "dispenser_cross_section", (0.2 * bigger_cyl_size[1], 0.2 * bigger_cyl_size[2]))
    dispenser_density = get(params, "dispenser_density", 0.325)

    # Extract fluid parameters
    len_cut = get(params, "len_cut", 0.3 * bigger_cyl_size[3])
    menisc_radius = get(params, "menisc_radius", bigger_cyl_size[2])  # Default to cylinder radius
    dual_phase_percentage = get(params, "dual_phase_percentage", 0.5)
    density_inside = get(params, "density_inside", 0.25)
    density_inside_b = get(params, "density_inside_b", 0.8 * density_inside)

    # Extract meniscus parameters
    menisc_radius_mult = get(params, "menisc_radius_mult", 1.5)
    menisc_cut_total_height = get(params, "menisc_cut_total_height", 0.2 * bigger_cyl_size[3])
    menisc_cut_height = get(params, "menisc_cut_height", 0.075 * bigger_cyl_size[3])

    # Extract feature flags
    rounded_bottom = get(params, "rounded_bottom", true)
    double_bottom_curvature = get(params, "double_bottom_curvature", true)
    add_pipe = get(params, "add_pipe", true)
    first_ball = get(params, "first_ball", false)
    second_ball = get(params, "second_ball", false)

    # Extract shape parameters
    curvature = get(params, "curvature", 0.04 * bigger_cyl_size[3])
    tube_ball_fraction = get(params, "tube_ball_fraction", 0.75)
    rel_dist = get(params, "rel_dist", 0.5)
    torus_height_factor = get(params, "torus_height_factor", 0.1)
    r_xy_small_factor = get(params, "r_xy_small_factor", 0.2)
    r_z_small_factor = get(params, "r_z_small_factor", 0.02)

    # Get spacing from JSON or use default
    spacing = get(params, "spacing", (0.05, 0.05, 0.05))

    # Initialize objects as empty arrays (same as in original function)
    ob4 = []
    ob5 = []
    ob4b = []
    ob5b = []
    ob_cut = []
    ob_menisc_cut = []
    ob5b_mask = []
    density_map = Dict{String,Float64}()
    name_type_list = []

    # Debug output similar to the original function
    if is_debug
        println("==== CYLINDER VOLUME DEBUG ====")

        # Calculate current parameters' volume
        curr_radius = bigger_cyl_size[1]
        curr_radius_inner = curr_radius - cylinder_wall_thickness
        curr_drop = curr_radius_inner * sqrt(2) * sin(max(abs(x_cut_angle), abs(y_cut_angle)))
        curr_cut_center = (bigger_cyl_size[3] / 2) - (len_cut / 2)
        curr_lowest_cut = curr_cut_center - curr_drop - (len_cut / 2)
        curr_cylinder_top = bigger_cyl_size[3] / 2
        curr_distance_from_top = curr_cylinder_top - curr_lowest_cut
        curr_height_below_cut = bigger_cyl_size[3] - curr_distance_from_top

        curr_calc_cyl = cylinder(
            0.0cm, 0.0cm,
            (-bigger_cyl_size[3] / 2 + curr_height_below_cut / 2)cm,
            curr_radius_inner * cm, curr_radius_inner * cm,
            curr_height_below_cut * cm,
            0.0, 0.0, 0.0, 1.0f0
        )
        curr_volume = IP.volume(curr_calc_cyl)
        curr_volume_scalar = Unitful.ustrip(cm^3, curr_volume)
        println("Current parameters: height=$(bigger_cyl_size[3]), radius=$(curr_radius), wall=$(cylinder_wall_thickness), len_cut=$(len_cut)")
        println("Current cylinder effective: radius=$(curr_radius_inner)cm, height=$(curr_height_below_cut)cm")
        println("Current cylinder volume: $(curr_volume_scalar) cm³")
        println("================================")
    end

    # Create the geometric object (same call as in the original function)
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
        torus_height_factor,
        r_xy_small_factor,
        r_z_small_factor
    )

    # Return the same structure as the original function
    return (vol, [("density_inside", density_inside),
            ("center_cylinder", center_cylinder),
            ("bigger_cyl_size", bigger_cyl_size),
            ("cylinder_wall_thickness", cylinder_wall_thickness),
            ("cylinder_bottom_curvature", cylinder_bottom_curvature),
            ("cylinder_top_curvature", cylinder_top_curvature),
            ("angle", angle),
            ("pipe_len", pipe_len),
            ("pipe_cross_section", pipe_cross_section),
            ("pipe_density", pipe_density),
            ("dispenser_len", dispenser_len),
            ("dispenser_cross_section", dispenser_cross_section),
            ("dispenser_density", dispenser_density / 2),
            ("len_cut", len_cut),
            ("dual_phase_percentage", dual_phase_percentage),
            ("x_cut_angle", x_cut_angle),
            ("y_cut_angle", y_cut_angle),
            ("density_inside_b", density_inside_b),
            ("cylinder_bottom_curvature_b", cylinder_bottom_curvature_b),
            ("rel_size_bottom_curvature", rel_size_bottom_curvature),
            ("double_bottom_curvature", double_bottom_curvature),
            ("menisc_radius_mult", menisc_radius_mult),
            ("menisc_cut_height", menisc_cut_height),
            ("menisc_cut_total_height", menisc_cut_total_height),
            ("rounded_bottom", rounded_bottom),
            ("curvature", curvature),
            ("tube_ball_fraction", tube_ball_fraction),
            ("rel_dist", rel_dist),
            ("add_pipe", add_pipe),
            ("torus_height_factor", torus_height_factor)],
        density_map, name_type_list, double_bottom_curvature, rounded_bottom, first_ball, second_ball, add_pipe)
end



# wandb.log(Dict("setting up wandb process"=>uuid_main ))

"""
    random_can(ig, image_size_cm, is_2d, seed, is_debug=false) -> Tuple

Generate a randomized cylindrical phantom container with fluid and other structures.

# Arguments
- `ig`: ImageGeom object for the phantom container
- `image_size_cm`: Size of the image in centimeters
- `is_2d`: Boolean indicating if 2D mode should be used
- `seed`: Seed for random number generation
- `is_debug`: Boolean to enable debugging output for cylinder volume

# Returns
A tuple containing:
- `vol`: Dictionary with volume information
- `args`: List of parameter tuples describing the phantom
- `density_map`: Dictionary mapping object names to density values
- `name_type_list`: List of objects and their types
- `double_bottom_curvature`: Boolean flag for bottom curvature type
- `rounded_bottom`: Boolean indicating if bottom is rounded
- `first_ball`: Boolean indicating if first ball is present
- `second_ball`: Boolean indicating if second ball is present
"""
function random_can(ig, image_size_cm, is_2d, seed, spacing, is_debug=false)
    # Define the exact minimum and maximum values from debugging calculations
    MIN_HEIGHT = 0.8 * max_len
    MIN_RADIUS_FACTOR = 0.25  # min_radius = min_height * 0.25



    MIN_LEN_CUT_FACTOR = 0.15 # min_len_cut = min_height * 0.45

    MAX_HEIGHT = 1.0 * max_len
    MAX_RADIUS_FACTOR = 0.3   # max_radius = max_height * 0.3
    MAX_LEN_CUT_FACTOR = 0.6 # max_len_cut = max_height * 0.25

    # Maximum tilt angle that affects volume calculations
    MAX_TILT_ANGLE = 0.05

    # Create a Mersenne Twister RNG with the specified seed for reproducibility
    rng = Random.MersenneTwister(seed)

    rounded_bottom = rand(rng) > 0.5
    two_fluids = true #rand(rng) > 0.5
    add_pipe = rand(rng) > 0.5
    #aluminium cans are thicker
    if (rounded_bottom)
        cylinder_wall_thickness = rand(rng, 0.028:0.001:0.044)
    else
        cylinder_wall_thickness = rand(rng, 0.019:0.001:0.025)

    end
    cylinder_wall_thickness = max(spacing[1], cylinder_wall_thickness)
    # Generate random parameters within the specified constraints
    center_cylinder = (0.0, 0.0, 0.0)

    # Use exact values from debug calculations for height
    bigger_cyl_size = [1.0, 1.0, 1.0]
    bigger_cyl_size[3] = MIN_HEIGHT + rand(rng) * (MAX_HEIGHT - MIN_HEIGHT)

    # Use exact radius factors
    radius_factor = MIN_RADIUS_FACTOR + rand(rng) * (MAX_RADIUS_FACTOR - MIN_RADIUS_FACTOR)
    bigger_cyl_size[1] = bigger_cyl_size[3] * radius_factor
    bigger_cyl_size[2] = bigger_cyl_size[1]

    # Use exact cut length factors
    cut_factor = MAX_LEN_CUT_FACTOR + rand(rng) * (MIN_LEN_CUT_FACTOR - MAX_LEN_CUT_FACTOR)
    len_cut = cut_factor * bigger_cyl_size[3]
    menisc_radius = bigger_cyl_size[2]



    cylinder_bottom_curvature = rand(rng, 0.05:0.1) * bigger_cyl_size[3]
    cylinder_bottom_curvature_b = cylinder_bottom_curvature * rand(rng, 1.8:0.1:2.3)
    rel_size_bottom_curvature = rand(rng, 0.4:0.01:0.6)

    cylinder_top_curvature = rand(rng, 0.05:0.01:0.25) * bigger_cyl_size[3]
    angle = 0.0
    pipe_len = rand(rng, 0.8:0.01:0.85) * (bigger_cyl_size[3] - cylinder_bottom_curvature + cylinder_top_curvature)

    pipe_cross = rand(rng, 0.1:0.2) * bigger_cyl_size[1]
    pipe_cross_section = (pipe_cross, pipe_cross)
    pipe_density = rand(rng, 0.2:0.01:0.4)
    dispenser_len = rand(rng, 0.15:0.01:0.25) * bigger_cyl_size[3]
    dispenser_cross_section = (rand(rng, 0.1:0.01:0.3) * bigger_cyl_size[1], rand(rng, 0.1:0.01:0.3) * bigger_cyl_size[2])
    dispenser_density = rand(rng, 0.25:0.4)

    dual_phase_percentage = rand(rng, 0.45:0.75)
    if (two_fluids)
        dual_phase_percentage = 1.0
    end

    density_inside = rand(rng, 0.15:0.35)
    density_inside_b = rand(rng, 0.7:0.9) * density_inside
    if (rounded_bottom)
        density_inside = min(density_inside, 0.6)
        density_inside_b = min(density_inside_b,0.6)
    end

    double_bottom_curvature = true
    first_ball = rand(rng) > 0.5
    second_ball = rand(rng) > 0.5

    curvature = rand(rng, 0.03:0.01:0.05) * bigger_cyl_size[3]
    tube_ball_fraction = rand(rng, 0.6:0.01:0.9)
    rel_dist = rand(rng, 0.3:0.01:0.7)
    torus_height_factor = rand(rng, 0.1:0.01:0.15)

    # Use maximum tilt angle from debug calculations
    x_cut_angle = rand(rng, -MAX_TILT_ANGLE:0.001:MAX_TILT_ANGLE)
    y_cut_angle = rand(rng, -MAX_TILT_ANGLE:0.001:MAX_TILT_ANGLE)
    if (is_2d)
        x_cut_angle = 0.0
        y_cut_angle = 0.0
    end

    menisc_radius_mult = 1.5
    menisc_cut_total_height = rand(rng, 0.1:0.01:0.3) * bigger_cyl_size[3]
    menisc_cut_height = rand(rng, 0.05:0.01:0.1) * bigger_cyl_size[3]

    # Initialize objects as empty arrays
    ob4 = []
    ob5 = []
    ob4b = []
    ob5b = []
    ob_cut = []
    ob_menisc_cut = []
    ob5b_mask = []
    density_map = Dict{String,Float64}()
    name_type_list = []
    double_bottom_curvature = true

    r_xy_small_factor = rand(rng, 0.19:0.01:0.22)
    r_z_small_factor = rand(rng, 0.019:0.001:0.022)

    # If debugging is enabled, calculate and print cylinder volumes for min/max parameters
    if is_debug
        println("==== CYLINDER VOLUME DEBUG ====")

        # Calculate volume with minimum parameters from our constants
        min_height = MIN_HEIGHT
        min_radius = min_height * MIN_RADIUS_FACTOR
        min_wall = 0.02
        min_len_cut = min_height * MIN_LEN_CUT_FACTOR
        min_radius_inner = min_radius - min_wall

        # Simulate max tilt impact (MAX_TILT_ANGLE rad)
        min_drop = min_radius_inner * sqrt(2) * sin(MAX_TILT_ANGLE)
        min_cut_center = (min_height / 2) - (min_len_cut / 2)
        min_lowest_cut = min_cut_center - min_drop - (min_len_cut / 2)
        min_cylinder_top = min_height / 2
        min_distance_from_top = min_cylinder_top - min_lowest_cut
        min_height_below_cut = min_height - min_distance_from_top

        # Min volume cylinder
        min_calc_cyl = cylinder(
            0.0cm, 0.0cm,
            (-min_height / 2 + min_height_below_cut / 2)cm,
            min_radius_inner * cm, min_radius_inner * cm,
            min_height_below_cut * cm,
            0.0, 0.0, 0.0, 1.0f0
        )
        min_volume = IP.volume(min_calc_cyl)
        min_volume_scalar = Unitful.ustrip(cm^3, min_volume)
        println("Minimum parameters: height=$(min_height), radius=$(min_radius), wall=$(min_wall), len_cut=$(min_len_cut)")
        println("Min cylinder effective: radius=$(min_radius_inner)cm, height=$(min_height_below_cut)cm")
        println("Minimum cylinder volume: $(min_volume_scalar) cm³")

        # Calculate volume with maximum parameters from our constants
        max_height = MAX_HEIGHT
        max_radius = max_height * MAX_RADIUS_FACTOR
        max_wall = 0.045
        max_len_cut = max_height * MAX_LEN_CUT_FACTOR
        max_radius_inner = max_radius - max_wall

        # Simulate min tilt (0.0 rad - no tilt)
        max_drop = 0.0
        max_cut_center = (max_height / 2) - (max_len_cut / 2)
        max_lowest_cut = max_cut_center - max_drop - (max_len_cut / 2)
        max_cylinder_top = max_height / 2
        max_distance_from_top = max_cylinder_top - max_lowest_cut
        max_height_below_cut = max_height - max_distance_from_top

        # Max volume cylinder
        max_calc_cyl = cylinder(
            0.0cm, 0.0cm,
            (-max_height / 2 + max_height_below_cut / 2)cm,
            max_radius_inner * cm, max_radius_inner * cm,
            max_height_below_cut * cm,
            0.0, 0.0, 0.0, 1.0f0
        )
        max_volume = IP.volume(max_calc_cyl)
        max_volume_scalar = Unitful.ustrip(cm^3, max_volume)
        println("Maximum parameters: height=$(max_height), radius=$(max_radius), wall=$(max_wall), len_cut=$(max_len_cut)")
        println("Max cylinder effective: radius=$(max_radius_inner)cm, height=$(max_height_below_cut)cm")
        println("Maximum cylinder volume: $(max_volume_scalar) cm³")
        println("Expected volume range: $(min_volume_scalar) to $(max_volume_scalar) cm³")

        # Calculate current parameters' volume
        curr_radius = bigger_cyl_size[1]
        curr_radius_inner = curr_radius - cylinder_wall_thickness
        curr_drop = curr_radius_inner * sqrt(2) * sin(max(abs(x_cut_angle), abs(y_cut_angle)))
        curr_cut_center = (bigger_cyl_size[3] / 2) - (len_cut / 2)
        curr_lowest_cut = curr_cut_center - curr_drop - (len_cut / 2)
        curr_cylinder_top = bigger_cyl_size[3] / 2
        curr_distance_from_top = curr_cylinder_top - curr_lowest_cut
        curr_height_below_cut = bigger_cyl_size[3] - curr_distance_from_top

        curr_calc_cyl = cylinder(
            0.0cm, 0.0cm,
            (-bigger_cyl_size[3] / 2 + curr_height_below_cut / 2)cm,
            curr_radius_inner * cm, curr_radius_inner * cm,
            curr_height_below_cut * cm,
            0.0, 0.0, 0.0, 1.0f0
        )
        curr_volume = IP.volume(curr_calc_cyl)
        curr_volume_scalar = Unitful.ustrip(cm^3, curr_volume)
        println("Current parameters: height=$(bigger_cyl_size[3]), radius=$(curr_radius), wall=$(cylinder_wall_thickness), len_cut=$(len_cut)")
        println("Current cylinder effective: radius=$(curr_radius_inner)cm, height=$(curr_height_below_cut)cm")
        println("Current cylinder volume: $(curr_volume_scalar) cm³")
        println("================================")
    end

    # Create the geometric object
    curr_time = Dates.now()
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

    return (vol, [("density_inside", density_inside),
            ("center_cylinder", center_cylinder),
            ("bigger_cyl_size", bigger_cyl_size),
            ("cylinder_wall_thickness", cylinder_wall_thickness),
            ("cylinder_bottom_curvature", cylinder_bottom_curvature),
            ("cylinder_top_curvature", cylinder_top_curvature),
            ("angle", angle),
            ("pipe_len", pipe_len),
            ("pipe_cross_section", pipe_cross_section),
            ("pipe_density", pipe_density),
            ("dispenser_len", dispenser_len),
            ("dispenser_cross_section", dispenser_cross_section),
            ("dispenser_density", dispenser_density / 2),
            ("len_cut", len_cut),
            ("dual_phase_percentage", dual_phase_percentage),
            ("x_cut_angle", x_cut_angle),
            ("y_cut_angle", y_cut_angle),
            ("density_inside_b", density_inside_b),
            ("cylinder_bottom_curvature_b", cylinder_bottom_curvature_b),
            ("rel_size_bottom_curvature", rel_size_bottom_curvature),
            ("double_bottom_curvature", double_bottom_curvature),
            ("menisc_radius_mult", menisc_radius_mult),
            ("menisc_cut_height", menisc_cut_height),
            ("menisc_cut_total_height", menisc_cut_total_height),
            ("rounded_bottom", rounded_bottom),
            ("curvature", curvature),
            ("tube_ball_fraction", tube_ball_fraction),
            ("rel_dist", rel_dist),
            ("add_pipe", add_pipe),
            ("torus_height_factor", torus_height_factor)],
        density_map, name_type_list, double_bottom_curvature, rounded_bottom, first_ball, second_ball,add_pipe)
end

"""
    get_cylinder_bool_mask(ob_el) -> BitArray{3}

Convert a cylinder object to a Boolean mask.

# Arguments
- `ob_el`: Cylinder object or array of objects

# Returns
A 3D boolean array representing the cylinder volume
"""
function get_cylinder_bool_mask(ob_el)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end

"""
    get_half_s_bool(ob_el) -> BitArray{3}

Convert a half sphere object to a Boolean mask.

# Arguments
- `ob_el`: Half sphere object or array of objects

# Returns
A 3D boolean array representing the half sphere volume
"""
function get_half_s_bool(ob_el)
    return (phantom(axes(ig)..., ob_el) .!= 0)

end

"""
    get_temp_folder() -> String

Get the path to the temporary directory for saving phantom data.

# Returns
Path to the temporary directory
"""
function get_temp_folder()
    temp_dir = mktempdir(; cleanup=false)
    return temp_dir
end

"""
save arguments to json file
"""
function save_args(args, vol, main_folder)
    json_path = "$(main_folder)/argss.json"
    # Convert the list of tuples to a dictionary
    args_dict = Dict(arg[1] => arg[2] for arg in args)
    args_dict["vol"] = vol
    # Save the dictionary to a JSON file
    open(json_path, "w") do io
        JSON.print(io, args_dict)
    end
end


"""
    create_boolean_density_dict(name_type_list, density_map) -> Dict{String, Tuple{Array{Bool}, Array{Float32}}}

Construct a dictionary where each key is an object name, and each value is a
3‐element tuple of:
1. A Boolean mask (true/false coverage).
2. A Float32 mask scaled by the corresponding density.
3. The original object.

# Arguments
- `name_type_list`: A list of tuples of the form (object_name, object, object_type).
- `density_map`: A dictionary mapping object names to their densities.

# Returns
A dictionary where each key is the object name and each value is a tuple of
(boolean_mask, float_mask × density, original_object).
"""
function create_boolean_density_dict(name_type_list, density_map)
    # Now store a third entry for the original object
    result_dict = Dict{String,Tuple{Array{Bool},Array{Float32},Any}}()

    for (obj_name, obj_data, obj_type) in name_type_list
        # Skip if obj_data is empty
        emptyy = false
        if isa(obj_data, Array)
            if isempty(obj_data)
                emptyy = true
            end
        end

        if (!emptyy)
            # Retrieve density
            if haskey(density_map, obj_name)
                dens = density_map[obj_name]
            else
                dens = 1.0
            end
            # print("\n name $(obj_name) \n")
            # Generate the Boolean mask depending on object type
            bool_mask = if obj_type == "cylinder"
                get_cylinder_bool_mask([obj_data])
            elseif obj_type == "array"
                obj_data
            else
                print("\n name $(obj_name) \n")
                get_half_s_bool([obj_data])
            end

            # Create a Float32 mask scaled by density
            float_mask = Float32.(bool_mask) .* Float32(dens)

            # Store tuple in the result dictionary
            result_dict[obj_name] = (bool_mask, float_mask, obj_data)
        end
    end
    print("\n ffffffffffffff finished dict \n")

    return result_dict
end

"""
    add_noise_at(phantoms_dict, key, density_map) -> Array{Float32,3}

Add random noise to a specific object in the phantom.

# Arguments
- `phantoms_dict`: Dictionary containing phantom boolean masks
- `key`: Object key to add noise to
- `density_map`: Dictionary mapping object names to density values

# Returns
A float array with added noise at the specified object location
"""
function add_noise_at(phantoms_dict, key, density_map)
    noise_mean = density_map[key]
    noise_std = 0.02
    noise_tensor = noise_mean .+ noise_std .* randn(size(phantoms_dict[key][1]))
    noise_tensor[(.!phantoms_dict[key][1])] .= 0.0
    return noise_tensor
end



function get_random_can_uploaded(is_2d, seed, uuid=nothing)
    # Use provided UUID or generate one if not provided
    uuid = isnothing(uuid) ? UUIDs.uuid4() : uuid

    main_folder = "$(get_temp_folder())/$(uuid)"
    # Create the main folder if it does not exist
    if !isdir(main_folder)
        mkpath(main_folder)
    end

    # Add the output folder to params so volume calculation can use it
    params = Dict{String,Any}()
    params["output_folder"] = main_folder

    # try
    #to monitor how long it take
    curr_time = Dates.now()
    #unique name
    #get random can

    vol, args, density_map, name_type_list, double_bottom_curvature, rounded_bottom, first_ball, second_ball,add_pipe = random_can(ig, (dims[1] * spacing[1] * 0.9, dims[2] * spacing[2] * 0.9, dims[3] * spacing[3] * 0.9), is_2d, seed, spacing, true)
    if args_json_path != " "

        print(" \n loading json\n")

        vol, args, density_map, name_type_list, double_bottom_curvature, rounded_bottom, first_ball, second_ball, add_pipe = json_based_can(ig, args_json_path, is_2d, true)
    else
        vol, args, density_map, name_type_list, double_bottom_curvature, rounded_bottom, first_ball, second_ball, add_pipe = random_can(ig, (dims[1] * spacing[1] * 0.9, dims[2] * spacing[2] * 0.9, dims[3] * spacing[3] * 0.9), is_2d, seed, spacing, true)
    end


    save_args(args, vol, main_folder)
    #get dictionary of phantoms
    phantoms_dict = create_boolean_density_dict(name_type_list, density_map)
    result_tensor = zeros(dims)
    if (!rounded_bottom)
        result_tensor = sum([phantoms_dict["ob4"][2], phantoms_dict["ob5"][2], phantoms_dict["ob4b"][2], phantoms_dict["ob5b"][2]
            # ,phantoms_dict["ob2_a"][2]
            # ,phantoms_dict["ob2_b"][2]
            , phantoms_dict["ob3"][2], phantoms_dict["pipe"][2], phantoms_dict["plastic_dispenser"][2]])
    else
        result_tensor = sum([phantoms_dict["ob4b"][2], phantoms_dict["ob5b"][2]
            # ,phantoms_dict["ob2_a"][2]
            # ,phantoms_dict["ob2_b"][2]
            , phantoms_dict["ob3"][2], phantoms_dict["pipe"][2], phantoms_dict["plastic_dispenser"][2]])
    end

    # ob2a_orig=copy(phantoms_dict["ob2_a"][1])

    ob_menisc_cut = phantoms_dict["ob_menisc_cut"][1] .& (phantoms_dict["ob_cut"][1].|phantoms_dict["ob_cut_orig"][1])
    # ob_menisc_cut=Float32.(ob_menisc_cut)
    # # Perform Gaussian smoothing on ob_menisc_cut
    # σ = (1.0, 1.0, 1.0)  # adjust as needed
    # gaussian_kernel = Kernel.gaussian(σ)
    # gaussian_kernel = Kernel.gaussian(σ)
    # ob_menisc_cut = imfilter(ob_menisc_cut, gaussian_kernel)
    # ob_menisc_cut=(ob_menisc_cut .> 0.02)

    # ob_menisc_cut = ob_menisc_cut.& phantoms_dict["ob3"][1]




    ## now we add the tilt to fluids of ob2_a and ob2_b
    # we add to ob2_a part that is not present in ob2_a and oblique_middle_cut but is in straight
    addded_ob2a = (.!(phantoms_dict["ob2_a"][1] .| phantoms_dict["oblique_middle_cut"][1])) .& phantoms_dict["straight_middle_cut"][1]
    #we add to ob2_b part that is not present in ob2_b but is in oblique and ob2_a
    addded_ob2b = ((.!(phantoms_dict["ob2_b"][1])) .& (phantoms_dict["oblique_middle_cut"][1])) .& phantoms_dict["straight_middle_cut"][1]
    #we need now to update entries in the dictionary so we add to the ob_2a addded_ob2a and subtract addded_ob2a from ob2_b
    #we need yo update a pair in dict fist is simple boolean mask second is multiplied by selected density
    new_ob2a_bool = ((phantoms_dict["ob2_a"][1]) .| addded_ob2a) .& (.!(addded_ob2b))
    #we do the same for ob2_b
    new_ob2b_bool = ((phantoms_dict["ob2_b"][1]) .| addded_ob2b) .& (.!(addded_ob2a))

    if (rounded_bottom)
        result_tensor[phantoms_dict["outer_sphere"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["outer_sphere"][2]
        result_tensor[phantoms_dict["inner_sphere"][1]] .= 0.0

        result_tensor[phantoms_dict["main_torus"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["main_torus"][2]

        result_tensor[phantoms_dict["inner_torus"][1]] .= 0.0

        new_ob2a_bool = new_ob2a_bool .| ((phantoms_dict["inner_torus"][1]))
        new_ob2a_bool = new_ob2a_bool .& (.!(phantoms_dict["outer_sphere"][1]))

        new_ob2b_bool = new_ob2b_bool .& (.!(phantoms_dict["outer_sphere"][1]))
    end


    new_ob2b_float32 = Float32.(new_ob2b_bool) .* Float32(density_map["ob2_b"])
    new_ob2a_float32 = Float32.(new_ob2a_bool) .* Float32(density_map["ob2_a"])


    #we update the dictionary
    local ob2a_orig = phantoms_dict["ob2_a"]
    local ob2b_orig = phantoms_dict["ob2_b"]
    phantoms_dict["ob2_a"] = (new_ob2a_bool, new_ob2a_float32, ob2a_orig[3])
    phantoms_dict["ob2_b"] = (new_ob2b_bool, new_ob2b_float32, ob2b_orig[3])

    #the cut above is a semicircle as we want to simulate miniscus
    ob_menisc_cut = ob_menisc_cut .& (phantoms_dict["ob2_a"][1] .| phantoms_dict["ob2_b"][1])


    result_tensor[phantoms_dict["ob2_a"][1]] .= 0.0
    result_tensor[phantoms_dict["ob2_b"][1]] .= 0.0
    result_tensor = result_tensor + add_noise_at(phantoms_dict, "ob2_a", density_map)
    result_tensor = result_tensor + add_noise_at(phantoms_dict, "ob2_b", density_map)

    # result_tensor= result_tensor+phantoms_dict["ob2_b"][2]
    # result_tensor= result_tensor+phantoms_dict["ob2_a"][2]


    result_tensor[ob_menisc_cut] .= 0.0

    result_tensor[phantoms_dict["ob4b"][1]] .= 0.0
    if (!rounded_bottom)
        if (double_bottom_curvature)
            result_tensor[phantoms_dict["center_bottom_out"][1]] .= 0.0
            result_tensor[phantoms_dict["center_bottom_in"][1]] .= 0.0
            result_tensor = result_tensor + phantoms_dict["center_bottom_out"][2]
            result_tensor = result_tensor + phantoms_dict["center_bottom_in"][2]
        end

        result_tensor[phantoms_dict["ob4"][1]] .= 0.0
        result_tensor[phantoms_dict["center_bottom_in"][1]] .= 0.0
    end



    result_tensor[phantoms_dict["ob_cut_orig"][1]] .= 0.0
    # result_tensor[(phantoms_dict["pipe"][1].&(.!phantoms_dict["pipe_in"][1]))] .= 0.0

    # result_tensor=result_tensor+Float32.(ob_menisc_cut.&phantoms_dict["ob_cut"][1])

    if (!rounded_bottom)
        result_tensor[phantoms_dict["ob5"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["ob5"][2]
        result_tensor[phantoms_dict["ob4"][1]] .= 0.0
        result_tensor[phantoms_dict["center_bottom_in"][1]] .= 0.0
    else
        result_tensor[phantoms_dict["outer_sphere"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["outer_sphere"][2]
        result_tensor[phantoms_dict["inner_sphere"][1]] .= 0.0

        # result_tensor[phantoms_dict["new_torus"][1]].=0.0
        # result_tensor=result_tensor+phantoms_dict["new_torus"][2]

        # result_tensor[phantoms_dict["inner_torus"][1]].=0.0

    end

    if (first_ball)
        result_tensor[phantoms_dict["ball1"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["ball1"][2]
    end

    if (second_ball)
        result_tensor[phantoms_dict["ball2"][1]] .= 0.0
        result_tensor = result_tensor + phantoms_dict["ball2"][2]
    end

    result_tensor[phantoms_dict["ob5b"][1]] .= 0.0
    result_tensor = result_tensor + phantoms_dict["ob5b"][2]
    result_tensor[phantoms_dict["ob4b"][1]] .= 0.0



    result_tensor[phantoms_dict["plastic_dispenser"][1]] .= 0.0
    result_tensor = result_tensor + phantoms_dict["plastic_dispenser"][2]

    if(rounded_bottom)
        result_tensor[phantoms_dict["inner_torus"][1]].=0.0
        result_tensor = result_tensor + phantoms_dict["inner_torus"][1].*Float32(density_map["ob2_a"])
    end


    if(add_pipe)
        result_tensor[phantoms_dict["pipe"][1].&(.!phantoms_dict["pipe_in"][1])] .= 0.0

        pipee = phantoms_dict["pipe"][2]
        pipee[phantoms_dict["pipe_in"][1]] .= 0.0
        result_tensor = result_tensor + pipee#+phantoms_dict["pipe_in"][2]
    end


    result_tensor = result_tensor ./ (maximum(result_tensor))




    # Rotate the image around the first axis by 45 degrees
    if (is_2d)

        for angle in [0, 45, 90, 135, 180, 225, 270, 315]
            rotated_image = rotation3d(sitk.GetImageFromArray(result_tensor), 1, angle)
            im_arr = sitk.GetArrayFromImage(rotated_image)

            # Create a linearly spaced vector for weighting along the first dimension
            incrr = range(0.5, stop=1.0, length=size(im_arr, 1))
            # Reshape to a (n,1,1) array so it broadcasts correctly
            incrr_reshaped = reshape(collect(incrr), size(im_arr,1), 1, 1)

            # Sum along the first axis after multiplying by the weights
            weighted_sum = sum(im_arr .* incrr_reshaped, dims=1)
            # Remove the singleton first dimension to yield a 2D result
            im_arr_weighted = dropdims(weighted_sum, dims=1)

            # Apply Gaussian smoothing to a 2D image (img_slice)
            gaussian_kernel = Kernel.gaussian((1.0, 1.0))
            im_arr_filtered = imfilter(im_arr_weighted, gaussian_kernel)
            im_arr_filtered=im_arr_filtered.-minimum(im_arr_filtered)
            im_arr_filtered=im_arr_filtered./maximum(im_arr_filtered)
            im_arr_filtered=(im_arr_filtered.*(-1)).*1000
            im_arr_filtered=im_arr_filtered.+minimum(im_arr_filtered)
            im_arr_filtered=Float32.(im_arr_filtered)

            # immm = sitk.GetImageFromArray(im_arr_filtered)
            # immm.SetSpacing((spacing[1] * 10, spacing[2] * 10, spacing[3] * 10))
            # immm.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
            # immm.SetOrigin((0.0, 0.0, 0.0))


            immm = sitk.GetImageFromArray(im_arr_filtered)
            immm.SetSpacing((spacing[1] * 10, spacing[2] * 10, spacing[3] * 10))


            sitk.WriteImage(immm, "$(main_folder)/rotated_example_can_$angle.nii.gz")
            dic_folder = "$(main_folder)/rotated_example_can_$angle"
            if !isdir(dic_folder)
                mkpath(dic_folder)
            end
            save_2d_sitk_image_as_dicom(immm, dic_folder)


        end

    end

    #saving raw files
    reference_nifti_path = "$(main_folder)/example_can_raw.nii.gz"
    save_nifti_with_meta(result_tensor, false, spacing, reference_nifti_path)
    save_sitk_image_as_dicom(sitk.ReadImage(reference_nifti_path), "$(main_folder)/example_can_raw")

    if(add_smooth)
        # Apply Gaussian smoothing to result_tensor
        σ = (1.0, 1.0, 1.0)  # adjust as needed
        gaussian_kernel = Kernel.gaussian(σ)
        result_tensor = imfilter(result_tensor, gaussian_kernel)
        result_tensor = imfilter(result_tensor, gaussian_kernel)
    end
    # Multiply result_tensor by the noise factor
    result_tensor = result_tensor .* (rand(size(result_tensor)))
    result_tensor = result_tensor .+ (rand(size(result_tensor)))

    result_tensor = result_tensor .+ (additive_noise .* randn(size(result_tensor)))

    reference_nifti_path = "$(main_folder)/example_can.nii.gz"
    save_nifti_with_meta(result_tensor, false, spacing, reference_nifti_path)

    save_sitk_image_as_dicom(sitk.ReadImage(reference_nifti_path), "$(main_folder)/example_can")



    if (add_radon)
        # Execute the Python script using the run function
        input_path = "$(main_folder)/example_can.nii.gz"
        output_path = "$(main_folder)/after_radon.nii.gz"
        output_path_b = "$(main_folder)/after_radon_plus_before.nii.gz"
        script_path = joinpath(@__DIR__, "get_approximate_radon_inverse.py")

        command = `python3 $script_path $input_path $output_path $output_path_b`

        run(command)

        save_sitk_image_as_dicom(sitk.ReadImage(output_path), "$(main_folder)/after_radon")
        save_sitk_image_as_dicom(sitk.ReadImage(output_path_b), "$(main_folder)/after_radon_plus_before")
    end


    fluid_mask = phantoms_dict["ob2_a"][1] .| phantoms_dict["ob2_b"][1]
    if first_ball
        fluid_mask .&= .!phantoms_dict["ball1"][1]
    end
    if second_ball
        fluid_mask .&= .!phantoms_dict["ball2"][1]
    end
    if(add_pipe)
        fluid_mask = fluid_mask .& (.!phantoms_dict["pipe"][1])

    end
    fluid_mask = fluid_mask .& (.!phantoms_dict["ob_cut_orig"][1])
    fluid_mask = fluid_mask .& (.!ob_menisc_cut)
    if (!rounded_bottom)
        fluid_mask = fluid_mask .& (.!phantoms_dict["center_bottom_out"][1])
        fluid_mask = fluid_mask .& (.!phantoms_dict["ob4"][1])
        fluid_mask = fluid_mask .& (.!phantoms_dict["ob4b"][1])
        fluid_mask = fluid_mask .& (.!phantoms_dict["ob5"][1])
    end
    #add inside pipe
    if(add_pipe)
        pipe_inn_=(phantoms_dict["calc_cyl"][1].& phantoms_dict["pipe_in"][1])
        fluid_mask = fluid_mask .| pipe_inn_
        voxel_volume = spacing[1] * spacing[2] * spacing[3]
        print("\n nnnnnnnnnnn pipe_inn_ numerical $(sum(pipe_inn_)*voxel_volume) \n")


    end

    fluid_mask_uint8 = UInt8.(fluid_mask)
    fluid_mask_img = sitk.GetImageFromArray(fluid_mask_uint8)
    save_nifti_with_meta(UInt8.(fluid_mask), false, spacing, joinpath(main_folder, "fluid_mask_debug.nii.gz"))

    # Save fluid mask as NIfTI
    fluid_mask_uint8 = UInt8.(fluid_mask)
    fluid_mask_img = sitk.GetImageFromArray(fluid_mask_uint8)
    fluid_mask_img.SetSpacing((spacing[1], spacing[2], spacing[3]))
    fluid_mask_img.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    fluid_mask_img.SetOrigin((0.0, 0.0, 0.0))

    fluid_mask_path = joinpath(main_folder, "fluid_mask.nii.gz")
    save_nifti_with_meta(UInt8.(fluid_mask), false, spacing, fluid_mask_path)

    # Reference DICOM - use the directory containing the DICOM slices
    reference_dicom_path = joinpath(main_folder, "example_can")

    # Make sure we have created the DICOM files before conversion
    if !isdir(reference_dicom_path)
        # If we haven't saved the example_can as DICOM yet, do it now
        if isfile(joinpath(main_folder, "example_can.nii.gz"))
            save_sitk_image_as_dicom(sitk.ReadImage(joinpath(main_folder, "example_can.nii.gz")), reference_dicom_path)
        else
            println("Warning: example_can.nii.gz not found, cannot create reference DICOM")
        end
    end

    # Convert fluid mask to DICOM-SEG
    seg_output_folder = joinpath(main_folder, "fluid_mask")
    convert_nifti_to_dicom_seg(fluid_mask_path, reference_dicom_path, seg_output_folder, reference_nifti_path)

    # Save additional masks as NIfTI files
    # 1. Ball masks
    if first_ball && haskey(phantoms_dict, "ball1")
        save_mask_as_nifti(
            phantoms_dict["ball1"][1],
            joinpath(main_folder, "ball1_mask.nii.gz"),
            spacing
        )
        convert_nifti_to_dicom_seg(joinpath(main_folder, "ball1_mask.nii.gz"), reference_dicom_path, joinpath(main_folder, "ball1_mask"))

    end

    if second_ball && haskey(phantoms_dict, "ball2")
        save_mask_as_nifti(
            phantoms_dict["ball2"][1],
            joinpath(main_folder, "ball2_mask.nii.gz"),
            spacing
        )
        convert_nifti_to_dicom_seg(joinpath(main_folder, "ball2_mask.nii.gz"), reference_dicom_path, joinpath(main_folder, "ball2_mask"))

    end
    save_mask_as_nifti(
        Array(ob_menisc_cut),
        joinpath(main_folder, "ob_menisc_cut.nii.gz"),
        spacing
    )


    # 2. Pipe mask
    if(add_pipe)
        if haskey(phantoms_dict, "pipe")
            println("Type of pipe mask: ", typeof(phantoms_dict["pipe"][1]))#Type of pipe mask: Array{Bool, 3}
            save_mask_as_nifti(
                phantoms_dict["pipe"][1],
                joinpath(main_folder, "pipe_mask.nii.gz"),
                spacing
            )
            convert_nifti_to_dicom_seg(joinpath(main_folder, "pipe_mask.nii.gz"), reference_dicom_path, joinpath(main_folder, "pipe_mask"))

        end
    end

    # 3. Plastic dispenser mask
    if haskey(phantoms_dict, "plastic_dispenser")

        # println("Size of arr: ", size(arr))
        save_mask_as_nifti(phantoms_dict["plastic_dispenser"][1]
            ,
            joinpath(main_folder, "dispenser_mask.nii.gz"),
            spacing
        )
        convert_nifti_to_dicom_seg(joinpath(main_folder, "dispenser_mask.nii.gz"), reference_dicom_path, joinpath(main_folder, "dispenser_mask"))

    end

    save_mask_as_nifti(
        phantoms_dict["ob_cut"][1],
        joinpath(main_folder, "ob_cut.nii.gz"),
        spacing
    )
    save_mask_as_nifti(
        phantoms_dict["ob_cut"][1],
        joinpath(main_folder, "ob_cut.nii.gz"),
        spacing
    )



    save_mask_as_nifti(
        phantoms_dict["ob_cut_orig"][1],
        joinpath(main_folder, "ob_cut_orig.nii.gz"),
        spacing
    )
    save_mask_as_nifti(
        phantoms_dict["calc_cyl"][1],
        joinpath(main_folder, "calc_cyl.nii.gz"),
        spacing
    )
    save_mask_as_nifti(
        phantoms_dict["inside_of_pipe_in_fluid"][1],
        joinpath(main_folder, "inside_of_pipe_in_fluid.nii.gz"),
        spacing
    )







    # Store original state for meniscus volume calculation
    original_ob2_a = copy(phantoms_dict["ob2_a"][1])
    original_ob2_b = copy(phantoms_dict["ob2_b"][1])
    original_ob_menisc_cut = copy(phantoms_dict["ob_menisc_cut"][1])
    original_ob_cut = copy(phantoms_dict["ob_cut"][1])

    # Extract parameters into a dict for easier function call and ensure all needed values are included
    params = Dict{String,Any}()
    for (key, value) in args
        params[key] = value
    end

    # Add essential parameters to ensure they're available
    params["spacing"] = spacing
    params["output_folder"] = main_folder
    params["fluid_mask"] = fluid_mask

    # Calculate meniscus volume before modifying the object arrays
    meniscus_mask = (original_ob_menisc_cut .& original_ob_cut) .&
                    (phantoms_dict["ob2_a"][1] .| phantoms_dict["ob2_b"][1])
    meniscus_voxel_count = sum(meniscus_mask)
    meniscus_volume_numerical = meniscus_voxel_count * (spacing[1] * spacing[2] * spacing[3])
    params["meniscus_volume_numerical"] = meniscus_volume_numerical
    println("Meniscus volume (numerical): $(meniscus_volume_numerical) cm³")

    # Initialize required values with defaults
    params["inner_torus_volume_analytical"] = 0.0
    params["outer_sphere_volume"] = 0.0

    # Extract object volumes based on geometry type
    if rounded_bottom && haskey(phantoms_dict, "outer_sphere") && haskey(phantoms_dict, "inner_torus")
        # For rounded bottom, calculate outer sphere and inner torus volumes
        outer_sphere_volume = sum(phantoms_dict["outer_sphere"][1]) * (spacing[1] * spacing[2] * spacing[3])
        inner_torus_volume = sum(phantoms_dict["inner_torus"][1]) * (spacing[1] * spacing[2] * spacing[3])

        params["outer_sphere_volume"] = outer_sphere_volume
        params["inner_torus_volume_analytical"] = inner_torus_volume

        println("Outer sphere volume: $(outer_sphere_volume) cm³")
        println("Inner torus volume: $(inner_torus_volume) cm³")
    end

    # Calculate fluid volume using the fixed function
    numerical_vol, analytical_vol = compute_accurate_fluid_volume_fixed(
        phantoms_dict,
        spacing,
        first_ball,
        second_ball,
        rounded_bottom,
        params,
        ig
    )

    # Add volumes to arguments for saving
    push!(args, ("fluid_volume_numerical", numerical_vol))
    fluid_volume_numerical_my = sum(fluid_mask) * (spacing[1] * spacing[2] * spacing[3])
    # push!(args, ("fluid_volume_numerical_my", fluid_volume_numerical_my))
    push!(args, ("fluid_volume_analytical", analytical_vol))
    push!(args, ("spacing", spacing))
    push!(args, ("uuid", uuid))
    push!(args, ("seed", seed))

    save_args(args, vol, main_folder)

    zip_path = "$(get_temp_folder())/$(uuid).zip"
    run(`zip -r $zip_path $main_folder`)
    # Ensure the Google Cloud Storage folder exists
    file_name = "can_is_2D_$(is_2d)_is_radon_$(add_radon)_$(dims[1])|$(dims[2])|$(dims[3])_$uuid"

    if !haskey(ENV, "SKIP_UPLOAD")
        command = `gcloud storage cp $zip_path gs://metro_tk_kplayground/cansx128/$file_name`
        # Execute the command
        run(command)
        rm(main_folder; force=true, recursive=true)
        rm(zip_path; force=true)
    else
        println("Skipping upload and deletion for debugging/testing.")
        println("Output stored in: $main_folder")
        println("Zip stored in: $zip_path")
    end

    return numerical_vol, analytical_vol, fluid_volume_numerical_my
end

"""
    calculate_and_save_fluid_volume(seed, is_2d=false, uuid=nothing) -> Tuple{Float64, Float64, Float64}

Calculate and save the fluid volume for a randomly generated can phantom.
Also saves various masks as NIfTI files and converts the fluid mask to DICOM-SEG.

# Arguments
- `seed`: Random seed for reproducibility
- `is_2d`: Boolean indicating if 2D mode should be used (default: false)
- `uuid`: Optional UUID to use for folder naming (default: nothing)

# Returns
A tuple containing numerical and analytical volumes in cubic centimeters
"""
function calculate_and_save_fluid_volume(seed, is_2d=false, uuid=nothing)
    numerical_vol, analytical_vol, fluid_volume_numerical_my = get_random_can_uploaded(is_2d, seed, uuid)
    return numerical_vol, analytical_vol, fluid_volume_numerical_my
end

# Main loop with UUID-based seeding
    # Generate a UUID and convert it to an integer for the seed
# while(true)

# uuid = UUIDs.uuid4()
seed = abs(hash(uuid))


    print("\n UUID: $uuid")
    print("\n seed: $seed \n")

    numerical_vol, analytical_vol, fluid_volume_numerical_my = calculate_and_save_fluid_volume(seed, is_2d, uuid)

    println("Main execution completed with calculated fluid volumes:")
    println("Numerical: $(numerical_vol) cm³")
    println("Analytical: $(analytical_vol) cm³")
    println("fluid_volume_numerical_my: $(fluid_volume_numerical_my) cm³")
# end


# julia in_docker_organized/main_create_phantom_can.jl 128x128x128 true true ddac5ff5-1cc4-4447-9f44-96c372134657 true true 0.1
# /workspaces/synthethic_tomo/.devcontainer# julia in_docker_organized/main_create_phantom_can.jl 256x256x256 true true ddac5ff5-1cc4-4447-9c44-96c372134657 true '/workspaces/synthethic_tomo/data/args.json'
#julia in_docker_organized/main_create_phantom_can.jl 256x256x256 false true b230290a-ec7f-413a-921f-b62dead16418 true true 0.1

# /workspaces/devcontainer/.devcontainer/in_docker_organized/main_create_phantom_can.jl