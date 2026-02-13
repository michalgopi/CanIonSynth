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

# __precompile__(false)  # Add this at the very top
using Pkg;
Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git");

using Pkg
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

isinteractive() ? jim(:prompt, true) : prompt(:draw)
includet.([
    "get_geometry_main.jl",
    "geometry_utils.jl"
])

# ---------------------------------------------------------------------
# Helper functions and docstrings


#=
This program creates an image composed of several cylinders that share the same long axis.
A "base" cylinder is defined with fixed center (c_x, c_y). All other cylinders are specified
by a user-provided dictionary (with keys for density, radius, length, and beginning offset along the z‐axis).
Each additional cylinder’s center z is computed relative to the base bottom, and its density is adjusted
to compensate for overlapping cylinders. When the z‐intervals (projections onto the long axis) of two cylinders
overlap, their densities add up. To display the desired density for each additional cylinder, the effective density
passed to the cylinder constructor is computed as:
  actual_density = desired_density − sum(overlapping densities from previously defined cylinders)
=#

# Function to compute the z-interval of a cylinder
"""
    cylinder_interval(c_z::Float64, l_z::Float64) -> (Float64, Float64)

Compute the z-interval for a cylinder given its center `c_z` and length `l_z`.

# Returns
A tuple `(z_start, z_end)`, where:
- `z_start = c_z - (l_z/2)`
- `z_end = c_z + (l_z/2)`
"""
function cylinder_interval(c_z::Float64, l_z::Float64)
    return (c_z - l_z / 2, c_z + l_z / 2)
end

# Function to check if two cylinders overlap along the z-axis
"""
    overlaps(c1::Float64, l1::Float64, c2::Float64, l2::Float64) -> Bool

Determine if two cylinders overlap along the z-axis by comparing their intervals.

# Arguments
- `c1`, `l1`: Center and length of the first cylinder.
- `c2`, `l2`: Center and length of the second cylinder.

# Returns
True if the intervals intersect, false otherwise.
"""
function overlaps(c1::Float64, l1::Float64, c2::Float64, l2::Float64)
    s1, e1 = cylinder_interval(c1, l1)
    s2, e2 = cylinder_interval(c2, l2)
    return max(s1, s2) < min(e1, e2)
end

# Function to compute adjusted density for a new cylinder
function compute_adjusted_density(new_desired::Float64, new_z::Float64, new_lz::Float64, existing::Vector{NamedTuple})
    total_overlap_density = 0.0
    for cyl in existing
        if overlaps(new_z, new_lz, cyl.c_z, cyl.l_z)
            total_overlap_density += cyl.d
        end
    end
    return new_desired - total_overlap_density
end

# Function to create a new cylinder given its specification dictionary
"""
    create_cylinder(base_x::Float64, base_y::Float64, base_bottom_z::Float64,
                    spec::Dict{String,Float64}) -> NamedTuple

Create a cylinder using the specification in `spec` and the base cylinder position.
`spec` must contain:
  - "d": desired density,
  - "radius": radius of the cylinder (l_x),
  - "l_z": length of the cylinder,
  - "offset": beginning offset relative to the base bottom (i.e., where z = 0 at base bottom).

The center z for the cylinder is computed as:
  c_z = base_bottom_z + offset + (l_z/2)

# Returns
A NamedTuple with the parameters used to call the cylinder constructor:
  (c_x, c_y, c_z, l_x, l_y, l_z, angle, d)
"""
function create_cylinder(base_x::Float64, base_y::Float64, base_bottom_z::Float64,
    spec, angle::Float64, name)
    # Extract specification values
    d_desired = spec["d"]
    radius = spec["radius"]
    l_z = spec["l_z"]
    offset = spec["offset"]
    # Determine is_halph_sphere and is_ellipsoid based on spec; default to false if not provided
    is_halph_sphere = haskey(spec, "is_halph_sphere") ? spec["is_halph_sphere"] : false
    is_ellipsoid = haskey(spec, "is_ellipsoid") ? spec["is_ellipsoid"] : false

    # Compute center z based on base bottom
    c_z = base_bottom_z + offset + (l_z / 2)
    if offset < 0
        c_z = 0.0
    end
    # Check if "theta" and "center_circle_radius" are present in spec and add them to the result
    if haskey(spec, "theta")
        if spec["theta"] > 0.0
            return (name=name, is_spiral=true, is_halph_sphere=is_halph_sphere, is_ellipsoid=is_ellipsoid, c_x=base_x, c_y=base_y, c_z=c_z,
                l_x=radius, l_y=radius, l_z=l_z, angle=angle, d=d_desired,
                theta=spec["theta"], center_circle_radius=spec["center_circle_radius"], orig_d=d_desired)
        end
    end

    # Return a NamedTuple representing this cylinder with original desired density
    return (name=name, is_spiral=false, is_halph_sphere=is_halph_sphere, is_ellipsoid=is_ellipsoid, c_x=base_x, c_y=base_y, c_z=c_z,
        l_x=radius, l_y=radius, l_z=l_z, angle=angle, d=d_desired)
end

"""
    adjust_densities!(cylinders::Vector{NamedTuple})

Adjust densities of all cylinders to account for overlaps.
Modifies the cylinders vector in-place.
"""
function adjust_densities!(cylinders::Vector{NamedTuple})
    n = length(cylinders)
    adjusted_cylinders = Vector{NamedTuple}(undef, n)
    adjusted_cylinders[1] = cylinders[1]
    for i in 2:n
        d_adjusted = compute_adjusted_density(cylinders[i].d, cylinders[i].c_z, cylinders[i].l_z, cylinders)
        adjusted_cylinders[i] = (name=cylinders[i].name, c_x=cylinders[i].c_x, c_y=cylinders[i].c_y, c_z=cylinders[i].c_z,
            l_x=cylinders[i].l_x, l_y=cylinders[i].l_y, l_z=cylinders[i].l_z,
            angle=cylinders[i].angle, d=d_adjusted, orig_d=cylinders[i].d)
    end
    return adjusted_cylinders
end

"""
    find_immediate_parent(child::NamedTuple, cylinders::Vector{NamedTuple}) -> Union{NamedTuple, Nothing}

Find the immediate parent cylinder that contains the given child cylinder.
A cylinder is considered a parent if:
1. It overlaps with the child on the z-axis
2. Its radius (l_x) is larger than the child's
Among valid parents, selects the one with the smallest radius difference.

Returns nothing if no valid parent is found.
"""
function find_immediate_parent(child::NamedTuple, cylinders::Vector{NamedTuple})
    smallest_diff = Inf
    immediate_parent = nothing

    for potential_parent in cylinders
        # Skip if same cylinder or smaller/equal radius
        if potential_parent.l_x <= child.l_x
            continue
        end

        # Check z-axis overlap
        if overlaps(child.c_z, child.l_z, potential_parent.c_z, potential_parent.l_z)
            radius_diff = potential_parent.l_x - child.l_x
            if radius_diff < smallest_diff
                smallest_diff = radius_diff
                immediate_parent = potential_parent
            end
        end
    end

    return immediate_parent
end

"""
    sort_cylinders_by_radius!(cylinders::Vector{NamedTuple})

Sort cylinders by radius in descending order, keeping base cylinder first.
This ensures we process larger cylinders before smaller ones that might be inside them.
"""
function sort_cylinders_by_radius!(cylinders::Vector{NamedTuple})
    base = cylinders[1]
    rest = sort(cylinders[2:end], by=x -> -x.l_x)  # Sort by descending radius
    return vcat([base], rest)
end

"""
    parent_index(child::NamedTuple, adjusted::Vector{NamedTuple}) -> Int

Given a child cylinder, find the index of its immediate parent in `adjusted`.
Returns -1 if no parent is found.
"""
function parent_index(child::NamedTuple, adjusted::Vector{NamedTuple})
    parent_indicies = []
    for (idx, cyl) in enumerate(adjusted)
        if cyl !== nothing && cyl.l_x > child.l_x && overlaps(child.c_z, child.l_z, cyl.c_z, cyl.l_z)
            push!(parent_indicies, (idx, abs((cyl.is_halph_sphere ? cyl.l_x * 0.5 : cyl.l_x) - (child.is_halph_sphere ? child.l_x * 0.5 : child.l_x))))
        end
    end

    if length(parent_indicies) > 0
        sorted = sort(parent_indicies, by=x -> x[2])
        return sorted[1][1]
    end

    return -1
end

"""
    compute_final_sum(child_d::Float64, parent_final::Float64) -> Float64

Compute the child's final sum in its region, ensuring the child ends up with density = child_d.
The 'partial' is child_d - parent_final,
so child_final = parent_final + partial = parent_final + (child_d - parent_final) = child_d.
"""
function compute_final_sum(child_d::Float64, parent_final::Float64)
    return child_d
end

"""
    compute_partial(child_d::Float64, parent_final::Float64) -> Float64

Compute the actual partial density for the child cylinder, ensuring final sum = child_d.
partial = child_d - parent_final
"""
function compute_partial(child_d::Float64, parent_final::Float64)
    return child_d - parent_final
end

"""
    compute_density_adjustment(child::NamedTuple, parent::NamedTuple, is_base::Bool) -> Float64

Compute how much density the child cylinder should add above the parent.
If the parent is the base, return (parent.d - child.d).
Otherwise, return (child.d - parent.d).
"""
function compute_density_adjustment(child::NamedTuple, parent::NamedTuple, is_base::Bool)
    if is_base
        return parent.d - child.d
    else
        return child.d - parent.d
    end
end

"""
    compute_effective_density(child::NamedTuple, parent_final::Float64) -> Float64

Compute the effective density for a child cylinder given its immediate parent’s displayed density.
This effective density is defined as:
  effective = parent_final - child.d

So that the final displayed density becomes (parent_final - effective) = child.d.
"""
function compute_effective_density(child::NamedTuple, parent_final)
    return child.d - parent_final
end

"""
    adjust_densities!(cylinders::Vector{NamedTuple}) -> Vector{NamedTuple}

Sort cylinders by descending radius (keeping the base first). For each cylinder with an immediate parent,
set its effective density to (parent's displayed density – desired density). If no parent is found, the
cylinder keeps its own density. In all cases, the final displayed density will equal the user-specified value.
"""
function adjust_densities!(cylinders::Vector{NamedTuple})::Vector{NamedTuple}
    sorted_cyls = sort_cylinders_by_radius!(cylinders)
    n = length(sorted_cyls)
    adjusted = Vector{NamedTuple}(undef, n)

    # Base cylinder: final displayed density equals its own dictionary value.
    adjusted[1] = sorted_cyls[1]

    for i in 2:n
        child = sorted_cyls[i]
        pinx = parent_index(child, sorted_cyls)
        if pinx == -1
            effective = child.d
        else
            p_final = sorted_cyls[pinx].d # parent's displayed density
            effective = compute_effective_density(child, p_final)
        end
        effective = effective + 0.000001
        adjusted[i] = (@set child.d = effective)
    end

    return adjusted
end

function get_beg_voxel(spacing, center_cylinder, cyl_size, axis)
    return ((center_cylinder[axis] - cyl_size[axis] / 2) / spacing[axis])
end

"""
    get_cyl_obj(new_cyl, ig, spacing)

If new_cyl contains "theta" (and optionally "center_circle_radius"), return the spiral_screw phantom.
Otherwise, create a standard cylinder phantom.
"""
function get_cyl_obj(new_cyl, ig, spacing)
    voll=-100000.0
    if new_cyl.is_halph_sphere
        ob_el = half_sphere_z((new_cyl.c_x)cm, (new_cyl.c_y)cm, (new_cyl.c_z - new_cyl.l_z / 2)cm,
            (new_cyl.l_x)cm, (new_cyl.l_x)cm, (new_cyl.l_z)cm,
            0, 0, 0, new_cyl.d)
        voll=IP.volume(ob_el)
        ob_el = (phantom(axes(ig)..., [ob_el]) .!= 0)
    elseif new_cyl.is_ellipsoid
        ob_el = ellipsoid((new_cyl.c_x)cm, (new_cyl.c_y)cm, (new_cyl.c_z)cm,
            (new_cyl.l_x)cm, (new_cyl.l_y)cm, (new_cyl.l_z)cm,
            0, 0, 0, new_cyl.d)
        voll=IP.volume(ob_el)
        ob_el = (phantom(axes(ig)..., [ob_el]) .!= 0)
    else
        ob_el = ImagePhantoms.cylinder((new_cyl.c_x)cm, (new_cyl.c_y)cm, (new_cyl.c_z)cm,
            (new_cyl.l_x)cm, (new_cyl.l_x)cm, (new_cyl.l_z)cm,
            0, 0, 0, new_cyl.d)

        voll=IP.volume(ob_el)
        ob_el = phantom(axes(ig)..., [ob_el])

        ob_el = (ob_el + reverse(ob_el, dims=1)) .!= 0

        if new_cyl.is_spiral
            bigger_cyl_size = (new_cyl.l_x, new_cyl.l_y, new_cyl.l_z)
            dims = ig.dims
            beginning_center = ((round(dims[1] / 2)), (round(dims[2] / 2)), 2.0)
            z_length = (new_cyl.l_z * 4) / spacing[3]
            theta = new_cyl.theta
            center_circle_radius = haskey(new_cyl, :center_circle_radius) ? new_cyl.center_circle_radius : 0.0
            spiral = spiral_screw(beginning_center, z_length, theta, bigger_cyl_size[1] / spacing[1] - center_circle_radius / spacing[1], center_circle_radius, ig.dims)
        end
    end

    res = Float32.(ob_el) * new_cyl.d
    return res,voll
end

# Main function to build the image from the base cylinder and additional cylinders dictionary
"""
    build_image(base_params::Dict{String, Float64}, specs::Dict{String, Dict{String, Float64}}, angle::Float64) -> Vector{NamedTuple}

Construct an image composed of cylinders. All cylinders share the same base (c_x, c_y).
`base_params` must include:
  - "c_x", "c_y", "c_z": center of the base cylinder,
  - "radius": radius for the base,
  - "l_z": length of the base,
  - "d": density of the base.

`specs` is a dictionary where each key is a cylinder name and each value is a dictionary with keys:
  - "d", "radius", "l_z", "offset" (offset from the bottom of the base cylinder).

The function returns a vector of NamedTuples, each representing a cylinder with final parameters.
"""
function build_image(base_params, specs, angle::Float64, spacing, dims)
    ig = ImageGeom(dims=(dims[1],dims[2],dims[3]), deltas=(spacing[1] * cm, spacing[2] * cm, spacing[3] * cm))

    # Base cylinder's parameters
    base_x = 0.0#(dims[1]/4)*spacing[1]
    base_y = 0.0#(dims[2]/4)*spacing[2]
    base_cz = base_params["c_z"]
    base_lz = base_params["l_z"]
    base_radius = base_params["radius"]
    base_d = base_params["d"]

    # Compute base bottom z
    base_bottom_z = base_cz - (base_lz / 2)

    # Create the base cylinder as a NamedTuple and add to the list
    cylinders = Vector{NamedTuple}(undef, 0)
    base_cyl = (name="base_cyl", is_halph_sphere=false, is_spiral=false, is_ellipsoid=false, c_x=base_x, c_y=base_y, c_z=base_cz, l_x=base_radius, l_y=base_radius,
        l_z=base_lz, angle=angle, d=base_d, orig_d=base_d)
    push!(cylinders, base_cyl)

    # Process each additional cylinder specification
    for (name, spec) in specs
        new_cyl = create_cylinder(base_x, base_y, base_bottom_z, spec, angle, name)
        push!(cylinders, new_cyl)
    end

    # Adjust densities after all cylinders are created
    cylinders = adjust_densities!(cylinders)

    cylinders_prim = deepcopy(cylinders)

    # Print cylinder information after density adjustment
    for (i, cyl) in enumerate(cylinders)
        name = i == 1 ? "base" : "cylinder_$(i-1)"
        println(" is_spiral : $(cyl.is_spiral)  ;$name: center = (", cyl.c_x, ", ", cyl.c_y, ", ", cyl.c_z,
            ") radius = ", cyl.l_x, " length = ", cyl.l_z, " density = ", cyl.d)
    end
    # Create a dictionary where the key is the cylinder name and the value is the cylinder density
    cylinder_densities = Dict{String,Float64}()
    for cyl in cylinders
        cylinder_densities[cyl.name] = cyl.d
    end
    res = map(new_cyl -> get_cyl_obj(new_cyl, ig, spacing), cylinders)

    # Create dictionary associating the cylinder "name" with its entry in "res"
    named_res = Dict{String,Any}()
    for (i, cyl) in enumerate(cylinders)
        named_res[cyl.name] = res[i][1]
    end

    vol_res = Dict{String,Any}()
    for (i, cyl) in enumerate(cylinders)
        vol_res[cyl.name] = res[i][2]
    end

    return cylinders_prim, named_res, cylinder_densities,vol_res
end


"""
    iterative_position_adjustment(base_image, rotated_image, axis, spacing)

Iteratively move rotated_image relative to base_image along specified axis until
finding the optimal position where they just touch without excessive overlap.

# Arguments
- `base_image`: Boolean array representing the base object
- `rotated_image`: Boolean array representing the object to be positioned
- `axis`: Axis along which to move (typically 3 for z-axis)
- `spacing`: Array of voxel spacing in each dimension

# Returns
- `distance_moved`: Total distance moved (steps × spacing[axis])
- `direction`: Final movement direction ("up" or "down")
- `adjusted_image`: The final positioned image
"""
function iterative_position_adjustment(base_image, rotated_image, axis, spacing)
    # Create a copy of the image to adjust
    adjusted_image = copy(rotated_image)
    steps_moved = 0
    max_iterations = 1000  # Safety limit

    # Check initial overlap
    has_overlap = any(base_image .& adjusted_image)

    if has_overlap
        # Case 1: Images overlap - move up until no overlap
        direction = "up"
        while has_overlap && steps_moved < max_iterations
            # Move one layer at a time
            adjusted_image = move_image(adjusted_image, axis, spacing[axis], direction, spacing)
            steps_moved += 1
            has_overlap = any(base_image .& adjusted_image)
        end

        # Move back one step to ensure minimal contact
        if steps_moved > 0
            adjusted_image = move_image(adjusted_image, axis, spacing[axis], "down", spacing)
            steps_moved -= 1
        end

        # Report net movement
        distance_moved = steps_moved * spacing[axis]
        final_direction = steps_moved > 0 ? "up" : "down"

    else
        # Case 2: No initial overlap - move down until overlap
        direction = "down"
        steps_down = 0
        while !has_overlap && steps_down < max_iterations
            # Move one layer at a time
            adjusted_image = move_image(adjusted_image, axis, spacing[axis], direction, spacing)
            steps_down += 1
            has_overlap = any(base_image .& adjusted_image)
        end

        # If we found an overlap, move back one step
        if has_overlap && steps_down > 0
            adjusted_image = move_image(adjusted_image, axis, spacing[axis], "up", spacing)
            steps_down -= 1
        end

        # Report net movement (negative because we're moving down)
        steps_moved = -steps_down
        distance_moved = abs(steps_moved) * spacing[axis]
        final_direction = "down"
    end

    # Clean up edge cases
    if steps_moved == 0
        final_direction = "none"  # No movement occurred
        distance_moved = 0.0
    end

    println("Position adjustment complete: moved $(distance_moved) cm $(final_direction)")
    return distance_moved, final_direction, adjusted_image
end


function reverse_and_add(array, axis)
    # arr=array.!=0.0
    # vall=mean(array[arr])
    return reverse(array, dims=axis) .| array
    # return arr.*vall

end



"""
    create_ionic_chamber_phantom(params)

Main function to create an ionic chamber phantom based on the provided parameters.
Takes a dictionary of parameters and produces a 3D phantom image.

# Chamber Type Selection Logic
The chamber type is determined by the following parameters:
- square_top: When true, creates a flat-topped chamber
- ball_like: When true, creates a spherical-topped chamber
- lolipop_like: When true, creates a lollipop-shaped chamber
- new_flat_sizes: When true, creates a standardized chamber with precise volume

If none of these are true, a default rounded-top chamber is created.

# Layer Structure
All ionic chamber types share a common layer structure:
- Central copper electrode: Innermost conductive component
- Inner insulator: Polyethylene layer surrounding the copper electrode
- Inner aluminum layer: First metallic shielding layer
- Outer insulator: Second polyethylene insulation layer
- Outer aluminum layer: Second metallic shielding layer
- Graphite housing: Outermost structural layer

The specific geometry of these layers varies by chamber type.
"""
function create_ionic_chamber_phantom(params)
    # Extract parameters from the dictionary
    graphite_density = params["graphite_density"]
    copper_density = params["copper_density"]
    aluminium_density = params["aluminium_density"]
    insulation_density = params["insulation_density"]
    square_top = params["square_top"]
    ball_like = params["ball_like"]
    lolipop_like = params["lolipop_like"]
    add_graphite_in_copper = params["add_graphite_in_copper"]
    add_spiral = params["add_spiral"]
    elongate_copper = params["elongate_copper"]
    total_len = params["total_len"]
    base_len = params["base_len"]
    main_radius = params["main_radius"]
    main_graphite_thickness = params["main_graphite_thickness"]
    side_cut_size = params["side_cut_size"]
    air_thickness = params["air_thickness"]
    graphite_electrode_radius = params["graphite_electrode_radius"]
    copper_radius = params["copper_radius"]
    inner_insluation_thickness = params["inner_insluation_thickness"]
    aluminium_inner_thicness = params["aluminium_inner_thicness"]
    outer_insluation_thickness = params["outer_insluation_thickness"]
    aluminium_outer_thicness = params["aluminium_outer_thicness"]
    theta_aluminium_outer = params["theta_aluminium_outer"]
    graphite_electrode_connection_radius = params["graphite_electrode_connection_radius"]
    graphite_electrode_thickness = params["graphite_electrode_thickness"]

    spacing = params["spacing"]
    dims = params["dims"]

    rounded_top=params["rounded_top"]

    # Apply spiral modifications if needed
    if !add_spiral
        theta_aluminium_outer = -1.0
        aluminium_outer_thicness = aluminium_outer_thicness * 0.4
    end
    center_circle_radius_aluminium_outer = 2.5

    total_radius = air_thickness + (copper_radius * 1.2) + main_graphite_thickness

    # Set base parameters based on shape type
    if ball_like || lolipop_like
        base_params = Dict("c_x" => main_radius, "c_y" => main_radius, "c_z" => (-1) * (((dims[3] / 2) * spacing[3]) * 0.9), "radius" => total_radius, "l_z" => base_len, "d" => 1.0)
    else
        base_params = Dict("c_x" => main_radius, "c_y" => main_radius, "c_z" => (-1) * (((dims[3] / 4) * spacing[3]) * 0.9), "radius" => total_radius, "l_z" => total_len, "d" => 1.0)
    end

    # Define specifications based on shape type
    if square_top
        specs = Dict(
            "copper_el" => Dict("d" => copper_density, "radius" => copper_radius, "l_z" => base_len, "offset" => 0.0),
            "inner_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius, "l_z" => base_len, "offset" => 0.0),
            "inner_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness, "l_z" => base_len, "offset" => 0.0),
            "outer_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness, "l_z" => base_len, "offset" => 0.0),
            "outer_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness + aluminium_outer_thicness, "l_z" => base_len, "offset" => 0.0, "theta" => theta_aluminium_outer, "center_circle_radius" => center_circle_radius_aluminium_outer),
            "base_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => base_len, "offset" => 0.0),
            "cut" => Dict("d" => 0.0, "radius" => total_radius - 0.00001, "l_z" => base_len, "offset" => 0.0),
            "bottom_main_graphite" => Dict("d" => graphite_density, "radius" => total_radius - 0.00001, "l_z" => main_graphite_thickness, "offset" => base_len),
            "graphite_electrode" => Dict("d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => total_len - base_len - (main_graphite_thickness * 2) - air_thickness, "offset" => base_len + main_graphite_thickness),
            "air" => Dict("d" => 0.000001, "radius" => total_radius - main_graphite_thickness, "l_z" => total_len - base_len - (main_graphite_thickness * 2), "offset" => base_len + main_graphite_thickness),
        )
    end
    if ball_like
        top_sphere_radius = total_len - base_len
        specs = Dict(
            "copper_el" => Dict("d" => copper_density, "radius" => copper_radius, "l_z" => base_len, "offset" => 0.0),
            "inner_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius, "l_z" => base_len, "offset" => 0.0),
            "inner_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness, "l_z" => base_len, "offset" => 0.0),
            "outer_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness, "l_z" => base_len, "offset" => 0.0),
            "outer_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness + aluminium_outer_thicness, "l_z" => base_len, "offset" => 0.0, "theta" => theta_aluminium_outer, "center_circle_radius" => center_circle_radius_aluminium_outer),
            "base_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => base_len, "offset" => 0.0),
            "cut" => Dict("d" => 0.0, "radius" => total_radius - 0.00001, "l_z" => base_len, "offset" => 0.0),
            "bottom_main_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => main_graphite_thickness, "offset" => base_len),
            "graphite_electrode" => Dict("d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => total_len - base_len - main_graphite_thickness, "offset" => base_len + main_graphite_thickness),
            "top_electrode" => Dict("is_halph_sphere" => true, "d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => graphite_electrode_radius, "offset" => (base_len) + total_len - base_len),
            "air" => Dict("is_ellipsoid" => true, "d" => 0.00001, "radius" => top_sphere_radius - (main_graphite_thickness), "l_z" => top_sphere_radius - (main_graphite_thickness), "offset" => (top_sphere_radius / 2) + base_len + (main_graphite_thickness / 2)),
            "outer_graphite" => Dict("is_ellipsoid" => true, "d" => graphite_density, "radius" => top_sphere_radius, "l_z" => top_sphere_radius, "offset" => base_len + (top_sphere_radius / 2)),
        )

    end
    if lolipop_like
        top_sphere_radius = total_len - base_len
        specs = Dict(
            "copper_el" => Dict("d" => copper_density, "radius" => copper_radius, "l_z" => base_len + main_graphite_thickness, "offset" => 0.0),
            "inner_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius, "l_z" => base_len + main_graphite_thickness, "offset" => 0.0),
            "inner_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness, "l_z" => base_len + main_graphite_thickness, "offset" => 0.0),
            "outer_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness, "l_z" => base_len + main_graphite_thickness, "offset" => 0.0),
            "outer_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness + aluminium_outer_thicness, "l_z" => base_len + main_graphite_thickness, "offset" => 0.0, "theta" => theta_aluminium_outer, "center_circle_radius" => center_circle_radius_aluminium_outer),
            "base_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => base_len, "offset" => 0.0),
            "cut" => Dict("d" => 0.0, "radius" => total_radius - 0.00001, "l_z" => base_len, "offset" => 0.0),
            "bottom_main_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => main_graphite_thickness, "offset" => base_len),

            "air" => Dict("d" => 0.00001, "radius" => top_sphere_radius-((air_thickness-graphite_electrode_thickness)/2), "l_z" => air_thickness, "offset" => -1.0),
            "outer_graphite" => Dict("d" => graphite_density, "radius" => top_sphere_radius, "l_z" => air_thickness+main_graphite_thickness * 2, "offset" => -1.0),
            "graphite_electrode" => Dict("d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => graphite_electrode_thickness, "offset" => -1.0),
            "graphite_electrode_connection" => Dict("d" => graphite_density, "radius" => graphite_electrode_connection_radius, "l_z" => ((air_thickness-graphite_electrode_thickness)/2)*1.25, "offset" => base_len + main_graphite_thickness),
        )
    end
    if(rounded_top)
        specs = Dict(
                "copper_el" => Dict("d" => copper_density, "radius" => copper_radius, "l_z" => base_len, "offset" => 0.0),
                "inner_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius, "l_z" => base_len, "offset" => 0.0),
                "inner_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness, "l_z" => base_len, "offset" => 0.0),
                "outer_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness, "l_z" => base_len, "offset" => 0.0),
                "outer_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness + aluminium_outer_thicness, "l_z" => base_len, "offset" => 0.0, "theta" => theta_aluminium_outer, "center_circle_radius" => center_circle_radius_aluminium_outer),
                "base_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => base_len, "offset" => 0.0),
                "cut" => Dict("d" => 0.0, "radius" => total_radius - 0.00001, "l_z" => base_len, "offset" => 0.0),
                "bottom_main_graphite" => Dict("d" => graphite_density, "radius" => total_radius - 0.00001, "l_z" => main_graphite_thickness, "offset" => base_len),
                "graphite_electrode" => Dict("d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => total_len - base_len - main_graphite_thickness, "offset" => base_len + main_graphite_thickness),
                "air" => Dict("d" => 00001, "radius" => total_radius - (main_graphite_thickness * 1.05), "l_z" => total_len - base_len - main_graphite_thickness, "offset" => base_len + main_graphite_thickness),
                "top_inner" => Dict("is_halph_sphere" => true, "d" => 0.0, "radius" => (total_radius - (main_graphite_thickness * 0.95)), "l_z" => (total_radius), "offset" => total_len - main_graphite_thickness),
                "top_outer" => Dict("is_halph_sphere" => true, "d" => graphite_density, "radius" => total_radius, "l_z" => total_radius, "offset" => total_len),
                "top_electrode" => Dict("is_halph_sphere" => true, "d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => graphite_electrode_radius, "offset" => (base_len) + total_len - base_len),
            )
    end

    # Add additional specifications based on parameters
    if add_graphite_in_copper
        specs["in_copper_el"] = Dict("d" => graphite_density, "radius" => copper_radius * 0.4, "l_z" => base_len, "offset" => 0.0)
    end

    if elongate_copper
        specs["electrode_through_graphite"] = Dict("d" => copper_density, "radius" => copper_radius, "l_z" => main_graphite_thickness, "offset" => base_len)

        if add_graphite_in_copper
            specs["in_copper_el_through_graphite"] = Dict("d" => graphite_density, "radius" => copper_radius * 0.4, "l_z" => main_graphite_thickness, "offset" => base_len)
        end
    end

    # Define common angle
    angle = 0.0

    # Build the phantom image
    cylinders_prim, named_res, cylinder_densities,vol_res = build_image(base_params, specs, angle, spacing, dims)

    # Create the final result tensor
    res = zeros(size(named_res["base_cyl"]))
    keyss = keys(named_res)

    print("\n  keyss $keyss \n")

    if !ball_like && !lolipop_like
        for key in keyss
            if (key != "top_outer") && (key != "top_inner")
                res .+= named_res[key]
            end
        end

        if !square_top
            res = res + ((named_res["top_outer"] + named_res["top_inner"]) .* Float64.(named_res["base_cyl"] .== 0))
        end
        res[(named_res["air"].!=0)] .= 0.0
        res = res + named_res["graphite_electrode"]
    end

    if ball_like
        res = zeros(size(named_res["base_cyl"]))
        for key in keyss
            if (key != "bottom_main_graphite") && (key != "graphite_electrode")
                res .+= named_res[key]
            end
        end

        bottom_main_graphite = named_res["bottom_main_graphite"] .!= 0.0
        outer_graphite = named_res["outer_graphite"] .!= 0.0
        graphite_electrode = named_res["graphite_electrode"] .!= 0.0
        air = named_res["air"] .!= 0.0

        res = res + ((bottom_main_graphite .& (.!outer_graphite)) .* graphite_density)
        res = res + ((graphite_electrode .& air) .* graphite_density)
    end

    if lolipop_like
        res = zeros(size(named_res["base_cyl"]))
        for key in keyss
            if (key != "air") && (key != "outer_graphite") && (key != "graphite_electrode") && (key != "bottom_main_graphite")  && (key != "graphite_electrode_connection")
                res .+= named_res[key]
            end
        end

        air = named_res["air"] .!= 0.0
        outer_graphite = named_res["outer_graphite"] .!= 0.0
        graphite_electrode = named_res["graphite_electrode"] .!= 0.0

        axis = 2
        theta = 90

        air_rotated = rotate_and_retrieve_array(air, axis, theta)
        outer_graphite_rotated = rotate_and_retrieve_array(outer_graphite, axis, theta)
        graphite_electrode_rotated = rotate_and_retrieve_array(graphite_electrode, axis, theta)

        air_rotated = reverse_and_add(air_rotated, 3)
        outer_graphite_rotated = reverse_and_add(outer_graphite_rotated, 3)
        graphite_electrode_rotated = reverse_and_add(graphite_electrode_rotated, 3)

        top_sphere_radius = total_len - base_len
        function get_distance_from_edge(dims, spacing)
            return (dims[3] * spacing[3]) / 2
        end

        distance_to_center = get_distance_from_edge(dims, spacing)
        base_end = base_len
        top_sphere_radius = total_len - base_len
        direction = "down"
        distance = abs(distance_to_center - top_sphere_radius - base_end)

        if (distance_to_center - top_sphere_radius) < base_end
            direction = "up"
        end

        axis = 3

        base = named_res["base_graphite"] .!= 0.0

        distance_moved, final_direction, positioned_outer_graphite = iterative_position_adjustment(
            base,
            outer_graphite_rotated,
            3,
            spacing
        )

        res .= res + (named_res["bottom_main_graphite"].!=0.0) .* graphite_density

        res_bool = res .!= 0.0


        if(final_direction!="none")
            air_rotated = move_image(air_rotated, axis, distance_moved, final_direction, spacing)
            outer_graphite_rotated = move_image(outer_graphite_rotated, axis, distance_moved, final_direction, spacing)
            graphite_electrode_rotated = move_image(graphite_electrode_rotated, axis, distance_moved, final_direction, spacing)
        end

        air_rotated[res_bool.&air_rotated] .= false
        outer_graphite_rotated[res_bool.&outer_graphite_rotated] .= false
        graphite_electrode_rotated[res_bool.&graphite_electrode_rotated] .= false

        combined_bool_array_electrode = copy(graphite_electrode_rotated)

        graphite_electrode_connection_bool=((named_res["graphite_electrode_connection"].!=0.0).&(.!combined_bool_array_electrode) )
        graphite_electrode_connection_bool=graphite_electrode_connection_bool .& (.!res_bool)
        res .= res + graphite_electrode_connection_bool.* graphite_density
        res_bool = res .!= 0.0

        air_bool= air_rotated.&(.!graphite_electrode_rotated)
        air_bool = air_bool .& (.!graphite_electrode_connection_bool)


        air_rotated = air_rotated .* cylinder_densities["air"]
        outer_graphite_rotated = outer_graphite_rotated .* cylinder_densities["outer_graphite"]
        graphite_electrode_rotated = graphite_electrode_rotated .* cylinder_densities["graphite_electrode"]

        numerical_graphite_electrode_connection_bool = sum(graphite_electrode_connection_bool.& (.!combined_bool_array_electrode))
        numerical_graphite_electrode_connection_bool = numerical_graphite_electrode_connection_bool * (spacing[1] * spacing[2] * spacing[3])
        res .= res + air_rotated .+ outer_graphite_rotated .+ graphite_electrode_rotated
    end


    copper_el_bool= named_res["copper_el"] .!= 0.0
    if("electrode_through_graphite" in keys(named_res))
        copper_el_bool= copper_el_bool .| (named_res["electrode_through_graphite"] .!= 0.0)
    end

    if !lolipop_like

        if "graphite_electrode" in keys(named_res) && "top_electrode" in keys(named_res)
            combined_bool_array_electrode = (named_res["graphite_electrode"] .!= 0.0) .| (named_res["top_electrode"] .!= 0.0)
        elseif "graphite_electrode" in keys(named_res)
            combined_bool_array_electrode = named_res["graphite_electrode"] .!= 0.0

        end

        air_bool= named_res["air"] .!= 0.0

        if("top_inner" in keys(named_res))
            air_bool = air_bool .| (named_res["top_inner"] .!= 0.0)
        end


    end

    air_bool = air_bool .& (.!combined_bool_array_electrode)

    ### calculate the volume of the air
    #in case of graphite_electrode it is just air minus graphite_electrode
    #in case of ball_like it is just air minus (graphite_electrode+top_electrode)
    #in case of lolipop_like it is just air minus (graphite_electrode+graphite_electrode_connection)
    # in rounded top it is just air plus top_inner minus (graphite_electrode+top_outer)
    air_vol=vol_res["air"]
    if("graphite_electrode" in keys(vol_res))
        air_vol=air_vol-vol_res["graphite_electrode"]
    end
    if("top_electrode" in keys(vol_res))
        air_vol=air_vol-vol_res["top_electrode"]
    end
    if("graphite_electrode_connection" in keys(vol_res))
        air_vol=air_vol-(numerical_graphite_electrode_connection_bool)cm^3#(vol_res["graphite_electrode_connection"]/1.25)
    end
    if("top_inner" in keys(vol_res))
        air_vol=air_vol+vol_res["top_inner"]
    end

    res_raw=copy(res)
    res = res .+ (randn(size(res)) .* additive_noise)
    # res = res .* ((rand(size(res))))
    if(add_smooth)
        σ = (0.5, 0.5, 0.5)
        gaussian_kernel = Kernel.gaussian(σ)
        res = imfilter(res, gaussian_kernel)
    end

    # Calculate the numerical volume of air
    air_volume_voxels = sum(air_bool)
    air_volume_numerical = air_volume_voxels * (spacing[1] * spacing[2] * spacing[3])


    return res,res_raw,air_bool,copper_el_bool,combined_bool_array_electrode,air_vol,air_volume_numerical,named_res
end

"""
    save_ionic_chamber_params(params, filename)

Save the ionic chamber parameters to a JSON file.
"""
function save_ionic_chamber_params(params, filename)
    open(filename, "w") do io
        JSON.print(io, params)
    end
end


"""
    json_params(args_json_path, uuid=nothing, dims=(500, 500, 500))

Generate a dictionary of parameters for the ionic chamber by loading them from a JSON file.

# Arguments
- `args_json_path`: Path to the JSON file containing chamber parameters
- `uuid`: Optional UUID (if not provided, uses the one from the JSON or generates a new one)
- `dims`: Optional dimensions (if not in JSON)

# Returns
A dictionary with all parameters needed for creating an ionic chamber phantom
"""
function json_params(args_json_path, uuid=nothing, dims=(128, 128, 128))
    # Check if file exists
    if !isfile(args_json_path)
        error("JSON file not found at path: $args_json_path")
    end

    # Load parameters from JSON file
    params = JSON.parsefile(args_json_path)

    # Use provided UUID or the one from JSON, or generate a new one
    if uuid === nothing
        if haskey(params, "uuid")
            uuid = params["uuid"]
        else
            uuid = string(UUIDs.uuid4())
        end
    end

    # Extract dimensions from JSON or use default
    dims = haskey(params, "dims") ? params["dims"] : dims

    # Extract or set default values for all parameters
    graphite_density = get(params, "graphite_density", 1.0)
    copper_density = get(params, "copper_density", 2.67)
    aluminium_density = get(params, "aluminium_density", 1.2)
    insulation_density = get(params, "insulation_density", 0.13)

    square_top = get(params, "square_top", false)
    ball_like = get(params, "ball_like", false)
    lolipop_like = get(params, "lolipop_like", false)
    rounded_top = get(params, "rounded_top", false)

    add_graphite_in_copper = get(params, "add_graphite_in_copper", false)
    add_spiral = get(params, "add_spiral", false)
    elongate_copper = get(params, "elongate_copper", false)
    new_flat_sizes = get(params, "new_flat_sizes", false)

    total_len = get(params, "total_len", 2.65)
    base_len = get(params, "base_len", 0.8)
    main_radius = get(params, "main_radius", 0.625)
    main_graphite_thickness = get(params, "main_graphite_thickness", 0.25)
    side_cut_size = get(params, "side_cut_size", 0.2)
    air_thickness = get(params, "air_thickness", 0.265)
    graphite_electrode_radius = get(params, "graphite_electrode_radius", 0.125)
    copper_radius = get(params, "copper_radius", 0.075)
    inner_insluation_thickness = get(params, "inner_insluation_thickness", 0.03)
    aluminium_inner_thicness = get(params, "aluminium_inner_thicness", 0.03)
    outer_insluation_thickness = get(params, "outer_insluation_thickness", 0.04)
    aluminium_outer_thicness = get(params, "aluminium_outer_thicness", 0.2)
    theta_aluminium_outer = get(params, "theta_aluminium_outer", -1.0)
    graphite_electrode_connection_radius = get(params, "graphite_electrode_connection_radius", 0.15)
    graphite_electrode_thickness = get(params, "graphite_electrode_thickness", 0.1)

    # Get spacing from JSON or calculate based on dimensions
    spacing = if haskey(params, "spacing")
        params["spacing"]
    else
        max_y_x = 2.0 * 1.2
        max_z = 7.0 * 1.2
        spacing_calc = (max_y_x / dims[1], max_y_x / dims[2], max_z / dims[3])
        (maximum(spacing_calc), maximum(spacing_calc), maximum(spacing_calc))
    end

    # Create and return the parameter dictionary
    result = Dict(
        "seed" => get(params, "seed", abs(hash(uuid))),
        "uuid" => string(uuid),
        "graphite_density" => graphite_density,
        "copper_density" => copper_density,
        "aluminium_density" => aluminium_density,
        "insulation_density" => insulation_density,
        "square_top" => square_top,
        "ball_like" => ball_like,
        "lolipop_like" => lolipop_like,
        "rounded_top" => rounded_top,
        "add_graphite_in_copper" => add_graphite_in_copper,
        "add_spiral" => add_spiral,
        "elongate_copper" => elongate_copper,
        "total_len" => total_len,
        "base_len" => base_len,
        "main_radius" => main_radius,
        "main_graphite_thickness" => main_graphite_thickness,
        "side_cut_size" => side_cut_size,
        "air_thickness" => air_thickness,
        "graphite_electrode_radius" => graphite_electrode_radius,
        "copper_radius" => copper_radius,
        "inner_insluation_thickness" => inner_insluation_thickness,
        "aluminium_inner_thicness" => aluminium_inner_thicness,
        "outer_insluation_thickness" => outer_insluation_thickness,
        "aluminium_outer_thicness" => aluminium_outer_thicness,
        "theta_aluminium_outer" => theta_aluminium_outer,
        "graphite_electrode_connection_radius" => graphite_electrode_connection_radius,
        "graphite_electrode_thickness" => graphite_electrode_thickness,
        "spacing" => spacing,
        "dims" => dims,
        "new_flat_sizes"=>new_flat_sizes,
        "rand_ver"=>get(params, "rand_ver", 1)
    )

    # Print parameters for debugging
    for (key, value) in result
        println("$key => $value")
    end

    return result
end

"""
    generate_random_params(uuid,dims::Tuple{Int,Int,Int}=(500, 500, 500), randomize=true,variable_spacing=false)

Generate a dictionary of random parameters for the ionic chamber.

# Chamber Types
The system can generate five major types of ionic chambers:

## 1. Square-Top Chamber (flat-topped)
Selected when rand_int == 1. Key parameters:
- Flat top surface instead of curved
- Main structure is cylindrical
- Air cavity follows the cylindrical shape
- Typical specs:
  * total_len: 1.15-1.25 cm  (base_len + top portion)
  * base_len: ~0.8 cm
  * air_thickness: ~0.265 cm

## 2. Ball-Shaped Chamber (spherical top)
Selected when rand_int == 2. Key parameters:
- Spherical top structure with various axis ratios
- Can be perfectly spherical or oval
- Base remains cylindrical
- Special parameters:
  * total_len: base_len + 1.3 cm (shorter than other types)
  * top_sphere_radius: total_len - base_len
  * main_graphite_thickness: 0.3 cm (thicker than default)
  * elongate_copper: automatically true

## 3. Lollipop-Shaped Chamber (thin stem with flat head)
Selected when rand_int == 3. Key parameters:
- Thin stem connecting to flat cylindrical head
- Extended top structure
- Special parameters:
  * air_thickness: 0.6 cm (larger than default)
  * graphite_electrode_radius: sized to match top diameter (up to 2.0 cm)
  * total_len: base_len + 2.5 cm (longer than default)
  * graphite_electrode_connection_radius: 0.15 cm

## 4. Standardized Chamber (precise volumes)
Selected when rand_int == 4. Enables new_flat_sizes=true with predefined volume standards.
- Creates chambers with exact volumes based on rand_ver (1-3):
  * Version 1 (rand_ver=1): 60,000 mm³, radius=20mm, height=48mm (largest)
  * Version 2 (rand_ver=2): 30,000 mm³, radius=15mm, height=42.9mm (medium)
  * Version 3 (rand_ver=3): 10,000 mm³, radius=10mm, height=32.5mm (smallest)
- Special parameters:
  * main_graphite_thickness: random between 0.01-0.05 cm
  * main_radius: depends on version (1.0, 1.5, or 2.0 cm) plus graphite thickness
  * main_part_len: depends on version (3.25, 4.29, or 4.8 cm)
  * graphite_electrode_radius: fixed at 0.15 cm
  * Always sets ball_like=false and lolipop_like=false
  * square_top is randomly set (50% probability)

## 5. Rounded-Top Chamber (default)
Selected when rand_int == 5 or others are false. Key parameters:
- Standard cylindrical base with half-ellipsoid top
- Balanced proportions between parts
- Standard parameter ranges apply

# Shared Parameters
Parameters common to all chamber types (some with varying defaults):

## Material Densities
- graphite_density: Base density unit (relative value: 1.0)
- copper_density: Copper electrode density (relative to graphite: ~2.67)
- aluminium_density: Aluminum layer density (relative to graphite: ~1.2)
- insulation_density: Insulation layer density (relative to graphite: ~0.13)

## Structural Components
- total_len: Total chamber length (typically 2.65 cm but varies by type)
- base_len: Length of base portion (typically 0.8 cm)
- main_radius: Main radius of the chamber (typically 0.625 cm)
- main_graphite_thickness: Thickness of graphite shell (typically 0.25 cm)
- side_cut_size: Size of side cuts (typically 0.2 cm)
- air_thickness: Thickness of the air cavity (typically 0.265 cm)
- copper_radius: Radius of the copper electrode (typically 0.075 cm)
- inner_insluation_thickness: Inner insulation layer thickness (typically 0.03 cm)
- aluminium_inner_thicness: Inner aluminum layer thickness (typically 0.03 cm)
- outer_insluation_thickness: Outer insulation layer thickness (typically 0.04 cm)
- aluminium_outer_thicness: Outer aluminum layer thickness (typically 0.2 cm)
- graphite_electrode_radius: Radius of the graphite electrode (typically 0.125 cm)
- graphite_electrode_connection_radius: Connection radius (typically 0.15 cm)
- graphite_electrode_thickness: Electrode thickness (typically 0.1 cm)

## Optional Features
- add_graphite_in_copper: Boolean to include graphite inside copper electrode
- add_spiral: Boolean to add spiral features to outer aluminum layer
- elongate_copper: Boolean to extend copper through graphite
- theta_aluminium_outer: Spiral angle for outer aluminum (-1.0 when spirals disabled)

## Control Parameters
- randomize: When true, all dimensional parameters are varied by factors between 0.9-1.1
- variable_spacing: When true, spacing is adapted to make chamber fill most of image

"""
function generate_random_params(uuid,dims::Tuple{Int,Int,Int}=(500, 500, 500), randomize=true,variable_spacing=false)


    if(!variable_spacing)
        max_y_x = (3.0 * 1.2)#/10
        max_z = (7.0 * 1.2)#/10

        spacing = (max_y_x / dims[1], max_y_x / dims[2], max_z / dims[3])
        spacing = (maximum(spacing), maximum(spacing), maximum(spacing))
    end
    seed = abs(hash(uuid))
    rng = Random.MersenneTwister(seed)

    graphite_density=15
    # graphite_density = 4
    copper_density = 40
    aluminium_density = 18
    insulation_density = 2

    copper_density = copper_density / graphite_density
    aluminium_density = aluminium_density / graphite_density
    insulation_density = insulation_density / graphite_density
    graphite_density = graphite_density / graphite_density


    rand_int=rand(rng, 1:5)


    # square_top = false
    # ball_like = false
    # lolipop_like = true

    square_top = round(rand_int) == 1
    ball_like = round(rand_int) == 2
    lolipop_like = round(rand_int) == 3
    rounded_top = (round(rand_int) == 5)
    new_flat_sizes = round(rand_int) == 4

    rand_ver=rand(rng, 1:3)

    # new_flat_sizes are abour flat end and rounded ones
    if(new_flat_sizes)
        ball_like=false
        lolipop_like=false
        square_top=rand() < 0.5
        rounded_top=!square_top
    end




    add_graphite_in_copper = rand(rng) < 0.5
    add_spiral = rand(rng) < 0.5
    elongate_copper = rand(rng) < 0.5

    if (ball_like)
        elongate_copper = false

    end




    total_len = 11.5 + 15
    # base_len = 11.5 # len at least 7
    base_len = 8.0 # len at least 7
    main_radius = 6.25
    # main_radius = 6.25

    main_graphite_thickness = 2.5
    side_cut_size = 2.0
    air_thickness = 2.65
    graphite_electrode_radius = 1.25

    copper_radius = 0.75
    inner_insluation_thickness = 0.3
    aluminium_inner_thicness = 0.3
    outer_insluation_thickness = 0.4
    aluminium_outer_thicness = 2.0
    theta_aluminium_outer = 0.45

    if (ball_like)
        total_len = base_len + 13
        main_graphite_thickness = 3.0

        elongate_copper = true
    end



    # odległość w środku między zewnętrznymi częściami komory wynosi 6mm
    #, szerokość wewnętrznej elektrody 1mm.
    #  zewnętrzny promień komory 25mm, wewnętrzny 22,25
    #, promień elektrody 20mm. elektroda jest ustawiona centrycznie
    # zamocowana na trzpieniu
    #  o promieniu 1,5mm natomiast sama komora zamocowana
    #jest do trzebienia o promieniu 3,5mm.
    graphite_electrode_connection_radius = 1.5
    graphite_electrode_thickness = 1.0

    if (lolipop_like)
        air_thickness = 6.0
        graphite_electrode_radius = 20.0

        total_len = base_len + 25.0
        # graphite_density = graphite_density * 0.5

    end




    #divided to get mm
    # total_len = total_len / 10.0
    # base_len = base_len / 10.0
    # main_radius = main_radius / 10.0
    # main_graphite_thickness = main_graphite_thickness / 10.0
    # side_cut_size = side_cut_size / 10.0
    # air_thickness = air_thickness / 10.0
    # graphite_electrode_radius = graphite_electrode_radius / 10.0
    # copper_radius = copper_radius / 10.0
    # inner_insluation_thickness = inner_insluation_thickness / 10.0
    # aluminium_inner_thicness = aluminium_inner_thicness / 10.0
    # outer_insluation_thickness = outer_insluation_thickness / 10.0
    # aluminium_outer_thicness = aluminium_outer_thicness / 10.0
    # graphite_electrode_connection_radius = graphite_electrode_connection_radius / 10.0
    # graphite_electrode_thickness = graphite_electrode_thickness / 10.0

    theta_aluminium_outer = theta_aluminium_outer * π / 180.0

    if (randomize)
        scale_total_len = rand(rng, 0.9:0.001:1.1)
        scale_base_len = rand(rng, 0.9:0.001:1.1)
        scale_main_radius = rand(rng, 0.9:0.001:1.1)
        scale_main_graphite_thickness = rand(rng, 0.8:0.001:1.1)
        scale_side_cut_size = rand(rng, 0.9:0.001:1.1)
        scale_air_thickness = rand(rng, 0.9:0.001:1.1)
        scale_graphite_electrode_radius = rand(rng, 0.9:0.001:1.1)
        scale_copper_radius = rand(rng, 0.9:0.001:1.1)
        scale_inner_insluation_thickness = rand(rng, 0.9:0.001:1.1)
        scale_aluminium_inner_thicness = rand(rng, 0.9:0.001:1.1)
        scale_outer_insluation_thickness = rand(rng, 0.9:0.001:1.1)
        scale_aluminium_outer_thicness = rand(rng, 0.9:0.001:1.1)
        scale_graphite_electrode_connection_radius = rand(rng, 0.9:0.001:1.1)
        scale_graphite_electrode_thickness = rand(rng, 0.9:0.001:1.1)

        # Apply the scale factors and divide by 10.0
        total_len = (total_len * scale_total_len) / 10.0
        base_len = (base_len * scale_base_len) / 10.0
        main_radius = (main_radius * scale_main_radius) / 10.0
        main_graphite_thickness = (main_graphite_thickness * scale_main_graphite_thickness) / 10.0
        side_cut_size = (side_cut_size * scale_side_cut_size) / 10.0
        air_thickness = (air_thickness * scale_air_thickness) / 10.0
        graphite_electrode_radius = (graphite_electrode_radius ) / 10.0 # we do not want it random
        copper_radius = (copper_radius * scale_copper_radius) / 10.0
        inner_insluation_thickness = (inner_insluation_thickness * scale_inner_insluation_thickness) / 10.0
        aluminium_inner_thicness = (aluminium_inner_thicness * scale_aluminium_inner_thicness) / 10.0
        outer_insluation_thickness = (outer_insluation_thickness * scale_outer_insluation_thickness) / 10.0
        aluminium_outer_thicness = (aluminium_outer_thicness * scale_aluminium_outer_thicness) / 10.0
        graphite_electrode_connection_radius = (graphite_electrode_connection_radius * scale_graphite_electrode_connection_radius) / 10.0
        graphite_electrode_thickness = (graphite_electrode_thickness * scale_graphite_electrode_thickness) / 10.0


        """
to 60 000 mm3 mają wewnętrzny walec o promieniu 20mm, mniejsze 30 000mm3 o promieniu 15 mm
, a 10 0000mm3 mają o promieniu 10mm
        """


        if (new_flat_sizes)
            """
            wzorce grafitowe i w nich wymiary wewnętrznych walców będą wynosiły - r 20mm, h 48mm, mniejsza r 15mm, h 42,9mm, i najmniejsza r 10mm, h 32,5mm
            """
            main_graphite_thickness = rand(rng, 0.1:0.5)
            if(rand_ver==1)
                main_radius =(2.0*scale_main_radius)+main_graphite_thickness*2
                main_part_len=(4.8*scale_total_len)
                graphite_electrode_radius=0.15
            end
            if(rand_ver==2)
                main_radius =(1.5*scale_main_radius)+main_graphite_thickness*2
                main_part_len=(4.29*scale_total_len)
                graphite_electrode_radius=0.15
            end
            if(rand_ver==3)
                main_radius =(1.0*scale_main_radius)+main_graphite_thickness*2
                main_part_len=(3.25*scale_total_len)
            end
            if(rounded_top)
                main_part_len=main_part_len-main_radius*0.9
                main_part_len = min(main_part_len, 3.25 * main_radius)
            else
                main_part_len = min(main_part_len, 3.25 * main_radius)
            end


            total_len=base_len+main_part_len

            graphite_electrode_radius=0.15
            main_radius = rand(rng, 1.0:0.001:2.0)
        end


    end

    if(lolipop_like)
        graphite_electrode_radius=   (total_len-base_len)-0.5
    end


    if (!add_spiral)
        theta_aluminium_outer = -1.0
        aluminium_outer_thicness = aluminium_outer_thicness * 0.4
    else

    end
    center_circle_radius_aluminium_outer = 2.5


    total_radius = air_thickness + (copper_radius * 1.2) + main_graphite_thickness

    if(variable_spacing)
        max_y_x = total_radius*3.0#/10
        max_z = total_len*2.0#/10

        spacing = (max_y_x / dims[1], max_y_x / dims[2], max_z / dims[3])
        spacing = (maximum(spacing), maximum(spacing), maximum(spacing))
    end


    # Return parameters in a dictionary
    res= Dict(
        "seed" => seed,
        "uuid" => string(uuid),
        "graphite_density" => graphite_density,
        "copper_density" => copper_density,
        "aluminium_density" => aluminium_density,
        "insulation_density" => insulation_density,
        "square_top" => square_top,
        "ball_like" => ball_like,
        "lolipop_like" => lolipop_like,
        "add_graphite_in_copper" => add_graphite_in_copper,
        "add_spiral" => add_spiral,
        "elongate_copper" => elongate_copper,
        "total_len" => total_len,
        "base_len" => base_len,
        "main_radius" => main_radius,
        "main_graphite_thickness" => main_graphite_thickness,
        "side_cut_size" => side_cut_size,
        "air_thickness" => air_thickness,
        "graphite_electrode_radius" => graphite_electrode_radius,
        "copper_radius" => copper_radius,
        "inner_insluation_thickness" => inner_insluation_thickness,
        "aluminium_inner_thicness" => aluminium_inner_thicness,
        "outer_insluation_thickness" => outer_insluation_thickness,
        "aluminium_outer_thicness" => aluminium_outer_thicness,
        "theta_aluminium_outer" => theta_aluminium_outer,
        "graphite_electrode_connection_radius" => graphite_electrode_connection_radius,
        "graphite_electrode_thickness" => graphite_electrode_thickness,
        "spacing" => spacing,
        "dims" => dims,
        "new_flat_sizes"=>new_flat_sizes,
        "rand_ver"=>rand_ver,
        "rounded_top"=>rounded_top
    )


    println("new_flat_sizes => ", new_flat_sizes)

    for (key, value) in res
        println("$key => $value")
    end

    return res
end


function get_random_chamber(dims,uuid,temp_fold,variable_spacing,randomize)

    if args_json_path != " "
        params = json_params(args_json_path)

    else
        params = generate_random_params(uuid,dims, variable_spacing,randomize)

    end

    add_str=""
    if(params["square_top"])
        add_str="_square_top"
    elseif (params["ball_like"])
        add_str="_ball_like"
    elseif (params["lolipop_like"])
        add_str="_lolipop_like"
    else
        add_str="_rounded_top"
    end
    add_str=add_str*"_variable_spacing_$(variable_spacing)_randomize_$(randomize)_add_radon_$(add_radon)_"

    main_folder="$temp_fold/$(add_str)$(uuid)"
    # Create base_folder if it does not exist
    if !isdir(main_folder)
        mkpath(main_folder)
    end

    phantom_image,res_raw ,air_bool, copper_el_bool, combined_bool_array_electrode,air_vol,air_volume_numerical,named_res = create_ionic_chamber_phantom(params)


    #### save raw data
    img = sitk.GetImageFromArray(Float32.(res_raw))
    img.SetSpacing((params["spacing"][1] * 10, params["spacing"][2] * 10, params["spacing"][3] * 10))
    sitk.WriteImage(img, "$main_folder/ionic_chamber_raw.nii.gz")

    reference_dicom_path = joinpath(main_folder, "ionic_chamber_raw")
    save_sitk_image_as_dicom(sitk.ReadImage("$main_folder/ionic_chamber_raw.nii.gz"), reference_dicom_path)




    # Save air volume information to params
    params["air_vol"] = Unitful.ustrip(cm^3,air_vol)
    params["air_volume_numerical"] = air_volume_numerical

    # Print air volume information
    println("Air volume (analytical): ", air_vol)
    println("Air volume (numerical): ", air_volume_numerical)

    # Compare analytical and numerical air volumes
    volume_difference = abs(Unitful.ustrip(cm^3,air_vol) - air_volume_numerical)
    volume_difference=volume_difference/(Unitful.ustrip(cm^3,air_vol)+air_volume_numerical)
    println("Difference between analytical and numerical air volumes: $volume_difference %")

    save_ionic_chamber_params(params, "$main_folder/ionic_chamber_params.json")

    # Save the phantom image
    spacing= params["spacing"]
    phantom_image=Float32.(phantom_image)
    img = sitk.GetImageFromArray(phantom_image)
    img.SetSpacing((params["spacing"][1] * 10, params["spacing"][2] * 10, params["spacing"][3] * 10))
    sitk.WriteImage(img, "$main_folder/ionic_chamber.nii.gz")

    print(" \n add_radon $add_radon \n")
    if (add_radon)
        print("\n execuuute radon \n")
        # Execute the Python script using the run function
        input_path = "$(main_folder)/ionic_chamber.nii.gz"
        output_path = "$(main_folder)/after_radon.nii.gz"
        output_path_b = "$(main_folder)/after_radon_plus_before.nii.gz"
        script_path = joinpath(@__DIR__, "get_approximate_radon_inverse.py")

        command = `python3 $script_path $input_path $output_path $output_path_b`

        run(command)

        save_sitk_image_as_dicom(sitk.ReadImage(output_path), "$(main_folder)/after_radon")
        save_sitk_image_as_dicom(sitk.ReadImage(output_path_b), "$(main_folder)/after_radon_plus_before")
        print("\n finish  radon \n")

    end




    reference_dicom_path = joinpath(main_folder, "ionic_chamber")
    save_sitk_image_as_dicom(sitk.ReadImage(joinpath(main_folder, "ionic_chamber.nii.gz")), reference_dicom_path)
    spacing = (params["spacing"][1], params["spacing"][2], params["spacing"][3])

    save_mask_as_nifti(
        Array(air_bool),
        joinpath(main_folder, "air_bool.nii.gz"),
        spacing)
    convert_nifti_to_dicom_seg(joinpath(main_folder, "air_bool.nii.gz"), reference_dicom_path, joinpath(main_folder, "air_bool"))


    save_mask_as_nifti(
        Array(copper_el_bool),
        joinpath(main_folder, "copper_el_bool.nii.gz"),
        spacing
    )

    convert_nifti_to_dicom_seg(joinpath(main_folder, "copper_el_bool.nii.gz"), reference_dicom_path, joinpath(main_folder, "copper_el_bool"))

    save_mask_as_nifti(
        Array(combined_bool_array_electrode),
        joinpath(main_folder, "combined_bool_array_electrode.nii.gz"),
        spacing
    )
    convert_nifti_to_dicom_seg(joinpath(main_folder, "combined_bool_array_electrode.nii.gz"), reference_dicom_path, joinpath(main_folder, "combined_bool_array_electrode"))

    if("graphite_electrode_connection" in keys(named_res))
        save_mask_as_nifti(
            Array(named_res["graphite_electrode_connection"].!=0.0),
            joinpath(main_folder, "graphite_electrode_connection.nii.gz"),
            spacing
        )
        # convert_nifti_to_dicom_seg(joinpath(main_folder, "graphite_electrode_connection.nii.gz"), reference_dicom_path, joinpath(main_folder, "graphite_electrode_connection"))
    end

    file_name = "ionic_ch_$(add_str)_$(dims[1])|$(dims[2])|$(dims[3])_add_radon_$(add_radon)_$uuid"

    zip_path = "$(temp_fold)/$(file_name).zip"
    run(`zip -r $zip_path $main_folder`)

    if !haskey(ENV, "SKIP_UPLOAD")
        command = `gcloud storage cp $zip_path gs://metro_tk_kplayground/ionic-chambersx128/$file_name`
        # Execute the command
        run(command)
        rm(main_folder; force=true, recursive=true)
        rm(zip_path; force=true)
    else
        println("Skipping upload and deletion for debugging/testing.")
        println("Output stored in: $main_folder")
        println("Zip stored in: $zip_path")
    end

    return add_str
end




# Example usage

# dims = (128, 128, 128)

# while true
    # uuid = UUIDs.uuid4()
    # temp_fold="/workspaces/synthethic_tomo/data/debug/"
    temp_fold = mktempdir(; cleanup=false)
    # variable_spacing=true
    get_random_chamber(dims,uuid,temp_fold,variable_spacing,randomize)



# end





# julia in_docker_organized/main_create_phantom_ionic_chamber.jl 256x256x256 true true ddac5ff5-1cc4-4447-9c44-96c372134657 true true 0.1












# using Printf

# volumes = [10000.0, 30000.0, 60000.0]  # mm^3
# for v in volumes
#     r = ((3v) / (4π))^(1/3)
#     @printf "Volume: %.1f mm^3 -> Radius: %.3f mm\n" v r
# end