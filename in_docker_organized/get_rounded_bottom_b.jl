__precompile__(false)

using Pkg
using Sinograms: SinoPar, rays, plan_fbp, fbp, fbp_sino_filter, CtFanArc, CtFanFlat, Window, Hamming, fdk, ct_geom_plot3, project_bdd, backproject_bdd
using ImageGeoms: ImageGeom, fovs, MaskCircle, axesf
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom, Object, spectrum, Cylinder, cylinder, ellipsoid_parameters, ellipsoid, Ellipsoid
using Unitful
using Unitful: mm, unit, °, cm, ustrip
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Plots: plot, plot!, scatter!, default, gui
using PyCall
# using Revise, Statistics # Remove Revise
using Statistics
using UUIDs
using JSON
using Dates, Logging
using Accessors
using ImageFiltering, Revise
import ImagePhantoms as IP

includet("get_geometry_main.jl")
includet("geometry_utils.jl")
sitk = pyimport_conda("SimpleITK", "simpleitk")



"""
    get_cylinder_bool_mask(ig, ob_el) -> BitArray{3}

Generate a Boolean mask for a cylinder object using the provided ImageGeom.

# Arguments
- `ig`: ImageGeom defining the sampling grid
- `ob_el`: Cylinder object or array of objects

# Returns
A 3D boolean array representing the cylinder volume
"""
function get_cylinder_bool_mask(ig, ob_el)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end




"""
    create_sphere(ig, center, r_xy, r_z)

Create a Boolean array of an ellipsoid using the provided center and radii.
Returns the boolean array and the volume (in cm³).
"""
function create_sphere(ig, center::Tuple{Float64,Float64,Float64}, r_xy::Float64, r_z::Float64,return_raw=false)
    r_xy = max(r_xy, 0.1)
    r_z = max(r_z, 0.1)
    test_sphere = ellipsoid(
        center[1] * cm,
        center[2] * cm,
        center[3] * cm,
        r_xy * cm, r_xy * cm, r_z * cm,
        0, 0, 0, 1.0
    )

    if return_raw
        return get_cylinder_bool_mask(ig, [test_sphere]), ustrip(IP.volume(test_sphere)),test_sphere
    end

    # Return the boolean array and the volume as a pure number (strip units)
    return get_cylinder_bool_mask(ig, [test_sphere]), ustrip(IP.volume(test_sphere))
end


"""
    create_torus(ig, center_torus, r_xy, r_z, base_cylinder_radius, spacing)

Create a torus-like shape by:
1) Iterating over angles around the x-axis (in small increments determined by the arc length).
2) For each angle, compute the center of a new sphere at distance (base_cylinder_radius - r_xy) in the YZ-plane relative to `center_torus`.
3) Combine each new sphere via logical OR to form the full torus.
"""
function create_torus(
    ig,
    center_torus::Vector{Float64},
    r_xy::Float64,
    r_z::Float64,
    base_cylinder_radius::Float64,
    spacing::Tuple{Float64,Float64,Float64}, dims
)
    # Calculate the effective radius in the YZ-plane
    orbit_radius = base_cylinder_radius - r_xy

    # Prepare an empty Boolean volume to accumulate the torus
    result = fill(false, dims)

    # Calculate angle step based on arc length
    min_spacing = minimum(spacing)
    angle_rad = min_spacing / orbit_radius
    angle_deg = rad2deg(angle_rad)
    num_steps = ceil(Int, 360 / angle_deg)
    angle_step = 360 / num_steps
    println("Angle step: $angle_step num_steps $num_steps")

    # Iterate over angles, place sphere at each angle, then combine
    # Initialize vector of vectors to store sphere positions
    sphere_positions = Vector{Vector{CartesianIndex{3}}}(undef, num_steps)
    total_sphere_volume = 0.0  # Initialize as dimensionless

    Threads.@threads for i in 0:(num_steps-1)
        angle = i * angle_step
        angle_rad2 = deg2rad(angle)
        print("*")

        # Rotate around x-axis in the YZ-plane
        # The center shifts in X and Z by orbit_radius
        x_offset = orbit_radius * cos(angle_rad2)
        y_offset = orbit_radius * sin(angle_rad2)
        new_center = (
            center_torus[1] + x_offset,
            center_torus[2] + y_offset,
            center_torus[3]
        )

        # Create a sphere at the new center and combine via OR
        # The volume returned by create_sphere is now dimensionless (pure number)
        new_sphere, sphere_volume = create_sphere(ig, new_center, r_xy, r_z)
        pos = findall(new_sphere)
        sphere_positions[i+1] = pos
        total_sphere_volume += sphere_volume  # Now safe to accumulate
    end

    sphere_positions = vcat(sphere_positions...)
    # Get unique positions
    unique_positions = unique(sphere_positions)

    # Set all entries in result to true where positions are present in unique_positions
    result[unique_positions] .= true

    return result, orbit_radius, r_xy, r_z, total_sphere_volume
end


"""
    calculate_torus_volumes(main_torus, inner_torus, main_torus_radius, main_torus_r_xy, main_torus_r_z,
                           inner_torus_radius, inner_torus_r_xy, inner_torus_r_z,
                           spacing, cylinder_wall_thickness)
                           -> Tuple{Float64, Float64, Float64, Float64}

Calculate volumes of main and inner torus structures using both analytical and numerical methods.

# Arguments
- `main_torus`: Boolean mask for the main torus
- `inner_torus`: Boolean mask for the inner torus
- `main_torus_radius`, `main_torus_r_xy`, `main_torus_r_z`: Dimensions of main torus
- `inner_torus_radius`, `inner_torus_r_xy`, `inner_torus_r_z`: Dimensions of inner torus
- `spacing`: Voxel spacing tuple
- `cylinder_wall_thickness`: Wall thickness of the cylinder

# Returns
A tuple containing (main_torus_volume_analytical, main_torus_volume_numerical,
                   inner_torus_volume_analytical, inner_torus_volume_numerical)
"""
function calculate_torus_volumes(
    main_torus, inner_torus,
    main_torus_radius, main_torus_r_xy, main_torus_r_z,
    inner_torus_radius, inner_torus_r_xy, inner_torus_r_z,
    spacing, cylinder_wall_thickness
)
    # Method 1: Analytical formula for torus volumes
    # For elliptical torus: V = 2π² × R × r_xy × r_z
    main_torus_volume_analytical = 2 * π^2 * main_torus_radius * main_torus_r_xy * main_torus_r_z
    inner_torus_volume_analytical = 2 * π^2 * (main_torus_radius + cylinder_wall_thickness) * (main_torus_r_xy - cylinder_wall_thickness) * (main_torus_r_z - cylinder_wall_thickness)
    # inner_torus_volume_analytical = 2 * π^2 * inner_torus_radius * inner_torus_r_xy * inner_torus_r_z

    # Method 2: Count voxels in the arrays and multiply by voxel volume
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    main_torus_voxel_count = sum(main_torus)
    inner_torus_voxel_count = sum(inner_torus)

    main_torus_volume_numerical = main_torus_voxel_count * voxel_volume
    inner_torus_volume_numerical = inner_torus_voxel_count * voxel_volume

    # Check overlap and containment between toruses
    overlap = sum(main_torus .& inner_torus)
    overlap_percentage = (overlap / inner_torus_voxel_count) * 100

    # Analyze inner torus shape consistency
    dims = size(inner_torus)
    center_x = div(dims[1], 2)
    center_y = div(dims[2], 2)
    center_z = div(dims[3], 2)

    # Check if inner torus has a hollow center (true torus should)
    center_region = inner_torus[
        max(1, center_x - 3):min(dims[1], center_x + 3),
        max(1, center_y - 3):min(dims[2], center_y + 3),
        center_z
    ]

    return main_torus_volume_analytical, main_torus_volume_numerical,
    inner_torus_volume_analytical, inner_torus_volume_numerical
end


"""
    get_rounded_bottom(
        spacing,
        dims,
        ig,
        cylinder_wall_thickness,
        r_xy_small,
        r_z_small,
        r_xy_torus,
        center_torus,
        base_cylinder_radius,
        curvature,
        tube_ball_fraction,
        rel_dist
    )

Create outer/inner torus and outer/inner sphere. Returns (main_torus, inner_torus, outer_sphere, inner_sphere, ball1, ball2,
main_torus_volume_analytical, inner_torus_volume_analytical, outer_sphere_volume, ball1_volume, ball2_volume).
Where ball1 and ball2 are positioned at the bottom of inner_torus, sized by tube_ball_fraction and spaced by rel_dist.
"""
function get_rounded_bottom(
    spacing::Tuple{Float64,Float64,Float64},
    dims::Tuple{Int,Int,Int},
    ig::ImageGeom,
    cylinder_wall_thickness::Float64,
    r_xy_small::Float64,
    r_z_small::Float64,
    r_xy_torus::Float64,
    center_torus::Vector{Float64},
    base_cylinder_radius::Float64,
    curvature::Float64,
    tube_ball_fraction::Float64=0.9,
    rel_dist::Float64=0.5
)
    # Create main torus and other elements
    main_torus, main_torus_radius, main_torus_r_xy, main_torus_r_z, main_torus_total_volume = create_torus(ig, center_torus, r_xy_small, r_z_small, base_cylinder_radius, spacing, dims)

    inner_r_xy = r_xy_small - cylinder_wall_thickness
    inner_r_z = r_z_small - cylinder_wall_thickness
    inner_torus, inner_torus_radius, inner_torus_r_xy, inner_torus_r_z, inner_torus_total_volume = create_torus(ig, center_torus, inner_r_xy, inner_r_z, base_cylinder_radius - cylinder_wall_thickness, spacing, dims)

    main_torus_volume_analytical, main_torus_volume_numerical, inner_torus_volume_analytical, inner_torus_volume_numerical = calculate_torus_volumes(
        main_torus, inner_torus,
        main_torus_radius, main_torus_r_xy, main_torus_r_z,
        inner_torus_radius, inner_torus_r_xy, inner_torus_r_z,
        spacing, cylinder_wall_thickness
    )

    # Create spheres
    sphere_center = (center_torus[1], center_torus[2], center_torus[3])

    outer_sphere_r = base_cylinder_radius - (r_xy_small * 2.0)
    outer_sphere, outer_sphere_volume = create_sphere(ig, sphere_center, outer_sphere_r, curvature)

    inner_sphere_r = outer_sphere_r - cylinder_wall_thickness
    inner_sphere_z = curvature
    inner_sphere, inner_sphere_volume = create_sphere(ig, (sphere_center[1], sphere_center[2], sphere_center[3] - cylinder_wall_thickness), inner_sphere_r, inner_sphere_z)

    # Calculate voxel closest to center_torus
    voxel_center = (
        round(Int, (center_torus[1] / spacing[1]) + dims[1] / 2),
        round(Int, (center_torus[2] / spacing[2]) + dims[2] / 2),
        round(Int, (center_torus[3] / spacing[3]) + dims[3] / 2)
    )

    # Apply masks to restrict shapes to appropriate regions
    outer_sphere[:, :, 1:voxel_center[3]] .= false
    main_torus[:, :, voxel_center[3]+1:end] .= false
    inner_torus[:, :, voxel_center[3]+1:end] .= false

    # --- Create ball1 and ball2 ---

    # Calculate ball radius based on tube_ball_fraction
    ball_radius = inner_r_xy * tube_ball_fraction * 0.6  # Reduce to 60% to ensure fit

    # Calculate exact position where balls should be placed
    orbit_radius = base_cylinder_radius - inner_r_xy - cylinder_wall_thickness

    # Calculate Z position - the bottom of the torus
    ball_z_pos = center_torus[3] - inner_r_z + ball_radius + cylinder_wall_thickness

    # Use rel_dist to determine angle separation between balls
    max_angle_separation = π * rel_dist  # at most π radians (180 degrees) apart

    # Calculate ball positions on the circular path
    angle1 = π / 2 - max_angle_separation / 2  # Start from bottom, offset by half the separation
    angle2 = π / 2 + max_angle_separation / 2

    # Position balls on the circular orbit with both X and Y coordinates determined by the angles
    ball1_center = (
        center_torus[1] + orbit_radius * cos(angle1),
        center_torus[2] + orbit_radius * sin(angle1),  # Use sine for Y coordinate
        ball_z_pos
    )

    ball2_center = (
        center_torus[1] + orbit_radius * cos(angle2),
        center_torus[2] + orbit_radius * sin(angle2),  # Use sine for Y coordinate
        ball_z_pos
    )

    # Create the ball boolean arrays
    ball1, ball1_volume,ball1_obj = create_sphere(ig, ball1_center, ball_radius, ball_radius,true)
    ball2, ball2_volume,ball2_obj = create_sphere(ig, ball2_center, ball_radius, ball_radius,true)

    # Ensure balls don't overlap
    overlap = ball1 .& ball2
    if any(overlap)
        # If they still overlap, adjust positions slightly
        adjusted_angle1 = π / 2 - (max_angle_separation / 2) * 1.1  # Increase angle separation by 10%
        adjusted_angle2 = π / 2 + (max_angle_separation / 2) * 1.1

        adjusted_ball1_center = (
            center_torus[1] + orbit_radius * cos(adjusted_angle1),
            center_torus[2] + orbit_radius * sin(adjusted_angle1),
            ball_z_pos
        )

        adjusted_ball2_center = (
            center_torus[1] + orbit_radius * cos(adjusted_angle2),
            center_torus[2] + orbit_radius * sin(adjusted_angle2),
            ball_z_pos
        )

        # Recreate balls with adjusted positions
        ball1, ball1_volume,ball1_obj = create_sphere(ig, adjusted_ball1_center, ball_radius * 0.9, ball_radius * 0.9,true)
        ball2, ball2_volume,ball2_obj = create_sphere(ig, adjusted_ball2_center, ball_radius * 0.9, ball_radius * 0.9,true)
    end

    # Final check and correction for any residual overlap
    overlap = ball1 .& ball2
    if any(overlap)
        println("Warning: Adjusted balls still overlap slightly. Partitioning overlap region.")
        overlap_coords = findall(overlap)
        for (i, coord) in enumerate(overlap_coords)
            if i % 2 == 0
                ball1[coord] = false
            else
                ball2[coord] = false
            end
        end
    end

    # For verification - print ball positions relative to torus
    # println("Inner torus radius: $(inner_r_xy), Ball radius: $(ball_radius)")
    # println("Ball1 center: $(ball1_center), Ball2 center: $(ball2_center)")
    # println("Distance from ball1 to center_torus: $(sqrt((ball1_center[1] - center_torus[1])^2 + (ball1_center[2] - center_torus[2])^2))")
    # println("Distance from ball2 to center_torus: $(sqrt((ball2_center[1] - center_torus[1])^2 + (ball2_center[2] - center_torus[2])^2))")
    # println("Distance between ball centers: $(sqrt((ball1_center[1] - ball2_center[1])^2 + (ball1_center[2] - ball2_center[2])^2 + (ball1_center[3] - ball2_center[3])^2))")
    # println("Ball1 volume: $(ball1_volume), Ball2 volume: $(ball2_volume), Outer sphere volume: $(outer_sphere_volume)")

    return main_torus, inner_torus, outer_sphere, inner_sphere, ball1_obj, ball2_obj,
    main_torus_volume_analytical, inner_torus_volume_analytical,
    outer_sphere_volume, ball1_volume, ball2_volume
end



"""
    save_nifti(image_array, spacing, out_path)

Save a 3D array as a NIfTI file with the specified spacing.

# Arguments
- `image_array`: 3D array to save
- `spacing`: Voxel spacing tuple in physical units
- `out_path`: Output path for the NIfTI file
"""
function save_nifti(image_array, spacing, out_path)
    if eltype(image_array) == Bool
        image_array = UInt8.(image_array)
    end
    img = sitk.GetImageFromArray(image_array)
    img.SetSpacing((
        Unitful.ustrip(spacing[1]),
        Unitful.ustrip(spacing[2]),
        Unitful.ustrip(spacing[3])
    ))
    sitk.WriteImage(img, out_path)
end

# function main()
#     # Use provided parameters
#     spacing = (0.1, 0.1, 0.1)
#     dims = (256, 256, 512)
#     ig = ImageGeom(dims=dims, deltas=(spacing[1]*cm, spacing[2]*cm, spacing[3]*cm))
#     cylinder_wall_thickness = 0.1
#     r_xy_small = 0.7
#     r_z_small = 0.7
#     r_xy_rorus = 3.4
#     center_torus = [0.0, 0.0, -3.85]
#     base_cylinder_radius = 2.0
#     curvature = 0.8
#     tube_ball_fraction = 0.75
#     rel_dist = 0.577

#     results = get_rounded_bottom(
#         spacing, dims, ig, cylinder_wall_thickness, r_xy_small, r_z_small, 
#         r_xy_rorus, center_torus, base_cylinder_radius, curvature, 
#         tube_ball_fraction, rel_dist
#     )

#     main_torus, inner_torus, outer_sphere, inner_sphere = results

#     # Save using SimpleITK
#     main_fold = "/workspaces/synthethic_tomo/data/torus"
#     if isdir(main_fold)
#         rm(main_fold; force=true, recursive=true)
#     end
#     mkpath(main_fold)

#     save_nifti(main_torus, spacing, "$main_fold/main_torus.nii.gz")
#     save_nifti(inner_torus, spacing, "$main_fold/inner_torus.nii.gz")
#     save_nifti(outer_sphere, spacing, "$main_fold/outer_sphere.nii.gz")
#     save_nifti(inner_sphere, spacing, "$main_fold/inner_sphere.nii.gz")

#     println("Torus and spheres saved")
# end

# main()



# Main Torus:
#   - Analytical: 0.898587 cubic units
#   - Numerical (voxel count): 5416 voxels = 0.885815 cubic units
#   - Difference: 0.012773 (1.42%)

# Inner Torus:
#   - Analytical: 0.342319 cubic units
#   - Numerical (voxel count): 3184 voxels = 0.52076 cubic units
#   - Difference: 0.178441 (52.13%)
#   - Inner torus hollow at center: Yes (correct)
# Warning: Adjusted balls still overlap slightly. Partitioning overlap region.

# Inner Torus:
#   - Analytical: 0.293463 cubic units
#   - Numerical (voxel count): 2888 voxels = 0.472347 cubic units
#   - Difference: 0.178885 (60.96%)
#   - Inner torus hollow at center: Yes (correct)
# Warning: Adjusted balls still overlap slightly. Partitioning overlap region.


# pydicom_seg.template.from_dcmqi_metainfo NotImplementedError()