using Sinograms: SinoPar, rays, plan_fbp, fbp, fbp_sino_filter, CtFanArc, CtFanFlat, Window, Hamming, fdk, ct_geom_plot3
using ImageGeoms: ImageGeom, fovs, MaskCircle, axesf
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom, Object, spectrum, Cylinder, cylinder, ellipsoid, ellipsoid_parameters, sphere, half_sphere_x, half_sphere_y, half_sphere_z, half_sphere_z_b, Ellipsoid, Cuboid, cuboid, cylinder_irr
import ImagePhantoms as IP
using Unitful: mm, unit, °, cm
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Plots: plot, plot!, scatter!, default, gui
using PyCall
using Statistics, LinearAlgebra

"""
    get_cylinder_bool_mask(ob_el, ig)

get a phantom cylinder and return a boolean mask of it
"""
function get_cylinder_bool_mask(ob_el, ig)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end

"""
    volume_of_elliptical_cylinder(pipe_cross_section, overlay_length) -> Float64

Calculate the volume of an elliptical cylinder.

# Arguments
- `pipe_cross_section`: Tuple containing the ellipse semi-axes (a, b)
- `overlay_length`: Length of the cylinder

# Returns
The volume of the elliptical cylinder in cubic units
"""
function volume_of_elliptical_cylinder(pipe_cross_section, overlay_length)
    a = pipe_cross_section[1] / 2
    b = pipe_cross_section[2] / 2
    h = overlay_length
    volume = π * a * b * h
    return volume
end

"""
    ellipsoid_segment_intersection(
        center_bottom_curvature::NTuple{3,Float64},
        halph_s_bigger_size::NTuple{3,Float64},
        p1::NTuple{3,Float64},
        p2::NTuple{3,Float64}
    ) -> Union{Nothing, NTuple{3,Float64}}

Finds intersection between line segment [p1,p2] and ellipsoid surface.
Returns point closest to p1 if multiple intersections exist.
Returns nothing if no valid intersection found.
"""
function ellipsoid_segment_intersection(
    center_bottom_curvature::NTuple{3,Float64},
    halph_s_bigger_size::NTuple{3,Float64},
    p1::NTuple{3,Float64},
    p2::NTuple{3,Float64}
)::Union{Nothing,NTuple{3,Float64}}
    # Constants
    EPS = 1e-12

    # Get radii
    Rx, Ry, Rz = halph_s_bigger_size

    # Convert to vectors
    C = collect(center_bottom_curvature)
    P1 = collect(p1)
    P2 = collect(p2)

    # Check for degenerate segment
    if norm(P2 .- P1) < EPS
        px = P1 .- C
        on_surface = abs((px[1]^2) / (Rx^2) + (px[2]^2) / (Ry^2) + (px[3]^2) / (Rz^2) - 1.0) < EPS
        return on_surface ? p1 : nothing

    end

    # Direction and offset vectors
    dV = P2 .- P1
    px = P1 .- C

    # Quadratic coefficients
    A = (dV[1]^2) / (Rx^2) + (dV[2]^2) / (Ry^2) + (dV[3]^2) / (Rz^2)
    B = 2.0 * ((px[1] * dV[1]) / (Rx^2) + (px[2] * dV[2]) / (Ry^2) + (px[3] * dV[3]) / (Rz^2))
    Cc = ((px[1]^2) / (Rx^2) + (px[2]^2) / (Ry^2) + (px[3]^2) / (Rz^2)) - 1.0

    # Check discriminant
    disc = B^2 - 4.0 * A * Cc
    if disc < 0
        return nothing
    end

    # Find intersection parameters
    sqrt_disc = sqrt(disc)
    t1 = (-B - sqrt_disc) / (2.0 * A)
    t2 = (-B + sqrt_disc) / (2.0 * A)

    # Find valid intersections (0 ≤ t ≤ 1.0)
    valid_t = filter(t -> 0.0 ≤ t ≤ 1.0, [t1, t2])

    if isempty(valid_t)
        return nothing
    end

    # Get closest intersection to p1
    t = minimum(valid_t)

    # Calculate intersection point
    P = P1 .+ t .* dV

    # Verify point lies on segment
    d1 = norm(P .- P1)
    d2 = norm(P2 .- P)
    seg_len = norm(P2 .- P1)

    if abs(d1 + d2 - seg_len) < EPS
        return (P[1], P[2], P[3])
    end

    return nothing
end

"""
    random_line_pairs_on_cylinder(
        c::NTuple{3,Float64},
        r::Float64,
        h::Float64,
        d::Float64,
        w1::Float64,
        w2::Float64
    ) -> Tuple{NTuple{3,Float64}, NTuple{3,Float64}, NTuple{3,Float64}, NTuple{3,Float64}}

Given a cylinder of radius `r`, centered at `c` with height `h`, this function
uses two random weights `w1` and `w2` in [0,1] to generate two lines through the
cylinder parallel to its long axis (the z-direction). Each line is determined by
a point on the bottom face circumference and the corresponding point on the top
face circumference.

We require that the minimal distance in the xy-plane between these two lines be
at least `d`.

Returns a tuple of four 3D points: (bottom1, top1, bottom2, top2). If the distance
requirement cannot be satisfied, the function throws an error.
"""
function random_line_pairs_on_cylinder(
    c,
    r::Float64,
    h::Float64,
    d::Float64,
    w1::Float64,
    w2::Float64
)
    @assert 0.0 ≤ w1 ≤ 1.0 "w1 must be in [0,1]"
    @assert 0.0 ≤ w2 ≤ 1.0 "w2 must be in [0,1]"

    # Cylinder center
    (cx, cy, cz) = c

    # First angle and line
    angle1 = 2π * w1
    bottom1 = (cx + r * cos(angle1), cy + r * sin(angle1), cz - h / 2)
    top1 = (cx + r * cos(angle1), cy + r * sin(angle1), cz + h / 2)

    # Second angle and line
    angle2 = 2π * w2
    bottom2 = (cx + r * cos(angle2), cy + r * sin(angle2), cz - h / 2)
    top2 = (cx + r * cos(angle2), cy + r * sin(angle2), cz + h / 2)

    # Distance between lines in XY is the distance between the two bottom points (or top).
    # They are parallel, so the distance between the lines is the same for any z.
    dx = (cx + r * cos(angle2)) - (cx + r * cos(angle1))
    dy = (cy + r * sin(angle2)) - (cy + r * sin(angle1))
    dist_xy = sqrt(dx * dx + dy * dy)

    #if balls are too close get some other random number
    if dist_xy < d
        return random_line_pairs_on_cylinder(
            c,
            r,
            h,
            d,
            w1,
            rand()
        )
    end

    return (bottom1, top1, bottom2, top2)
end

"""
    point_on_line_by_distance(p1, p2, d)

Given a line in 3D defined by two distinct points `p1` and `p2`, returns a point
on this line that is at distance `d` from `p1` and lies in the direction from `p1`
towards `p2`. If `p1 == p2`, the function throws an error.

Example:
```julia
p1 = (0.0, 0.0, 0.0)
p2 = (1.0, 1.0, 1.0)
d  = 2.0
p  = point_on_line_by_distance(p1, p2, d)

"""
function point_on_line_by_distance(p1, p2, d::Float64)
    if p1 == p2
        error("Degenerate line: p1 and p2 must not be the same.")
    end
    # Convert tuples to vectors
    v1 = collect(p1)
    v2 = collect(p2)
    direction = v2 .- v1    # p2 - p1
    dist_dir = norm(direction)

    # Normalize the direction and scale by d
    unit_dir = direction ./ dist_dir
    p_target = v1 .+ d .* unit_dir

    # Convert vector back to tuple
    return (p_target[1], p_target[2], p_target[3])
end


"""
    find_ball_positions(
        center_cylinder::NTuple{3,Float64},
        bigger_cyl_size::NTuple{3,Float64},
        ball_radius::Float64,
        center_bottom_curvature::NTuple{3,Float64},
        cylinder_bottom_curvature::Float64,
        w1::Float64,
        w2::Float64
    ) -> Tuple{NTuple{3,Float64}, NTuple{3,Float64}}

This function demonstrates how to compose calls to:
1) `random_line_pairs_on_cylinder` to obtain two parallel lines in a cylinder.
2) `point_on_segment_at_distance` to find points on these lines at a specified distance
   from `center_bottom_curvature`.
3) `point_on_line_by_distance` to shift each point further by `ball_radius`.

Steps:
1) We use `random_line_pairs_on_cylinder(center_cylinder, r, h, ball_radius, w1, w2)` where:
   - `r` = `bigger_cyl_size[1]`
   - `h` = `bigger_cyl_size[3]`
   - `ball_radius` is passed as distance requirement
   - `w1`, `w2` are random weights in [0,1]
2) Extract the two pairs of points: `(bottom1, top1, bottom2, top2)`.
3) Call `point_on_segment_at_distance(center_bottom_curvature, bottom1, top1, (ball_radius*2) + cylinder_bottom_curvature)`
   to get `rp1`.
4) Call `point_on_segment_at_distance(center_bottom_curvature, bottom2, top2, (ball_radius*2) + cylinder_bottom_curvature)`
   to get `rp2`.
5) From each result (`rp1` and `rp2`), call `point_on_line_by_distance(rp, center_bottom_curvature, ball_radius)` to
   shift the point by `ball_radius` closer to `center_bottom_curvature`.
6) Return a tuple of these two final points, ensuring they are each exactly `ball_radius` away from the bottom half-sphere
   and from the cylinder walls in the geometry context of the problem.
"""
function find_ball_positions(
    center_cylinder,
    bigger_cyl_size,
    ball_radius::Float64,
    center_bottom_curvature,
    cylinder_bottom_curvature::Float64,
    w1::Float64,
    w2::Float64,
    halph_s_bigger_size
)::Tuple{NTuple{3,Float64},NTuple{3,Float64}}

    # 1) Cylinder parameters
    r = bigger_cyl_size[1]        # radius
    h = bigger_cyl_size[3]        # height

    # 2) Get two lines from random_line_pairs_on_cylinder
    (bottom1, top1, bottom2, top2) = random_line_pairs_on_cylinder(
        center_cylinder,
        r,
        h,
        ball_radius,
        w1,
        w2
    )

    # 4) Find rp1 on the segment (bottom1, top1) at distance dist_needed from center_bottom_curvature
    elipsoid_radiuses = (halph_s_bigger_size[1] + ball_radius * 2, halph_s_bigger_size[2] + ball_radius * 2, halph_s_bigger_size[3] + ball_radius * 2)
    rp1 = ellipsoid_segment_intersection(
        center_bottom_curvature,
        elipsoid_radiuses,
        bottom1,
        top1,
    )

    # 5) Find rp2 on the segment (bottom2, top2) at distance dist_needed from center_bottom_curvature
    rp2 = ellipsoid_segment_intersection(
        center_bottom_curvature,
        elipsoid_radiuses,
        bottom2,
        top2,
    )

    # 6) Shift rp1 by ball_radius along line from rp1 to center_bottom_curvature
    rp1_ball = point_on_line_by_distance(
        rp1,
        center_bottom_curvature,
        ball_radius
    )

    # 7) Shift rp2 by ball_radius along line from rp2 to center_bottom_curvature
    rp2_ball = point_on_line_by_distance(
        rp2,
        center_bottom_curvature,
        ball_radius
    )

    return (rp1_ball, rp2_ball)
end


"""
    get_sharp_bottom(
        center_cylinder,
        bigger_cyl_size,
        halph_s_bigger_size,
        halph_s_bigger_inner,
        cylinder_wall_thickness,
        center_bottom_curvature,
        angle,
        density_inside
    ) -> Tuple

Creates half-sphere bottoms with inside and outside geometry based on provided parameters.

# Arguments
- `center_cylinder`: Center coordinates of the cylinder
- `bigger_cyl_size`: Size of the main cylinder
- `halph_s_bigger_size`: Size parameters for outer half-sphere
- `halph_s_bigger_inner`: Size parameters for inner half-sphere
- `cylinder_wall_thickness`: Wall thickness of the cylinder
- `center_bottom_curvature`: Center point of bottom curvature
- `angle`: Rotation angle
- `density_inside`: Density of the interior material

# Returns
A tuple `(ob4, ob5, center_bottom_in, center_bottom_out)` containing the geometrical objects
"""
function get_sharp_bottom(
    center_cylinder,
    bigger_cyl_size,
    halph_s_bigger_size,
    halph_s_bigger_inner,
    cylinder_wall_thickness,
    center_bottom_curvature,
    angle,
    density_inside
)
    ob4 = half_sphere_z(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] - (bigger_cyl_size[3] / 2)) - cylinder_wall_thickness)cm,
        (halph_s_bigger_size[1] - cylinder_wall_thickness)cm,
        (halph_s_bigger_size[2] - cylinder_wall_thickness)cm,
        (halph_s_bigger_size[3] + (cylinder_wall_thickness / 2))cm,
        angle, 0.0, 0, -2 + density_inside
    )

    ob5 = half_sphere_z(
        center_bottom_curvature[1]cm, center_bottom_curvature[2]cm,
        center_bottom_curvature[3]cm,
        halph_s_bigger_size[1]cm, halph_s_bigger_size[2]cm,
        (halph_s_bigger_size[3] + cylinder_wall_thickness)cm,
        angle, 0.0, 0, 1.0f0
    )

    center_bottom_in = half_sphere_z(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (center_cylinder[3] - bigger_cyl_size[3] / 2)cm,
        (halph_s_bigger_inner[1] - cylinder_wall_thickness)cm,
        (halph_s_bigger_inner[2] - cylinder_wall_thickness)cm,
        (halph_s_bigger_inner[3] - cylinder_wall_thickness)cm,
        angle, 0.0, 0, -2 + density_inside
    )

    center_bottom_out = half_sphere_z(
        center_bottom_curvature[1]cm, center_bottom_curvature[2]cm,
        center_bottom_curvature[3]cm,
        halph_s_bigger_inner[1]cm, halph_s_bigger_inner[2]cm,
        halph_s_bigger_inner[3]cm, angle, 0.0, 0, 1.0f0
    )

    return ob4, ob5, center_bottom_in, center_bottom_out
end

"""
    create_top_cut_objects(
        center_cylinder::NTuple{3,Float64},
        bigger_cyl_size::NTuple{3,Float64},
        cylinder_wall_thickness::Float64,
        len_cut::Float64,
        menisc_radius::Float64,
        menisc_radius_mult::Float64,
        menisc_cut_height::Float64,
        angle::Float64,
        x_cut_angle::Float64,
        y_cut_angle::Float64,
        density_inside::Float64
    ) -> Tuple{Any, Any, Any}

Create the top cut objects for a cylindrical phantom with a meniscus.
The function calculates the optimal position for the half-sphere to ensure
it intersects the bottom of the cut cylinder at a precise distance from the wall.

# Arguments
- `center_cylinder`: Center coordinates of the cylinder (x,y,z)
- `bigger_cyl_size`: Size of the cylinder (x,y,z)
- `cylinder_wall_thickness`: Thickness of the cylinder wall
- `len_cut`: Length of the cut section
- `menisc_radius`: Z-radius of the meniscus half-sphere
- `menisc_radius_mult`: Multiplier for meniscus X/Y radius
- `menisc_cut_height`: Height parameter for meniscus positioning
- `angle`: Rotation angle around the Z axis
- `x_cut_angle`: Tilt angle around X axis
- `y_cut_angle`: Tilt angle around Y axis
- `density_inside`: Density value for the cylinder interior

# Returns
A tuple (ob_cut, ob_cut_orig, ob_menisc_cut) containing the three geometric objects
"""
function create_top_cut_objects(
    center_cylinder,
    bigger_cyl_size,
    cylinder_wall_thickness::Float64,
    len_cut::Float64,
    menisc_radius::Float64,
    menisc_radius_mult::Float64,
    menisc_cut_height::Float64,
    angle::Float64,
    x_cut_angle::Float64,
    y_cut_angle::Float64,
    density_inside::Float64
)

    for_menisc_translation=0.9
    # Calculate the meniscus ellipsoid's X/Y radius
    half_sphere_radius_xy = bigger_cyl_size[1] * menisc_radius_mult
    half_sphere_base_z = ((center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut)))+for_menisc_translation )

    # Create the tilted cut cylinder (ob_cut)
    ob_cut = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))))cm,
        ((bigger_cyl_size[1] - cylinder_wall_thickness) * 1.5)cm,
        ((bigger_cyl_size[2] - cylinder_wall_thickness) * 1.5)cm,
        (len_cut)cm, angle, x_cut_angle,
        y_cut_angle, (-density_inside)
    )

    # Create the non-tilted cut cylinder (ob_cut_orig)
    ob_cut_orig = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 4)))+for_menisc_translation)cm,
        ((bigger_cyl_size[1] - cylinder_wall_thickness))cm,
        ((bigger_cyl_size[2] - cylinder_wall_thickness))cm,
        ((len_cut*1.5) )cm, angle, 0,
        0, (-density_inside)
    )

    # Create the meniscus half-sphere with the optimized position
    ob_menisc_cut = half_sphere_z_b(
        center_cylinder[1]cm, center_cylinder[2]cm,
        half_sphere_base_z * cm,
        (bigger_cyl_size[1]-cylinder_wall_thickness) * cm, (bigger_cyl_size[1]-cylinder_wall_thickness) * cm,
        (bigger_cyl_size[1]-cylinder_wall_thickness) * cm, angle, x_cut_angle,
        y_cut_angle, 1.0
    )

    # Return all three objects
    return ob_cut, ob_cut_orig, ob_menisc_cut
end

"""
    empty_cylinder_with_half_sphere_bottom_p(
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
        torus_height_factor,
        r_xy_small_factor, r_z_small_factor
    ) -> (objects_list, volume_dict, ob4, ob5, ob4b, ob5b, ob_cut, ob_menisc_cut, ob_cyl_mask, ob5b_mask, density_map, name_type_list)

Creates a cylindrical phantom with optional bottom/top half spheres, internal pipes, and dispensers.
Returns:
1. A list of primary objects (`ob`).
2. A dictionary containing "can_inside" volume.
3. Individual objects (ob4, ob5, ob4b, ob5b, ob_cut, ob_menisc_cut, ob5b_mask).
4. A dictionary of densities for each object.
5. A list of tuples (`name_type_list`), where each tuple is of the form
   (object_name, object_itself, object_type).

# Notes:
- Objects are created via `cylinder` or `half_sphere_z`/`half_sphere_z_b`.
- The function calculates volumes and overlaps to derive final volumes.
- The returned `name_type_list` helps track each object's creation type.
"""
function empty_cylinder_with_half_sphere_bottom_p(
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
    torus_height_factor,
    r_xy_small_factor, r_z_small_factor
)

    # Initialize objects as empty arrays
    ob4 = []
    ob5 = []
    ob4b = []
    ob5b = []
    ob_cut = []
    center_bottom_in = []
    center_bottom_out = []
    ob_menisc_cut = []
    ob5b_mask = []
    density_map = Dict{String,Float64}()
    name_type_list = []
    vol = Dict{String,Float64}()

    # Initialize name_type_list and vol
    name_type_list = []
    vol = Dict{String,Float64}()

    # Prepare dispenser cross section
    dispenser_cross_section = (pipe_cross_section[1] + dispenser_cross_section[1],
        pipe_cross_section[2] + dispenser_cross_section[2])

    # Half sphere sizes
    halph_s_bigger_size = (bigger_cyl_size[1] - cylinder_wall_thickness,
        bigger_cyl_size[2] - cylinder_wall_thickness,
        cylinder_bottom_curvature)

    halph_s_bigger_inner = ((bigger_cyl_size[1] - cylinder_wall_thickness * 1.5) * rel_size_bottom_curvature,
        (bigger_cyl_size[2] - cylinder_wall_thickness * 1.5) * rel_size_bottom_curvature,
        cylinder_bottom_curvature_b)

    halph_s_bigger_size_top = (bigger_cyl_size[1] - cylinder_wall_thickness * 1.5,
        bigger_cyl_size[2] - cylinder_wall_thickness * 1.5,
        cylinder_top_curvature)

    # Liquid phases (dual_phase_percentage)
    if (dual_phase_percentage == 1.0)
        ob2_a = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (center_cylinder[3])cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[3])cm, angle, 0, 0, (-(1 - density_inside)) / 2)
        ob2_b = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (center_cylinder[3])cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[3])cm, angle, 0, 0, (-(1 - density_inside)) / 2)

        #just dummies
        straight_middle_cut = ob2_a
        oblique_middle_cut = ob2_b

    else
        len1 = bigger_cyl_size[3] * dual_phase_percentage
        len2 = bigger_cyl_size[3] - len1
        L = len1 + len2
        z_center = center_cylinder[3]
        z1 = z_center - (L / 2) + (len1 / 2)
        z2 = z_center + (L / 2) - (len2 / 2)

        ob2_a = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (z1)cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            (len1)cm, angle, 0, 0, -(1 - density_inside))

        ob2_b = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (z2)cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            (len2)cm, angle, 0, 0, -(1 - density_inside_b))

        #we want to define two cylinders that would be around the division between fluids
        # one will be straight as the long axis of can and second will be oblique as top fluid is

        len_cyls_for_cut = L * 0.4
        z_base = center_cylinder[3] - (bigger_cyl_size[3] / 2)
        z_base = z_base + len1
        straight_middle_cut = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (z_base)cm,#+(len_cyls_for_cut/2)
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            (len_cyls_for_cut * 0.8)cm, angle, 0, 0, 1)

        oblique_middle_cut = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (z_base + len_cyls_for_cut / 2)cm,#+(len_cyls_for_cut/2)
            (bigger_cyl_size[1] * 1.1)cm,
            (bigger_cyl_size[2] * 1.1)cm,
            (len_cyls_for_cut)cm, angle, x_cut_angle, y_cut_angle, 1)
    end

    # Add initial objects to name_type_list
    push!(name_type_list, ("ob2_a", ob2_a, "cylinder"))
    push!(name_type_list, ("ob2_b", ob2_b, "cylinder"))

    # A cut cylinder on top
    ob_cut, ob_cut_orig, ob_menisc_cut = create_top_cut_objects(
        center_cylinder,
        bigger_cyl_size,
        cylinder_wall_thickness,
        len_cut,
        menisc_radius,
        menisc_radius_mult,
        menisc_cut_height,
        angle,
        x_cut_angle,
        y_cut_angle,
        density_inside
    )

    # Outer cylinder
    ob3 = cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (center_cylinder[3] - (cylinder_wall_thickness / 2))cm,
        bigger_cyl_size[1]cm, bigger_cyl_size[2]cm, (bigger_cyl_size[3] + (cylinder_wall_thickness / 1.5))cm,
        angle, 0, 0, 1.0f0)

    # Add more objects to name_type_list
    push!(name_type_list, ("ob3", ob3, "cylinder"))
    push!(name_type_list, ("ob_cut", ob_cut, "cylinder"))
    push!(name_type_list, ("ob_cut_orig", ob_cut_orig, "cylinder"))
    push!(name_type_list, ("ob_menisc_cut", ob_menisc_cut, "half_sphere"))

    center_bottom_curvature = (center_cylinder[1], center_cylinder[2],
        ((center_cylinder[3] - bigger_cyl_size[3] / 2)))
    bottom_cyl_start = ((center_cylinder[3] - bigger_cyl_size[3] / 2)) + cylinder_wall_thickness


    if rounded_bottom
        # Calculate parameters for get_rounded_bottom
        r_xy_small = bigger_cyl_size[1] * r_xy_small_factor #0.21
        r_z_small = bigger_cyl_size[3] * r_z_small_factor#0.02
        r_xy_torus = bigger_cyl_size[1]        # Same as cylinder radius
        base_cylinder_radius = bigger_cyl_size[1]
        center_torus = [
            center_cylinder[1],
            center_cylinder[2],
            center_cylinder[3] - bigger_cyl_size[3] / 2
        ]

        main_torus, inner_torus, outer_sphere, inner_sphere, ball1, ball2, main_torus_volume_analytical, inner_torus_volume_analytical, outer_sphere_volume, ball1_volume, ball2_volume = get_rounded_bottom(
            spacing,                # spacing tuple
            dims,                  # dimensions tuple
            ig,                    # ImageGeom
            cylinder_wall_thickness,
            r_xy_small,           # small radius xy
            r_z_small,            # small radius z
            r_xy_torus,           # torus radius
            center_torus,         # center position vector
            base_cylinder_radius,  # base cylinder radius
            curvature,            # curvature parameter
            tube_ball_fraction,   # tube ball fraction
            rel_dist             # relative distance
        )

        # Add objects to name_type_list with correct types
        push!(name_type_list, ("main_torus", main_torus, "array"))
        push!(name_type_list, ("inner_torus", inner_torus, "array"))
        push!(name_type_list, ("outer_sphere", outer_sphere, "array"))
        push!(name_type_list, ("inner_sphere", inner_sphere, "array"))
        push!(name_type_list, ("ball1", ball1, "half_sphere"))
        push!(name_type_list, ("ball2", ball2, "half_sphere"))

    else
        ob4, ob5, center_bottom_in, center_bottom_out = get_sharp_bottom(
            center_cylinder,
            bigger_cyl_size,
            halph_s_bigger_size,
            halph_s_bigger_inner,
            cylinder_wall_thickness,
            center_bottom_curvature,
            angle,
            density_inside
        )

        # Add sharp bottom objects to name_type_list
        push!(name_type_list, ("ob4", ob4, "half_sphere"))
        push!(name_type_list, ("ob5", ob5, "half_sphere"))
        push!(name_type_list, ("center_bottom_in", center_bottom_in, "half_sphere"))
        push!(name_type_list, ("center_bottom_out", center_bottom_out, "half_sphere"))
    end

    # Top half sphere inside
    ob4b = half_sphere_z(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] + bigger_cyl_size[3] / 2))cm,
        (halph_s_bigger_size_top[1] - cylinder_wall_thickness * 0.3)cm,
        (halph_s_bigger_size_top[2] - cylinder_wall_thickness * 0.3)cm,
        (halph_s_bigger_size_top[3])cm,
        angle, 0.0, 0, -(1)
    )

    # Top half sphere outside
    ob5b = half_sphere_z(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] + bigger_cyl_size[3] / 2))cm,
        (halph_s_bigger_size_top[1] + cylinder_wall_thickness * 1.4)cm,
        (halph_s_bigger_size_top[2] + cylinder_wall_thickness * 1.4)cm,
        (halph_s_bigger_size_top[3] + (cylinder_wall_thickness))cm,
        angle, 0.0, 0, 1.0f0
    )

    ob5b_mask = half_sphere_z(
        center_cylinder[1]cm, center_cylinder[2]cm,
        ((center_cylinder[3] + bigger_cyl_size[3] / 2) - cylinder_wall_thickness * 0.05)cm,
        ((halph_s_bigger_size_top[1] + cylinder_wall_thickness * 0.95))cm,
        ((halph_s_bigger_size_top[2] + cylinder_wall_thickness * 0.95))cm,
        ((halph_s_bigger_size_top[3] + cylinder_wall_thickness))cm,
        angle, 0.0, 0, 1.0f0
    )

    # Add top half sphere objects to name_type_list
    push!(name_type_list, ("ob4b", ob4b, "half_sphere"))
    push!(name_type_list, ("ob5b", ob5b, "half_sphere"))
    push!(name_type_list, ("ob5b_mask", ob5b_mask, "half_sphere"))

    # Pipe geometry
    pipe_start = center_cylinder[3] + (bigger_cyl_size[3] / 2) + cylinder_top_curvature - (pipe_len / 2)
    pipe = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm, (pipe_start)cm,
        pipe_cross_section[1]cm, pipe_cross_section[2]cm, (pipe_len)cm,
        angle, 0, 0, pipe_density
    )
    pipe_in = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm, (pipe_start)cm,
        pipe_cross_section[1]cm / 2, pipe_cross_section[2]cm / 2, (pipe_len)cm,
        angle, 0, 0, -pipe_density
    )

    dispenser_start = center_cylinder[3] + (bigger_cyl_size[3] / 2) + cylinder_top_curvature
    plastic_dispenser = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm, (dispenser_start)cm,
        dispenser_cross_section[1]cm, dispenser_cross_section[2]cm,
        (dispenser_len)cm, angle, 0, 0, dispenser_density
    )

    # Add pipe and dispenser objects to name_type_list
    push!(name_type_list, ("pipe", pipe, "cylinder"))
    push!(name_type_list, ("pipe_in", pipe_in, "cylinder"))
    push!(name_type_list, ("plastic_dispenser", plastic_dispenser, "cylinder"))
    if (!rounded_bottom)
        ball_radius = pipe_cross_section[1]
        bc1, bc2 = find_ball_positions(
            center_cylinder,
            (bigger_cyl_size[1] - cylinder_wall_thickness, bigger_cyl_size[2] - cylinder_wall_thickness, bigger_cyl_size[3] - cylinder_wall_thickness),
            ball_radius,
            center_bottom_curvature,
            cylinder_bottom_curvature,
            rand(),
            rand(),
            halph_s_bigger_size
        )

        ball1 = ellipsoid(bc1[1]cm, bc1[2]cm, bc1[3]cm,
            (ball_radius)cm, (ball_radius)cm,
            (ball_radius)cm, angle, 0, 0, pipe_density - density_inside)

        ball2 = ellipsoid(bc2[1]cm, bc2[2]cm, bc2[3]cm,
            (ball_radius)cm, (ball_radius)cm,
            (ball_radius)cm, angle, 0, 0, pipe_density - density_inside)

        # Add ball objects to name_type_list
        push!(name_type_list, ("ball1", ball1, "half_sphere"))
        push!(name_type_list, ("ball2", ball2, "half_sphere"))
    end

    pipe_center_z = pipe_start
    pipe_beg_z = pipe_start - (pipe_len / 2)
    pipe_end_z = pipe_start + (pipe_len / 2)
    ob2_beg_z = center_cylinder[3] - (bigger_cyl_size[3] / 2)
    ob2_end_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)

    overlay = 0
    if pipe_beg_z <= ob2_end_z && pipe_end_z >= ob2_beg_z
        overlay = abs(min(pipe_end_z, ob2_end_z) - max(pipe_beg_z, ob2_beg_z))
    end

    whole_pipe = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (pipe_center_z - overlay / 2)cm,
        pipe_cross_section[1]cm, pipe_cross_section[2]cm,
        (overlay)cm, angle, 0, 0, pipe_density
    )
    whole_pipe_volume = IP.volume(whole_pipe)

    inside_of_pipe_volume = IP.volume(
        cylinder(center_cylinder[1]cm, center_cylinder[2]cm, (pipe_end_z)cm,
            pipe_cross_section[1]cm / 2, pipe_cross_section[2]cm / 2,
            (overlay)cm, angle, 0, 0, -pipe_density)
    )

       res_vol = 0.0#vol1 - vol2 - whole_pipe_volume + inside_of_pipe_volume - IP.volume(ob_cut)

    # Create a dictionary of densities
    density_map = Dict{String,Float64}()

    density_inside = density_inside
    density_inside_b = density_inside_b
    density_wall = 1.0
    #rounded bottom is for aluminium can
    if (rounded_bottom)
        density_wall = 0.8
    end
    density_inside = min(density_inside, density_wall - 0.2)
    density_inside_b = min(density_inside_b, density_wall - 0.2)

    if dual_phase_percentage == 1.0
        density_map["ob2_a"] = density_inside / 2#(-(1 - density_inside)) / 2
        density_map["ob2_b"] = density_inside / 2#(-(1 - density_inside)) / 2
    else
        density_map["ob2_a"] = density_inside#-(1 - density_inside)
        density_map["ob2_b"] = density_inside_b#-(1 - density_inside_b)
    end
    density_map["ob3"] = density_wall
    density_map["pipe"] = pipe_density
    density_map["pipe_in"] = -pipe_density
    density_map["plastic_dispenser"] = dispenser_density * 1.3
    density_map["ob_cut"] = -density_inside
    density_map["ob4"] = density_wall
    density_map["ob5"] = density_wall
    density_map["ob4b"] = -1.0
    density_map["ob5b"] = density_wall
    density_map["ob5b_mask"] = density_wall
    density_map["ob_menisc_cut"] = density_wall
    # density_map["ob_cyl_mask"]  = 1.0
    density_map["ob_cut_orig"] = density_wall
    density_map["straight_middle_cut"] = density_wall
    density_map["oblique_middle_cut"] = density_wall
    density_map["ball1"] = density_wall - density_inside
    density_map["ball2"] = density_wall - density_inside
    density_map["center_bottom_in"] = density_wall
    density_map["center_bottom_out"] = density_wall

    # Calculate the cylinder representing fluid below the cut (calc_cyl)
    # Get the top of the cylinder
    cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)

    # Calculate the lowest point of the cut plane intersection with ob3
    cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))

    # Due to tilt, calculate the maximum drop
    radius = bigger_cyl_size[1] - cylinder_wall_thickness

    max_drop_x = radius * sin(abs(x_cut_angle))
    max_drop_y = radius * sin(abs(y_cut_angle))
    max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

    # Calculate the z-coordinate of the lowest point of the cut plane
    lowest_cut_z = cut_center_z - max_drop - len_cut / 2

    # Distance from cylinder top to lowest cut point
    distance_from_top = cylinder_top_z - lowest_cut_z

    # Height of the cylindrical portion below the cut
    cylinder_height_below_cut = bigger_cyl_size[3] - distance_from_top

    # Create the analytical cylinder representing fluid below cut
    calc_cyl = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2)cm,
        (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
        (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
        cylinder_height_below_cut * cm,
        angle, 0, 0, 1.0f0
    )

    top_fluid=(center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut)
    bottom_pipe=pipe_start-pipe_len/2
    # bottom_pipe=top_fluid-1.0
    len_pipe_in_fluid=top_fluid-bottom_pipe

    inside_of_pipe_in_fluid = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (bottom_pipe + len_pipe_in_fluid/2)cm,
        (pipe_cross_section[1]/2)cm, (pipe_cross_section[2]/2)cm,
        (len_pipe_in_fluid)cm, angle, 0, 0, pipe_density
    )

    outside_of_pipe_in_fluid = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (bottom_pipe + len_pipe_in_fluid/2)cm,
        (pipe_cross_section[1])cm, (pipe_cross_section[2])cm,
        (len_pipe_in_fluid)cm, angle, 0, 0, pipe_density
    )

    # Add calc_cyl to density map and name_type_list
    density_map["calc_cyl"] = 1.0
    density_map["inside_of_pipe_in_fluid"] = 1.0
    density_map["outside_of_pipe_in_fluid"] = 1.0

    push!(name_type_list, ("calc_cyl", calc_cyl, "cylinder"))

    # ----------------------------------------------------------
    # Create a list of (name, object, type)
    # ----------------------------------------------------------

    # Push each entry into the name_type_list
    push!(name_type_list, ("ob2_a", ob2_a, "cylinder"))
    push!(name_type_list, ("ob2_b", ob2_b, "cylinder"))
    push!(name_type_list, ("ob3", ob3, "cylinder"))
    push!(name_type_list, ("pipe", pipe, "cylinder"))
    push!(name_type_list, ("pipe_in", pipe_in, "cylinder"))
    push!(name_type_list, ("plastic_dispenser", plastic_dispenser, "cylinder"))
    push!(name_type_list, ("ob_cut", ob_cut, "cylinder"))
    push!(name_type_list, ("ob4", ob4, "half_sphere"))
    push!(name_type_list, ("ob5", ob5, "half_sphere"))
    push!(name_type_list, ("ob4b", ob4b, "half_sphere"))
    push!(name_type_list, ("ob5b", ob5b, "half_sphere"))
    # push!(name_type_list, ("ob5b_mask", ob5b_mask, "half_sphere"))  # Uncomment if needed
    push!(name_type_list, ("ob_menisc_cut", ob_menisc_cut, "half_sphere"))
    # push!(name_type_list, ("ob_cyl_mask", ob_cyl_mask, "cylinder"))  # Uncomment if needed
    push!(name_type_list, ("ob_cut_orig", ob_cut_orig, "cylinder"))
    # push!(name_type_list, ("ball1", ball1, "half_sphere"))
    # push!(name_type_list, ("ball2", ball2, "half_sphere"))
    push!(name_type_list, ("straight_middle_cut", straight_middle_cut, "cylinder"))
    push!(name_type_list, ("oblique_middle_cut", oblique_middle_cut, "cylinder"))
    push!(name_type_list, ("center_bottom_in", center_bottom_in, "half_sphere"))
    push!(name_type_list, ("center_bottom_out", center_bottom_out, "half_sphere"))
    push!(name_type_list, ("inside_of_pipe_in_fluid", inside_of_pipe_in_fluid, "cylinder"))
    push!(name_type_list, ("outside_of_pipe_in_fluid", outside_of_pipe_in_fluid, "cylinder"))

    return vol, density_map, name_type_list
end

"""
    get_geometric_onjects() -> Tuple

Get default geometric objects for phantom simulation.

# Returns
A tuple containing the phantom container objects and their properties
"""
function get_geometric_onjects()
    # return ionic_chamber()
    return empty_cylinder_with_half_sphere_bottom()

end


"""
    distance2d(x1, y1, x2, y2)

Return the Euclidean distance between two points in 2D.
"""
function distance2d(x1, y1, x2, y2)::Float64
    return sqrt((x2 - x1)^2 + (y1 - y2)^2)
end

"""
    set_circle!(shape, cx, cy, z, radius)

For a given 3D Boolean array `shape`, set to true all points (x, y, z)
within distance `radius` of center `(cx, cy)` on slice `z`.
"""
function set_circle!(shape::Array{Bool,3},
    cx::Float64, cy::Float64, z::Int,
    radius::Float64)
    nx, ny, nz = size(shape)
    # Bounding box around the circle to limit iteration
    x_min = max(floor(Int, cx - radius), 1)
    x_max = min(ceil(Int, cx + radius), nx)
    y_min = max(floor(Int, cy - radius), 1)
    y_max = min(ceil(Int, cy + radius), ny)

    for xx in x_min:x_max
        for yy in y_min:y_max
            if distance2d(xx, yy, cx, cy) <= radius
                shape[xx, yy, z] = true
            end
        end
    end
end

"""
    move_center(cx, cy, center_circle_radius, theta_step)

Given a current center `(cx, cy)`, return `(nx, ny)` which is
moved by angle `theta_step` on a circle of radius `center_circle_radius`
around the *origin* `(0,0)`.
"""
function move_center(cx::Float64, cy::Float64,
    center_circle_radius::Float64,
    theta_step::Float64)
    # Convert current (cx, cy) to polar w.r.t. (0,0)
    r = sqrt(cx^2 + cy^2)
    θ = atan(cy, cx)

    # We want to keep the radius to the origin at `center_circle_radius`
    # and increment the angle by `theta_step`.
    new_angle = θ + theta_step
    nx = center_circle_radius * cos(new_angle)
    ny = center_circle_radius * sin(new_angle)
    return nx, ny
end

"""
    spiral_screw(
        beginning_center::Tuple{Float64,Float64,Float64},
        z_length::Int,
        theta::Float64,
        main_radius::Float64,
        center_circle_radius::Float64,
        dims::Tuple{Int,Int,Int}
    )::Array{Bool,3}

Create a 3D Boolean array `shape` of size `dims` that resembles a spiral "screw."
Works layer by layer in the z-direction. For each layer z, it:

1. Places the "draw center" on a circle of radius `center_circle_radius`
   around the main axis (starting from `beginning_center`).
2. Draws a circle of radius `main_radius` in that layer.
3. Moves the center by `theta` on the circle for the next layer.
4. Continues until z_length layers have been processed.

Returns the 3D Boolean mask with the spiral shape.
"""
function spiral_screw(
    beginning_center::Tuple{Float64,Float64,Float64},
    z_length,
    theta::Float64,
    main_radius::Float64,
    center_circle_radius::Float64,
    dims::Tuple{Int,Int,Int}
)::Array{Bool,3}

    shape = fill(false, dims)
    # Unpack beginning_center
    (cx0, cy0, cz0) = beginning_center
    # We assume (cx0, cy0) is the center in x-y plane at top; then we pick an initial offset
    # so that "new_center" is center_circle_radius away from (cx0, cy0).
    # For simplicity, shift to an origin-based approach:
    # We'll treat (cx0, cy0) as the axis center in the x-y plane,
    # so the "moving center" will revolve around that axis.

    # Initial angle or offset on the center_circle
    # e.g. place it at angle = 0 relative to (cx0, cy0)
    # That means initial new_center in x-y is (cx0 + center_circle_radius, cy0).
    new_cx = cx0 + center_circle_radius
    new_cy = cy0

    for layer_idx in 0:(z_length-1)
        z_pos = round(Int, cz0 + layer_idx)
        if z_pos < 1 || z_pos > dims[3]
            # If we've gone out of the bounding box in z, stop
            break
        end
        # Draw circle for this layer
        set_circle!(shape, new_cx, new_cy, z_pos, main_radius)
        # Move the center by "theta" around (cx0, cy0)
        # Translate new_cx/new_cy to origin, apply angle step, translate back
        rel_x = new_cx - cx0
        rel_y = new_cy - cy0
        nx, ny = move_center(rel_x, rel_y, center_circle_radius, theta)
        new_cx = nx + cx0
        new_cy = ny + cy0
    end

    return shape
end

"""
    get_beg_voxel(spacing, center_cylinder, cyl_size, axis, dims) -> Float64

Calculate the beginning voxel position for a cylinder along the specified axis.

# Arguments
- `spacing`: Voxel spacing tuple
- `center_cylinder`: Center coordinates of the cylinder
- `cyl_size`: Size of the cylinder
- `axis`: Axis index (1=x, 2=y, 3=z)
- `dims`: Dimensions of the volume

# Returns
The starting voxel position for the cylinder along the specified axis
"""
function get_beg_voxel(spacing, center_cylinder, cyl_size, axis, dims)
    return ((center_cylinder[axis] - cyl_size[axis] / 2) / spacing[axis])
end