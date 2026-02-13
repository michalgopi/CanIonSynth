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
const IP = ImagePhantoms
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

wandb.init(project="synth")
isinteractive() ? jim(:prompt, true) : prompt(:draw)
include.([
    "../in_docker_organized/get_geometry_main.jl",
    "../in_docker_organized/geometry_utils.jl"
])

# Test Configuration
dims = (32, 32, 32)
spacing = (0.1, 0.1, 0.1)

# Helper function mocks or copies from main_create_phantom_ionic_chamber.jl
# ... (Include necessary helper functions like create_cylinder, adjust_densities!, etc. if they are not in included files)
# Wait, main_create_phantom_ionic_chamber.jl defines these functions locally.
# I need to copy them or include the file.
# Including the file `main_create_phantom_ionic_chamber.jl` runs the code at the bottom.
# So I must redefine the necessary functions here.

# Copying essential functions from main_create_phantom_ionic_chamber.jl
# (Ideally, these should be in a module)

function cylinder_interval(c_z::Float64, l_z::Float64)
    return (c_z - l_z / 2, c_z + l_z / 2)
end

function overlaps(c1::Float64, l1::Float64, c2::Float64, l2::Float64)
    s1, e1 = cylinder_interval(c1, l1)
    s2, e2 = cylinder_interval(c2, l2)
    return max(s1, s2) < min(e1, e2)
end

function compute_adjusted_density(new_desired::Float64, new_z::Float64, new_lz::Float64, existing::Vector{NamedTuple})
    total_overlap_density = 0.0
    for cyl in existing
        if overlaps(new_z, new_lz, cyl.c_z, cyl.l_z)
            total_overlap_density += cyl.d
        end
    end
    return new_desired - total_overlap_density
end

function create_cylinder(base_x::Float64, base_y::Float64, base_bottom_z::Float64,
    spec, angle::Float64, name)
    d_desired = spec["d"]
    radius = spec["radius"]
    l_z = spec["l_z"]
    offset = spec["offset"]
    is_halph_sphere = haskey(spec, "is_halph_sphere") ? spec["is_halph_sphere"] : false
    is_ellipsoid = haskey(spec, "is_ellipsoid") ? spec["is_ellipsoid"] : false

    c_z = base_bottom_z + offset + (l_z / 2)
    if offset < 0
        c_z = 0.0
    end
    if haskey(spec, "theta")
        if spec["theta"] > 0.0
            return (name=name, is_spiral=true, is_halph_sphere=is_halph_sphere, is_ellipsoid=is_ellipsoid, c_x=base_x, c_y=base_y, c_z=c_z,
                l_x=radius, l_y=radius, l_z=l_z, angle=angle, d=d_desired,
                theta=spec["theta"], center_circle_radius=spec["center_circle_radius"], orig_d=d_desired)
        end
    end

    return (name=name, is_spiral=false, is_halph_sphere=is_halph_sphere, is_ellipsoid=is_ellipsoid, c_x=base_x, c_y=base_y, c_z=c_z,
        l_x=radius, l_y=radius, l_z=l_z, angle=angle, d=d_desired)
end

function sort_cylinders_by_radius!(cylinders::Vector{NamedTuple})
    base = cylinders[1]
    rest = sort(cylinders[2:end], by=x -> -x.l_x)
    return vcat([base], rest)
end

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

function compute_effective_density(child::NamedTuple, parent_final)
    return child.d - parent_final
end

function adjust_densities!(cylinders::Vector{NamedTuple})::Vector{NamedTuple}
    sorted_cyls = sort_cylinders_by_radius!(cylinders)
    n = length(sorted_cyls)
    adjusted = Vector{NamedTuple}(undef, n)
    adjusted[1] = sorted_cyls[1]

    for i in 2:n
        child = sorted_cyls[i]
        pinx = parent_index(child, sorted_cyls)
        if pinx == -1
            effective = child.d
        else
            p_final = sorted_cyls[pinx].d
            effective = compute_effective_density(child, p_final)
        end
        effective = effective + 0.000001
        adjusted[i] = (@set child.d = effective)
    end

    return adjusted
end

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
    end

    res = Float32.(ob_el) * new_cyl.d
    return res,voll
end

function build_image(base_params, specs, angle::Float64, spacing, dims)
    ig = ImageGeom(dims=(dims[1],dims[2],dims[3]), deltas=(spacing[1] * cm, spacing[2] * cm, spacing[3] * cm))

    base_x = 0.0
    base_y = 0.0
    base_cz = base_params["c_z"]
    base_lz = base_params["l_z"]
    base_radius = base_params["radius"]
    base_d = base_params["d"]

    base_bottom_z = base_cz - (base_lz / 2)

    cylinders = Vector{NamedTuple}(undef, 0)
    base_cyl = (name="base_cyl", is_halph_sphere=false, is_spiral=false, is_ellipsoid=false, c_x=base_x, c_y=base_y, c_z=base_cz, l_x=base_radius, l_y=base_radius,
        l_z=base_lz, angle=angle, d=base_d, orig_d=base_d)
    push!(cylinders, base_cyl)

    for (name, spec) in specs
        new_cyl = create_cylinder(base_x, base_y, base_bottom_z, spec, angle, name)
        push!(cylinders, new_cyl)
    end

    cylinders = adjust_densities!(cylinders)
    cylinders_prim = deepcopy(cylinders)

    cylinder_densities = Dict{String,Float64}()
    for cyl in cylinders
        cylinder_densities[cyl.name] = cyl.d
    end
    res = map(new_cyl -> get_cyl_obj(new_cyl, ig, spacing), cylinders)

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

function create_ionic_chamber_phantom(params)
    graphite_density = params["graphite_density"]
    copper_density = params["copper_density"]
    aluminium_density = params["aluminium_density"]
    insulation_density = params["insulation_density"]

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

    spacing = params["spacing"]
    dims = params["dims"]

    total_radius = air_thickness + (copper_radius * 1.2) + main_graphite_thickness
    base_params = Dict("c_x" => main_radius, "c_y" => main_radius, "c_z" => (-1) * (((dims[3] / 4) * spacing[3]) * 0.9), "radius" => total_radius, "l_z" => total_len, "d" => 1.0)

    specs = Dict(
        "copper_el" => Dict("d" => copper_density, "radius" => copper_radius, "l_z" => base_len, "offset" => 0.0),
        "inner_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius, "l_z" => base_len, "offset" => 0.0),
        "inner_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness, "l_z" => base_len, "offset" => 0.0),
        "outer_insulator" => Dict("d" => insulation_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness, "l_z" => base_len, "offset" => 0.0),
        "outer_aluminium" => Dict("d" => aluminium_density, "radius" => inner_insluation_thickness + copper_radius + aluminium_inner_thicness + outer_insluation_thickness + aluminium_outer_thicness, "l_z" => base_len, "offset" => 0.0, "theta" => -1.0, "center_circle_radius" => 0.0),
        "base_graphite" => Dict("d" => graphite_density, "radius" => total_radius - side_cut_size, "l_z" => base_len, "offset" => 0.0),
        "cut" => Dict("d" => 0.0, "radius" => total_radius - 0.00001, "l_z" => base_len, "offset" => 0.0),
        "bottom_main_graphite" => Dict("d" => graphite_density, "radius" => total_radius - 0.00001, "l_z" => main_graphite_thickness, "offset" => base_len),
        "graphite_electrode" => Dict("d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => total_len - base_len - main_graphite_thickness, "offset" => base_len + main_graphite_thickness),
        "air" => Dict("d" => 00001, "radius" => total_radius - (main_graphite_thickness * 1.05), "l_z" => total_len - base_len - main_graphite_thickness, "offset" => base_len + main_graphite_thickness),
        "top_inner" => Dict("is_halph_sphere" => true, "d" => 0.0, "radius" => (total_radius - (main_graphite_thickness * 0.95)), "l_z" => (total_radius), "offset" => total_len - main_graphite_thickness),
        "top_outer" => Dict("is_halph_sphere" => true, "d" => graphite_density, "radius" => total_radius, "l_z" => total_radius, "offset" => total_len),
        "top_electrode" => Dict("is_halph_sphere" => true, "d" => graphite_density, "radius" => graphite_electrode_radius, "l_z" => graphite_electrode_radius, "offset" => (base_len) + total_len - base_len),
    )

    angle = 0.0
    cylinders_prim, named_res, cylinder_densities,vol_res = build_image(base_params, specs, angle, spacing, dims)

    res = zeros(size(named_res["base_cyl"]))
    keyss = keys(named_res)

    for key in keyss
        if (key != "top_outer") && (key != "top_inner")
            res .+= named_res[key]
        end
    end

    res = res + ((named_res["top_outer"] + named_res["top_inner"]) .* Float64.(named_res["base_cyl"] .== 0))
    res[(named_res["air"].!=0)] .= 0.0
    res = res + named_res["graphite_electrode"]

    return res
end

function run_test()
    println("Starting Ionic Chamber Phantom Test...")

    params = Dict(
        "graphite_density" => 1.0,
        "copper_density" => 2.67,
        "aluminium_density" => 1.2,
        "insulation_density" => 0.13,
        "square_top" => false,
        "ball_like" => false,
        "lolipop_like" => false,
        "rounded_top" => true,
        "add_graphite_in_copper" => false,
        "add_spiral" => false,
        "elongate_copper" => false,
        "total_len" => 2.65,
        "base_len" => 0.8,
        "main_radius" => 0.625,
        "main_graphite_thickness" => 0.25,
        "side_cut_size" => 0.2,
        "air_thickness" => 0.265,
        "graphite_electrode_radius" => 0.125,
        "copper_radius" => 0.075,
        "inner_insluation_thickness" => 0.03,
        "aluminium_inner_thicness" => 0.03,
        "outer_insluation_thickness" => 0.04,
        "aluminium_outer_thicness" => 0.2,
        "theta_aluminium_outer" => -1.0,
        "graphite_electrode_connection_radius" => 0.15,
        "graphite_electrode_thickness" => 0.1,
        "spacing" => spacing,
        "dims" => dims
    )

    res = create_ionic_chamber_phantom(params)

    println("Phantom generated with size: ", size(res))

    # Save output
    output_dir = "tests/output_ionic"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    output_path = joinpath(output_dir, "test_ionic.nii.gz")

    img = sitk.GetImageFromArray(Float32.(res))
    img.SetSpacing((params["spacing"][1] * 10, params["spacing"][2] * 10, params["spacing"][3] * 10))
    sitk.WriteImage(img, output_path)

    println("Saved NIfTI to ", output_path)

    if isfile(output_path)
        println("Test PASSED: Output file exists.")
    else
        println("Test FAILED: Output file not found.")
    end
end

run_test()
