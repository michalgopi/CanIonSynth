#remove ImagePhantoms
#add https://github.com/jakubMitura14/ImagePhantoms.jl.git
using Pkg
# using Wandb

Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git")
# Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git")
# Pkg.add(url="https://github.com/jakubMitura14/Sinograms.jl.git")


#!!! we need to use a fork of image phantoms  and sinograms  
using Sinograms: SinoPar, rays, plan_fbp, fbp, fbp_sino_filter
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using Unitful: mm
using MIRTjim: jim, prompt
using ImagePhantoms: Object, phantom, radon, spectrum
using ImagePhantoms: Cylinder, cylinder
import ImagePhantoms as IP
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Unitful: mm, unit, °
using Plots: plot, plot!, scatter!, default
using Plots # gif @animate
using PyCall
import Unitful
using Unitful
using Plots: plot, gui
using Unitful: cm
using Sinograms: CtFanArc, CtFanFlat # CtPar
using Sinograms: rays, plan_fbp, Window, Hamming, fdk, ct_geom_plot3
using ImageGeoms: ImageGeom, MaskCircle, fovs
using ImagePhantoms: ellipsoid_parameters, ellipsoid
using ImagePhantoms: radon, phantom
using MIRTjim: jim, prompt
using ImagePhantoms: Ellipsoid, ellipsoid
using ImagePhantoms: Cuboid, cuboid
import ImagePhantoms
using Revise, Statistics
using Sinograms: project_bdd, backproject_bdd
using UUIDs
using JSON
using Dates, Logging

wandb_api_key="5a6c2e2caccb19a9c3cbfde388b61c7104eab632"

# Initialize the project
# lg = WandbLogger(; project = "synth", name = nothing)
# using RadonKA
sitk = pyimport_conda("SimpleITK", "simpleitk")
np = pyimport_conda("numpy", "numpy")
wandb = pyimport_conda("wandb", "wandb")



uuid_main = UUIDs.uuid4()

wandb.init(
    # set the wandb project where this run will be logged
    project="synth",

)



isinteractive() ? jim(:prompt, true) : prompt(:draw);
# includet("/workspaces/synthethic_tomo/.devcontainer/get_geometry_main.jl") #TODO remove
# includet("/workspaces/synthethic_tomo/.devcontainer/get_geometry_main.jl") #TODO remove

includet("/root/.devcontainer/get_geometry_main.jl") 
includet("/root/.devcontainer/CtFanArc_params.jl") 

#controls how irregular the cylinders are
ImagePhantoms.IRREGULARITY = 0.3

#get size of the generated image
spacing = (0.1, 0.1, 0.1)
# spacing = (0.3,0.3,0.3)
# dims=(66,64,50)
low_res = false
dims = (800, 800, 500)
# dims = (300, 300, 300)#works

if (low_res)
    dims = (200, 200, 200)

end
# dims=(2000,2000,2000)
ig = ImageGeom(dims=dims, deltas=(spacing[1]cm, spacing[2]cm, spacing[3]cm))
# global_logger(lg)


wandb.log(Dict("setting up wandb process"=>uuid_main ))


function random_can(ig, image_size_cm)
    # Generate random parameters within the specified constraints
    center_cylinder = (
        0.0,
        0.0,
        0.0
    )
    # center_cylinder = (
    #     rand(-image_size_cm[1]/2:image_size_cm[1]/2),
    #     rand(-image_size_cm[2]/2:image_size_cm[2]/2),
    #     rand(-image_size_cm[3]/2:image_size_cm[3]/2)
    # )


    bigger_cyl_size = [
        max(1.0 + 0.25 * 2, rand(0.05:0.01:0.5) * image_size_cm[1] / 2),
        max(1.0 + 0.25 * 2, rand(0.05:0.01:0.5) * image_size_cm[2] / 2),
        1.0
    ]
    bigger_cyl_size[3] = rand(1.0:0.01:8.0) * maximum(collect(bigger_cyl_size))
    bigger_cyl_size[3] = min(bigger_cyl_size[3], image_size_cm[3] * 0.6)


    bigger_cyl_size[2] = bigger_cyl_size[1]

    len_cut = rand(0.1:0.3) * bigger_cyl_size[3]

    menisc_radius = 1.15 * len_cut

    cylinder_wall_thickness = rand(0.02:0.1) * minimum(bigger_cyl_size)

    cylinder_bottom_curvature = rand(0.05:0.01:0.1) * bigger_cyl_size[3]
    cylinder_top_curvature = rand(0.05:0.01:0.2) * bigger_cyl_size[3]
    angle = 0.0#rand(0:0.01:2π)
    density_inside = rand(0.05:0.3:0.7)
    pipe_len = rand(0.8:1.5) * (bigger_cyl_size[3] - cylinder_bottom_curvature + cylinder_top_curvature) #0.2 to 0.8 of the cylinder size[3]

    pipe_cross = rand(0.1:0.3) * bigger_cyl_size[1]
    pipe_cross_section = (pipe_cross, pipe_cross) #0.02 to 0.1 of the cylinder size[1] and cylinder size[2]
    # pipe_density=rand(0.1:0.9) #0.1 to 0.9
    pipe_density = rand(0.1:0.9) #0.1 to 0.9
    dispenser_len = rand(0.2:0.4) * bigger_cyl_size[3] #0.1 to 0.3 of the cylinder size[3]
    dispenser_cross_section = (rand(0.1:0.3) * bigger_cyl_size[1], rand(0.1:0.3) * bigger_cyl_size[2]) #0.02 to 0.1 of the cylinder size[1] and cylinder size[2]
    dispenser_density = rand(0.1:0.5)
    
    dual_phase_percentage = rand(0.6:1.0)

    density_inside_b = rand(0.3:0.01:0.9)*density_inside

    # Create the geometric object
    curr_time = Dates.now()
    ob, vol, ob4, ob5, ob4b, ob5b, ob_cut, ob_menisc_cut, ob_cyl_mask, ob5b_mask = empty_cylinder_with_half_sphere_bottom_p(
        center_cylinder, # needs to be within the max(bigger cylinder size) from edges of the image
        bigger_cyl_size, # its max need to be less than 0.75 of image size
        cylinder_wall_thickness, # less than 1/4 of the bigger_cyl_size more than 0.1
        cylinder_bottom_curvature, # min 0 max half of the cylinder size[3]
        cylinder_top_curvature, # min 0 max half of the cylinder size[3]
        angle, # between 0 and 2 pi
        density_inside, # between 0.05 and 0.9
        pipe_len,
        pipe_cross_section,
        pipe_density,
        dispenser_len,
        dispenser_cross_section,
        dispenser_density,
        len_cut,
        menisc_radius,
        dual_phase_percentage, density_inside_b)
    return (ob, vol["can_inside"], [("density_inside", density_inside),
            ("center_cylinder", center_cylinder), # needs to be within the max(bigger cylinder size) from edges of the image
            ("bigger_cyl_size", bigger_cyl_size), # its max need to be less than 0.75 of image size
            ("cylinder_wall_thickness", cylinder_wall_thickness), # less than 1/4 of the bigger_cyl_size more than 0.1
            ("cylinder_bottom_curvature", cylinder_bottom_curvature), # min 0 max half of the cylinder size[3]
            ("cylinder_top_curvature", cylinder_top_curvature), # min 0 max half of the cylinder size[3]
            ("angle", angle), # between 0 and 2 pi
            ("pipe_len", pipe_len),
            ("pipe_cross_section", pipe_cross_section),
            ("pipe_density", pipe_density),
            ("dispenser_len", dispenser_len),
            ("dispenser_cross_section", dispenser_cross_section),
            ("dispenser_density", dispenser_density), ("len_cut", len_cut), ("dual_phase_percentage", dual_phase_percentage), ("density_inside_b", density_inside_b)], ob4, ob5, ob4b, ob5b, ob_cut, ob_menisc_cut, ob_cyl_mask, ob5b_mask)


end

"""
get a phantom cylinder and return a boolean mask of it
"""
function get_cylinder_bool_mask(ob_el)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end

function get_half_s_bool(ob_el)
    return (phantom(axes(ig)..., ob_el) .!= 0)

end


function get_temp_folder()
    temp_dir = mktempdir()
    return temp_dir
end

function get_random_can_uploaded()
    while(true)
        uuid = UUIDs.uuid4()
        wandb.log(Dict("start main uuid $uuid_main; loc uuid $uuid "=>uuid_main ))

        main_foloder = "$(get_temp_folder())/$(uuid)"
        # Create the main folder if it does not exist
        if !isdir(main_foloder)
            mkpath(main_foloder)
        end


        try
            #to monitor how long it take
            curr_time = Dates.now()
            #unique name
            #get random can 
            ob, vol, args, ob4, ob5, ob4b, ob5b, ob_cut, ob_menisc_cut, ob_cyl_mask, ob5b_mask = random_can(ig, (dims[1] * spacing[1] * 0.9, dims[2] * spacing[2] * 0.9, dims[3] * spacing[3] * 0.9))
            ob2_a, ob2_b, ob3, pipe, pipe_in, plastic_dispenser, ob_cut = ob
            wandb.log(Dict("generated objects $uuid_main; loc uuid $uuid " =>uuid_main ))

            
            density_inside = args[1][2]
            json_path = "$(main_foloder)/argss.json"
            hdf5_path = "$(main_foloder)/data.h5"
            # Convert the list of tuples to a dictionary
            args_dict = Dict(arg[1] => arg[2] for arg in args)
            f = h5open(hdf5_path, "w")
            for arg in args
                write(f, arg[1], collect(arg[2]))
            end
            close(f)

            # Save the dictionary to a JSON file
            open(json_path, "w") do io
                JSON.print(io, args_dict)
            end
            #get projection first of a main cylinders and a pipe then the rest of the elements
            proj_arc_b, rg = get_CTFAN_proj(ob, !low_res)
            proj_arc, rg = get_CTFAN_proj([ob4, ob5, ob4b, ob5b], !low_res)

            @info "generated ctfan_proj $uuid_main " 
            wandb.log(Dict("generated ctfan_proj $uuid_main "  =>uuid_main ))

            #get masks of diffrent elements
            cyl_inner_bool = get_cylinder_bool_mask([ob[1], ob[2]])
            cyl_inner_bool_a = get_cylinder_bool_mask([ob[1]])
            cyl_inner_bool_b = get_cylinder_bool_mask([ob[2]])
            plastic_dispenser_bool = get_cylinder_bool_mask([plastic_dispenser])



            cyl_main_bool = get_cylinder_bool_mask([ob_cyl_mask])

            pipe_bool = get_cylinder_bool_mask([ob[4]])
            pipe_in_bool = get_cylinder_bool_mask([pipe_in])
            pipe_only_bool = pipe_bool .& .!pipe_in_bool

            ob_cut = get_cylinder_bool_mask([ob_cut])
            ob_menisc_cut = get_half_s_bool([ob_menisc_cut])
            bottom_h_sph_bigger = get_half_s_bool([ob5])
            upper_h_sph_bigger = get_half_s_bool([ob5b_mask])


            sitk.WriteImage(sitk.GetImageFromArray(UInt8.(ob_menisc_cut)), "$(main_foloder)/menisc_cut_a.nii.gz")

            ob_menisc_cut = ob_menisc_cut .& ob_cut


            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(ob_menisc_cut)), "$(main_foloder)/menisc_cut_b.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_inner_bool_a)), "$(main_foloder)/cyl_inner_bool_a.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_inner_bool_b)), "$(main_foloder)/cyl_inner_bool_b.nii.gz")


            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_inner_bool)), "$(main_foloder)/cyl_inner_bool.nii.gz")
            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_main_bool)), "$(main_foloder)/cyl_main_bool.nii.gz")
            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(pipe_bool)), "$(main_foloder)/pipe_bool.nii.gz")


            # proj_arc=reshape(proj_arc,(dims[1],dims[2],dims[3],10))
            # proj_arc_fl=Float32.(ustrip.(proj_arc))
            sitk.WriteImage(sitk.GetImageFromArray(proj_arc), "$(main_foloder)/example_cylinder_proj_arc.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(proj_arc_b), "$(main_foloder)/example_cylinder_proj_arc_halph_spheres.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(proj_arc_b + proj_arc), "$(main_foloder)/example_cylinder_added.nii.gz")

            f = h5open(hdf5_path, "r+")
            write(f, "cyl_inner_bool_a", UInt8.(cyl_inner_bool_a))
            write(f, "cyl_inner_bool_b", UInt8.(cyl_inner_bool_b))

            write(f, "example_cylinder_proj_arc", proj_arc)
            write(f, "example_cylinder_proj_arc_halph_spheres", proj_arc_b)
            write(f, "example_cylinder_added", proj_arc_b + proj_arc)

            write(f, "ob_cut", UInt8.(ob_cut))
            write(f, "ob_menisc_cut", UInt8.(ob_menisc_cut))
            close(f)

            #half spheres
            plan = plan_fbp(rg, ig; window=Window(Hamming(), 1.0))
            fdk_arc = fdk(plan, proj_arc)
            projs_views = Float32.(ustrip.(fdk_arc[:, :, :, 1]))
            phhh = phantom(axes(ig)..., [ob4, ob5, ob4b, ob5b])
            #setting to 0 all values outside the phantom of halph spheres
            projs_views[phhh.!=0.0] .= 0


            #cylinders
            fdk_arc = fdk(plan, proj_arc_b)
            projs_viewsb = Float32.(ustrip.(fdk_arc[:, :, :, 1]))

            #copy pipe
            pipee = copy(projs_viewsb)
            pipe_bool_b = .!((pipe_bool .& ob_menisc_cut) .| ((pipe_bool .& upper_h_sph_bigger)))
            pipee[pipe_bool_b] .= 0.0

            # #setting all outside cylinders to zeros       
            # projs_viewsb[cyl_main_bool].=0.0
            # #removing pipe to add it later
            mask_a = cyl_main_bool .| plastic_dispenser_bool
            projs_viewsb[.!mask_a] .= 0.0

            # sitk.WriteImage(sitk.GetImageFromArray(projs_views), "$(main_foloder)/example_cylinder_prim.nii.gz")

            projs_views_c = projs_views + projs_viewsb
            # projs_views_c[ob_menisc_cut].=0.0
            projs_views_c[phantom(axes(ig)..., [ob4]).!=0.0] .= 0

            # msmm=mean( projs_views_c[ob_cut] )
            projs_views_c[ob_menisc_cut] .= 0.0

            mask = (.!cyl_main_bool) .& (.!upper_h_sph_bigger) .& (.!plastic_dispenser_bool)
            projs_views_c[mask] .= 0.0
            projs_views_c = projs_views_c + (pipee)

            # projs_views_c[(.!(cyl_main_bool).|(.!upper_h_sph_bigger).|(.!plastic_dispenser_bool))].=0.0

            result_tensor = cyl_inner_bool .& .!ob_menisc_cut
            result_tensor = result_tensor .& .!bottom_h_sph_bigger
            result_tensor = result_tensor .& .!pipe_only_bool


            # Generate Gaussian noise tensor with the same shape as projs_views_c
            noise_mean = 0.0
            noise_std = 0.01
            noise_tensor = noise_mean .+ noise_std .* randn(size(projs_views_c))
            noise_tensor[(.!result_tensor)] .= 0.0

            # Add the noise tensor to the projection views
            projs_views_c_noisy = projs_views_c .+ noise_tensor
            # Save the noisy projection views
            sitk.WriteImage(sitk.GetImageFromArray((projs_views_c_noisy .* 1000)), "$(main_foloder)/example_cylinder_fulll_noisy.nii.gz")

            f = h5open(hdf5_path, "r+")
            write(f, "example_cylinder_fulll_noisy", projs_views_c_noisy .* 1000)
            write(f, "example_cylinder_fulll_disc", projs_views_c .* 1000)
            write(f, "main_segm", UInt8.(result_tensor))

            close(f)


            sitk.WriteImage(sitk.GetImageFromArray(projs_views_c), "$(main_foloder)/example_cylinder_fulll.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(projs_views_c .* 1000), "$(main_foloder)/example_cylinder_fulll_disc.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(UInt8.(ob_cut)), "$(main_foloder)/ob_cut.nii.gz")
            sitk.WriteImage(sitk.GetImageFromArray(UInt8.(result_tensor)), "$(main_foloder)/main_segm.nii.gz")
            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_inner_bool)), "/workspaces/synthethic_tomo/data/cylinder_inner.nii.gz")
            # Compress the folder to a zip archive
            zip_path = "$(main_foloder).zip"
            run(`zip -r $zip_path $main_foloder`)
            command = `gcloud storage cp $zip_path gs://bucket_metro_tk/`
            # Execute the command
            run(command)
            rm(main_foloder; force=true, recursive=true)
            rm(zip_path; force=true)

            println("Time to get single image: ", Dates.now() - curr_time)
            # @info "precess $uuid_main ; Time to get single image: $( Dates.now() - curr_time)" 
            wandb.log(Dict("precess $uuid_main ; Time to get single image: $( Dates.now() - curr_time)" =>uuid_main ))


        catch e
            println(e)
            # @debug "precess $uuid_main ; error $e" 
            wandb.log(Dict("precess $uuid_main ; error $e"  =>uuid_main ))

        end
    end    
end

# get_random_can_uploaded()






####################################################3
#
#  Ionic Chamber 
#
#
####################################################3

# Define material densities
const GRAPHITE_DENSITY = 1.0f0
const COPPER_DENSITY = 10.0f0
const ALUMINUM_DENSITY = 3.7f0
const POLYETHYLENE_DENSITY = 0.94f0
const AIR_DENSITY = -1.0f0

function get_cylinder_bool_mask(ob_el,ig)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end


function get_CTFAN_proj(ob,is_high_res=false)
    p = (ns = 300, ds = 0.15cm, nt = 300, dt = 0.15cm, na = 300, dsd = 250cm, dod = 100cm)
    # p = (ns = 800, ds = 0.15cm, nt = 800, dt = 0.15cm, na = 800, dsd = 200cm, dod = 40cm)
    # p = (ns = 800, ds = 0.15cm, nt = 800, dt = 0.2cm, na = 800, dsd = 200cm, dod = 40cm)
    # p = (ns = 500, ds = 0.2cm, nt = 500, dt = 0.2cm, na = 300, dsd = 40cm, dod = 40cm)
    if(is_high_res)
        p = (ns = 650, ds = 0.15cm, nt = 650, dt = 0.15cm, na = 650, dsd = 250cm, dod = 100cm)
    end
    
    rg = CtFanArc( ; p...)
    rayss=rays(rg)
    curr_time = Dates.now()
    proj_arc=zeros(Float32, p.ns*p.nt*p.na) #initialization
    radon(rayss, ob,proj_arc)
    # proj_arc=RadonKA.radon(collect(rayss), ob)
    proj_arc=reshape(proj_arc, (p.ns,p.nt, p.na))
    println("Time to get radon: ", Dates.now() - curr_time)
    return proj_arc,rg
end    


function generate_image_pure(ob, ob_long, ig)
    # Use ob_long for projection
    proj_arc, rg = get_CTFAN_proj([ob_long], false)
    plan = plan_fbp(rg, ig; window=Window(Hamming(), 1.0))

    # Use ob for mask
    mask = get_cylinder_bool_mask([ob], ig)

    fdk_arc = fdk(plan, proj_arc)
    projs_views = Float32.(ustrip.(fdk_arc[:, :, :, 1]))
    projs_views[.!mask] .= 0.0
    return projs_views, mask
end

function create_two_cyl(x, y, z, rx, ry, length, angle_x, angle_y, angle_z, density)
    normal = cylinder(x, y, z, rx, ry, length, angle_x, angle_y, angle_z, density)
    # Create elongated version with 2x length and adjusted z center
    elongated = cylinder(x, y, z , rx*5, ry*5, length*4, angle_x, angle_y, angle_z, density)
    return normal, elongated
end


function ionic_chamber_p(
    center_main_cylinder,
    main_cylinder_size,
    graphite_top_thickness,
    air_top_thickness,
    air_thicknes,
    copper_el_size,
    polyethylene_size_dif,
    insulator_thicknes,
    central_electrode_size,
    stem_height,
    angle,
    thick_part_length,
    thin_part_length,
    graphite_insertion_len,
    thick_radius
)
    # Convert input values to have units
    main_cylinder_size = [x*cm for x in main_cylinder_size]
    copper_el_size = [x*cm for x in copper_el_size]
    polyethylene_size_dif = [x*cm for x in polyethylene_size_dif]
    central_electrode_size = [x*cm for x in central_electrode_size]
    stem_height = stem_height*cm
    graphite_top_thickness = graphite_top_thickness*cm
    air_top_thickness = air_top_thickness*cm
    air_thicknes = air_thicknes*cm
    insulator_thicknes = insulator_thicknes*cm
    thick_part_length = thick_part_length*cm
    thin_part_length = thin_part_length*cm
    graphite_insertion_len = graphite_insertion_len*cm
    thick_radius = thick_radius*cm

    # Calculate positions with units
    total_height = main_cylinder_size[3]
    stem_bottom_z = (center_main_cylinder[3]*cm) - (total_height/2)
    stem_z = stem_bottom_z + stem_height/2        # <-- Define stem_z here
    # All electrode cylinders start here
    electrode_bottom_z = stem_bottom_z

    # Radii hierarchy
    copper_radius = thick_radius * 0.2
    # Inner polyethylene == thin electrode
    inner_poly_radius = thick_radius * 0.8
    # Aluminium strip
    aluminium_strip_radius = inner_poly_radius * 1.3
    # Outer polyethylene
    outer_poly_radius = aluminium_strip_radius * 1.2

    # Thick part cylinders: same starting Z, full thick_part_length
    copper_wire, copper_wire_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        electrode_bottom_z + (thick_part_length/2),
        copper_radius, copper_radius, thick_part_length,
        angle, 0., 0, COPPER_DENSITY
    )

    polyethylene, polyethylene_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        electrode_bottom_z + (thick_part_length/2),
        outer_poly_radius, outer_poly_radius, thick_part_length,
        angle, 0., 0, -0.9
    )

    correction, correction_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        electrode_bottom_z + (thick_part_length/2),
        outer_poly_radius, outer_poly_radius, thick_part_length,
        angle, 0., 0, ((-1)*(ALUMINUM_DENSITY))
    )

    aluminium_strip, aluminium_strip_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        electrode_bottom_z + (thick_part_length/2),
        aluminium_strip_radius, aluminium_strip_radius, thick_part_length,
        angle, 0., 0, ALUMINUM_DENSITY-1.0
    )

    # Inner polyethylene shorter by graphite_insertion_len
    inner_polyethylene_length = thick_part_length - graphite_insertion_len
    inner_polyethylene, inner_polyethylene_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        electrode_bottom_z + (inner_polyethylene_length/2),
        inner_poly_radius, inner_poly_radius, inner_polyethylene_length,
        angle, 0., 0, -(ALUMINUM_DENSITY-1.0)
    )

    # Thin part specs: same radius as inner_poly, no overlap, minus air_top_thickness
    thin_total = graphite_insertion_len + thin_part_length
    thin_electrode_len = thin_total - air_top_thickness
    thin_electrode_start = electrode_bottom_z + inner_polyethylene_length
    thin_part_electrode, thin_part_electrode_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        thin_electrode_start + (thin_electrode_len/2),
        inner_poly_radius, inner_poly_radius, thin_electrode_len,
        angle, 0., 0, GRAPHITE_DENSITY
    )

    # Calculate air cylinder radius based on thin part
    air_radius = inner_poly_radius + air_thicknes

    # Shift air cylinder to start above thick part (avoiding overlap)
    air_start = electrode_bottom_z + thick_part_length
    air, air_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm,
        air_start + (thin_total/2),
        air_radius, air_radius,  # Use calculated air_radius instead of main_cylinder_size
        thin_total,
        angle, 0., 0, AIR_DENSITY
    )

    # Stem and widest remain unchanged
    # Create cylinders with proper units and positions
    widest, widest_long = create_two_cyl(
        (center_main_cylinder[1]*cm), (center_main_cylinder[2]*cm), (center_main_cylinder[3]*cm),
        main_cylinder_size[1], main_cylinder_size[2], main_cylinder_size[3],
        angle, 0., 0, GRAPHITE_DENSITY
    )

    aluminium_stem, aluminium_stem_long = create_two_cyl(
        center_main_cylinder[1]*cm, center_main_cylinder[2]*cm, stem_z,
        main_cylinder_size[1], main_cylinder_size[2], stem_height,
        angle, 0., 0, ALUMINUM_DENSITY
    )

    ob = [
        widest, aluminium_stem, copper_wire, polyethylene, aluminium_strip,
        inner_polyethylene, thin_part_electrode, air, correction
    ]

    ob_long = [
        widest_long, aluminium_stem_long, copper_wire_long, polyethylene_long,
        aluminium_strip_long, inner_polyethylene_long, thin_part_electrode_long,
        air_long, correction_long
    ]

    return ob, ob_long
end

function get_random_ionic_ch_uploaded()
    while(true)
        uuid = UUIDs.uuid4()
        wandb.log(Dict("start main uuid $uuid_main; loc uuid $uuid "=>uuid_main ))

        image_size_cm = (dims[1] * spacing[1], dims[2] * spacing[2], dims[3] * spacing[3])

        main_foloder = "$(get_temp_folder())/$(uuid)"
        # Create the main folder if it does not exist
        if !isdir(main_foloder)
            mkpath(main_foloder)
        end


        # try
            #to monitor how long it take
            curr_time = Dates.now()
            #unique name
            #get random can 
         
            center_main_cylinder = (0.0, 0.0, 0.0)
            center_cylinder = center_main_cylinder
        
            bigger_cyl_size = [
                max(1.0+0.25*2,rand(0.3:0.01:0.5) * image_size_cm[1]/2),
                max(1.0+0.25*2,rand(0.3:0.01:0.5) * image_size_cm[2]/2),
                1.0
            ]
            bigger_cyl_size[3] = rand(1.0:0.01:8.0) * maximum(collect(bigger_cyl_size))
            bigger_cyl_size[3] = min(bigger_cyl_size[3],  image_size_cm[3]*0.6)
        
            main_cylinder_size = bigger_cyl_size
        
        
            graphite_top_thickness = rand(0.05:0.01:0.35)*bigger_cyl_size[3]
            air_top_thickness = rand(0.15:0.01:0.35)*bigger_cyl_size[3]
            air_thicknes = rand(0.15:0.01:0.4) * min(bigger_cyl_size[2], bigger_cyl_size[1])
            graphite_top_thickness = air_thicknes * rand(0.8:0.01:1.2)  # Make proportional to air_thicknes
        
            copper_el_size = [
                rand(0.05:0.19)  * bigger_cyl_size[1],
                rand(0.05:0.19)  * bigger_cyl_size[2],
                rand(0.25:0.6)  * bigger_cyl_size[3]
            ]    
            polyethylene_size_dif = [
                rand(0.03:0.01:0.3) * bigger_cyl_size[1],
                rand(0.03:0.01:0.3) * bigger_cyl_size[2],
                rand(0.03:0.01:0.3) * bigger_cyl_size[3]
            ]    
            insulator_thicknes = rand(0.1:0.01:1.0)
            central_electrode_size = [
                rand(0.2:0.01:0.6) * bigger_cyl_size[1],
                rand(0.2:0.01:0.6) * bigger_cyl_size[2],
                rand(0.4:0.01:0.6) * bigger_cyl_size[3]
            ]
        
            stem_height = rand(0.1:0.4) * bigger_cyl_size[3]
            angle = 0.0 #rand(0:0.01:2π)
        
        
        
            # Example usage (outside ionic_chamber_p):
            electrode_length = (stem_height + main_cylinder_size[3]) - (graphite_top_thickness + air_top_thickness)
            thick_part_length = rand(0.2:0.01:0.5) * electrode_length
            thin_part_length = electrode_length - thick_part_length
            graphite_insertion_len = rand(0.1:0.01:0.3) * thick_part_length
            thick_radius = clamp(rand(0.2:0.01:0.4) * main_cylinder_size[1], 0.01, main_cylinder_size[1])
        
            ob, ob_long = ionic_chamber_p(
                center_main_cylinder,
                main_cylinder_size,
                graphite_top_thickness,
                air_top_thickness,
                air_thicknes,
                copper_el_size,
                polyethylene_size_dif,
                insulator_thicknes,
                central_electrode_size,
                stem_height,
                angle,
                thick_part_length,
                thin_part_length,
                graphite_insertion_len,
                thick_radius
            )
            widest, aluminium_stem, copper_wire, polyethylene, aluminium_strip,
                inner_polyethylene, thin_part_electrode, air = ob
            widest_long, aluminium_stem_long, copper_wire_long, polyethylene_long, aluminium_strip_long,
                inner_polyethylene_long, thin_part_electrode_long, air_long=ob_long


        #parameters for CT reconstruction
        widest_im, widest_mask = generate_image_pure(widest, widest_long, ig)
        aluminium_stem_im, aluminium_stem_mask = generate_image_pure(aluminium_stem, aluminium_stem_long, ig)
        air_im, air_mask = generate_image_pure(air, air_long, ig)
        copper_wire_im, copper_wire_mask = generate_image_pure(copper_wire, copper_wire_long, ig)
        polyethylene_im, polyethylene_el_mask = generate_image_pure(polyethylene, polyethylene_long, ig)
        inner_polyethylene_im, inner_polyethylene_mask = generate_image_pure(inner_polyethylene, inner_polyethylene_long, ig)
        aluminium_strip_im, aluminium_strip_mask = generate_image_pure(aluminium_strip, aluminium_strip_long, ig)
        thin_part_electrode_im, thin_part_electrode_mask = generate_image_pure(thin_part_electrode, thin_part_electrode_long, ig)

        to_subtract = Float32.((.!air_mask).&thin_part_electrode_mask)*(-(ALUMINUM_DENSITY-1))

        all_cyls_gen = [to_subtract,widest_im, aluminium_stem_im, copper_wire_im, air_im, polyethylene_im, inner_polyethylene_im, aluminium_strip_im, thin_part_electrode_im]
        projs_views = sum(all_cyls_gen)
        sitk.WriteImage(sitk.GetImageFromArray(projs_views), "$(main_foloder)/projs_views_$(uuid).nii.gz")

        # Add small Gaussian noise to projs_views
        noise_level = 0.015  # Adjust the noise level as needed
        projs_views = projs_views .+ (noise_level * randn(size(projs_views)))

        sitk.WriteImage(sitk.GetImageFromArray(projs_views), "$(main_foloder)/projs_views_with_noise_$(uuid).nii.gz")

         
        comm=(air_mask.&thin_part_electrode_mask)
        vol=air_mask.&(.!comm)
        vol=sum(vol)
        vol=spacing[1]*spacing[2]*spacing[3]*vol

        json_path="$(main_foloder)/example_$(uuid).json"

        variables = Dict(
            "center_main_cylinder" => center_main_cylinder,
            "main_cylinder_size" => main_cylinder_size,
            "graphite_top_thickness" => graphite_top_thickness,
            "air_top_thickness" => air_top_thickness,
            "air_thicknes" => air_thicknes,
            "copper_el_size" => copper_el_size,
            "polyethylene_size_dif" => polyethylene_size_dif,
            "insulator_thicknes" => insulator_thicknes,
            "central_electrode_size" => central_electrode_size,
            "stem_height" => stem_height,
            "angle" => angle,
            "vol"=>vol
        )

        # write(json_path, JSON.json(variables, 4)) # 4 is for pretty printing with indent

        # Save the dictionary to a JSON file
        open(json_path, "w") do io
            JSON.print(io, variables)
        end
            # sitk.WriteImage(sitk.GetImageFromArray(UInt8.(cyl_inner_bool)), "/workspaces/synthethic_tomo/data/cylinder_inner.nii.gz")
            # Compress the folder to a zip archive
            zip_path = "$(main_foloder).zip"
            run(`zip -r $zip_path $main_foloder`)
            command = `gcloud storage cp $zip_path gs://bucket_metro_ionic_ch/high_res`
            # Execute the command
            run(command)
            rm(main_foloder; force=true, recursive=true)
            rm(zip_path; force=true)

            println("Time to get single image: ", Dates.now() - curr_time)
            # @info "precess $uuid_main ; Time to get single image: $( Dates.now() - curr_time)" 
            wandb.log(Dict("precess $uuid_main ; Time to get single image: $( Dates.now() - curr_time)" =>uuid_main ))


        # catch e
        #     println(e)
        #     # @debug "precess $uuid_main ; error $e" 
        #     wandb.log(Dict("precess $uuid_main ; error $e"  =>uuid_main ))

        # end
    end    
end

get_random_ionic_ch_uploaded()


#docker build -t my-container .
#docker run -d -p 8080:80 my-container
