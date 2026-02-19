using Unitful: cm
using ImagePhantoms
import ImagePhantoms as IP
using Unitful
using LinearAlgebra
using PyCall

# Include the ellipsoid intersection functions

"""
    compute_fluid_volume_in_can(
        phantoms_dict::Dict,
        spacing::Tuple{Float64,Float64,Float64},
        first_ball::Bool,
        second_ball::Bool,
        rounded_bottom::Bool,
        params::Dict{String,Any}
    ) -> Tuple{Float64,Float64}

Calculate fluid volume in a phantom can using both numerical and analytical approaches.

# Arguments
- `phantoms_dict`: Dictionary containing all phantom boolean masks and their density values
- `spacing`: Voxel spacing tuple (dx, dy, dz) in cm
- `first_ball`: Boolean flag indicating if first ball is present
- `second_ball`: Boolean flag indicating if second ball is present  
- `rounded_bottom`: Boolean flag indicating if the bottom is rounded
- `params`: Dictionary containing geometry parameters and calculated values

# Returns
A tuple containing (numerical_volume, analytical_volume) in cubic centimeters
"""
function compute_fluid_volume_in_can(
    phantoms_dict::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    first_ball::Bool,
    second_ball::Bool,
    rounded_bottom::Bool,
    params::Dict{String,Any}
)
    # Safe parameter extraction with defaults
    function safe_get(dict, key, default=0.0)
        return haskey(dict, key) ? dict[key] : default
    end

    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    dual_phase_percentage = safe_get(params, "dual_phase_percentage", 1.0)

    # ================= NUMERICAL VOLUME CALCULATION =================
    # Simply use the fluid_mask already calculated during phantom creation
    fluid_mask = params["fluid_mask"]
    numerical_vol = sum(fluid_mask) * voxel_volume

    # ================= ANALYTICAL VOLUME CALCULATION =================
    # Base fluid volume: Different handling based on dual phase percentage
    local main_volume
    if dual_phase_percentage == 1.0
        # Single fluid phase
        main_volume = IP.volume(phantoms_dict["ob2_a"][3])
        println("Base fluid volume (single phase): $(main_volume)")
    else
        # Dual fluid phase
        main_volume = IP.volume(phantoms_dict["ob2_a"][3]) + IP.volume(phantoms_dict["ob2_b"][3])
        println("Base fluid volume (dual phase): $(main_volume)")
    end

    # STEP 1: Calculate the top cut plane and cylinder geometry
    # Extract parameters needed for cut calculations
    center_cylinder = params["center_cylinder"]
    bigger_cyl_size = params["bigger_cyl_size"]
    cylinder_wall_thickness = params["cylinder_wall_thickness"]
    len_cut = params["len_cut"]
    x_cut_angle = params["x_cut_angle"]
    y_cut_angle = params["y_cut_angle"]
    output_folder = params["output_folder"]

    # Get the top of the cylinder (from ob3)
    cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)

    # Calculate the lowest point of the cut plane intersection with ob3
    # The cut is at: center_cylinder[3] + ((bigger_cyl_size[3]/2) - (len_cut/2))
    cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))

    # Due to tilt, the lowest point depends on angle and radius
    radius = (bigger_cyl_size[1] - cylinder_wall_thickness)

    # Calculate the maximum drop due to tilt
    max_drop_x = radius * sin(abs(x_cut_angle))
    max_drop_y = radius * sin(abs(y_cut_angle))
    max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

    # Calculate the z-coordinate of the lowest point of the cut
    lowest_cut_z = cut_center_z - max_drop - len_cut / 2

    # Calculate the distance from cylinder top to lowest cut point
    distance_from_top = cylinder_top_z - lowest_cut_z

    # Calculate the height of the cylindrical portion below the cut
    cylinder_height_below_cut = bigger_cyl_size[3] - distance_from_top

    # Create the cylinder object representing the fluid below the cut
    calc_cyl = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2)cm,
        (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
        (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
        cylinder_height_below_cut * cm,
        0, 0, 0, 1.0f0
    )

    # Calculate the analytic volume of this cylinder
    cylinder_below_cut_vol = IP.volume(calc_cyl)

    # Generate a boolean mask for this cylinder
    calc_cyl_mask = phantom(axes(ig)..., [calc_cyl]) .!= 0

    # Save the mask for debugging
    if haskey(params, "output_folder")
        sitk = pyimport("SimpleITK")
        calc_cyl_uint8 = UInt8.(calc_cyl_mask) 
        calc_cyl_img = sitk.GetImageFromArray(calc_cyl_uint8)
        sitk.WriteImage(calc_cyl_img, joinpath(output_folder, "analytical_bottom_cyl.nii.gz"))
    end

    # Calculate difference between fluid_mask and calculated cylinder
    # We're only interested in the part above the cylinder center
    cylinder_center_z = Int(round(center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2))
    diff_mask = fluid_mask .& .!calc_cyl_mask

    # Apply z-filter to keep only the part above cylinder center
    for i in 1:size(diff_mask, 1)
        for j in 1:size(diff_mask, 2)
            for k in 1:size(diff_mask, 3)
                if k < cylinder_center_z && diff_mask[i, j, k]
                    diff_mask[i, j, k] = false
                end
            end
        end
    end

    # Calculate the volume of this difference
    diff_vol = sum(diff_mask) * voxel_volume

    # Save the difference mask for debugging
    if haskey(params, "output_folder")
        diff_uint8 = UInt8.(diff_mask) .* 255
        diff_img = sitk.GetImageFromArray(diff_uint8)
        sitk.WriteImage(diff_img, joinpath(output_folder, "fluid_cylinder_diff.nii.gz"))
    end

    # Update the main volume calculation
    main_volume = cylinder_below_cut_vol + (diff_vol * cm^3)
    println("Cylinder below cut: $(cylinder_below_cut_vol)")
    println("Difference above cylinder: $(diff_vol) cm³")

    # STEP 2: Handle bottom geometry
    if rounded_bottom
        # For rounded bottom: add half inner torus, subtract half outer sphere
        inner_torus_volume_analytical = safe_get(params, "inner_torus_volume_analytical")
        outer_sphere_volume = safe_get(params, "outer_sphere_volume")

        if inner_torus_volume_analytical > 0.0
            inner_torus_cm3 = inner_torus_volume_analytical * cm^3
            println("Adding half inner torus: $(inner_torus_cm3/2)")
            main_volume += inner_torus_cm3 / 2
        end

        if outer_sphere_volume > 0.0
            outer_sphere_cm3 = outer_sphere_volume * cm^3
            println("Subtracting half outer sphere: $(outer_sphere_cm3/2)")
            main_volume -= outer_sphere_cm3 / 2
        end
    else
        # For non-rounded bottom
        if haskey(phantoms_dict, "ob5") && haskey(phantoms_dict, "center_bottom_out")
            ob5_obj = phantoms_dict["ob5"][3]  # Get the actual object, not the mask
            center_bottom_obj = phantoms_dict["center_bottom_out"][3]

            # Calculate the volume of ob5
            ob5_vol = IP.volume(ob5_obj)

            # Calculate the volume of center_bottom_out
            center_bottom_vol = IP.volume(center_bottom_obj)

            # Calculate the volume of intersection between ob5 and center_bottom_out
            # Using boolean operations on masks
            ob5_mask = phantoms_dict["ob5"][1]
            center_bottom_mask = phantoms_dict["center_bottom_out"][1]
            intersection_mask = ob5_mask .& center_bottom_mask
            intersection_vol = sum(intersection_mask) * voxel_volume

            # Calculate the part of center_bottom_out not in ob5
            center_bottom_unique_vol = center_bottom_vol - (intersection_vol * cm^3)

            println("Volume of ob5: $(ob5_vol)")
            println("Volume of center_bottom_out: $(center_bottom_vol)")
            println("Volume of intersection: $(intersection_vol) cm³")
            println("Unique volume of center_bottom_out: $(center_bottom_unique_vol)")

            # Subtract these volumes from main_volume
            main_volume = main_volume - ob5_vol - center_bottom_unique_vol
        end
    end

    # STEP 3: Subtract ball volumes if present
    if first_ball && haskey(phantoms_dict, "ball1")
        ball1_obj = phantoms_dict["ball1"][3]  # Get the actual object, not mask
        ball1_volume = IP.volume(ball1_obj)
        println("Subtracting ball1: $(ball1_volume)")
        main_volume -= ball1_volume
    end

    if second_ball && haskey(phantoms_dict, "ball2")
        ball2_obj = phantoms_dict["ball2"][3]  # Get the actual object, not mask
        ball2_volume = IP.volume(ball2_obj)
        println("Subtracting ball2: $(ball2_volume)")
        main_volume -= ball2_volume
    end

    # STEP 4: Calculate pipe overlap with calculated cylinder
    if haskey(phantoms_dict, "pipe")
        # Get pipe parameters
        pipe_obj = phantoms_dict["pipe"][3]
        pipe_mask = phantoms_dict["pipe"][1]

        # Calculate intersection volume between pipe and calc_cyl
        pipe_calc_cyl_intersection = pipe_mask .& calc_cyl_mask
        pipe_overlap_vol = sum(pipe_calc_cyl_intersection) * voxel_volume

        println("Pipe overlap with fluid: $(pipe_overlap_vol) cm³")
        main_volume -= pipe_overlap_vol * cm^3
    end

    # Convert analytic volume to scalar (removing units)
    analytical_vol = Unitful.ustrip(cm^3, main_volume)

    # Print final comparison
    println("Final analytical volume: $(main_volume) ($(analytical_vol) cm³)")
    println("Final numerical volume: $(numerical_vol) cm³")
    println("Difference: $(abs(analytical_vol - numerical_vol)) cm³")

    if numerical_vol > 0
        rel_diff = 100.0 * abs(analytical_vol - numerical_vol) / numerical_vol
        println("Relative difference: $(rel_diff)%")
    end

    return numerical_vol, analytical_vol
end

"""
    calculate_fluid_volume_for_phantom(
        phantoms_dict::Dict,
        params::Dict,
        spacing::Tuple{Float64,Float64,Float64},
        rounded_bottom::Bool,
        first_ball::Bool,
        second_ball::Bool,
        dual_phase_percentage::Float64
    ) -> Tuple{Float64, Float64}

Legacy wrapper function for backward compatibility
"""
function calculate_fluid_volume_for_phantom(
    phantoms_dict::Dict,
    params::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    rounded_bottom::Bool,
    first_ball::Bool,
    second_ball::Bool,
    dual_phase_percentage::Float64
)
    # Set the dual_phase_percentage in the params
    params["dual_phase_percentage"] = dual_phase_percentage

    # Call the main function
    return compute_fluid_volume_in_can(
        phantoms_dict,
        spacing,
        first_ball,
        second_ball,
        rounded_bottom,
        params
    )
end

function compute_fluid_volume_in_can_v2(
    phantoms_dict::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    params::Dict{String,Any},
    first_ball::Bool,
    second_ball::Bool,
    rounded_bottom::Bool
)
    # 1) Numerical solution
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    fluid_mask = params["fluid_mask"]
    numerical_vol = sum(fluid_mask) * voxel_volume

    # 2) Analytical base volume
    if params["dual_phase_percentage"] == 1.0
        main_volume = IP.volume(phantoms_dict["ob2_a"][3])
    else
        main_volume = IP.volume(phantoms_dict["ob2_a"][3]) + IP.volume(phantoms_dict["ob2_b"][3])
    end

    # 3) Handle top cut geometry and partial difference above cylinder
    # ...existing code for tilt & new cylinder creation...
    # let "cyl_vol" be IP.volume(...) of the newly created cylinder object
    # let "diff_above" be the fluid above cylinder center from fluid_mask
    # Add these volumes to main_volume:
    # main_volume = cyl_vol + diff_above

    # 4) Handle bottom geometry
    if rounded_bottom
        main_volume += (params["inner_torus_volume_analytical"] / 2)
        main_volume -= (params["outer_sphere_volume"] / 2)
    else
        # ob5 volume
        ob5_vol = IP.volume(phantoms_dict["ob5"][3])
        # center_bottom volume outside ob5
        # ...existing code for boolean difference...
        # main_volume -= ob5_vol + center_bottom_vol
    end

    # 5) Subtract ball volumes
    if first_ball
        ball1_vol = IP.volume(phantoms_dict["ball1"][3])
        main_volume -= ball1_vol
    end
    if second_ball
        ball2_vol = IP.volume(phantoms_dict["ball2"][3])
        main_volume -= ball2_vol
    end

    # 6) Adjust with outer_sphere and inner_torus for both A/B (if required)
    main_volume -= (params["outer_sphere_volume"] / 2)
    main_volume += (params["inner_torus_volume_analytical"] / 2)

    # 7) Subtract pipe volume inside the new cylinder
    # ...existing code for pipe object & overlap with cylinder...
    # main_volume -= pipe_overlap

    return numerical_vol, main_volume
end

"""
    compute_fluid_volume_in_can_v3(
        phantoms_dict::Dict,
        spacing::Tuple{Float64,Float64,Float64},
        first_ball::Bool,
        second_ball::Bool,
        rounded_bottom::Bool,
        params::Dict{String,Any},
        ig::ImageGeom
    ) -> Tuple{Float64,Float64}

Precisely calculate fluid volume in a phantom can, accounting for top cut tilt
and bottom geometry.

# Arguments
- `phantoms_dict`: Dictionary containing phantom objects and masks
- `spacing`: Voxel spacing (dx, dy, dz) in cm
- `first_ball`: Whether first ball is present
- `second_ball`: Whether second ball is present
- `rounded_bottom`: Whether can has rounded bottom
- `params`: Dictionary of parameters
- `ig`: ImageGeom object for phantom generation

# Returns
A tuple containing (numerical_volume, analytical_volume) in cubic centimeters
"""
function compute_fluid_volume_in_can_v3(
    phantoms_dict::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    first_ball::Bool,
    second_ball::Bool,
    rounded_bottom::Bool,
    params::Dict{String,Any},
    ig::ImageGeom
)
    # Safe parameter extraction
    function safe_get(dict, key, default=0.0)
        return haskey(dict, key) ? dict[key] : default
    end

    # Extract needed parameters
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    center_cylinder = params["center_cylinder"]
    bigger_cyl_size = params["bigger_cyl_size"]
    cylinder_wall_thickness = params["cylinder_wall_thickness"]
    dual_phase_percentage = params["dual_phase_percentage"]
    output_folder = params["output_folder"]
    angle = safe_get(params, "angle", 0.0)

    # NUMERICAL CALCULATION (simple count of fluid mask voxels)
    fluid_mask = params["fluid_mask"]
    numerical_vol = sum(fluid_mask) * voxel_volume

    # ANALYTICAL CALCULATION
    # 1. Base fluid volume based on dual phase
    if dual_phase_percentage == 1.0
        main_volume = IP.volume(phantoms_dict["ob2_a"][3])
    else
        main_volume = IP.volume(phantoms_dict["ob2_a"][3]) + IP.volume(phantoms_dict["ob2_b"][3])
    end
    println("Base fluid volume: $(main_volume)")

    # 2. Account for top cut with tilt
    # Calculate lowest point of tilted cut plane
    len_cut = params["len_cut"]
    x_cut_angle = params["x_cut_angle"]
    y_cut_angle = params["y_cut_angle"]

    # Get cylinder top z
    cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)

    # Calculate cut center and lowest point due to tilt
    cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))
    radius = bigger_cyl_size[1] - cylinder_wall_thickness

    # Calculate drop due to tilt angles
    max_drop_x = radius * sin(abs(x_cut_angle))
    max_drop_y = radius * sin(abs(y_cut_angle))
    max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

    # Calculate lowest point of cut
    lowest_cut_z = cut_center_z - max_drop - len_cut / 2

    # Calculate cylinder height below the cut
    cylinder_height_below_cut = bigger_cyl_size[3] - (cylinder_top_z - lowest_cut_z)

    # Create cylinder object representing fluid below cut
    calc_cyl = cylinder(
        center_cylinder[1]cm, center_cylinder[2]cm,
        (center_cylinder[3] - (bigger_cyl_size / 2) + cylinder_height_below_cut / 2)cm,
        (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
        (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
        cylinder_height_below_cut * cm,
        angle, 0, 0, 1.0f0
    )

    # Calculate analytical volume of cylinder
    cylinder_vol = IP.volume(calc_cyl)
    println("Volume of fluid cylinder below cut: $(cylinder_vol)")

    # Create boolean mask for the cylinder
    calc_cyl_mask = phantom(axes(ig)..., [calc_cyl]) .!= 0

    # Save the cylinder mask as a nifti file
    sitk = pyimport("SimpleITK")
    calc_cyl_uint8 = UInt8.(calc_cyl_mask) .* 255
    calc_cyl_img = sitk.GetImageFromArray(calc_cyl_uint8)
    sitk.WriteImage(calc_cyl_img, joinpath(output_folder, "analytical_bottom_cyl.nii.gz"))

    # Find difference between fluid mask and cylinder (only above center)
    cylinder_center_z = Int(round(center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2))
    diff_mask = fluid_mask .& .!calc_cyl_mask

    # Filter diff_mask to keep only parts above cylinder center
    for i in 1:size(diff_mask, 1)
        for j in 1:size(diff_mask, 2)
            for k in 1:size(diff_mask, 3)
                if k < cylinder_center_z && diff_mask[i, j, k]
                    diff_mask[i, j, k] = false
                end
            end
        end
    end

    # Calculate volume of the difference
    diff_vol = sum(diff_mask) * voxel_volume
    println("Volume of fluid above cylinder center: $(diff_vol) cm³")

    # Save the difference mask
    diff_uint8 = UInt8.(diff_mask) .* 255
    diff_img = sitk.GetImageFromArray(diff_uint8)
    sitk.WriteImage(diff_img, joinpath(output_folder, "fluid_cylinder_diff.nii.gz"))

    # Set main_volume to cylinder volume plus fluid above
    main_volume = cylinder_vol + (diff_vol * cm^3)

    # 3. Account for bottom geometry
    inner_torus_volume_analytical = safe_get(params, "inner_torus_volume_analytical", 0.0)
    outer_sphere_volume = safe_get(params, "outer_sphere_volume", 0.0)

    if rounded_bottom
        # A) For rounded bottom: add half inner torus, subtract half outer sphere
        if inner_torus_volume_analytical > 0.0
            inner_torus_vol_cm3 = inner_torus_volume_analytical * cm^3
            println("Adding half inner torus: $(inner_torus_vol_cm3/2)")
            main_volume += inner_torus_vol_cm3 / 2
        end

        if outer_sphere_volume > 0.0
            outer_sphere_vol_cm3 = outer_sphere_volume * cm^3
            println("Subtracting half outer sphere: $(outer_sphere_vol_cm3/2)")
            main_volume -= outer_sphere_vol_cm3 / 2
        end
    else
        # B) For flat bottom: handle ob5 and center_bottom_out
        if haskey(phantoms_dict, "ob5") && haskey(phantoms_dict, "center_bottom_out")
            # Get the actual geometric objects (third item in tuple)
            ob5_obj = phantoms_dict["ob5"][3]
            center_bottom_obj = phantoms_dict["center_bottom_out"][3]

            # Calculate analytical volumes
            ob5_vol = IP.volume(ob5_obj)
            center_bottom_vol_full = IP.volume(center_bottom_obj)

            # Get boolean masks
            ob5_mask = phantoms_dict["ob5"][1]
            center_bottom_mask = phantoms_dict["center_bottom_out"][1]

            # Calculate intersection numerically using masks
            intersection = ob5_mask .& center_bottom_mask
            intersection_vol = sum(intersection) * voxel_volume * cm^3

            # Calculate part of center_bottom_out not in ob5
            center_bottom_unique = center_bottom_vol_full - intersection_vol

            println("ob5 volume: $(ob5_vol)")
            println("center_bottom_out volume: $(center_bottom_vol_full)")
            println("Intersection volume: $(intersection_vol)")
            println("Unique center_bottom_out volume: $(center_bottom_unique)")

            # Save center_bottom_out as nifti
            center_bottom_uint8 = UInt8.(center_bottom_mask) .* 255
            center_bottom_img = sitk.GetImageFromArray(center_bottom_uint8)
            sitk.WriteImage(center_bottom_img, joinpath(output_folder, "center_bottom_out.nii.gz"))

            # Subtract both volumes (full ob5 and unique center_bottom)
            main_volume -= ob5_vol
            main_volume -= center_bottom_unique
        end
    end

    # 4. For both A and B: subtract balls, add/subtract torus/sphere, subtract pipe
    # 4.1 Subtract ball volumes if present
    if first_ball && haskey(phantoms_dict, "ball1")
        ball1_vol = IP.volume(phantoms_dict["ball1"][3])
        println("Subtracting ball1: $(ball1_vol)")
        main_volume -= ball1_vol
    end

    if second_ball && haskey(phantoms_dict, "ball2")
        ball2_vol = IP.volume(phantoms_dict["ball2"][3])
        println("Subtracting ball2: $(ball2_vol)")
        main_volume -= ball2_vol
    end

    # 4.2 Adjust for rounded bottom geometric correction for both A and B
    # We always subtract outer_sphere/2 and add inner_torus/2 as a final adjustment
    if outer_sphere_volume > 0.0 && !rounded_bottom
        outer_sphere_vol_cm3 = outer_sphere_volume * cm^3
        println("Final adjustment: subtracting half outer sphere: $(outer_sphere_vol_cm3/2)")
        main_volume -= outer_sphere_vol_cm3 / 2
    end

    if inner_torus_volume_analytical > 0.0 && !rounded_bottom
        inner_torus_vol_cm3 = inner_torus_volume_analytical * cm^3
        println("Final adjustment: adding half inner torus: $(inner_torus_vol_cm3/2)")
        main_volume += inner_torus_vol_cm3 / 2
    end

    # 4.3 Calculate pipe overlap with calc_cyl
    if haskey(phantoms_dict, "pipe") && haskey(params, "pipe_len")
        # Get pipe parameters for calculation
        pipe_cross_section = safe_get(params, "pipe_cross_section", (0.0, 0.0))
        pipe_len = params["pipe_len"]
        pipe_start = cylinder_top_z + safe_get(params, "cylinder_top_curvature", 0.0) - (pipe_len / 2)

        # Create pipe object similar to original definition
        pipe_obj = phantoms_dict["pipe"][3]

        # Calculate intersection between pipe and calc_cyl masks
        pipe_mask = phantoms_dict["pipe"][1]
        pipe_cyl_intersection = pipe_mask .& calc_cyl_mask
        pipe_overlap_vol = sum(pipe_cyl_intersection) * voxel_volume

        println("Pipe overlap with fluid: $(pipe_overlap_vol) cm³")
        main_volume -= pipe_overlap_vol * cm^3
    end

    # Convert main_volume to scalar (removing Unitful units)
    analytical_vol = Unitful.ustrip(cm^3, main_volume)

    # Print final results
    println("Final analytical volume: $(main_volume) ($(analytical_vol) cm³)")
    println("Final numerical volume: $(numerical_vol) cm³")
    println("Difference: $(abs(analytical_vol - numerical_vol)) cm³")
    if numerical_vol > 0
        rel_diff = 100.0 * abs(analytical_vol - numerical_vol) / numerical_vol
        println("Relative difference: $(rel_diff)%")
    end

    return numerical_vol, analytical_vol
end

"""
    compute_accurate_fluid_volume(
        phantoms_dict::Dict,
        spacing::Tuple{Float64,Float64,Float64},
        first_ball::Bool,
        second_ball::Bool,
        rounded_bottom::Bool,
        params::Dict{String,Any},
        ig::ImageGeom
    ) -> Tuple{Float64,Float64}

Calculate fluid volume in a phantom can with precise handling of complex geometries
including the tilted top cut, bottom shape, balls, and pipe overlaps.

# Arguments
- `phantoms_dict`: Dictionary containing phantom objects and their masks
- `spacing`: Voxel spacing (dx, dy, dz) in cm
- `first_ball`: Boolean flag indicating if first ball is present
- `second_ball`: Boolean flag indicating if second ball is present
- `rounded_bottom`: Boolean flag indicating if the bottom is rounded
- `params`: Dictionary containing geometry parameters and calculated values
- `ig`: ImageGeom object for phantom generation

# Returns
A tuple containing (numerical_volume, analytical_volume) in cubic centimeters
"""
function compute_accurate_fluid_volume(
    phantoms_dict::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    first_ball::Bool,
    second_ball::Bool,
    rounded_bottom::Bool,
    params::Dict{String,Any},
    ig::ImageGeom
)
    # Safe parameter extraction with defaults
    function safe_get(dict, key, default=0.0)
        return haskey(dict, key) ? dict[key] : default
    end

    # Extract base parameters
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    output_folder = params["output_folder"]
    dual_phase_percentage = safe_get(params, "dual_phase_percentage", 1.0)
    center_cylinder = params["center_cylinder"]
    bigger_cyl_size = params["bigger_cyl_size"]
    cylinder_wall_thickness = params["cylinder_wall_thickness"]
    angle = safe_get(params, "angle", 0.0)

    # STEP 1: NUMERICAL VOLUME CALCULATION
    # Simply use the fluid_mask that has been prepared
    fluid_mask = params["fluid_mask"]
    numerical_vol = sum(fluid_mask) * voxel_volume
    println("Numerical volume from fluid_mask: $(numerical_vol) cm³")

    # STEP 2: ANALYTICAL VOLUME CALCULATION - Base volume
    # Determine base volume from ob2_a/ob2_b based on dual phase percentage
    local main_volume
    if dual_phase_percentage == 1.0
        main_volume = IP.volume(phantoms_dict["ob2_a"][3])
        println("Base fluid volume (single phase): $(main_volume)")
    else
        main_volume = IP.volume(phantoms_dict["ob2_a"][3]) + IP.volume(phantoms_dict["ob2_b"][3])
        println("Base fluid volume (dual phases): $(main_volume)")
    end

    # STEP 3: Handle the top cut with tilt
    # Get the calc_cyl object that was created in empty_cylinder_with_half_sphere_bottom_p
    if haskey(phantoms_dict, "calc_cyl")
        calc_cyl_mask = phantoms_dict["calc_cyl"][1]

        # Calculate the analytical volume of this cylinder
        cylinder_below_cut_vol = IP.volume(phantoms_dict["calc_cyl"][3])
        println("Volume of fluid cylinder below cut: $(cylinder_below_cut_vol)")
    else
        # If calc_cyl wasn't added to phantoms_dict, we need to create it here
        x_cut_angle = params["x_cut_angle"]
        y_cut_angle = params["y_cut_angle"]
        len_cut = params["len_cut"]

        # Calculate the tilted cut plane parameters
        cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)
        cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))
        radius = bigger_cyl_size[1] - cylinder_wall_thickness

        # Calculate maximum drop due to tilt
        max_drop_x = radius * sin(abs(x_cut_angle))
        max_drop_y = radius * sin(abs(y_cut_angle))
        max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

        # Calculate lowest point of cut plane
        lowest_cut_z = cut_center_z - max_drop - len_cut / 2
        distance_from_top = cylinder_top_z - lowest_cut_z
        cylinder_height_below_cut = bigger_cyl_size[3] - distance_from_top

        # Create the cylinder object representing fluid below the cut
        calc_cyl = cylinder(
            center_cylinder[1]cm, center_cylinder[2]cm,
            (center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2)cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            cylinder_height_below_cut * cm,
            angle, 0, 0, 1.0f0
        )

        # Calculate the analytical volume
        cylinder_below_cut_vol = IP.volume(calc_cyl)
        println("Volume of fluid cylinder below cut (created): $(cylinder_below_cut_vol)")

        # Create the boolean mask
        calc_cyl_mask = phantom(axes(ig)..., [calc_cyl]) .!= 0
    end

    # Save calc_cyl mask as a NIFTI file for debugging
    sitk = pyimport("SimpleITK")
    calc_cyl_uint8 = UInt8.(calc_cyl_mask) .* 255
    calc_cyl_img = sitk.GetImageFromArray(calc_cyl_uint8)
    sitk.WriteImage(calc_cyl_img, joinpath(output_folder, "analytical_bottom_cyl.nii.gz"))

    # Calculate difference between fluid_mask and calc_cyl
    # Only consider the part above the cylinder center
    cylinder_center_z = Int(round(center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2))

    # Start with the difference mask
    diff_mask = fluid_mask .& .!calc_cyl_mask

    # Also subtract inner_torus if present as requested
    if rounded_bottom && haskey(phantoms_dict, "inner_torus")
        diff_mask .&= .!phantoms_dict["inner_torus"][1]
    end

    # Keep only parts above the cylinder center
    for i in 1:size(diff_mask, 1)
        for j in 1:size(diff_mask, 2)
            for k in 1:size(diff_mask, 3)
                if k < cylinder_center_z && diff_mask[i, j, k]
                    diff_mask[i, j, k] = false
                end
            end
        end
    end

    # Calculate the volume of this difference and convert to cubic centimeters
    diff_vol = sum(diff_mask) * voxel_volume
    println("Volume of fluid above cylinder center: $(diff_vol) cm³")

    # Save the difference mask for debugging
    diff_uint8 = UInt8.(diff_mask) .* 255
    diff_img = sitk.GetImageFromArray(diff_uint8)
    sitk.WriteImage(diff_img, joinpath(output_folder, "fluid_cylinder_diff.nii.gz"))

    # Update main_volume with cylinder and difference volumes
    main_volume = cylinder_below_cut_vol + (diff_vol * cm^3)
    println("Adjusted fluid volume after top cut: $(main_volume)")

    # STEP 4: Handle bottom geometry based on rounded_bottom
    if rounded_bottom
        # Case A: Rounded bottom
        inner_torus_volume_analytical = safe_get(params, "inner_torus_volume_analytical", 0.0)
        outer_sphere_volume = safe_get(params, "outer_sphere_volume", 0.0)

        # Add half inner torus
        if inner_torus_volume_analytical > 0
            inner_torus_vol_cm3 = inner_torus_volume_analytical * cm^3
            println("Adding half inner torus: $(inner_torus_vol_cm3/2)")
            main_volume += inner_torus_vol_cm3 / 2
        end

        # Subtract half outer sphere
        if outer_sphere_volume > 0
            outer_sphere_vol_cm3 = outer_sphere_volume * cm^3
            println("Subtracting half outer sphere: $(outer_sphere_vol_cm3/2)")
            main_volume -= outer_sphere_vol_cm3 / 2
        end
    else
        # Case B: Non-rounded bottom
        if haskey(phantoms_dict, "ob5") && haskey(phantoms_dict, "center_bottom_out")
            # Get the actual objects
            ob5_obj = phantoms_dict["ob5"][3]
            center_bottom_obj = phantoms_dict["center_bottom_out"][3]

            # Calculate analytical volumes
            ob5_vol = IP.volume(ob5_obj)
            center_bottom_vol_full = IP.volume(center_bottom_obj)

            # Get boolean masks
            ob5_mask = phantoms_dict["ob5"][1]
            center_bottom_mask = phantoms_dict["center_bottom_out"][1]

            # Save center_bottom_out mask for debugging
            center_bottom_uint8 = UInt8.(center_bottom_mask) .* 255
            center_bottom_img = sitk.GetImageFromArray(center_bottom_uint8)
            sitk.WriteImage(center_bottom_img, joinpath(output_folder, "center_bottom_out.nii.gz"))

            # Calculate intersection volume using boolean operations
            intersection_mask = ob5_mask .& center_bottom_mask
            intersection_vol = sum(intersection_mask) * voxel_volume * cm^3

            # Calculate the unique part of center_bottom_out
            center_bottom_unique_vol = center_bottom_vol_full - intersection_vol

            println("Volume of ob5: $(ob5_vol)")
            println("Volume of center_bottom_out: $(center_bottom_vol_full)")
            println("Volume of intersection: $(intersection_vol)")
            println("Unique volume of center_bottom_out: $(center_bottom_unique_vol)")

            # Subtract both volumes
            main_volume -= ob5_vol
            main_volume -= center_bottom_unique_vol

            println("Adjusted fluid volume after bottom subtract: $(main_volume)")
        end
    end

    # STEP 5: Common adjustments for both bottom types
    # Subtract ball volumes if present
    if first_ball && haskey(phantoms_dict, "ball1")
        ball1_volume = IP.volume(phantoms_dict["ball1"][3])
        println("Subtracting ball1: $(ball1_volume)")
        main_volume -= ball1_volume
    end

    if second_ball && haskey(phantoms_dict, "ball2")
        ball2_volume = IP.volume(phantoms_dict["ball2"][3])
        println("Subtracting ball2: $(ball2_volume)")
        main_volume -= ball2_volume
    end

    # Apply additional adjustments for both A and B cases
    # Subtract outer_sphere_volume/2 and add inner_torus_volume_analytical/2
    if !rounded_bottom
        # These should only be applied if they weren't already
        inner_torus_volume_analytical = safe_get(params, "inner_torus_volume_analytical", 0.0)
        outer_sphere_volume = safe_get(params, "outer_sphere_volume", 0.0)

        if outer_sphere_volume > 0
            outer_sphere_vol_cm3 = outer_sphere_volume * cm^3
            println("Additional: subtracting half outer sphere: $(outer_sphere_vol_cm3/2)")
            main_volume -= outer_sphere_vol_cm3 / 2
        end

        if inner_torus_volume_analytical > 0
            inner_torus_vol_cm3 = inner_torus_volume_analytical * cm^3
            println("Additional: adding half inner torus: $(inner_torus_vol_cm3/2)")
            main_volume += inner_torus_vol_cm3 / 2
        end
    end

    # STEP 6: Calculate pipe overlap with calc_cyl
    if haskey(phantoms_dict, "pipe")
        pipe_mask = phantoms_dict["pipe"][1]
        pipe_calc_cyl_overlap = pipe_mask .& calc_cyl_mask
        pipe_overlap_vol = sum(pipe_calc_cyl_overlap) * voxel_volume * cm^3

        println("Pipe overlap with fluid: $(pipe_overlap_vol)")
        main_volume -= pipe_overlap_vol
    end

    # Convert to scalar value (strip units)
    analytical_vol = Unitful.ustrip(cm^3, main_volume)

    # Print final comparison
    println("Final analytical volume: $(main_volume) ($(analytical_vol) cm³)")
    println("Final numerical volume: $(numerical_vol) cm³")
    println("Absolute difference: $(abs(analytical_vol - numerical_vol)) cm³")

    if numerical_vol > 0
        rel_diff = 100.0 * abs(analytical_vol - numerical_vol) / numerical_vol
        println("Relative difference: $(rel_diff)%")
    end

    return numerical_vol, analytical_vol
end

"""
    get_cylinder_bool_mask(ob_el, ig::ImageGeom) -> BitArray{3}

Create a boolean mask from a cylinder object using the proper masking technique.

# Arguments
- `ob_el`: Cylinder object or array of objects
- `ig`: ImageGeom defining the phantom dimensions and spacing

# Returns
A 3D boolean array representing the cylinder volume
"""
function get_cylinder_bool_mask(ob_el, ig::ImageGeom)
    pp = phantom(axes(ig)..., ob_el)
    pp_reversed_1 = reverse(pp, dims=1)
    cyl_inner_bool = (pp + pp_reversed_1) .!= 0
    return cyl_inner_bool
end

"""
    compute_accurate_fluid_volume_fixed(
        phantoms_dict::Dict,
        spacing::Tuple{Float64,Float64,Float64},
        first_ball::Bool,
        second_ball::Bool,
        rounded_bottom::Bool,
        params::Dict{String,Any},
        ig::ImageGeom
    ) -> Tuple{Float64,Float64}

Calculate fluid volume in a phantom can with robust handling of complex geometries.
Fixed version that addresses issues with cylinder_height_below_cut variable.

# Arguments
- `phantoms_dict`: Dictionary containing phantom objects and their masks
- `spacing`: Voxel spacing (dx, dy, dz) in cm
- `first_ball`: Whether first ball is present
- `second_ball`: Whether second ball is present
- `rounded_bottom`: Whether can has rounded bottom
- `params`: Dictionary of parameters
- `ig`: ImageGeom object for phantom generation

# Returns
A tuple containing (numerical_volume, analytical_volume) in cubic centimeters
"""
function compute_accurate_fluid_volume_fixed(
    phantoms_dict::Dict,
    spacing::Tuple{Float64,Float64,Float64},
    first_ball::Bool,
    second_ball::Bool,
    rounded_bottom::Bool,
    params::Dict{String,Any},
    ig::ImageGeom
)
    # Safe parameter extraction with defaults
    function safe_get(dict, key, default=0.0)
        return haskey(dict, key) ? dict[key] : default
    end

    # Extract needed parameters
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    output_folder = safe_get(params, "output_folder", ".")
    dual_phase_percentage = safe_get(params, "dual_phase_percentage", 1.0)
    center_cylinder = params["center_cylinder"]
    bigger_cyl_size = params["bigger_cyl_size"]
    cylinder_wall_thickness = params["cylinder_wall_thickness"]
    len_cut = params["len_cut"]
    x_cut_angle = params["x_cut_angle"]
    y_cut_angle = params["y_cut_angle"]
    angle = safe_get(params, "angle", 0.0)

    # NUMERICAL VOLUME CALCULATION - Use fluid_mask directly
    fluid_mask = params["fluid_mask"]
    numerical_vol = sum(fluid_mask) * voxel_volume
    println("Numerical volume from fluid_mask: $(numerical_vol) cm³")

    # ANALYTICAL VOLUME - Base calculation depending on fluid phases
    if dual_phase_percentage == 1.0
        main_volume = IP.volume(phantoms_dict["ob2_a"][3])
        println("Base fluid volume (single phase): $(main_volume)")
    else
        main_volume = IP.volume(phantoms_dict["ob2_a"][3]) + IP.volume(phantoms_dict["ob2_b"][3])
        println("Base fluid volume (dual phases): $(main_volume)")
    end

    # Handle the tilted top cut geometry
    local cylinder_height_below_cut, calc_cyl_mask, calc_cyl

    # Retrieve or calculate the calc_cyl cylinder that represents fluid below the cut
    if haskey(phantoms_dict, "calc_cyl")
        # If calc_cyl was created in empty_cylinder_with_half_sphere_bottom_p
        println("Using pre-calculated calc_cyl from phantoms_dict")
        calc_cyl = phantoms_dict["calc_cyl"][3]

        # Calculate cylinder parameters based on object properties
        cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)
        cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))
        radius = bigger_cyl_size[1] - cylinder_wall_thickness

        # Calculate max drop due to tilt angles
        max_drop_x = radius * sin(abs(x_cut_angle))
        max_drop_y = radius * sin(abs(y_cut_angle))
        max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

        # Calculate lowest point and height
        lowest_cut_z = cut_center_z - max_drop - len_cut / 2
        distance_from_top = cylinder_top_z - lowest_cut_z
        cylinder_height_below_cut = bigger_cyl_size[3] - distance_from_top

        # Get boolean mask using the correct function
        calc_cyl_mask = get_cylinder_bool_mask([calc_cyl], ig)
    else
        # If we need to create the calc_cyl here
        println("Creating calc_cyl cylinder")

        # Calculate geometry for the tilted cut
        cylinder_top_z = center_cylinder[3] + (bigger_cyl_size[3] / 2)
        cut_center_z = center_cylinder[3] + ((bigger_cyl_size[3] / 2) - (len_cut / 2))
        radius = bigger_cyl_size[1] - cylinder_wall_thickness

        # Calculate drop due to tilt angles
        max_drop_x = radius * sin(abs(x_cut_angle))
        max_drop_y = radius * sin(abs(y_cut_angle))
        max_drop = sqrt(max_drop_x^2 + max_drop_y^2)

        # Calculate lowest point of cut and cylinder height
        lowest_cut_z = cut_center_z - max_drop - len_cut / 2
        distance_from_top = cylinder_top_z - lowest_cut_z
        cylinder_height_below_cut = bigger_cyl_size[3] - distance_from_top

        # Create the cylinder representing fluid below cut
        calc_cyl = cylinder(
            center_cylinder[1]cm, center_cylinder[2]cm,
            (center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2)cm,
            (bigger_cyl_size[1] - cylinder_wall_thickness)cm,
            (bigger_cyl_size[2] - cylinder_wall_thickness)cm,
            cylinder_height_below_cut * cm,
            angle, 0, 0, 1.0f0
        )

        # Get boolean mask directly from phantom function
        calc_cyl_mask = phantom(axes(ig)..., [calc_cyl]) .!= 0
    end

    # Calculate volume of the cylinder below cut
    cylinder_below_cut_vol = IP.volume(calc_cyl)
    println("Volume of fluid cylinder below cut: $(cylinder_below_cut_vol)")

    # Save the calculated cylinder mask as a NIFTI file
    sitk = pyimport("SimpleITK")
    calc_cyl_uint8 = UInt8.(calc_cyl_mask) .* 255
    calc_cyl_img = sitk.GetImageFromArray(calc_cyl_uint8)
    sitk.WriteImage(calc_cyl_img, joinpath(output_folder, "analytical_bottom_cyl.nii.gz"))

    # Calculate the cylinder center z-position for filtering
    cylinder_center_z = Int(round(center_cylinder[3] - (bigger_cyl_size[3] / 2) + cylinder_height_below_cut / 2))

    # Calculate difference between fluid mask and calc_cyl mask
    diff_mask = fluid_mask .& .!calc_cyl_mask

    # Remove inner_torus from the difference if present (as requested)
    if rounded_bottom && haskey(phantoms_dict, "inner_torus")
        diff_mask .&= .!phantoms_dict["inner_torus"][1]
    end

    # Keep only the part of the difference above the cylinder center
    for i in 1:size(diff_mask, 1)
        for j in 1:size(diff_mask, 2)
            for k in 1:size(diff_mask, 3)
                if k < cylinder_center_z && diff_mask[i, j, k]
                    diff_mask[i, j, k] = false
                end
            end
        end
    end

    # Calculate volume of the difference
    diff_vol = sum(diff_mask) * voxel_volume
    println("Volume of fluid above cylinder center: $(diff_vol) cm³")

    # Save the difference mask for debugging
    diff_uint8 = UInt8.(diff_mask) .* 255
    diff_img = sitk.GetImageFromArray(diff_uint8)
    sitk.WriteImage(diff_img, joinpath(output_folder, "fluid_cylinder_diff.nii.gz"))

    # Update main volume with cylinder and difference volumes
    main_volume = cylinder_below_cut_vol + (diff_vol * cm^3)
    println("Volume after adjusting for top cut: $(main_volume)")

    # Handle bottom geometry based on rounded_bottom parameter
    inner_torus_volume_analytical = safe_get(params, "inner_torus_volume_analytical", 0.0)
    outer_sphere_volume = safe_get(params, "outer_sphere_volume", 0.0)

    if rounded_bottom
        # Case A: Rounded bottom
        # Add half inner torus volume
        if inner_torus_volume_analytical > 0.0
            inner_torus_cm3 = inner_torus_volume_analytical * cm^3
            println("Adding half inner torus: $(inner_torus_cm3/2)")
            main_volume += inner_torus_cm3 / 2
        end

        # Subtract half outer sphere volume
        if outer_sphere_volume > 0.0
            outer_sphere_cm3 = outer_sphere_volume * cm^3
            println("Subtracting half outer sphere: $(outer_sphere_cm3/2)")
            main_volume -= outer_sphere_cm3 / 2
        end
    else
        # Case B: Non-rounded bottom with ob5 and center_bottom_out
        if haskey(phantoms_dict, "ob5") && haskey(phantoms_dict, "center_bottom_out")
            # Get the geometric objects
            ob5_obj = phantoms_dict["ob5"][3]
            center_bottom_obj = phantoms_dict["center_bottom_out"][3]

            # Calculate volumes using IP.volume()
            ob5_vol = IP.volume(ob5_obj)
            center_bottom_vol_full = IP.volume(center_bottom_obj)

            # Get boolean masks
            ob5_mask = phantoms_dict["ob5"][1]
            center_bottom_mask = phantoms_dict["center_bottom_out"][1]


            # Calculate intersection volume using mask
            intersection_mask = (.!ob5_mask) .& center_bottom_mask

            intersection_vol = sum(intersection_mask) * voxel_volume * cm^3


            # Save center_bottom_out mask as requested
            center_bottom_uint8 = UInt8.(intersection_mask)
            center_bottom_img = sitk.GetImageFromArray(center_bottom_uint8)
            sitk.WriteImage(center_bottom_img, joinpath(output_folder, "intersection_bottom_mask.nii.gz"))


            # Calculate unique part of center_bottom_out not in ob5
            center_bottom_unique_vol = center_bottom_vol_full - intersection_vol

            println("Volume of ob5: $(ob5_vol)")
            println("Volume of center_bottom_out: $(center_bottom_vol_full)")
            println("Intersection volume: $(intersection_vol)")
            println("Unique center_bottom volume: $(center_bottom_unique_vol)")

            # Subtract both volumes as required
            main_volume -= ob5_vol
            main_volume -= center_bottom_unique_vol
        end
    end

    # For both A and B cases: Handle ball volumes and pipe overlap

    # Subtract ball volumes if present
    if first_ball && haskey(phantoms_dict, "ball1")
        ball1_vol = IP.volume(phantoms_dict["ball1"][3])
        println("Subtracting ball1: $(ball1_vol)")
        main_volume -= ball1_vol
    end

    if second_ball && haskey(phantoms_dict, "ball2")
        ball2_vol = IP.volume(phantoms_dict["ball2"][3])
        println("Subtracting ball2: $(ball2_vol)")
        main_volume -= ball2_vol
    end

    # For both A and B: subtract outer_sphere/2 and add inner_torus/2 if not already done
    if !rounded_bottom
        # Apply the additional adjustments
        if outer_sphere_volume > 0.0
            outer_sphere_cm3 = outer_sphere_volume * cm^3
            println("Additional: subtracting half outer sphere: $(outer_sphere_cm3/2)")
            main_volume -= outer_sphere_cm3 / 2
        end

        if inner_torus_volume_analytical > 0.0
            inner_torus_cm3 = inner_torus_volume_analytical * cm^3
            println("Additional: adding half inner torus: $(inner_torus_cm3/2)")
            main_volume += inner_torus_cm3 / 2
        end
    end

    # Calculate overlap between pipe and calc_cyl
    if haskey(phantoms_dict, "pipe")
        # pipe_mask = phantoms_dict["pipe"][1]
        # pipe_calc_cyl_intersection = pipe_mask .& calc_cyl_mask
        # pipe_overlap_vol = sum(pipe_calc_cyl_intersection) * voxel_volume

        # println("Pipe overlap with fluid: $(pipe_overlap_vol) cm³")
        # main_volume -= pipe_overlap_vol * cm^3
        main_volume -= IP.volume(phantoms_dict["outside_of_pipe_in_fluid"][3])
        main_volume += IP.volume(phantoms_dict["inside_of_pipe_in_fluid"][3])

    end
    

    # Convert main_volume to scalar (removing units)
    analytical_vol = Unitful.ustrip(cm^3, main_volume)

    # Print final results
    println("Final analytical volume: $(analytical_vol) cm³")
    println("Final numerical volume: $(numerical_vol) cm³")
    println("Difference: $(abs(analytical_vol - numerical_vol)) cm³")

    if numerical_vol > 0
        rel_diff = 100.0 * abs(analytical_vol - numerical_vol) / numerical_vol
        println("Relative difference: $(rel_diff)%")
    end

    return numerical_vol, analytical_vol
end


