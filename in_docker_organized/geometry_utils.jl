
# python implementation taken from https://stackoverflow.com/questions/56171643/simpleitk-rotation-of-volumetric-data-e-g-mri
function matrix_from_axis_angle(a)
    """ Compute rotation matrix from axis-angle.
    This is called exponential map or Rodrigues' formula.
    Parameters
    ----------
    a : array-like, shape (4,)
        Axis of rotation and rotation angle: (x, y, z, angle)
    Returns
    -------
    R : array-like, shape (3, 3)
        Rotation matrix
    """
    ux, uy, uz, theta = a
    c = cos(theta)
    s = sin(theta)
    ci = 1.0 - c
    R = [[ci * ux * ux + c,
            ci * ux * uy - uz * s,
            ci * ux * uz + uy * s],
        [ci * uy * ux + uz * s,
            ci * uy * uy + c,
            ci * uy * uz - ux * s],
        [ci * uz * ux - uy * s,
            ci * uz * uy + ux * s,
            ci * uz * uz + c],
    ]

    # This is equivalent to
    # R = (np.eye(3) * np.cos(theta) +
    #      (1.0 - np.cos(theta)) * a[:3, np.newaxis].dot(a[np.newaxis, :3]) +
    #      cross_product_matrix(a[:3]) * np.sin(theta))

    return R
end #matrix_from_axis_angle

function resample(image, transform)
    """
    This function resamples (updates) an image using a specified transform
    :param image: The sitk image we are trying to transform
    :param transform: An sitk transform (ex. resizing, rotation, etc.
    :return: The transformed sitk image
    """
    sitk = pyimport("SimpleITK")
    np = pyimport("numpy")

    reference_image = image
    interpolator = sitk.sitkLinear
    default_value = 0
    return sitk.Resample(image, reference_image, transform,
        interpolator, default_value)
end#resample

function get_center(img)
    """
    from python to test
    """
    width, height, depth = img.GetSize()
    centt = (Int(ceil(width / 2)), Int(ceil(height / 2)), Int(ceil(depth / 2)))
    # return img.TransformIndexToPhysicalPoint((np.ceil(width/2), np.ceil(height/2), np.ceil(depth/2)))
    return img.TransformIndexToPhysicalPoint(centt)
end #get_center


function rotation3d(image, axis, theta)
    """
    This function rotates an image across each of the x, y, z axes by theta_x, theta_y, and theta_z degrees
    respectively
    :return: The rotated image
    """
    sitk = pyimport("SimpleITK")
    np = pyimport("numpy")

    theta = np.deg2rad(theta)
    euler_transform = sitk.Euler3DTransform()
    image_center = get_center(image)
    euler_transform.SetCenter(image_center)

    direction = image.GetDirection()

    if (axis == 3)
        axis_angle = (direction[3], direction[6], direction[9], theta)
    elseif (axis == 2)
        axis_angle = (direction[2], direction[5], direction[8], theta)
    elseif (axis == 1)
        axis_angle = (direction[1], direction[4], direction[7], theta)
    end
    np_rot_mat = matrix_from_axis_angle(axis_angle)
    euler_transform.SetMatrix([np_rot_mat[1][1], np_rot_mat[1][2], np_rot_mat[1][3], np_rot_mat[2][1], np_rot_mat[2][2], np_rot_mat[2][3], np_rot_mat[3][1], np_rot_mat[3][2], np_rot_mat[3][3]])

    # if(axis==3)
    #     axis_angle = (direction[2], direction[5], direction[8], theta)
    # elseif (axis==2)
    #     axis_angle = (direction[1], direction[4], direction[7], theta)
    # elseif (axis==1)
    #     axis_angle = (direction[0], direction[3], direction[6], theta)
    # end
    # np_rot_mat = matrix_from_axis_angle(axis_angle)
    # euler_transform.SetMatrix([np_rot_mat[0][0],np_rot_mat[0][1],np_rot_mat[0][2]
    #                             ,np_rot_mat[1][0],np_rot_mat[1][1],np_rot_mat[1][2]
    #                             ,np_rot_mat[2][0],np_rot_mat[2][1],np_rot_mat[2][2] ])
    resampled_image = resample(image, euler_transform)
    return resampled_image
end #rotation3d


"""
    save_sitk_image_as_dicom(img::SimpleITK.Image, output_folder::String)

Saves a given SimpleITK image to a temporary NIfTI file, creates the specified
output folder if it does not exist, then invokes the Python DICOM converter script
(e.g. save_dicom.py) via the command line, and finally removes the temporary file.
"""
function save_sitk_image_as_dicom(img, output_folder::String)

    arr = sitk.GetArrayFromImage(img)
    arr = arr .- minimum(arr)
    arr = arr ./ maximum(arr)
    arr = arr .* 1000
    arr = UInt16.(round.(arr))
    print("\n savinggggg $(size(arr)) \n")
    # arr=UInt32.(round.(arr))
    is_two_dim=0
    if ndims(arr) == 2
        is_two_dim=1
        # arr=reshape(arr, (size(arr)...,1))
        # arr=cat(arr,arr,arr,dims=3)
        print("\n two dimm  $(size(arr))\n ")
    end

    # copy spatial metadata
    img_origin = img.GetOrigin()
    img_spacing = img.GetSpacing()
    img_direction = img.GetDirection()

    new_img = sitk.GetImageFromArray(arr)
    new_img.SetOrigin(img_origin)
    new_img.SetSpacing(img_spacing)
    new_img.SetDirection(img_direction)

    # 1) Create a unique temporary file
    timestamp = Dates.format(now(), "yyyyMMdd_HHMMSS")
    rnd_str = string(rand(1000:9999))
    tmp_nifti_path = "/tmp/temp_dicom_" * timestamp * "_" * rnd_str * ".nii"

    # 2) Write the image as a NIfTI file
    sitk.WriteImage(new_img, tmp_nifti_path)

    # 3) Create the output folder if needed
    isdir(output_folder) || mkpath(output_folder)

    # 4) Construct and run the Python command
    cmd = `nii2dcm $tmp_nifti_path $output_folder -d MR` # -d MR

    try
        if Sys.which("nii2dcm") !== nothing
            run(cmd)
        else
            println("Warning: nii2dcm not found in path. Skipping DICOM conversion.")
        end
    catch e
        println("Warning: Failed to run nii2dcm: $e")
    end

    # 5) Remove the temporary NIfTI
    rm(tmp_nifti_path, force=true)
end


"""
    get_per_slice_reconstruction_parallel(arr)::Array{Float32,3}

Reconstructs each slice of `arr` in parallel using radon/iradon transforms and
stores the result in a preallocated 3D array of Float32 with the same shape as `arr`.
"""
function get_per_slice_reconstruction(arr)
    # Dimensions
    nr, nc, nz = size(arr)
    # Preallocate result array
    recon_arr = Array{Float32}(undef, nr, nc, nz)

    # The Python radon/iradon functions
    # are presumably imported as skimage.transform
    Threads.@threads for z in 1:nr
        # for z in 1:nr
        # slice_2d = arr[z, :, :]

        theta = collect(1:256)
        radon_slice = skimage.transform.radon(arr[z, :, :], theta=theta, circle=true)
        slice_recon = skimage.transform.iradon(radon_slice, theta=theta, circle=true)

        # Ensure the result slice has the same 2D shape as the original slice
        # Convert to Float32
        recon_arr[z, :, :] = Float32.(slice_recon)
    end

    return recon_arr
end


"""
    save_2d_sitk_image_as_dicom(immm, output_path)

Save a SimpleITK 2D image as a DICOM file with proper intensity scaling.

# Arguments
- `immm`: SimpleITK Image to save
- `output_path`: Output path for the DICOM file (can be a file or directory)
"""
function save_2d_sitk_image_as_dicom(immm, output_path::String)
    """
    Given a 2D SimpleITK.Image (immm), save it into a DICOM file at `output_path`.
    """
    # Check if the output_path is a directory
    if isdir(output_path)
        # If it's a directory, append a default filename
        filename = joinpath(output_path, "slice0000.dcm")
    else
        filename = output_path
    end

    # Handle floating point data by converting to int16
    pixel_id = immm.GetPixelID()
    if pixel_id == sitk.sitkFloat32 || pixel_id == sitk.sitkFloat64
        # Convert float to int16 with rescaling for better contrast
        stats = sitk.StatisticsImageFilter()
        stats.Execute(immm)
        min_val = stats.GetMinimum()
        max_val = stats.GetMaximum()

        # Rescale to use the full int16 range
        output_min = 0
        output_max = 32767  # Max value for int16

        # Create the rescale filter
        rescale_filter = sitk.RescaleIntensityImageFilter()
        rescale_filter.SetOutputMinimum(output_min)
        rescale_filter.SetOutputMaximum(output_max)

        # Apply the filter
        rescaled_img = rescale_filter.Execute(immm)

        # Cast to int16
        cast_filter = sitk.CastImageFilter()
        cast_filter.SetOutputPixelType(sitk.sitkInt16)
        img_to_save = cast_filter.Execute(rescaled_img)

        # Set DICOM-specific metadata for proper display
        img_to_save.SetMetaData("0028|1052", "0")  # Rescale Intercept
        img_to_save.SetMetaData("0028|1053", "1")  # Rescale Slope
    else
        img_to_save = immm
    end

    writer = sitk.ImageFileWriter()
    writer.SetFileName(filename)
    writer.SetImageIO("GDCMImageIO")
    writer.Execute(img_to_save)
end

function save_nifti_with_meta(arr, cast_to_uint8::Bool, spacing::Tuple{Float64,Float64,Float64}, output_path::String)
    # Convert to SimpleITK image, cast if requested
    img = cast_to_uint8 ? sitk.GetImageFromArray(UInt8.(arr)) : sitk.GetImageFromArray(arr)

    # Set spacing in mm (from cm)
    img.SetSpacing((spacing[1] * 10, spacing[2] * 10, spacing[3] * 10))
    # Always fix direction & origin
    img.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    img.SetOrigin((0.0, 0.0, 0.0))

    # Write the image
    sitk.WriteImage(img, output_path)
end


"""
    save_mask_as_nifti(mask::Array{Bool,3}, output_path::String, spacing::Tuple{Float64,Float64,Float64})

Save a boolean mask as a NIfTI file with the specified spacing.

# Arguments
- `mask`: Boolean 3D array to save
- `output_path`: Path to save the NIfTI file
- `spacing`: Voxel spacing (dx, dy, dz) in cm
"""
function save_mask_as_nifti(mask::Array{Bool,3}, output_path::String, spacing::Tuple{Float64,Float64,Float64})
    # Convert boolean mask to UInt8
    save_nifti_with_meta(UInt8.(mask), true, spacing, output_path)
    println("Saved mask to: $(output_path)")
end

"""
    convert_nifti_to_dicom_seg(nifti_path::String, reference_dicom_path::String, output_folder::String, reference_nifti_path::String="")

Convert a NIfTI binary mask to DICOM-SEG format using the Python script.

# Arguments
- `nifti_path`: Path to the NIfTI mask file
- `reference_dicom_path`: Path to a reference DICOM file or folder containing DICOM files
- `output_folder`: Folder to save the DICOM-SEG file
- `reference_nifti_path`: Optional path to a reference NIfTI file with matching geometry
"""
function convert_nifti_to_dicom_seg(nifti_path::String, reference_dicom_path::String, output_folder::String, reference_nifti_path::String="")

    if !isdir(output_folder)
        mkpath(output_folder)
    end

    # Ensure nifti_to_dicom_seg.py is in the correct path
    # script_path = "/workspaces/synthethic_tomo/src/organised/nifti_to_dicom_seg.py"
    script_path = joinpath(@__DIR__, "nifti_to_dicom_seg.py")

    if !isfile(script_path)
        println("Warning: nifti_to_dicom_seg.py script not found at $(script_path)")
        return false
    end

    # Run the Python script with additional reference_nifti argument if provided
    cmd = if isempty(reference_nifti_path)
        `python3 $(script_path) $(nifti_path) $(reference_dicom_path) $(output_folder)`
    else
        `python3 $(script_path) $(nifti_path) $(reference_dicom_path) $(output_folder) --reference_nifti $(reference_nifti_path)`
    end

    println("Running command: $(cmd)")

    try
        run(cmd)
        println("Successfully converted NIfTI to DICOM-SEG")
        return true
    catch e
        println("Error converting NIfTI to DICOM-SEG: $(e)")
        return false
    end
end

function rotate_and_retrieve_array(array, axis, theta)
    img = sitk.GetImageFromArray(UInt8.(array))
    rotated_img = rotation3d(img, axis, theta)
    return sitk.GetArrayFromImage(rotated_img).!=0
end

"""
Function to move the image represented as boolean array in given axis up or down
The axis will be given as an argument as well as distance and direction ("up" or "down")
distance will be given in centimiters, spacing in centimiters need also to be taken into account
Goal is to move the image up or down in the given axis do it by removing layers of voxels from one side and adding them to the other side of the given axis
How many layers you will calculate from given distance and spacing for given axis.

 move_image(
        array::Array{Bool,3},
        axis::Int,
        distance::Float64,
        direction::String,
        spacing::Tuple{Float64,Float64,Float64}
    ) -> Array{Bool,3}

Move a binary 3D image along a specified axis by a given distance.

# Arguments
- `array`: Boolean 3D array representing the binary image
- `axis`: Axis along which to move (1=x, 2=y, 3=z)
- `distance`: Distance to move in centimeters
- `direction`: Direction to move, either "up" (positive direction) or "down" (negative direction)
- `spacing`: Voxel spacing in centimeters for each dimension (dx, dy, dz)

# Returns
A new boolean array with the image moved to the new position

# Examples
```julia
# Move image 5mm up along z-axis
moved_array = move_image(array, 3, 0.5, "up", (0.1, 0.1, 0.1))

"""
function move_image( array, axis::Int, distance::Float64, direction::String, spacing::Tuple{Float64,Float64,Float64} )::Array{Bool,3}
    # Check if axis is valid
if !(axis in 1:3) error("Axis must be 1 (x), 2 (y), or 3 (z)") end
# Check if direction is valid
if !(direction in ["up", "down"])
    error("Direction must be either 'up' or 'down'")
end

# Calculate number of voxel layers to shift
layers = round(Int, distance / spacing[axis])
if layers == 0
    return copy(array)
end

# Get array dimensions
dims = size(array)

# Create a new array of the same size, filled with false values
result = fill(false, dims)

# For all axes: "up" shifts in the positive direction, "down" in the negative
if direction == "up"
    # Moving in positive direction
    if axis == 1
        result[(layers+1):end, :, :] = array[1:(end-layers), :, :]
    elseif axis == 2
        result[:, (layers+1):end, :] = array[:, 1:(end-layers), :]
    else  # axis == 3
        result[:, :, (layers+1):end] = array[:, :, 1:(end-layers)]
    end
else  # direction == "down"
    # Moving in negative direction
    if axis == 1
        result[1:(end-layers), :, :] = array[(layers+1):end, :, :]
    elseif axis == 2
        result[:, 1:(end-layers), :] = array[:, (layers+1):end, :]
    else  # axis == 3
        result[:, :, 1:(end-layers)] = array[:, :, (layers+1):end]
    end
end

return result
end