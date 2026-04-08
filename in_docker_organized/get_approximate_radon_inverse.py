import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
import adrt
import SimpleITK as sitk
import sys
from skimage.transform import radon
from skimage.transform import iradon
from scipy import interpolate


class ADRTNormalOperator(LinearOperator):
    def __init__(self, img_size, dtype=None):
        super().__init__(dtype=dtype, shape=(img_size**2, img_size**2))
        self._img_size = img_size

    def _matmat(self, x):
        # Use batch dimensions to handle columns of matrix x
        n_batch = x.shape[-1]
        batch_img = np.moveaxis(x, -1, 0).reshape(
            (n_batch, self._img_size, self._img_size)
        )
        ret = adrt.utils.truncate(adrt.bdrt(adrt.adrt(batch_img))).mean(axis=1)
        return np.moveaxis(ret, 0, -1).reshape((self._img_size**2, n_batch))

    def _adjoint(self):
        return self


def iadrt_cg(b, /, *, op_cls=ADRTNormalOperator, **kwargs):
    if b.ndim > 3:
        raise ValueError("batch dimension not supported for iadrt_cg")
    img_size = b.shape[-1]
    linop = op_cls(img_size=img_size, dtype=b.dtype)
    tb = adrt.utils.truncate(adrt.bdrt(b)).mean(axis=0).ravel()
    x, info = cg(linop, tb, **kwargs)
    if info != 0:
        raise ValueError(f"convergence failed (cg status {info})")
    return x.reshape((img_size, img_size))


def next_power_of_two(n):
    """Return the next power of 2 greater than or equal to n."""
    return 1 if n == 0 else 2**int(np.ceil(np.log2(n)))

def pad_to_power_of_two(arr):
    """Pad each dimension of a 3D array to the next power of two."""
    d1, d2, d3 = arr.shape
    new_d1 = next_power_of_two(d1)
    new_d2 = next_power_of_two(d2)
    new_d3 = next_power_of_two(d3)
    
    if (d1, d2, d3) == (new_d1, new_d2, new_d3):
        return arr  # Array is already padded to powers of 2
    corners = [
        arr[0, 0, 0],
        arr[0, 0, -1],
        arr[0, -1, 0],
        arr[0, -1, -1],
        arr[-1, 0, 0],
        arr[-1, 0, -1],
        arr[-1, -1, 0],
        arr[-1, -1, -1],
    ]
    mean_corners = np.mean(corners)
    padded_arr = np.zeros((new_d1, new_d2, new_d3), dtype=arr.dtype) + mean_corners
    
    padded_arr[:d1, :d2, :d3] = arr
    return padded_arr

# New main block to allow command-line execution with two arguments
if __name__ == "__main__":
    # if len(sys.argv) != 3:
    #     print("Usage: python get_approximate_radon_inverse.py <input_nifti> <output_nifti>")
    #     sys.exit(1)
    input_path, output_path, output_path_b  = sys.argv[1], sys.argv[2], sys.argv[3]
    # Load input NIfTI as array
    input_img = sitk.ReadImage(input_path)
    arr = sitk.GetArrayFromImage(input_img)
    arr_orig=np.copy(arr)

    # Ensure that dimensions are powers of 2
    sizz = arr.shape
    arr=pad_to_power_of_two(arr)
    # Preallocate output with the same shape
    out = np.empty_like(arr)
    
    # Process each slice (first axis)
    for z in range(arr.shape[0]):
        img_slice = arr[z, :, :]
        
        n = img_slice.shape[1]
        # img_slice = adrt.adrt(img_slice)
        # img_slice = iadrt_cg(img_slice)
        # noise = np.random.normal(loc=0.0, scale=1.0, size=img_slice.shape)
        # img_slice += noise
        # img_slice = iadrt_cg(img_slice)
        
        th_array1 = np.unique(adrt.utils.coord_adrt(n).angle)
        theta = 90.0 + np.rad2deg(th_array1.squeeze())
        sinogram = radon(img_slice, theta=theta,circle=False)
        # t_array = np.linspace(-0.5, 0.5, n)
        # spline = interpolate.RectBivariateSpline(t_array, th_array1, sinogram)
        # s_array, th_array = adrt.utils.coord_adrt(n)
        # adrt_data = spline(s_array, th_array, grid=False)

        
        # # img_slice=adrt.iadrt_fmg(adrt_data)
        # img_slice=iadrt_cg(adrt_data)
        img_slice=iradon(sinogram, theta=theta, circle=False)

        # ...any processing on img_slice...
        out[z, :, :] = img_slice  # overwrite corresponding slice in the output
    # Create output image and copy metadata
    if out.shape != sizz:
        out = out[:sizz[0], :sizz[1], :sizz[2]]
        
    print(f"\n out {out.shape} sizz {sizz}  \n")
    out_img = sitk.GetImageFromArray(out)
    out_img.CopyInformation(input_img)
    
    sitk.WriteImage(out_img, output_path)
    
    
    combined_im=(arr_orig+out)/2
    combined_im = sitk.GetImageFromArray(combined_im)
    combined_im.CopyInformation(combined_im)
    sitk.WriteImage(combined_im, output_path_b)
        
    
    
    
# python3 /workspaces/CanIonSynth/src/organised/get_approximate_radon_inverse.py /workspaces/CanIonSynth/data/debug/25de96dc-bfb6-4022-be00-5df965b49283/example_can.nii.gz /workspaces/CanIonSynth/data/debug/25de96dc-bfb6-4022-be00-5df965b49283/after_radon.nii.gz /workspaces/CanIonSynth/data/debug/25de96dc-bfb6-4022-be00-5df965b49283/after_radon_plus_before.nii.gz

#/workspaces/CanIonSynth# python3 /workspaces/CanIonSynth/src/organised/get_approximate_radon_inverse.py "/workspaces/CanIonSynth/data/debug/0ae4e1be-5220-4c57-85a5-aaa0afe177c9/example_can.nii.gz" "/workspaces/CanIonSynth/data/debug/iadrt_test.nii.gz"
