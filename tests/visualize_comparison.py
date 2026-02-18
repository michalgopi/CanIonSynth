import os
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

REF_DIR = "tests/vis_out_orig/Can_Phantom___Default__Random_"
NEW_DIR = "tests/vis_out_ref/Can_Phantom___Default__Random_"
OUTPUT_DIR = "tests/comparisons"
FILENAME = "rotated_example_can_180.nii.gz"

os.makedirs(OUTPUT_DIR, exist_ok=True)

ref_path = os.path.join(REF_DIR, FILENAME)
new_path = os.path.join(NEW_DIR, FILENAME)

print(f"Loading reference: {ref_path}")
img_ref = nib.load(ref_path).get_fdata()

print(f"Loading new: {new_path}")
img_new = nib.load(new_path).get_fdata()

# Select middle slice
if len(img_ref.shape) == 3:
    slice_idx = img_ref.shape[2] // 2
    slice_ref = img_ref[:, :, slice_idx]
    slice_new = img_new[:, :, slice_idx]
else:
    # 2D image
    slice_ref = img_ref
    slice_new = img_new
    slice_idx = "2D"

# Calculate absolute difference
diff = np.abs(slice_ref - slice_new)

print(f"Max difference: {np.max(diff)}")
print(f"Mean difference: {np.mean(diff)}")

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im0 = axes[0].imshow(slice_ref, cmap='gray')
axes[0].set_title("Reference (Original)")
plt.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(slice_new, cmap='gray')
axes[1].set_title("Refactored")
plt.colorbar(im1, ax=axes[1])

im2 = axes[2].imshow(diff, cmap='hot')
axes[2].set_title(f"Difference (Max: {np.max(diff):.6f})")
plt.colorbar(im2, ax=axes[2])

plt.suptitle(f"Comparison of {FILENAME} (Slice {slice_idx})")
plt.tight_layout()

out_file = os.path.join(OUTPUT_DIR, "comparison_rotated_180.png")
plt.savefig(out_file)
print(f"Saved visualization to {out_file}")
