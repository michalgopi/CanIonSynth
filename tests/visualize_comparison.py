import os
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

REF_BASE = "tests/vis_out_orig"
NEW_BASE = "tests/vis_out_ref"
OUTPUT_DIR = "tests/comparisons"
FILENAME = "ionic_chamber.nii.gz"

os.makedirs(OUTPUT_DIR, exist_ok=True)

test_cases = [
    "Ionic_Chamber___Lolipop",
    "Ionic_Chamber___Ball",
    "Ionic_Chamber___Square"
]

for test_case in test_cases:
    print(f"Processing {test_case}...")
    ref_path = os.path.join(REF_BASE, test_case, FILENAME)
    new_path = os.path.join(NEW_BASE, test_case, FILENAME)

    if not os.path.exists(ref_path) or not os.path.exists(new_path):
        print(f"Skipping {test_case}: Files not found.")
        continue

    img_ref = nib.load(ref_path).get_fdata()
    img_new = nib.load(new_path).get_fdata()

    # Select middle slice
    slice_idx = img_ref.shape[2] // 2
    slice_ref = img_ref[:, :, slice_idx]
    slice_new = img_new[:, :, slice_idx]

    # Calculate absolute difference
    diff = np.abs(slice_ref - slice_new)

    print(f"  Max difference: {np.max(diff)}")
    print(f"  Mean difference: {np.mean(diff)}")

    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    im0 = axes[0].imshow(slice_ref, cmap='gray')
    axes[0].set_title(f"Reference ({test_case})")
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(slice_new, cmap='gray')
    axes[1].set_title("Refactored")
    plt.colorbar(im1, ax=axes[1])

    im2 = axes[2].imshow(diff, cmap='hot')
    axes[2].set_title(f"Difference (Max: {np.max(diff):.6f})")
    plt.colorbar(im2, ax=axes[2])

    plt.suptitle(f"Comparison of {test_case} (Slice {slice_idx})")
    plt.tight_layout()

    out_file = os.path.join(OUTPUT_DIR, f"comparison_{test_case}.png")
    plt.savefig(out_file)
    print(f"  Saved visualization to {out_file}")
    plt.close(fig)
