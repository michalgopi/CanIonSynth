#!/usr/bin/env python3
import argparse
import os
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

def visualize_nifti(nifti_path, output_path=None):
    """
    Loads a NIfTI file and saves a visualization of orthogonal slices.
    """
    if not os.path.exists(nifti_path):
        print(f"Error: File not found: {nifti_path}")
        return

    try:
        img = nib.load(nifti_path)
        data = img.get_fdata()
    except Exception as e:
        print(f"Error loading NIfTI: {e}")
        return

    # Handle 4D data (time series) by taking first volume
    if len(data.shape) > 3:
        print(f"Warning: Data is {len(data.shape)}D, taking first volume.")
        data = data[..., 0]

    # Normalize for display
    data_min = np.min(data)
    data_max = np.max(data)

    # Avoid division by zero
    if data_max > data_min:
        norm_data = (data - data_min) / (data_max - data_min)
    else:
        norm_data = data

    # Find center slices
    x, y, z = data.shape
    slice_x = data[x // 2, :, :]
    slice_y = data[:, y // 2, :]
    slice_z = data[:, :, z // 2]

    # Create plot
    fig = plt.figure(figsize=(15, 10))

    # 1. Sagittal View (fixed x)
    ax1 = fig.add_subplot(2, 2, 1)
    # Transpose might be needed depending on orientation, typically Z is up in these plots
    im1 = ax1.imshow(np.rot90(slice_x), cmap='gray', aspect='auto')
    ax1.set_title(f'Sagittal View (X={x//2})')
    plt.colorbar(im1, ax=ax1)
    ax1.axis('off')

    # 2. Coronal View (fixed y)
    ax2 = fig.add_subplot(2, 2, 2)
    im2 = ax2.imshow(np.rot90(slice_y), cmap='gray', aspect='auto')
    ax2.set_title(f'Coronal View (Y={y//2})')
    plt.colorbar(im2, ax=ax2)
    ax2.axis('off')

    # 3. Axial View (fixed z)
    ax3 = fig.add_subplot(2, 2, 3)
    im3 = ax3.imshow(np.rot90(slice_z), cmap='gray', aspect='auto')
    ax3.set_title(f'Axial View (Z={z//2})')
    plt.colorbar(im3, ax=ax3)
    ax3.axis('off')

    # 4. Histogram
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.hist(data.flatten(), bins=50, color='c', alpha=0.7)
    ax4.set_title('Intensity Histogram')
    ax4.set_xlabel('Intensity')
    ax4.set_ylabel('Count')
    ax4.grid(True, alpha=0.3)

    plt.suptitle(f"Visualization: {os.path.basename(nifti_path)}\nMin: {data_min:.2f}, Max: {data_max:.2f}, Shape: {data.shape}", fontsize=16)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path)
        print(f"Visualization saved to: {output_path}")
    else:
        # If no output specified, replace extension
        base, _ = os.path.splitext(nifti_path)
        if base.endswith('.nii'):
            base = base[:-4]
        out_png = base + "_vis.png"
        plt.savefig(out_png)
        print(f"Visualization saved to: {out_png}")

    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize NIfTI file as orthogonal slices.")
    parser.add_argument("file", help="Path to input NIfTI file")
    parser.add_argument("-o", "--output", help="Path to output PNG file (optional)")

    args = parser.parse_args()

    visualize_nifti(args.file, args.output)
