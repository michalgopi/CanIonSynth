#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import json
import numpy as np
import SimpleITK as sitk
import pydicom
import pydicom_seg
from pathlib import Path
from pydicom import config
config.replace_un_with_known_vr = False

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Path to metadata template for DICOM-SEG creation
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
METAINFO_PATH = os.path.join(SCRIPT_DIR, "for_dicom_seg.json")

def parse_args():
    """
    Parse command line arguments for the NIfTI to DICOM-SEG conversion.

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(description='Convert 3D NIfTI binary mask to DICOM-SEG')
    parser.add_argument('nifti_path', type=str, help='Path to NIfTI binary mask file (uint8)')
    parser.add_argument('reference_dicom', type=str, help='Path to reference DICOM directory or single file')
    parser.add_argument('output_folder', type=str, help='Path to output folder for DICOM-SEG')
    parser.add_argument('--reference_nifti', type=str, help='Path to reference NIfTI file with matching geometry', default=None)

    return parser.parse_args()

def validate_inputs(nifti_path, reference_dicom, output_folder, reference_nifti=None):
    """Validate the input files and output directory."""
    if not os.path.exists(nifti_path):
        logger.error(f"NIfTI file not found: {nifti_path}")
        return False

    if not os.path.exists(reference_dicom):
        logger.error(f"Reference DICOM path not found: {reference_dicom}")
        return False

    if reference_nifti and not os.path.exists(reference_nifti):
        logger.error(f"Reference NIfTI file not found: {reference_nifti}")
        return False

    if not os.path.exists(output_folder):
        logger.info(f"Output folder does not exist, creating: {output_folder}")
        try:
            os.makedirs(output_folder)
        except Exception as e:
            logger.error(f"Failed to create output folder: {e}")
            return False

    return True

def read_nifti_mask(nifti_path):
    """
    Read a NIfTI binary mask file and ensure it contains only 0 and 1 values.

    Args:
        nifti_path (str): Path to the NIfTI file

    Returns:
        SimpleITK.Image: The loaded NIfTI image as binary mask
    """
    try:
        logger.info(f"Reading NIfTI file: {nifti_path}")
        nifti_image = sitk.ReadImage(nifti_path)

        # Get array data to check values
        pixel_array = sitk.GetArrayFromImage(nifti_image)
        unique_values = np.unique(pixel_array)

        logger.info(f"Unique values in mask: {unique_values}")
        logger.info(f"Mask size: {nifti_image.GetSize()}")
        logger.info(f"Mask spacing: {nifti_image.GetSpacing()}")
        logger.info(f"Mask direction: {nifti_image.GetDirection()}")

        # Ensure the mask is binary (0 and 1)
        if not np.all(np.isin(unique_values, [0, 1])):
            logger.warning(f"NIfTI mask contains values other than 0 and 1: {unique_values}")

            # Convert to binary (0 and 1)
            binary_array = (pixel_array > 0).astype(np.uint8)
            logger.info("Converting to binary (0 and 1)")

            # Create a new SimpleITK image with binary values
            binary_image = sitk.GetImageFromArray(binary_array)
            binary_image.CopyInformation(nifti_image)  # Copy metadata

            return binary_image

        return nifti_image

    except Exception as e:
        logger.error(f"Failed to read NIfTI file: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def load_reference_dicoms(dicom_path):
    """
    Load DICOM files from a directory or a single file.

    Args:
        dicom_path (str): Path to DICOM directory or single DICOM file

    Returns:
        tuple: (list of DICOM file paths, SimpleITK image of the reference)
    """
    try:
        if os.path.isdir(dicom_path):
            # It's a directory, try to read series
            logger.info(f"Reading DICOM series from directory: {dicom_path}")

            # Use SimpleITK to find DICOM series
            series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(dicom_path)
            if not series_IDs:
                logger.error("No DICOM series found in the directory")
                return None, None

            # Use the first series
            series_ID = series_IDs[0]
            dicom_filenames = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(dicom_path, series_ID)

            # Read the DICOM series
            reader = sitk.ImageSeriesReader()
            reader.SetFileNames(dicom_filenames)
            reference_image = reader.Execute()

            logger.info(f"Loaded reference DICOM series with size: {reference_image.GetSize()}")
            logger.info(f"Reference DICOM spacing: {reference_image.GetSpacing()}")
            logger.info(f"Reference DICOM direction: {reference_image.GetDirection()}")

            return dicom_filenames, reference_image

        else:
            # It's a single file
            logger.info(f"Reading single DICOM file: {dicom_path}")

            # Read the DICOM file for reference
            reader = sitk.ImageFileReader()
            reader.SetFileName(dicom_path)
            reference_image = reader.Execute()

            logger.info(f"Loaded single DICOM with size: {reference_image.GetSize()}")

            return [dicom_path], reference_image

    except Exception as e:
        logger.error(f"Failed to load reference DICOMs: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None, None

def ensure_spacing_consistency(nifti_mask, dicom_reference, reference_nifti=None):
    """
    Ensure the NIfTI mask and DICOM reference have consistent spacing and orientation.
    Reorients or resamples the mask if needed.

    Args:
        nifti_mask (SimpleITK.Image): Mask as SimpleITK Image
        dicom_reference (SimpleITK.Image): Reference DICOM as SimpleITK Image
        reference_nifti (SimpleITK.Image, optional): Reference NIfTI image

    Returns:
        SimpleITK.Image: Reoriented/resampled mask that matches DICOM geometry
    """
    try:
        # If reference_nifti is provided, use it for orientation information
        if reference_nifti is not None:
            logger.info("Using reference NIfTI for orientation")

            # Check if the mask needs to be reoriented to match reference NIfTI
            if nifti_mask.GetSize() != reference_nifti.GetSize():
                logger.warning(f"Mask size {nifti_mask.GetSize()} differs from reference NIfTI {reference_nifti.GetSize()}")
                logger.info("Will attempt to match orientations...")

                # Try to reorient by permuting axes if needed
                nifti_array = sitk.GetArrayFromImage(nifti_mask)

                # Test different permutations to find a matching one
                for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
                    # Skip identity permutation
                    if perm == (0,1,2):
                        continue

                    # Transpose the array
                    permuted_array = np.transpose(nifti_array, perm)
                    if permuted_array.shape == sitk.GetArrayFromImage(reference_nifti).shape:
                        logger.info(f"Found matching permutation: {perm}")

                        # Create new image with permuted array
                        permuted_mask = sitk.GetImageFromArray(permuted_array)
                        permuted_mask.CopyInformation(reference_nifti)

                        # Now permuted_mask should match reference_nifti's orientation
                        nifti_mask = permuted_mask
                        break

        # Ensure mask image has matching physical space as the DICOM reference
        if nifti_mask.GetSize() != dicom_reference.GetSize():
            logger.warning(f"Mask size {nifti_mask.GetSize()} differs from DICOM reference {dicom_reference.GetSize()}")
            logger.info("Resampling mask to match DICOM reference geometry...")

            # Resample mask to match DICOM reference
            resampled_mask = sitk.Resample(
                nifti_mask,
                dicom_reference,
                sitk.Transform(),
                sitk.sitkNearestNeighbor,  # Use nearest neighbor to preserve binary labels
                0,  # Default pixel value
                dicom_reference.GetPixelID()
            )

            logger.info(f"Resampled mask size: {resampled_mask.GetSize()}")
            return resampled_mask
        else:
            # Sizes match, just copy the geometry
            result = sitk.Image(nifti_mask)
            result.CopyInformation(dicom_reference)
            return result

    except Exception as e:
        logger.error(f"Failed to ensure spacing consistency: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return nifti_mask  # Return original as fallback

def convert_nifti_to_dicom_seg(nifti_mask, dicom_filenames, dicom_reference, output_path, seg_name, reference_nifti=None):
    """
    Convert a NIfTI binary mask to DICOM-SEG using pydicom_seg.

    Args:
        nifti_mask (SimpleITK.Image): Binary mask as SimpleITK Image
        dicom_filenames (list): List of DICOM file paths
        dicom_reference (SimpleITK.Image): Reference DICOM as SimpleITK Image
        output_path (str): Path to save the DICOM-SEG file
        seg_name (str): Name of the segmentation
        reference_nifti (SimpleITK.Image, optional): Reference NIfTI image

    Returns:
        bool: True if conversion is successful
    """
    logger.info(f"Converting NIfTI mask to DICOM-SEG: {seg_name}")

    # Create template from metadata file
    if os.path.exists(METAINFO_PATH):
        template = pydicom_seg.template.from_dcmqi_metainfo(METAINFO_PATH)
        logger.info("Created DICOM-SEG template from metadata file")
    else:
        # Create a basic template if metadata file doesn't exist
        logger.warning(f"Metadata file not found: {METAINFO_PATH}")
        logger.info("Creating basic DICOM-SEG template")

        template = pydicom_seg.template.SegmentationTemplate()
        template.add_segment(
            segment_number=1,
            segment_label=seg_name,
            segmented_property_category="T-D0050",
            segmented_property_type="T-D0050",
            algorithm_type="AUTOMATIC",
            algorithm_name="NIfTI to DICOM-SEG"
        )


    # Create the writer with options
    logger.info("Creating DICOM-SEG writer")
    writer = pydicom_seg.MultiClassWriter(
        template=template,
        inplane_cropping=True,     # Crop to minimum bounding box
        skip_empty_slices=True,    # Skip slices with no segmentation
        skip_missing_segment=False  # Raise error for missing segments
    )

    try:
        # Ensure mask has proper orientation and spacing
        aligned_mask = ensure_spacing_consistency(nifti_mask, dicom_reference, reference_nifti)

        # Load the DICOM files for reference
        source_images = []
        for dcm_file in dicom_filenames:
            try:
                source_images.append(pydicom.dcmread(dcm_file, stop_before_pixels=True))
            except Exception as e:
                logger.warning(f"Failed to read DICOM file {dcm_file}: {e}")

        if not source_images:
            logger.error("No valid DICOM reference files found")
            return False

        # Write the DICOM-SEG
        logger.info(f"Writing DICOM-SEG with {len(source_images)} reference files")

        # Convert to uint8 binary mask if needed
        binary_mask = sitk.Cast(aligned_mask, sitk.sitkUInt8)

        dcm = writer.write(binary_mask, source_images)

        # Add custom series description
        dcm.SeriesDescription = f"{seg_name} Segmentation"

        # Save the file
        logger.info(f"Saving DICOM-SEG to {output_path}")
        dcm.save_as(output_path)
        logger.info(f"DICOM-SEG saved successfully to {output_path}")

        return True

    except Exception as e:
        logger.error(f"Failed to convert to DICOM-SEG: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def main():
    """Main function to orchestrate the NIfTI to DICOM-SEG conversion."""
    args = parse_args()

    if not validate_inputs(args.nifti_path, args.reference_dicom, args.output_folder, args.reference_nifti):
        return 1

    # Read the NIfTI mask
    nifti_mask = read_nifti_mask(args.nifti_path)
    if nifti_mask is None:
        return 1

    # Read reference NIfTI if provided
    reference_nifti = None
    if args.reference_nifti:
        try:
            logger.info(f"Reading reference NIfTI: {args.reference_nifti}")
            reference_nifti = sitk.ReadImage(args.reference_nifti)
            logger.info(f"Reference NIfTI size: {reference_nifti.GetSize()}")
            logger.info(f"Reference NIfTI spacing: {reference_nifti.GetSpacing()}")
        except Exception as e:
            logger.error(f"Failed to read reference NIfTI: {e}")
            reference_nifti = None

    # Load reference DICOMs
    dicom_filenames, dicom_reference = load_reference_dicoms(args.reference_dicom)
    if dicom_filenames is None or dicom_reference is None:
        return 1

    # Extract segmentation name from NIfTI filename
    nifti_filename = os.path.basename(args.nifti_path)
    seg_name = os.path.splitext(nifti_filename)[0]
    if seg_name.endswith('.nii'):
        seg_name = seg_name[:-4]

    # Remove suffix like '_mask' if present
    for suffix in ['_mask', '_seg', '_segmentation']:
        seg_name = seg_name.replace(suffix, '')

    # Create output file path
    output_path = os.path.join(args.output_folder, f'{seg_name}_segmentation.dcm')

    # Convert and save
    success = convert_nifti_to_dicom_seg(
        nifti_mask, dicom_filenames, dicom_reference,
        output_path, seg_name, reference_nifti
    )

    if success:
        logger.info(f"Successfully created DICOM-SEG file: {output_path}")
        return 0
    else:
        logger.error("DICOM-SEG conversion failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
