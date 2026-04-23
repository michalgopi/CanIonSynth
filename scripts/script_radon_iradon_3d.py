import click
import sys
import logging
from joblib import Parallel, delayed
from pathlib import Path
from time import perf_counter

import numpy as np
import SimpleITK as sitk
from skimage.transform import radon, iradon

log = logging.getLogger(__name__)


def basic_log_config() -> None:
    """Configure the logging format and level."""
    logging.basicConfig(
        format='[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stdout,
        level=logging.INFO,
    )


def add_noise_to_sinogram(sinogram: np.ndarray, noise_level: float) -> np.ndarray:
    """Add Gaussian noise to the sinogram. Set value from 0-no noise to 1-large noise."""
    sino_2d = sinogram.reshape(sinogram.shape[0], -1)
    noise = np.random.normal(0, noise_level * 100, sino_2d.shape)
    noisy_sinogram = sino_2d + noise
    return noisy_sinogram.reshape(sinogram.shape)


def post_processing(org_img: sitk.Image, reco_array: np.ndarray) -> sitk.Image:
    """ Post-process the reconstructed image to match the original image size and spacing"""
    org_array = sitk.GetArrayFromImage(org_img)
    org_shape = np.array(org_array.shape)
    reco_shape = np.array(reco_array.shape)
    diff = reco_shape - org_shape
    start = np.maximum(diff // 2, 0)
    end = start + org_shape
    cropped_reco = reco_array[start[0]:end[0], start[1]:end[1], start[2]:end[2]]
    pad_width = np.maximum(-diff // 2, 0)
    pad_end = org_shape - cropped_reco.shape - pad_width
    padded_reco = np.pad(cropped_reco,
                         pad_width=[(pad_width[0], pad_end[0]), (pad_width[1], pad_end[1]), (pad_width[2], pad_end[2])],
                         mode='constant')
    reco_img = sitk.GetImageFromArray(padded_reco)
    reco_img.CopyInformation(org_img)
    return reco_img


@click.command()
@click.argument('phantom-file', type=click.Path(exists=True))
@click.option('--n-theta', default=10, help='Number of angles in geometry for Radon transform')
@click.option('--output-filename', default=None, help='Output file name for the reco image')
@click.option('--noise-level', type=float, default=0.0, help='Add noise to the sinogram 0 - no noise, 1 - large noise')
def main(phantom_file: str, n_theta: int, output_filename: str, noise_level: float) -> None:
    """Process a 3D .nii.gz file using Radon and iRadon transforms and save as new .nii.gz file."""

    start_time = perf_counter()
    phantom_img = sitk.ReadImage(phantom_file)
    object_array = sitk.GetArrayFromImage(phantom_img)

    # Simplest geometry ever, identical for azimuthal and polar angles (theta)
    theta = np.linspace(90.0 / n_theta / 2, 180 - (90.0 / n_theta) / 2, n_theta)

    log.info(
        f'Loaded {Path(phantom_file).name}. '
        f'Starting with geometry: alpha: {n_theta}, beta: {n_theta} ({n_theta ** 2} projections). '
        f'Entering Radon stage 1 ...')

    # Radon for azimuthal angles (stage 1)
    r1 = np.stack(Parallel(n_jobs=-1)(
        delayed(radon)(object_array[:, :, k], theta=theta, circle=False)
        for k in range(object_array.shape[2])
    ), axis=2)

    log.info(f'Radon stage 1 finished. Entering Radon stage 2 ...')

    # Radon for polar angles (stage 2)
    r2 = np.stack(Parallel(n_jobs=-1)(
        delayed(radon)(r1[:, l, :], theta=theta, circle=False)
        for l in range(n_theta)
    ), axis=1)
    log.info(f'Radon stage 2 (sinogram) finished. Entering iRadon stage 1 ...')

    # Add noise to the sinogram

    if noise_level > 0:
        log.info(f'Adding noise to the sinogram ...')
        r2 = add_noise_to_sinogram(r2, noise_level=noise_level)

    # iRadon for polar angles (stage 1)
    i1 = np.stack(Parallel(n_jobs=-1)(
        delayed(iradon)(r2[:, m, :], theta=theta, circle=False)
        for m in range(n_theta)
    ), axis=1)
    log.info(f'iRadon stage 1 finished. Entering iRadon stage 2 ...')

    # iRadon for azimuthal angles (stage 2)
    reco = np.stack(Parallel(n_jobs=-1)(
        delayed(iradon)(i1[:, :, n], theta=theta, circle=False)
        for n in range(i1.shape[2])
    ), axis=2)
    log.info(f'iRadon stage 2 (reconstruction) finished. Starting postprocessing...')

    # Postprocessing
    reco_final = post_processing(phantom_img, reco)
    log.info(f'Postprocessing finished. Saving results...')

    # Save the reconstructed image
    if output_filename is None:
        source_path = Path(phantom_file).parent
        filename = Path(phantom_file).stem.split('.')[0]
        output_filename = f'{source_path / filename}_reco_{n_theta}x{n_theta}.nii.gz'

    sitk.WriteImage(reco_final, output_filename)
    log.info(f'Finished with {Path(phantom_file).name} in: {round(perf_counter() - start_time, 2)} s')


if __name__ == '__main__':
    basic_log_config()
    main()
