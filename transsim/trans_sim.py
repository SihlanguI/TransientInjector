import numpy as np
import xarray as xr
from astropy.io import fits
from typing import Optional, Tuple
from scipy.signal import fftconvolve


# Dictionary mapping Stokes parameter names to indices
stokes_mapping = {
    'I': 0,
    'Q': 1,
    'U': 2,
    'V': 3
}

def get_cube(path: str, stokes_par: str, varname: Optional[str] = None,
             slice_indices: Optional[Tuple[slice]] = None):
    """
    Load a Zarr file, extract data for a specific Stokes parameter, and optionally slice it.

    Parameters:
    ------------
    path : str
        Path to the Zarr file containing the dataset.
    stokes_par : str
        The Stokes parameter to select (e.g., 'I', 'Q', 'U', 'V').
    varname : str, optional
        Name of the variable to extract. Defaults to 'cube' if not provided.
    slice_indices : tuple of slice, optional
        Tuple specifying slices to extract from the Zarr array.

    Returns:
    --------
    ds : xarray.core.dataset.Dataset
        The complete dataset from the Zarr file.
    data : numpy.ndarray
        The extracted and optionally sliced data array.
    mjr_axis : numpy.ndarray
        Array of major axis values per time step.
    mnr_axis : numpy.ndarray
        Array of minor axis values per time step.
    pa : numpy.ndarray
        Array of position angle values per time step.
    header : fits.Header
        The FITS header associated with the data.
    """
    # Open the Zarr file
    ds = xr.open_zarr(path)
    mjr_axis = ds.psf_maj.values
    mnr_axis = ds.psf_min.values
    pa = ds.psf_pa.values
    
    # Extract variable (default to 'cube' if varname is not specified)
    cube = getattr(ds, varname) if varname else ds.cube
    
    # Extract FITS header from dataset attributes
    header = fits.Header(dict(ds.attrs['fits_header']))
    nt, nx, ny = cube.sizes['TIME'], cube.sizes['X'], cube.sizes['Y']

    # Check validity of Stokes parameter
    if stokes_par not in stokes_mapping:
        raise ValueError(f"Invalid Stokes parameter: {stokes_par}. Must be one of {list(stokes_mapping.keys())}")
    
    # Get the corresponding Stokes index
    stokes_index = stokes_mapping[stokes_par]
    
    # Slice data if indices are provided
    if slice_indices:
        data = cube[stokes_index][slice_indices]
    else:
        data = cube[stokes_index, slice(0, nt), slice(0, nx), slice(0, ny)]

    # Update header to match sliced data shape
    header['NAXIS'] = len(data.shape)
    header['NAXIS1'] = data.shape[1]
    header['NAXIS2'] = data.shape[2]
    header['NAXIS3'] = data.shape[0]
    header.pop('NAXIS4', None)  # Remove NAXIS4 if present

    return ds, data, mjr_axis, mnr_axis, pa, header


def gaussian_beam(header, x_size, y_size, fwhm_x, fwhm_y , position_angle_deg, peak_flux):
    """
    Create a 2D elliptical Gaussian beam aligned with FITS convention.
    
    Parameters:
    -----------
    header : dict
        FITS header containing pixel scale ('CDELT1/2') for converting sigma from degrees to pixels.
    x_size : int
        The width of the 2D grid in pixels.
    y_size : int
        The height of the 2D grid in pixels.
    sigma_major : float
        Major axis standard deviation of the Gaussian in degrees.
    sigma_minor : float
        Minor axis standard deviation of the Gaussian in degrees.
    position_angle_deg : float
        Position angle of the Gaussian's major axis in degrees.
    peak flux : float
        peak flux

    Returns:
    --------
    image : numpy.ndarray
        A 2D array representing the transient source.
    """
    x = np.linspace(-x_size // 2, x_size // 2, x_size)
    y = np.linspace(-y_size // 2, y_size // 2, y_size)
    X, Y = np.meshgrid(x, y)

    # Convert position angle: FITS PA from north to east -> rotation from x-axis
    position_angle_rad = np.radians(90 - position_angle_deg)

    # Rotate coordinates by position angle
    X_rot = X * np.cos(position_angle_rad) - Y * np.sin(position_angle_rad)
    Y_rot = X * np.sin(position_angle_rad) + Y * np.cos(position_angle_rad)

    # Convert sigma from degrees to pixels
    hwhm = 2 * np.sqrt(2 * np.log(2)) # Half-Width at half max
    sigma_major = fwhm_x / hwhm
    sigma_minor = fwhm_y / hwhm
    sigma_x_pix = sigma_major / header['CDELT1']  # x-direction
    sigma_y_pix = sigma_minor / header['CDELT2']  # y-direction

    # 2D Gaussian
    beam = np.exp(-0.5 * ((X_rot / sigma_x_pix) ** 2 + (Y_rot / sigma_y_pix) ** 2))
    image = peak_flux * beam

    return image


def insert_transient_into_zarr(cube, t_index, transient, pos_x, pos_y):
    """
    Insert a transient source into a 3D Zarr cube at a specified position and time index.

    Parameters:
    -----------
    cube : xarray.core.dataarray.DataArray
        The 3D data array where dimensions are typically (time, height, width).
    t_index : int
        The time index where the transient will be inserted.
    transient : numpy.ndarray
        The 2D array representing the transient source.
    pos_x : int
        The x-coordinate where the transient will be centered.
    pos_y : int
        The y-coordinate where the transient will be centered.

    Notes:
    ------
    - The transient is inserted by adding its values to the existing cube at the specified coordinates.
    - The function handles cases where the transient exceeds image boundaries.
    
    Returns:
    --------
    None
    """
    transient_size = transient.shape[1]
    img_height, img_width = cube.shape[1:]
    
    half_size = transient_size // 2
    x_start = max(pos_x - half_size, 0)
    y_start = max(pos_y - half_size, 0)
    x_end = min(pos_x + half_size, img_width)
    y_end = min(pos_y + half_size, img_height)
    
    transient_x_start = half_size - min(pos_x, half_size)
    transient_y_start = half_size - min(pos_y, half_size)
    
    cube[t_index, y_start:y_end, x_start:x_end] += transient[transient_y_start:transient_y_start + (y_end - y_start), 
                                                       transient_x_start:transient_x_start + (x_end - x_start)]


def update_zarr(out_zarr, ds, cube, varname, stokes_par, overwrite_cube):
    """
    Update the Zarr file with the modified cube data.

    Parameters:
    -----------
    zarr_path : str
        Path to the Zarr file to be updated.
    ds : xarray.core.dataset.Dataset
        The dataset object that contains the Zarr data.
    cube : xarray.core.dataarray.DataArray
        The modified cube data.
    varname : str
        The name of the variable within the dataset to be updated.
    stokes_par : str
        The Stokes parameter (e.g., 'I', 'Q', 'U', 'V') being updated.

    Returns:
    --------
    None
    """
    # Update the cube for the selected Stokes parameter
    stokes_index = stokes_mapping[stokes_par]
    ds[varname][stokes_index, :, :, :] = cube
    
    # Save the updated data back to the Zarr file
    if overwrite_cube:
        ds.to_zarr(out_zarr, mode='r+')
    else:
        ds.to_zarr(out_zarr, mode='w')