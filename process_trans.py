#!/usr/bin/env python

import numpy as np
import sys
import yaml
from trans_sim import get_cube, create_transient_source, insert_transient_into_zarr, update_zarr

def main(config):
    """
    Main function to insert transient sources at multiple positions and time indices into a Zarr cube.

    Parameters:
    -----------
    config : dict
        Dictionary containing configuration parameters.

    Returns:
    --------
    None
    """

    # Step 1: Load the cube and its associated parameters from the Zarr file
    zarr_path = config['zarr_path']
    varname = config['varname']
    stokes_par = config['stokes_par']
    start_index = config['start_index']
    dur = config['duration']
    inte = config['integration']
    positions = config['positions']
    amplitude = config['amplitude']

    ds, cube, mjr_axis, mnr_axis, pa, header = get_cube(zarr_path, stokes_par, varname=varname)
    print(f"Loaded data cube with shape: {cube.shape}")
    no_timestamps = dur//2
    t_index = [start_index+i for i in range(no_timestamps)]
    # Step 2: Insert the transient source at multiple time indices and positions
    for t_index, (pos_x, pos_y) in zip(t_indices, positions):
        # Get the corresponding major axis, minor axis, and position angle for the current time index
        sigma_major = mjr_axis[t_index]
        sigma_minor = mnr_axis[t_index]
        position_angle_deg = pa[t_index]

        # Create a transient source for the current time index
        transient_source = create_transient_source(
            header=header,
            x_size=cube.shape[1],
            y_size=cube.shape[2],
            sigma_major=sigma_major,
            sigma_minor=sigma_minor,
            position_angle_deg=position_angle_deg,
            amplitude=amplitude
        )
        print(f"Created transient source at time index {t_index} with major axis: {sigma_major}, minor axis: {sigma_minor}, PA: {position_angle_deg}")

        # Insert the transient source into the cube
        insert_transient_into_zarr(cube, t_index=t_index, transient=transient_source, pos_x=pos_x, pos_y=pos_y)
        print(f"Inserted transient into the cube at time index {t_index}, position ({pos_x}, {pos_y})")

    # Step 3: Save the modified cube back to the same Zarr file
    update_zarr(zarr_path, ds, cube, varname, stokes_par)
    print(f"Saved the modified cube to {zarr_path}")


if __name__ == "__main__":
    # Load the YAML configuration file
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config_file>")
        sys.exit(1)

    config_file = sys.argv[1]

    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)

    # Call the main function with the loaded configuration
    main(config)
