#!/usr/bin/env python
import argparse
import logging
import os
import sys
import yaml
from typing import List

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from transsim.trans_sim import get_cube, conv_source, insert_transient_into_zarr, update_zarr


def create_parser():
    parser = argparse.ArgumentParser(description='MeerKAT Transient Simulation.')
    print('MEERKAT Transient Simulation')
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='config file')
    
    return parser

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

    # Step 1: Get parmeters from the config file
    in_zarr = config['in_zarr']
    out_zarr = config['out_zarr']
    overwrite_cube = config['overwrite_cube']
    source = config['source_type']
    varname = config['varname']
    stokes_par = config['stokes_par']
    start_index = config['start_index']
    dur = config['duration']
    int_time= config['integration']
    positions = config['positions']
    amplitude = config['amplitude']

    if not any(isinstance(i, list) for i in positions):
            positions = [positions]
    
    print("Reading in the zarr cube")
    ds, cube, mjr_axis, mnr_axis, pa, header = get_cube(in_zarr, stokes_par, varname=varname)
    no_timestamps = dur//int_time
    t_indices = [start_index+i for i in range(no_timestamps)]
    # Step 2: Insert the transient source at multiple time indices and positions
    for t_index, amp, (pos_x, pos_y) in zip(t_indices, amplitude, positions):
        # Get the corresponding major axis, minor axis, and position angle for the current time index
        sigma_major = mjr_axis[t_index]
        sigma_minor = mnr_axis[t_index]
        position_angle_deg = pa[t_index]

        # Create a transient source for the current time index
        print("Inserting a transient into the cube at time index {}, position {}, and amplitude {}".format(t_index, (pos_x, pos_y), amp))
        convolved_source = conv_source(
            header=header,
            x_size=cube.shape[1],
            y_size=cube.shape[2],
            sigma_major=sigma_major,
            sigma_minor=sigma_minor,
            position_angle_deg=position_angle_deg,
            amplitude=amp,
            source_type=source
        )
        # Insert the transient source into the cube
        insert_transient_into_zarr(cube, t_index=t_index, transient=convolved_source, pos_x=pos_x, pos_y=pos_y)
        print("Completed the insertion of the transient source")

    # Step 3: Save the modified cube to disk
    if overwrite_cube:
        out_zarr = in_zarr
    else:
        path, filename = os.path.split(in_zarr)
        out_zarr = os.path.join(out_zarr, filename)
    print("Saving the modified cube to {}".format(out_zarr))
    update_zarr(out_zarr, ds, cube, varname, stokes_par, overwrite_cube)
    print("Cube is saved")


if __name__ == "__main__":
    # Load the YAML configuration file
    parser = create_parser()
    args = parser.parse_args()
    config_file = args.config

    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)

    # Call the main function with the loaded configuration
    main(config)
