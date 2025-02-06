# TransientInjector
This project simulates and injects transient astronomical sources into high-time-cadence Zarr image cubes derived from MeerKAT radio telescope data. By modelling transient sources and incorporating them into existing observational data, the project aims to test and refine transient detection algorithms. The process includes accounting for the time-varying point spread function (PSF) for each time step, ensuring realistic integration of transient sources.

###  Run the code:
     python3 process_trans.py config.yaml

