# TransientInjector
This project simulates and injects transient astronomical sources into high-time-cadence Zarr image cubes derived from MeerKAT radio telescope data. By modelling transient sources and incorporating them into existing observational data, the project aims to test and refine transient detection algorithms. The process includes accounting for the time-varying point spread function (PSF) for each time step, ensuring realistic integration of transient sources.

###  Key Components:
#### 1. Simulation of Transient Sources:
    The project involves generating synthetic transient sources that mimic real astrophysical phenomena, such as pulsars, fast radio bursts (FRBs), or supernova remnants. These sources are characterized by their temporary, short-lived signals.

#### 2. Injection into MeerKAT Data:
    These simulated transients are injected into actual MeerKAT observational data, represented in high-time-cadence Zarr cubes. This allows for the testing of detection and analysis techniques within realistic conditions.

#### 3. Consideration of Dirty PSF:
    The dirty PSF (Point Spread Function), which represents the telescope array's response at a given time, is applied for each time step. This ensures the injected sources are impacted by the same observational effects as the real data.

#### Goal:
    The primary objective is to enhance methods of detecting transient sources in radio astronomy, improving the accuracy and robustness of signal detection in MeerKAT data. This work is crucial for advancing transient astronomy and understanding short-lived cosmic events.
