# TransientInjector
This project simulates and injects transient astronomical sources into high-time-cadence Zarr image cubes derived from MeerKAT radio telescope data. By modelling transient sources and incorporating them into existing observational data, the project aims to test and refine transient detection algorithms. The process includes accounting for the time-varying point spread function (PSF) for each time step, ensuring realistic integration of transient sources.

## Installation

```
pip install TransientInjector
```
Alternatively, install from source

```bash
git clone https://github.com/SihlanguI/TransientInjector
cd TransientInjector
pip install .
```

---

## Usage

The package insert realistic transients in the data:

Run the package
```
python3 -m transsim.process_trans --help
python3 -m transsim.process_trans config.yaml

```
- ```<confif.yaml>```: This the config file required to run the package
     

NB: Edit the config file accordingly

