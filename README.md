# TransientInjector

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)

---

**TransientInjector** is a Python package for simulating and injecting astronomical transient sources into high-time-cadence Zarr image cubes derived from MeerKAT radio telescope data and interactively visualising radio-frequency interference (RFI) statistics for MeerKAT telescope observations.  By modelling transient sources and incorporating them into existing observational data, the package aims to test and refine transient detection algorithms. The process includes accounting for the time-varying point spread function (PSF) for each time step, ensuring realistic integration of transient sources.

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

