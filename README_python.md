# D-TOPS Signal Analysis (Python Version)

Python translation of the MATLAB codebase for analyzing D-TOPS (frequency-domain photothermal beam deflection) measurements on a metal-coated bulk sample using an isotropic free thermal expansion model.

It fits the bulk material's thermal conductivity ($\lambda$) and coefficient of thermal expansion (CTE, $\alpha_T$) from measured in-phase / out-of-phase beam deflection signals.

## Requirements & Installation

1. Python 3.8 or higher.
2. Install dependencies:
```bash
pip install -r requirements.txt
```

Dependencies:
- `numpy`
- `scipy`
- `matplotlib`

## Modules Overview

- `main.py`: Entry point script. Configures parameters, loads raw lock-in data, performs leaking frequency-response correction, runs fitting routines, outputs numerical results, and generates plots.
- `data_io.py`: Functionality for reading 4-column lock-in text files (`get_data_out_in_ratio_f_vsum`) and correcting raw complex voltage data (`datacorrection_complex_leaking`).
- `heat_diffusion.py`: Multi-layer bidirectional thermal transfer model (`bifdtr_bo_temp`) and steady-state heating calculation (`ss_heat`).
- `deflection_model.py`: Core physics model (`theta_iso_free_thermal_expansion_model`) summing thermoelastic deformation, air lens mirage effect, and Fresnel reflection contributions using `scipy.integrate.quad`.
- `fitting.py`: Nonlinear optimization fitting functions (`fit_inout`, `fit_ratio`) using `scipy.optimize.least_squares`.

## Usage

Run `main.py` with default sample data:
```bash
python main.py
```

Pass a custom data file:
```bash
python main.py --file DTOPS_AlCaF2_10x_1p0-1p0_100k-100_313p2mV.txt
```

Save fitted model predictions and experimental data to `test.dat`:
```bash
python main.py --save
```

Run in headless mode without displaying plots:
```bash
python main.py --no-plot
```
