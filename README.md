# D-TOPS Signal Analysis (Isotropic Free Thermal Expansion Model)

MATLAB code for analyzing D-TOPS (frequency-domain photothermal beam deflection) measurements on a metal-coated bulk sample, using an isotropic free thermal expansion model. Fits the bulk material's thermal conductivity and (optionally) its coefficient of thermal expansion (CTE) from the measured in-phase / out-of-phase deflection signal.

The model sums up to three contributions to the deflection signal:

1. **Thermoelastic** — surface deformation from the CTE of the bulk material (uses `T_bs`, the temperature at the top of the substrate, corrected from `T_s` via the interface conductance and coating transfer matrix)
2. **Air lens** (optional) — probe deflection from `dn_air/dT` in the air above the heated surface
3. **Fresnel reflection** (optional) — phase shift of the reflected probe beam from the temperature dependence of the refractive index

## Files

| File | Description |
|---|---|
| `Main_AlCaF2.m` | Main script: loads data, applies frequency-response correction, runs the fit, and plots results. Edit this file's parameters for your sample/setup. |
| `theta_iso_free_thermal_expansion_model.m` | Core physics model. Computes the predicted complex deflection angle `theta(f)` by summing the thermoelastic, air-lens, and Fresnel-reflection contributions. |
| `BiFDTR_BO_TEMP.m` | Computes the frequency- and wavevector-domain surface temperature `T_s` for a multilayer sample (bidirectional FDTR heat-diffusion solution). |
| `ss_heat.m` | Estimates steady-state (DC) surface heating from the multilayer thermal model. |
| `FIT_inout.m` | Cost function for fitting thermal conductivity **and** CTE simultaneously to the in-phase + out-of-phase data (`FDPBD_fitting1`). |
| `FIT_ratio.m` | Cost function for fitting thermal conductivity only, to the ratio signal (`FDPBD_fitting2`). |
| `GetData_out_in_ratio_f_VSUM.m` | Reads a 4-column lock-in data file (`V_in, V_out, f, V_SUM`). |
| `datacorrection_complex_leaking.m` | Corrects raw voltage data for pump-modulation/detector frequency-response ("leaking") imperfections. |
| `rombint.m` | Romberg numerical integration utility (used by `ss_heat.m`). |
| `*.txt` | Example raw lock-in data file. |

## Requirements

- MATLAB (uses `lsqcurvefit` / `optimoptions`, part of Optimization Toolbox; also uses the built-in `integral`).
- All `.m` files must be in the same folder, with the data file (`.txt`) alongside them.

## Usage

1. Place your raw lock-in data file (4 columns: `V_in`, `V_out`, `f [Hz]`, `V_SUM`, no header) in the same folder.
2. Open `Main_AlCaF2.m` and edit:
   - `FileNames_data` — data filename, without the `.txt` extension
   - Sample parameters (`lambda_down`, `eta_down`, `C_down`, `h_down`, `niu`, `alpha_T`) — thermal conductivity, anisotropy, heat capacity, thickness, and Poisson's ratio/CTE of each layer (coating, interface, bulk)
   - Experimental parameters (objective, beam radius/offset, laser power, metal optical constants)
   - `air_lens.enable` / `reflection.enable` — toggle the optional air-lens and Fresnel-reflection contributions
   - `FDPBD_fitting1` / `FDPBD_fitting2` — choose which fit to run (in/out-of-phase + CTE, or ratio-only)
3. Run the script. Raw data is automatically corrected for the pump/detector frequency response, trimmed to a frequency range around the out-of-phase peak, converted to deflection angle, and fit. Results are printed to the console and plotted (in-phase/out-of-phase and ratio vs. frequency, data vs. fitted model).

Set `flag_save = 1` to also write the data and fitted model curves to `test.dat`.
