import numpy as np
from scipy.special import jv
from scipy.integrate import quad
from heat_diffusion import bifdtr_bo_temp


def complex_quad(func, a, b, **kwargs):
    """Integrates a complex-valued function of a real variable over [a, b]."""
    real_res, _ = quad(lambda x: np.real(func(x)), a, b, **kwargs)
    imag_res, _ = quad(lambda x: np.imag(func(x)), a, b, **kwargs)
    return real_res + 1j * imag_res


def thermoelastic_integrand(
    kvect, freq, niu, alpha_T, q2_bulk,
    eta_bulk, lambda_metal, q2_metal, G_int, L_metal,
    r_pump, r_probe, A_pump, xoffset, C_probe,
    lambda_up, C_up, h_up, eta_up, lambda_down, C_down, h_down, eta_down
):
    k = 2.0 * np.pi * kvect
    zeta_metal = np.sqrt(q2_metal + k**2)
    zeta_metal_L = zeta_metal * L_metal
    gamma_metal = lambda_metal * zeta_metal

    conv = np.cosh(zeta_metal_L) + (gamma_metal / G_int) * np.sinh(zeta_metal_L)
    corr = np.sinh(zeta_metal_L) / gamma_metal + np.cosh(zeta_metal_L) / G_int

    T_s_Sprobe = bifdtr_bo_temp(
        kvect, freq, lambda_up, C_up, h_up, eta_up,
        lambda_down, C_down, h_down, eta_down,
        r_pump, r_probe, A_pump
    )

    flx = A_pump * np.exp(-(np.pi**2) * (r_pump**2) / 2.0 * (kvect**2))
    S_probe = np.exp(-(np.pi**2) * (r_probe**2) / 2.0 * (kvect**2))

    T_bs_Sprobe = conv * T_s_Sprobe - corr * flx * S_probe

    val = (-C_probe * 8.0 * (np.pi**2) * (kvect**2) *
           (-jv(1, 2.0 * np.pi * kvect * xoffset)) *
           (2.0 * (1.0 + niu) * alpha_T /
            (np.sqrt(4.0 * (np.pi**2) * eta_bulk * (kvect**2) + q2_bulk) + 2.0 * np.pi * kvect)) *
           T_bs_Sprobe)
    return val


def air_lens_integrand(
    kvect, freq, q2_air,
    lambda_up, C_up, h_up, eta_up,
    lambda_down, C_down, h_down, eta_down,
    r_pump, r_probe, A_pump, xoffset, C_probe, coef_air
):
    k = 2.0 * np.pi * kvect
    zeta_air = np.sqrt(q2_air + k**2)

    T_s_Sprobe = bifdtr_bo_temp(
        kvect, freq, lambda_up, C_up, h_up, eta_up,
        lambda_down, C_down, h_down, eta_down,
        r_pump, r_probe, A_pump
    )

    Z_air_Sprobe = (coef_air / zeta_air) * T_s_Sprobe

    val = (-C_probe * 8.0 * (np.pi**2) * (kvect**2) *
           (-jv(1, 2.0 * np.pi * kvect * xoffset)) * Z_air_Sprobe)
    return val


def reflection_integrand(
    kvect, freq, dphidT,
    lambda_up, C_up, h_up, eta_up,
    lambda_down, C_down, h_down, eta_down,
    r_pump, r_probe, A_pump, xoffset, C_probe, wavelength
):
    T_s_Sprobe = bifdtr_bo_temp(
        kvect, freq, lambda_up, C_up, h_up, eta_up,
        lambda_down, C_down, h_down, eta_down,
        r_pump, r_probe, A_pump
    )

    Z_r_Sprobe = -wavelength * dphidT * T_s_Sprobe / (4.0 * np.pi)

    val = (-C_probe * 8.0 * (np.pi**2) * (kvect**2) *
           (-jv(1, 2.0 * np.pi * kvect * xoffset)) * Z_r_Sprobe)
    return val


def theta_iso_free_thermal_expansion_model(
    niu, alpha_T, f,
    lambda_down, C_down, h_down, eta_down,
    lambda_up, C_up, h_up, eta_up,
    r_rms, C_probe, A_pump, xoffset,
    air_lens=None, reflection=None,
    n_k=100
):
    """
    Computes the complex deflection angle theta(f) using the isotropic free thermal expansion model.
    f: array-like or scalar frequencies (Hz)
    air_lens: dict with keys 'enable' (bool) and 'coef_air' (float)
    reflection: dict with keys 'enable' (bool), 'dphidT' (float), 'wavelength' (float)
    n_k: number of Gauss-Legendre quadrature points for k-space integration
    """
    if air_lens is None:
        air_lens = {'enable': False, 'coef_air': 0.0}
    if reflection is None:
        reflection = {'enable': False, 'dphidT': 0.0, 'wavelength': 0.0}

    f_arr = np.atleast_1d(f)
    r_pump = r_rms
    r_probe = r_rms
    kmax = 2.0 / np.sqrt(r_pump**2 + r_probe**2)

    nodes, weights = np.polynomial.legendre.leggauss(n_k)
    kvect = 0.5 * kmax * (nodes + 1.0)
    w_scaled = 0.5 * kmax * weights

    kvect_2d = kvect[None, :]  # shape (1, n_k)
    f_2d = f_arr[:, None]      # shape (N_f, 1)

    k = 2.0 * np.pi * kvect_2d
    q2_bulk = 1j * 2.0 * np.pi * f_2d * C_down[2] / lambda_down[2]
    q2_metal = 1j * 2.0 * np.pi * f_2d * C_down[0] / lambda_down[0]

    G_int = lambda_down[1] / h_down[1]
    L_metal = h_down[0]

    zeta_metal = np.sqrt(q2_metal + k**2)
    zeta_metal_L = zeta_metal * L_metal
    gamma_metal = lambda_down[0] * zeta_metal

    conv = np.cosh(zeta_metal_L) + (gamma_metal / G_int) * np.sinh(zeta_metal_L)
    corr = np.sinh(zeta_metal_L) / gamma_metal + np.cosh(zeta_metal_L) / G_int

    T_s_Sprobe = bifdtr_bo_temp(
        kvect_2d, f_2d, lambda_up, C_up, h_up, eta_up,
        lambda_down, C_down, h_down, eta_down,
        r_pump, r_probe, A_pump
    )

    flx = A_pump * np.exp(-(np.pi**2) * (r_pump**2) / 2.0 * (kvect_2d**2))
    S_probe = np.exp(-(np.pi**2) * (r_probe**2) / 2.0 * (kvect_2d**2))

    T_bs_Sprobe = conv * T_s_Sprobe - corr * flx * S_probe

    integrand = (-C_probe * 8.0 * (np.pi**2) * (kvect_2d**2) *
                 (-jv(1, 2.0 * np.pi * kvect_2d * xoffset)) *
                 (2.0 * (1.0 + niu) * alpha_T /
                  (np.sqrt(4.0 * (np.pi**2) * eta_down[2] * (kvect_2d**2) + q2_bulk) + 2.0 * np.pi * kvect_2d)) *
                 T_bs_Sprobe)

    if air_lens.get('enable', False):
        q2_air = 1j * 2.0 * np.pi * f_2d * C_up / lambda_up
        zeta_air = np.sqrt(q2_air + k**2)
        Z_air_Sprobe = (air_lens['coef_air'] / zeta_air) * T_s_Sprobe
        integrand_al = (-C_probe * 8.0 * (np.pi**2) * (kvect_2d**2) *
                        (-jv(1, 2.0 * np.pi * kvect_2d * xoffset)) * Z_air_Sprobe)
        integrand += integrand_al

    if reflection.get('enable', False):
        Z_r_Sprobe = -reflection['wavelength'] * reflection['dphidT'] * T_s_Sprobe / (4.0 * np.pi)
        integrand_ref = (-C_probe * 8.0 * (np.pi**2) * (kvect_2d**2) *
                         (-jv(1, 2.0 * np.pi * kvect_2d * xoffset)) * Z_r_Sprobe)
        integrand += integrand_ref

    theta = np.dot(integrand, w_scaled)

    if np.ndim(f) == 0:
        return theta[0]
    return theta
