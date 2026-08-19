import numpy as np
from scipy.special import jv
from scipy.integrate import quad


def bifdtr_bo_temp(
    kvect_in, freq,
    lambda_up, C_up, h_up, eta_up,
    lambda_down, C_down, h_down, eta_down,
    r_pump, r_probe, A_pump
):
    """
    Computes the frequency- and wavevector-domain surface temperature T_s for a multilayer sample.
    Works for scalar or array kvect_in.
    """
    kvect = np.asarray(kvect_in, dtype=complex)
    kvect2 = kvect**2

    # Liquid/upper layer
    alpha_up = lambda_up / C_up
    omega = 2.0 * np.pi * freq
    q2_up = 1j * omega / alpha_up

    un_up = np.sqrt(4.0 * (np.pi**2) * eta_up * kvect2 + q2_up)
    gamman_up = lambda_up * un_up
    g_up = 1.0 / gamman_up

    # Substrate/down layers
    n_layers = len(lambda_down)
    alpha_down = np.asarray(lambda_down, dtype=float) / np.asarray(C_down, dtype=float)

    q2_sub = 1j * omega / alpha_down[-1]
    un = np.sqrt(4.0 * (np.pi**2) * eta_down[-1] * kvect2 + q2_sub)
    gamman = lambda_down[-1] * un

    shape = un.shape
    b_plus = np.zeros(shape, dtype=complex)
    b_minus = np.ones(shape, dtype=complex)

    if n_layers > 1:
        for n in range(n_layers - 1, 0, -1):
            q2_prev = 1j * omega / alpha_down[n - 1]
            un_minus = np.sqrt(eta_down[n - 1] * 4.0 * (np.pi**2) * kvect2 + q2_prev)
            gamman_minus = lambda_down[n - 1] * un_minus

            aa = gamman_minus + gamman
            bb = gamman_minus - gamman

            temp1 = aa * b_plus + bb * b_minus
            temp2 = bb * b_plus + aa * b_minus

            expterm = np.exp(un_minus * h_down[n - 1])

            b_plus = (0.5 / (gamman_minus * expterm)) * temp1
            b_minus = (0.5 / gamman_minus) * expterm * temp2

            # Numerical stability check for thick or resistive layers
            penetration_logic = (h_down[n - 1] * np.abs(un_minus)) > 100
            if np.any(penetration_logic):
                if np.isscalar(penetration_logic):
                    if penetration_logic:
                        b_plus = 0.0
                        b_minus = 1.0
                else:
                    b_plus[penetration_logic] = 0.0
                    b_minus[penetration_logic] = 1.0

            un = un_minus
            gamman = gamman_minus

    g_down = (b_plus + b_minus) / (b_minus - b_plus) / gamman
    g = g_up * g_down / (g_up + g_down)

    s = np.exp(-(np.pi**2) * (r_probe**2) / 2.0 * kvect2)
    p = A_pump * np.exp(-(np.pi**2) * (r_pump**2) / 2.0 * kvect2)
    kernel = s * p
    integrand = g * kernel

    if np.ndim(kvect_in) == 0:
        return integrand.item() if np.iscomplexobj(integrand) else float(integrand)
    return integrand


def ss_heat(
    lambda_down, C_down, h_down, eta_down,
    lambda_up, C_up, h_up, eta_up,
    r_rms, A_pump, xoffset,
    n_k=100
):
    """
    Estimates steady-state (DC) surface heating from the multilayer thermal model.
    Uses Gauss-Legendre quadrature over k-space for fast vectorized integration.
    """
    r_pump = r_rms
    r_probe = r_rms
    kmax = 2.0 / np.sqrt(r_pump**2 + r_probe**2)
    kmin = 1.0 / (10000.0 * max(r_pump, r_probe))

    freq = 0.0
    nodes, weights = np.polynomial.legendre.leggauss(n_k)
    kvect = kmin + 0.5 * (kmax - kmin) * (nodes + 1.0)
    w_scaled = 0.5 * (kmax - kmin) * weights

    vals = kvect * jv(0, 2.0 * np.pi * kvect * xoffset) * bifdtr_bo_temp(
        kvect, freq,
        lambda_up, C_up, h_up, eta_up,
        lambda_down, C_down, h_down, eta_down,
        r_pump, r_probe, A_pump
    )
    return float(np.real(np.dot(vals, w_scaled)))
