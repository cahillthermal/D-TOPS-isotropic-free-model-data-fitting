import numpy as np
from scipy.optimize import least_squares
from deflection_model import theta_iso_free_thermal_expansion_model


def fit_inout(
    f, theta_exp_in, theta_exp_out,
    x_guess, bounds,
    niu, lambda_down, C_down, h_down, eta_down,
    lambda_up, C_up, h_up, eta_up,
    r_rms, C_probe, A_pump, xoffset,
    air_lens=None, reflection=None
):
    """
    Fits bulk thermal conductivity lambda_down[2] and CTE alpha_T simultaneously
    using in-phase and out-of-phase theta signals (FDPBD_fitting1).
    """
    lambda_down = list(lambda_down)
    max_theta_exp_in = np.max(np.abs(theta_exp_in))

    theta_exp_all = np.concatenate([theta_exp_out, theta_exp_in])
    scale_weights = np.concatenate([
        np.full_like(theta_exp_out, max_theta_exp_in),
        np.full_like(theta_exp_in, 3.0 * max_theta_exp_in)
    ])
    y_data_norm = theta_exp_all / scale_weights

    def residual_func(x):
        ld_copy = list(lambda_down)
        ld_copy[2] = x[0]
        alpha_t = x[1]

        theta_model = theta_iso_free_thermal_expansion_model(
            niu, alpha_t, f,
            ld_copy, C_down, h_down, eta_down,
            lambda_up, C_up, h_up, eta_up,
            r_rms, C_probe, A_pump, xoffset,
            air_lens=air_lens, reflection=reflection
        )
        theta_model_out = np.imag(theta_model)
        theta_model_in = np.real(theta_model)
        theta_model_all = np.concatenate([theta_model_out, theta_model_in])
        theta_model_norm = theta_model_all / scale_weights

        return theta_model_norm - y_data_norm

    lb, ub = bounds
    res = least_squares(residual_func, x0=x_guess, bounds=(lb, ub), verbose=0) # changed verbose from 1 to 0 to suppress printing

    x_sol = res.x
    # Approximate covariance and confidence interval calculation
    cov = None
    confidence_interval = None
    perr = np.array([np.nan, np.nan])
    if res.jac is not None and res.jac.size > 0:
        try:
            # J^T J
            jtj = res.jac.T @ res.jac
            inv_jtj = np.linalg.inv(jtj)
            s_sq = np.sum(res.fun**2) / (len(res.fun) - len(x_sol))
            cov = inv_jtj * s_sq
            perr = np.sqrt(np.diag(cov))
            confidence_interval = np.column_stack([x_sol - 1.96 * perr, x_sol + 1.96 * perr])
        except np.linalg.LinAlgError:
            pass

    return x_sol, res, confidence_interval, perr


def fit_ratio(
    f, v_corrected_ratio,
    x_guess,
    niu, alpha_t, lambda_down, C_down, h_down, eta_down,
    lambda_up, C_up, h_up, eta_up,
    r_rms, C_probe, A_pump, xoffset,
    air_lens=None, reflection=None
):
    """
    Fits bulk thermal conductivity lambda_down[2] using ratio signal (FDPBD_fitting2).
    """
    lambda_down = list(lambda_down)

    def residual_func(x):
        ld_copy = list(lambda_down)
        ld_copy[2] = x[0]

        theta_model = theta_iso_free_thermal_expansion_model(
            niu, alpha_t, f,
            ld_copy, C_down, h_down, eta_down,
            lambda_up, C_up, h_up, eta_up,
            r_rms, C_probe, A_pump, xoffset,
            air_lens=air_lens, reflection=reflection
        )
        theta_model_out = np.imag(theta_model)
        theta_model_in = np.real(theta_model)
        theta_model_ratio = -theta_model_in / theta_model_out

        return theta_model_ratio - v_corrected_ratio

    res = least_squares(residual_func, x0=[x_guess], verbose=1)
    x_sol = res.x[0]

    cov = None
    perr = np.nan
    if res.jac is not None and res.jac.size > 0:
        try:
            jtj = res.jac.T @ res.jac
            inv_jtj = np.linalg.inv(jtj)
            s_sq = np.sum(res.fun**2) / max(len(res.fun) - 1, 1)
            cov = inv_jtj * s_sq
            perr = np.sqrt(np.diag(cov))[0]
        except np.linalg.LinAlgError:
            pass

    return x_sol, res, perr
