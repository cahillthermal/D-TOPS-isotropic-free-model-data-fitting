import numpy as np


def get_data_out_in_ratio_f_vsum(file_path: str):
    """
    Reads a 4-column lock-in data file (V_in, V_out, f, V_SUM).
    Handles filenames provided with or without the '.txt' extension.
    Returns:
        V_out, V_in, V_ratio, V_SUM, f
    """
    if not file_path.endswith('.txt'):
        file_path = f"{file_path}.txt"

    data = np.loadtxt(file_path)
    v_in = data[:, 0]
    v_out = data[:, 1]
    f = data[:, 2]
    v_sum = data[:, 3]
    v_ratio = -v_in / v_out

    return v_out, v_in, v_ratio, v_sum, f


def datacorrection_complex_leaking(v_out_data, v_in_data, complex_leaking):
    """
    Corrects raw lock-in voltage data for pump-modulation/detector frequency-response ("leaking").
    Returns:
        Vcorrected_in, Vcorrected_out, Vcorrected_ratio
    """
    v_complex_data = v_in_data + 1j * v_out_data
    v_corrected_complex = v_complex_data / complex_leaking
    v_corrected_in = np.real(v_corrected_complex)
    v_corrected_out = np.imag(v_corrected_complex)
    v_corrected_ratio = -v_corrected_in / v_corrected_out

    return v_corrected_in, v_corrected_out, v_corrected_ratio
