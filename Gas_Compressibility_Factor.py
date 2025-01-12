import numpy as np


def carnahan_starling_hs_eos(pressure: float, temp: float, sg_gas: float) -> float:
    """
    Calculates the compressibility of the gas at a given pressure and temp_rerature using the Hall-Yarborough method

    The Carnahan-Starling Hard Sphere equation of state (EoS) is a model that describes the behavior of fluids by
    considering the molecules as hard spheres with finite volume. When solved using the Hall-Yarborough method,
    a numerical technique typically applied to cubic equations of state, it provides an efficient and accurate way to determine fluid properties

    Parameters:
    ----------
    pressure: float
        Pressure of the gas in psia
    temp: float
        Temperature of the gas in degree Rankine
    sg_gas: float
        Specific gravity of the gas relative to air

    Returns:
    -------
    float
        Compressibility of the gas

    Notes:
    ------
    The method is not recommended for application if the pseudo-reduced temperature is less than one

    Example:
    --------
    >>>carnahan_starling_hs_eos(pressure=1000, temp=500, sg_gas=0.7)
    0.7384
    """

    # Standing's equation to Calculate critical Pressure and temperature
    sg_vector: np.ndarray = np.array([[1], [sg_gas], [np.square(sg_gas)]])
    dry_gas_coeff_vector: np.ndarray = np.array([[168, 325, -12.5], [667, 15, -37.5]])
    wet_gas_coeff_vector: np.ndarray = np.array([[187, 330, -71.5], [706, -51.7, -11.1]])
    result = np.dot(dry_gas_coeff_vector, sg_vector) if sg_gas < 0.75 else np.dot(wet_gas_coeff_vector, sg_vector)
    temp_pc: float = result[0][0]
    pressure_pc: float = result[1][0]

    # Calculate the Pseudo-reduced Pressure and temperature
    temp_pr: float = temp / temp_pc
    pressure_pr: float = pressure / pressure_pc

    # Calculate the intermediate values
    theta: float = 1 / temp_pr
    alpha: float = 0.06125 * theta * np.exp(-1.2 * np.square(1 - theta))

    # Compute the reduced density parameter Rho Hat (rho_h) using Newton-Raphson Method
    rho_h_conv: float = 1.0
    rho_h_i: float = 0.8
    tol = 0.001  # Tolerance
    rho_h: float = rho_h_i
    for i in range(100):
        # Calculate the reduced density functions using the Hall-Yarborough Method
        f1_rho_h: float = -alpha * pressure_pr + (
                rho_h + np.square(rho_h) + np.power(rho_h, 3) - np.power(rho_h, 4)) / np.power(1 - rho_h, 3)
        f2_rho_h: float = -(14.76 * theta - 9.76 * np.square(theta) + 4.58 * np.power(theta, 3)) * np.square(rho_h)
        f3_rho_h: float = (90.7 * theta - 242.2 * np.square(theta) + 42.4 * np.power(theta, 3)) * np.power(rho_h,
                                                                                                           2.18 + 2.82 * theta)
        # Calculate the derivative of the reduced density functions
        f1_dash_rho_h: float = (1 + 4 * rho_h + 4 * np.square(rho_h) - 4 * np.power(rho_h, 3) + np.power(rho_h,
                                                                                                         4)) / np.power(
            1 - rho_h, 4)
        f2_dash_rho_h: float = 2 * f2_rho_h / rho_h
        f3_dash_rho_h: float = (2.18 + 2.82 * theta) * f3_rho_h / rho_h

        # Add up the functions to calculate the reduced density function and its derivative
        f_rho_h: float = f1_rho_h + f2_rho_h + f3_rho_h
        f_dash_rho_h: float = f1_dash_rho_h + f2_dash_rho_h + f3_dash_rho_h

        # Calculate the next value of the reduced density parameter
        rho_h_next = rho_h - f_rho_h / f_dash_rho_h
        if abs(rho_h_next - rho_h) < tol:
            rho_h_conv = rho_h_next
            break
        rho_h = rho_h_next

    compressibility: float = alpha * pressure_pr / rho_h_conv
    return compressibility
