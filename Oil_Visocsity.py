import numpy as np
import Gas_Solubility as Gas_Sol
from typing import Optional


def beggs_robinson(pressure: float, temp: float, oil_api: float, sg_gas: float, p_bubble: float,
                   gas_sol: Optional[float] = None, p_sep: Optional[float] = None,
                   t_sep: Optional[float] = None) -> float:
    """
    Calculates the viscosity of oil using Beggs-Robinson Correlation, viscosity above bubble point is calculated using
    Vasquez-Beggs Correlation



    Parameters:
    ----------
    pressure : float
        Pressure of the oil in psia
    temp : float
        Temperature of the oil in degrees Rankine
    oil_api : float
        API gravity of the oil, indicating its density (°API)
    p_bubble: float
        Bubble point pressure of the oil in psia
    sg_gas : float
        Specific gravity of the gas relative to air
    gas_comp_factor: float, optional
        Compressibility Factor of the gas
    p_sep : float, optional
         Actual pressure of the separator in psia
    t_sep : float, optional
        Actual temperature of the separator in degrees Rankine

    Returns:
    -------
    float
        Viscosity of the oil at the given conditions, expressed in centipoise (cp)

    Notes:
    ------
    - The author realised that the value of the specific gravity of the gas depends on the conditions under which
    it is separated from the oil, thus an adjustment is made in the value of the specific gravity

    Example:
    --------
    >>> beggs_robinson(pressure=1000, temp=620, oil_api=35, sg_gas=0.7, p_bubble=1500, p_sep=114.7, t_sep=520)

        """

    # Viscosity of dead oil
    y_param: float = np.power(10.0, 3.0324 - 0.02023 * oil_api)
    x_param: float = y_param * np.power(temp - 460, -1.163)
    oil_visc_dead: float = np.power(10.0, x_param) - 1

    # Use Beggs-Robinson Correlation for calculating viscosity of saturated oil
    if pressure <= p_bubble:
        gas_sol = gas_sol if gas_sol is not None else Gas_Sol.vasquez_beggs_gas_solubility(pressure, temp, oil_api,
                                                                                           sg_gas, p_bubble, p_sep,
                                                                                           t_sep)
        a_param: float = 10.715 * np.power(gas_sol + 100, -0.515)
        b_param: float = 5.44 * np.power(gas_sol + 150, -0.338)
        oil_visc: float = a_param * np.power(oil_visc_dead, b_param)
    # Use Vasquez-Beggs Correlation for calculating viscosity of undersaturated oil
    else:
        # First calculate viscosity of oil at bubble point using Beggs-Robinson Correlation
        gas_sol_pb: float = gas_sol if gas_sol is not None else Gas_Sol.vasquez_beggs_gas_solubility(p_bubble, temp,
                                                                                                     oil_api, sg_gas,
                                                                                                     p_bubble, p_sep,
                                                                                                     t_sep)
        oil_visc_pb: float = beggs_robinson(p_bubble, temp, oil_api, sg_gas, p_bubble, gas_sol=gas_sol_pb, p_sep=p_sep,
                                            t_sep=t_sep)
        # Vasquez-Beggs Correlation
        n_param: float = -3.9 * 0.00001 * pressure - 5
        m_param: float = 2.6 * np.power(pressure, 1.187) * np.power(10.0, n_param)
        oil_visc: float = oil_visc_pb * np.power(pressure / p_bubble, m_param)

    return oil_visc
