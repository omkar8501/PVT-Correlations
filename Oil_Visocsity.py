import numpy as np
import Gas_Solubility as Gas_Sol
from typing import Optional


def beggs_robinson(pressure: float, temp: float, oil_api: float, p_bubble: float, sg_gas: Optional[float] = None,
                   gas_comp_factor: Optional[float] = None, p_sep: Optional[float] = None,
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
    sg_gas : float, optional
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
