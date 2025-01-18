import Gas_Solubility as Gas_Sol
import numpy as np
from typing import List


def standings_oil_fvf(gas_sol: float, temp: float, oil_api: float, sg_gas: float) -> float:
    """
    Calculates the formation volume factor of the oil using Standing's Correlation.

    Standing's correlation is a widely used empirical relationship to estimate the oil formation volume factor (Bo)
    in petroleum engineering. It relates Bo to reservoir pressure, temperature, solution gas-oil ratio (Rs),
    gas-specific gravity, and oil gravity. This correlation is particularly effective for light to medium crude oils
    within specific pressure and temperature ranges.

    Parameters:
    ----------
    gas_sol : float
        Solubility of natural gas in oil at a given pressure in scf/bbl
    temp : float
        Temperature of the oil in degrees Rankine
    oil_api : float
        API gravity of the oil, indicating its density (°API)
    sp_gr_gas : float
        Specific gravity of the gas relative to air

    Returns:
    -------
    float
        Formation Volume Factor of the oil at the given conditions, expressed in reservoir barrel per
        barrel of oil at standard conditions (rb/stb)

    Notes:
    ------
    - This correlation is valid for light to medium crude oils and may not be accurate for heavy oils or unconventional resources.
    - It is typically used for pressures below the bubble point and within a specific temperature range (typically 50°C to 150°C).
    - The accuracy of this correlation decreases for extreme conditions, such as very high or very low reservoir temperatures or pressures.

    Example:
    --------
    >>> standings_oil_fvf(gas_sol=100, temp=617.67, oil_api=35, sg_gas=0.8)
    1.0855244491103584
    """

    # Convert the density of oil from deg API to Sp. Gr.
    sg_oil: float = 141.5 / (oil_api + 131.5)

    # Calculate the value of Oil Formation Volume Factor using Standing's Correlation
    oil_fvf = 0.9759 + 0.00012 * np.power(gas_sol * np.sqrt(sg_gas / sg_oil) + 1.25 * (temp - 460.67), 1.2)

    return oil_fvf


def vasquez_beggs_oil_fvf(pressure: float, temp: float, oil_api: float, sg_gas: float, p_bubble: float,
                          p_sep: float, t_sep: float) -> float:
    """
       Calculates the gas solubility (Rs) in oil using Vasquez-Begg's Correlation.

       The Vasquez-Begg correlation is a widely used empirical method to estimate the solubility of natural gas in
       crude oil under reservoir conditions. It relates gas solubility to temperature, pressure, and oil properties,
       such as API gravity and solution gas-oil ratio, providing reliable predictions for reservoir engineering with
       an average absolute error of 12.7%.

       Parameters:
        ----------
        pressure : float
            Pressure of the oil in psia
        temp : float
            Temperature of the oil in degrees Rankine
        oil_api : float
            API gravity of the oil, indicating its density (°API)
        sg_gas : float
            Specific gravity of the gas relative to air
        p_bubble: float
            Bubble point pressure of the oil in psia
        p_sep : float
             Actual pressure of the separator in psia
        t_sep : float
            Actual temperature of the separator in degrees Rankine

       Returns:
       -------
       float
           Formation Volume Factor of the oil at the given conditions, expressed in reservoir barrel per barrel of oil
           at standard conditiosn (rb/stb)

       Notes:
       ------
       - The author realised that the value of the specific gravity of the gas depends on the conditions under which
       it is separated from the oil, thus an adjustment is made in the value of the specific gravity

       Example:
       --------
       >>> vasquez_beggs_oil_fvf(pressure=1000, temp=620, oil_api=35, sg_gas=0.7, p_bubble=1500, p_sep=114.7, t_sep=520)

    """

    # Coefficients used in Vasquez-Beggs Correlation
    coeff_heavy: List[float] = [4.677 * np.power(10.0, -4), 1.751 * np.power(10.0, -5), -1.811 * np.power(10.0, -8)]
    coeff_light: List[float] = [4.67 * np.power(10.0, -4), 1.1 * np.power(10.0, -5), 1.337 * np.power(10.0, -9)]
    coeff: List[float] = coeff_heavy if oil_api <= 30 else coeff_light

    # Adjust the gas gravity for separator conditions
    sg_gas_sep: float = sg_gas * (1 + 5.912 * 0.00001 * oil_api * (t_sep - 460.67) * np.log10(p_sep / 114.7))

    # Calculate the solubility of gas at the given pressure temperature
    gas_sol: float = Gas_Sol.vasquez_beggs_gas_solubility(pressure, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep)
    # Calculate the solubility of gas in the oil and oil FVF at bubble point
    gas_sol_bp: float = Gas_Sol.vasquez_beggs_gas_solubility(p_bubble, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep)
    oil_fvf_pb: float = 1 + coeff[0] * gas_sol_bp + (temp - 520) * (oil_api / sg_gas_sep) * (
            coeff[1] + coeff[2] * gas_sol_bp)

    if pressure < p_bubble:
        oil_fvf: float = 1 + coeff[0] * gas_sol + (temp - 520) * (oil_api / sg_gas_sep) * (
                coeff[1] + coeff[2] * gas_sol)
    else:
        a_coeff: float = np.power(10.0, -5) * (
                -1433 + 5 * gas_sol_bp + 17.2 * (temp - 460) - 1180 * sg_gas_sep + 12.61 * oil_api)
        oil_fvf: float = oil_fvf_pb * np.exp(-a_coeff * np.log(pressure / p_bubble))

    return oil_fvf
