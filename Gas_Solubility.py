import numpy as np
from typing import List


def standings_gas_solubility(pressure: float, temp: float, oil_api: float, sg_gas: float) -> float:
    """
    Calculates the gas solubility (Rs) in oil using Standing's Correlation.

    Standing's Correlation is an empirical relationship widely used in petroleum engineering
    to estimate the solubility of natural gas in crude oil under given pressure and temperature
    conditions.

    Parameters:
    ----------
    pressure : float
        Pressure of the fluid in psia
    temp : float
        Temperature of the fluid in degrees Rankine
    oil_api : float
        API gravity of the oil, indicating its density (°API)
    sg_gas : float
        Specific gravity of the gas relative to air

    Returns:
    -------
    float
        Gas solubility (Rs) in the oil at the given conditions, expressed in standard cubic
        feet of gas per barrel of oil (scf/bbl)

    Notes:
    ------
    - This correlation assumes the oil is undersaturated (i.e., pressure below bubble-point pressure)
    - The temperature is converted to Rankine for internal calculations
    - Ensure input values are within the correlation's validity range to maintain accuracy

    Example:
    --------
    >>> standings_gas_solubility(pressure=1000, temp=600, oil_api=35, sg_gas=0.7)
    213.2417
    """

    exponent: float = 0.0125 * oil_api - 0.00091 * (temp - 460.67)
    solubility = sg_gas * np.power((pressure / 18.2 + 1.4) * np.power(10, exponent), 1.2048)

    return solubility


def vasquez_beggs_gas_solubility(pressure: float, temp: float, oil_api: float, sg_gas: float, p_bubble: float,
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
            Pressure of the fluid in psia
        temp : float
            Temperature of the fluid in degrees Rankine
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
            Gas solubility (Rs) in the oil at the given conditions, expressed in standard cubic
            feet of gas per barrel of oil (scf/bbl)

        Notes:
        ------
        - The author realised that the value of the specific gravity of the gas depends on the conditions under which
        it is separated from the oil, thus an adjustment is made in the value of the specific gravity

        Example:
        --------
        >>> vasquez_beggs_gas_solubility(pressure=1000, temp=620, oil_api=35, sg_gas=0.7, p_bubble=1500, p_sep=114.7, t_sep=520)

        """

    # Coefficients used in Vasquez-Beggs Correlation
    coeff_heavy: List[float] = [0.0362, 1.0937, 25.724]
    coeff_light: List[float] = [0.0178, 1.187, 23.931]
    coeff: List[float] = coeff_heavy if oil_api <= 30 else coeff_light

    # Adjust the gas gravity for separator conditions
    sg_gas_sep: float = sg_gas * (1 + 5.912 * 0.00001 * oil_api * (t_sep - 460.67) * np.log10(p_sep / 114.7))

    # Calculate the gas solubility using the Vasquez-Beggs Correlation
    # Use given pressure for saturated oil and bubble point pressure for under-saturated oil
    gas_sol: float = coeff[0] * sg_gas_sep * np.power(min(pressure, p_bubble), coeff[1]) * np.exp(
        coeff[2] * oil_api / temp)

    return gas_sol
