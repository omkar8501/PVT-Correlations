import numpy as np


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
    >>> standings_oil_fvf(gas_sol=100, temp=617.67, oil_api=35, sp_gr_gas=0.8)
    1.0891
    """

    # Convert the density of oil from deg API to Sp. Gr.
    sg_oil: float = 141.5 / (oil_api + 131.5)

    # Calculate the value of Oil Formation Volume Factor using Standing's Correlation
    oil_fvf = 0.9759 + 0.00012 * np.power(gas_sol * np.sqrt(sg_gas / sg_oil) + 1.25 * (temp - 460.67), 1.2)

    return oil_fvf
