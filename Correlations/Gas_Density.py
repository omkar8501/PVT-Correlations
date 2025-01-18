import Constants as const


def gas_density(pressure: float, temp: float, sg_gas: float) -> float:
    """
        Calculates the density of the gas at a given pressure temperature using the Ideal Gas Law

        Parameters:
        ----------
        pressure : float
            Pressure of the gas in psia
        temp : float
            Temperature of the gas in degree Rankine
        sg_gas : float
            Specific gravity of the gas relative to air

        Returns:
        -------
        float
            Density of the gas in lbm/ft3

        Example:
        --------
        >>> gas_density(pressure=2000, temp=600, sg_gas=0.78)
        7.0162948826781495
    """

    # Molecular weight of the gas
    mw_gas: float = sg_gas * const.mw_air

    # Ideal Gas Law
    gas_den: float = pressure * mw_gas / (const.gas_const * temp)

    return gas_den
