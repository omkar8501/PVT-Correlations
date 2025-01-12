def gas_formation_volume_factor(pressure: float, temp: float, gas_comp_factor: float) -> float:
    """
        Calculates the formation volume factor of the gas at the given conditions.

        Parameters:
        ----------
        pressure : float
            Pressure of the gas in psia
        temp : float
            Temperature of the gas in degrees Rankine
        gas_comp_factor : float
            Compressibility factor of the gas

        Returns:
        -------
        float
            Formation Volume Factor of the gas at the given conditions, expressed in reservoir cubic feet per
            cubic feet of gas at standard conditions (rcf/scf)

        Notes:
        ------
        - Value of Bg obtained is dependent on the correlation used to calculate compressibility factor

        Example:
        --------
        >>> gas_formation_volume_factor(pressure=1000, temp=617.67, gas_comp_factor=0.9)

        """

    gas_fvf: float = 0.02827 * gas_comp_factor * temp / pressure
    return gas_fvf
