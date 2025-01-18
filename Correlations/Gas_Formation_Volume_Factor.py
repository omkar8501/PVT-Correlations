import Gas_Compressibility_Factor as Gas_Comp
from typing import Optional


def gas_formation_volume_factor(pressure: float, temp: float, sg_gas: float, gas_comp_factor: Optional[float] = None) -> float:
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
        gas_comp_factor : float, optional
            Compressibility Factor of the gas

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
        >>> gas_formation_volume_factor(pressure=1000, temp=617.67, sg_gas=0.7)
        np.float64(0.015449190149601931)
        """
    # Calculate compressibility factor of the gas if not provided by the user
    gas_comp_factor = gas_comp_factor if gas_comp_factor is not None else Gas_Comp.carnahan_starling_hs_eos(pressure, temp, sg_gas)

    gas_fvf: float = 0.02827 * gas_comp_factor * temp / pressure
    return gas_fvf
