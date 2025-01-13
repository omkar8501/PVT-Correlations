import numpy as np
import Gas_Compressibility_Factor as Gas_Comp
from typing import Optional


def lee_gonzalez_eakin(pressure: float, temp: float, sg_gas: float, gas_comp_factor: Optional[float] = None) -> float:
    """
        Calculates theviscosity of the gas at a given pressure and temperature using the Lee-Gonzalez-Eakin method

        The Lee-Gonzalez-Eakin Method is a semi-empirical correlation used to estimate the viscosity of gas mixtures.
        It incorporates reduced temperature and molecular weight to predict viscosity, offering reliable results across a range of conditions.

        Parameters:
        ----------
        pressure : float
            Pressure of the gas in psia
        temp : float
            Temperature of the gas in degree Rankine
        sg_gas : float
            Specific gravity of the gas relative to air
        gas_comp_factor : float, optional
            Compressibility Factor of the gas

        Returns:
        -------
        float
            Viscosity of the gas in cp

        Notes:
        ------
        - The correlation is less accurate for gases with higher specific gravities
        - The method cannot be used for sour gases

        Example:
        --------
        >>> lee_gonzalez_eakin(pressure=2000, temp=600, sg_gas=0.78)
        np.float64(0.017746575283174415)
        """

    # Constants
    mw_air: float = 28.96  # Molecular weight of air in lb/lb-mol
    gas_const: float = 10.73159  # Gas Constant in psia-ft3/R-mol

    # Apparent Molecular Weight of the gas mixture
    mw_gas: float = sg_gas * mw_air

    # Gas Compressibility
    gas_comp_factor = gas_comp_factor if gas_comp_factor is not None else Gas_Comp.carnahan_starling_hs_eos(pressure,
                                                                                                            temp,
                                                                                                            sg_gas)

    # Density of the gas mixture calculated using the Real Gas Equation
    gas_density: float = pressure * mw_gas / (gas_comp_factor * gas_const * temp)

    # Parameters used in the equation
    k_param: float = (9.4 + 0.02 * mw_gas) * np.power(temp, 1.5) / (209 + 19 * mw_gas + temp)
    x_param: float = 3.5 + 986 / temp + 0.01 * mw_gas
    y_param: float = 2.4 - 0.2 * x_param

    # Viscosity equation
    gas_visc: float = 0.0001 * k_param * np.exp(x_param * np.power(gas_density / 62.4, y_param))

    return gas_visc
