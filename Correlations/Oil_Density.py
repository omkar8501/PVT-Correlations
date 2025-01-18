import Gas_Solubility as Gas_Sol
import Oil_Formation_Volume_Factor as Oil_FVF
from typing import Optional


def oil_density(pressure: float, temp: float, oil_api: float, sg_gas: float, p_bubble: float,
                gas_sol: Optional[float] = None, oil_fvf: Optional[float] = None, p_sep: Optional[float] = None,
                t_sep: Optional[float] = None) -> float:
    """
        Calculates the density of the oil at the given conditions.

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
        gas_sol: float, optional
            Solubility of gas in the oil in scf/bbl
        oil_fvf: float, optional
            Formation Volume Factor of the oil in rb/stb
        p_sep : float, optional
             Actual pressure of the separator in psia
        t_sep : float, optional
            Actual temperature of the separator in degrees Rankine

        Returns:
        -------
        float
            Density of the oil in lbm/ft3

        Example:
        --------
        >>> oil_density(pressure=1000, temp=620, oil_api=35, sg_gas=0.7, p_bubble=1500, p_sep=114.7, t_sep=520)
        np.float64(48.06738440645158)
    """

    # Gas Solubility and Formation Volume Factor of oil
    gas_sol = gas_sol if gas_sol is not None else Gas_Sol.vasquez_beggs_gas_solubility(pressure, temp, oil_api, sg_gas,
                                                                                       p_bubble, p_sep, t_sep)
    oil_fvf = oil_fvf if oil_fvf is not None else Oil_FVF.vasquez_beggs_oil_fvf(pressure, temp, oil_api, sg_gas,
                                                                                p_bubble, p_sep, t_sep)

    # Specific Gravity of the oil
    sg_oil: float = 141.5 / (131.5 + oil_api)

    # Density of the oil
    oil_den: float = (62.4 * sg_oil + 0.0136 * gas_sol * sg_gas) / oil_fvf

    return oil_den
