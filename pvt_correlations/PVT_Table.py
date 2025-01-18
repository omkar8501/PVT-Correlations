import numpy as np
import pandas as pd
from typing import Dict, Optional

import Wrappers
import Gas_Density as Gas_Den


def pvt_table(pressure_max: float, p_bubble: float, temp: float, oil_api: float, sg_gas: float, p_sep: float,
              t_sep: float, num_points: int, gas_sol_corr: Optional[str] = "Vasquez Beggs",
              gas_comp_corr: Optional[str] = "Carnahan Starling", gas_visc_corr: Optional[str] = "Lee Gonzalez Eakin",
              oil_fvf_corr: Optional[str] = "Vasquez Beggs",
              oil_visc_corr: Optional[str] = "Beggs Robinson", **kwargs) -> pd.DataFrame:
    gas_sol_correlations: Dict = {
        "Standing": Wrappers.standings_gas_solubility_wrapper,
        "Vasquez Beggs": Wrappers.vasquez_beggs_gas_solubility_wrapper
    }

    gas_comp_correlations: Dict = {
        "Carnahan Starling": Wrappers.carnahan_starling_hs_eos_wrapper
    }

    gas_visc_correlations: Dict = {
        "Lee Gonzalez Eakin": Wrappers.lee_gonzalez_eakin_wrapper
    }

    oil_fvf_correlations: Dict = {
        "Standing": Wrappers.standings_oil_fvf_wrapper,
        "Vasquez Beggs": Wrappers.vasquez_beggs_oil_fvf_wrapper
    }

    oil_visc_correlations: Dict = {
        "Beggs Robinson": Wrappers.beggs_robinson_wrapper
    }

    # Create a range of pressures to calculate parameters at
    p_range: np.ndarray = np.linspace(0, pressure_max, num_points)

    # Create an empty table with just the pressures
    pvt_tab: pd.DataFrame = pd.DataFrame(p_range, columns=['Pressure'])

    # Get all the required correlations
    rs_corr = gas_sol_correlations.get(gas_sol_corr)
    z_corr = gas_comp_correlations.get(gas_comp_corr)
    mu_g_corr = gas_visc_correlations.get(gas_visc_corr)
    bo_corr = oil_fvf_correlations.get(oil_fvf_corr)
    mu_o_corr = oil_visc_correlations.get(oil_visc_corr)

    # Add all the fluid properties
    pvt_tab['Gas Solubility'] = pvt_tab['Pressure'].apply(
        lambda p: rs_corr(pressure=p, temp=temp, oil_api=oil_api, sg_gas=sg_gas, p_bubble=p_bubble, p_sep=p_sep,
                          t_sep=t_sep))
    pvt_tab['Gas Compressibility Factor'] = pvt_tab['Pressure'].apply(
        lambda p: z_corr(pressure=p, temp=temp, sg_gas=sg_gas))
    pvt_tab['Gas Density'] = pvt_tab['Pressure'].apply(
        lambda p: Gas_Den.gas_density(pressure=p, temp=temp, sg_gas=sg_gas))
    pvt_tab['Gas FVF'] = pvt_tab.apply(
        lambda row: Wrappers.gas_formation_volume_factor_wrapper(pressure=row['Pressure'], temp=temp, sg_gas=sg_gas,
                                                                 gas_comp_factor=row['Gas Compressibility Factor']),
        axis=1)
    pvt_tab['Gas Viscosity'] = pvt_tab.apply(lambda row: mu_g_corr(pressure=row['Pressure'], temp=temp, sg_gas=sg_gas,
                                                                   gas_comp_factor=row['Gas Compressibility Factor']),
                                             axis=1)
    pvt_tab['Oil FVF'] = pvt_tab['Pressure'].apply(
        lambda p: bo_corr(pressure=p, temp=temp, oil_api=oil_api, sg_gas=sg_gas, p_bubble=p_bubble, p_sep=p_sep,
                          t_sep=t_sep))
    pvt_tab['Oil Viscosity'] = pvt_tab.apply(
        lambda row: mu_o_corr(pressure=row['Pressure'], temp=temp, oil_api=oil_api, sg_gas=sg_gas, p_bubble=p_bubble,
                              gas_sol=row['Gas Solubility']), axis=1)
    pvt_tab['Oil Density'] = pvt_tab.apply(
        lambda row: Wrappers.oil_density_wrapper(pressure=row['Pressure'], temp=temp, oil_api=oil_api, sg_gas=sg_gas,
                                                 p_bubble=p_bubble, gas_sol=row['Gas Solubility'],
                                                 oil_fvf=row['Oil FVF']), axis=1)
    return pvt_tab
