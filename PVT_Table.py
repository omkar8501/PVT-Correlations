import pandas as pd
from typing import Dict

import Gas_Solubility as Gas_Sol
import Gas_Compressibility_Factor as Gas_Comp
import Gas_Formation_Volume_Factor as Gas_FVF
import Gas_Viscosity as Gas_Visc
import Gas_Density as Gas_Den
import Oil_Formation_Volume_Factor as Oil_FVF
import Oil_Visocsity as Oil_Visc
import Oil_Density as Oil_Den


def pvt_table(pressure_max: float, p_bubble: float, temp: float, oil_api: float, sg_gas: float, p_sep: float,
              t_sep: float) -> pd.DataFrame:
    gas_sol_correlations: Dict = {
        "Standing": Gas_Sol.standings_gas_solubility,
        "Vasquez Beggs": Gas_Sol.vasquez_beggs_gas_solubility
    }

    gas_comp_correlations: Dict = {
        "Carnahan Starling": Gas_Comp.carnahan_starling_hs_eos
    }

    gas_visc_correlations: Dict = {
        "Lee Gonzalez Eakin": Gas_Visc.lee_gonzalez_eakin
    }

    oil_fvf_correlations: Dict = {
        "Standing": Oil_FVF.standings_oil_fvf,
        "Vasquez Beggs": Oil_FVF.vasquez_beggs_oil_fvf
    }

    oil_visc_correlations: Dict = {
        "Beggs Robinson": Oil_Visc.beggs_robinson
    }
