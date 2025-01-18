import inspect
from typing import Optional

import Gas_Solubility as Gas_Sol
import Gas_Compressibility_Factor as Gas_Comp
import Gas_Viscosity as Gas_Visc
import Gas_Formation_Volume_Factor as Gas_FVF
import Oil_Formation_Volume_Factor as Oil_FVF
import Oil_Visocsity as Oil_Visc
import Oil_Density as Oil_Den


# Gas Solubility Correlations
def standings_gas_solubility_wrapper(pressure, temp, oil_api, sg_gas, **kwargs):
    return Gas_Sol.standings_gas_solubility(pressure, temp, oil_api, sg_gas)


def vasquez_beggs_gas_solubility_wrapper(pressure, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep, **kwargs):
    return Gas_Sol.vasquez_beggs_gas_solubility(pressure, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep)


# Gas Compressibility Factor Correlations
def carnahan_starling_hs_eos_wrapper(pressure, temp, sg_gas, **kwargs):
    return Gas_Comp.carnahan_starling_hs_eos(pressure, temp, sg_gas)


# Gas Formation Volume Factor Function
def gas_formation_volume_factor_wrapper(pressure, temp, sg_gas, **kwargs):
    # Optional Argument
    gas_comp_factor = kwargs.get("gas_comp_factor", None)

    return Gas_FVF.gas_formation_volume_factor(pressure, temp, sg_gas, gas_comp_factor)


# Gas Viscosity Correlations
def lee_gonzalez_eakin_wrapper(pressure, temp, sg_gas, **kwargs):
    # Optional Argument
    gas_comp_factor = kwargs.get("gas_comp_factor", None)
    return Gas_Visc.lee_gonzalez_eakin(pressure, temp, sg_gas, gas_comp_factor)


# Oil Formation Volume Factor Correlations
def standings_oil_fvf_wrapper(gas_sol, temp, oil_api, sg_gas, **kwargs):
    return Oil_FVF.standings_oil_fvf(gas_sol, temp, oil_api, sg_gas)


def vasquez_beggs_oil_fvf_wrapper(pressure, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep, **kwargs):
    return Oil_FVF.vasquez_beggs_oil_fvf(pressure, temp, oil_api, sg_gas, p_bubble, p_sep, t_sep)


# Oil Viscosity Correlations
def beggs_robinson_wrapper(pressure, temp, oil_api, sg_gas, p_bubble, **kwargs):
    # Optional Arguments
    gas_sol = kwargs.get("gas_sol", None)
    p_sep = kwargs.get("p_sep", None)
    t_sep = kwargs.get("t_sep", None)

    return Oil_Visc.beggs_robinson(pressure, temp, oil_api, sg_gas, p_bubble, gas_sol, p_sep, t_sep)


# Function to calculate Oil Density
def oil_density_wrapper(pressure, temp, oil_api, sg_gas, p_bubble, **kwargs):
    # Optional Arguments
    gas_sol = kwargs.get("gas_sol", None)
    oil_fvf = kwargs.get("oil_fvf", None)
    p_sep = kwargs.get("p_sep", None)
    t_sep = kwargs.get("t_sep", None)

    return Oil_Den.oil_density(pressure, temp, oil_api, sg_gas, p_bubble, gas_sol, oil_fvf, p_sep, t_sep)


def dynamic_wrapper(function, **kwargs):
    signature = inspect.signature(function)

    # Create a dictionary of the required arguments
    valid_args = {key: val for key, val in kwargs.items() if key in signature.parameters}

    return function(**valid_args)
