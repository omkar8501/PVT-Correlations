import Fluid
import Gas_Solubility as Gas_Sol
import Oil_Formation_Volume_Factor as Oil_FVF
import Gas_Compressibility_Factor as Gas_Comp
import Gas_Formation_Volume_Factor as Gas_FVF

if __name__ == '__main__':
    oil = Fluid.Fluid(37, 0.7)

    pressure = 1000
    temp = 617.67
    oil_api = oil.oil_api
    sg_gas = oil.sg_gas

    Rs = Gas_Sol.standings_gas_solubility(pressure, temp, oil_api, sg_gas)
    Bo = Oil_FVF.standings_oil_fvf(Rs, temp, oil_api, sg_gas)
    Z = Gas_Comp.carnahan_starling_hs_eos(pressure, temp, sg_gas)
    Bg = Gas_FVF.gas_formation_volume_factor(pressure, temp, Z)

    print(Bg)
