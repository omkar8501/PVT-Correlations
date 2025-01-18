import Fluid
import PVT_Table

if __name__ == '__main__':
    oil = Fluid.Fluid(35, 0.68)

    pressure_max = 3000
    p_bubble = 1700
    temp = 620
    oil_api = oil.oil_api
    sg_gas = oil.sg_gas
    p_sep = 120
    t_sep = 520

    table = PVT_Table.pvt_table(pressure_max, p_bubble, temp, oil_api, sg_gas, p_sep, t_sep, 301)
    print(table)


