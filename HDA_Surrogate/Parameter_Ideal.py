import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from idaes.core.util.misc import extract_data

def BT_Params(b):

    b.Phase = Set(initialize=['Liq','Vap'])
    b.Comp = Set(initialize=['benzene','toluene','methane','hydrogen'])
    
# 安托因常数
    pressure_sat_coeff_data = {('benzene', 'A'): 4.202,
                               ('benzene', 'B'): 1322,
                               ('benzene', 'C'): -38.56,
                               ('toluene', 'A'): 4.216,
                               ('toluene', 'B'): 1435,
                               ('toluene', 'C'): -43.33,
                               ('hydrogen', 'A'): 3.543,
                               ('hydrogen', 'B'): 99.40,
                               ('hydrogen', 'C'): 7.726,
                               ('methane', 'A'): 3.990,
                               ('methane', 'B'): 443.0,
                               ('methane', 'C'): -0.49}

    b.pressure_sat_coeff_A = Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'A'},
        doc='Parameter A to compute saturated pressure'
        #units=dimensionless
    )

    b.pressure_sat_coeff_B = Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'B'},
        doc='Parameter B to compute saturated pressure'
        #units=K
    )

    b.pressure_sat_coeff_C = Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'C'},
        doc='Parameter C to compute saturated pressure'
        #units=K
    )

# 计算焓的常数
    # Sources: The Properties of Gases and Liquids (1987)
    #         4th edition, Chemical Engineering Series - Robert C. Reid
    #         Perry's Chemical Engineers Handbook
    #         - Robert H. Perry (Cp_liq)
    cp_ig_data = {('Liq', 'benzene', '1'): 1.29E5,
                ('Liq', 'benzene', '2'): -1.7E2,
                ('Liq', 'benzene', '3'): 6.48E-1,
                ('Liq', 'benzene', '4'): 0,
                ('Liq', 'benzene', '5'): 0,
                ('Vap', 'benzene', '1'): -3.392E1,
                ('Vap', 'benzene', '2'): 4.739E-1,
                ('Vap', 'benzene', '3'): -3.017E-4,
                ('Vap', 'benzene', '4'): 7.130E-8,
                ('Vap', 'benzene', '5'): 0,
                ('Liq', 'toluene', '1'): 1.40E5,
                ('Liq', 'toluene', '2'): -1.52E2,
                ('Liq', 'toluene', '3'): 6.95E-1,
                ('Liq', 'toluene', '4'): 0,
                ('Liq', 'toluene', '5'): 0,
                ('Vap', 'toluene', '1'): -2.435E1,
                ('Vap', 'toluene', '2'): 5.125E-1,
                ('Vap', 'toluene', '3'): -2.765E-4,
                ('Vap', 'toluene', '4'): 4.911E-8,
                ('Vap', 'toluene', '5'): 0,
                ('Liq', 'hydrogen', '1'): 1e-8,  # 6.6653e1,
                ('Liq', 'hydrogen', '2'): 1e-8,  # 6.7659e3,
                ('Liq', 'hydrogen', '3'): 1e-8,  # -1.2363e2,
                ('Liq', 'hydrogen', '4'): 1e-8,  # 4.7827e2, # Eqn 2
                ('Liq', 'hydrogen', '5'): 0,
                ('Vap', 'hydrogen', '1'): 2.714e1,
                ('Vap', 'hydrogen', '2'): 9.274e-3,
                ('Vap', 'hydrogen', '3'): -1.381e-5,
                ('Vap', 'hydrogen', '4'): 7.645e-9,
                ('Vap', 'hydrogen', '5'): 0,
                ('Liq', 'methane', '1'): 1e-8,  # 6.5708e1,
                ('Liq', 'methane', '2'): 1e-8,  # 3.8883e4,
                ('Liq', 'methane', '3'): 1e-8,  # -2.5795e2,
                ('Liq', 'methane', '4'): 1e-8,  # 6.1407e2, # Eqn 2
                ('Liq', 'methane', '5'): 0,
                ('Vap', 'methane', '1'): 1.925e1,
                ('Vap', 'methane', '2'): 5.213e-2,
                ('Vap', 'methane', '3'): 1.197e-5,
                ('Vap', 'methane', '4'): -1.132e-8,
                ('Vap', 'methane', '5'): 0}
    
    b.cp_ig_1 = Param(b.Phase,
                      b.Comp,
                      mutable=False,
                      initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == '1'},
                      doc="Parameter 1 to compute Cp_comp"
                      #units=pyunits.J/pyunits.mol/pyunits.K
                     )
    
    b.cp_ig_2 = Param(b.Phase,
                      b.Comp,
                      mutable=False,
                      initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == '2'},
                      doc="Parameter 2 to compute Cp_comp"
                      #units=pyunits.J/pyunits.mol/pyunits.K**2
                     )
    
    b.cp_ig_3 = Param(b.Phase,
                      b.Comp,
                      mutable=False,
                      initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == '3'},
                      doc="Parameter 3 to compute Cp_comp"
                      #units=pyunits.J/pyunits.mol/pyunits.K**3
                     )
    
    b.cp_ig_4 = Param(b.Phase,
                      b.Comp,
                      mutable=False,
                      initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == '4'},
                      doc="Parameter 4 to compute Cp_comp"
                      #units=pyunits.J/pyunits.mol/pyunits.K**4
                     )
    
    b.cp_ig_5 = Param(b.Phase,
                      b.Comp,
                      mutable=False,
                      initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == '5'},
                      doc="Parameter 5 to compute Cp_comp"
                      #units=pyunits.J/pyunits.mol/pyunits.K**5
                     )

    b.pressure_ref = Param(mutable=True,
                          default=101325,
                          #units=pyunits.Pa,
                          doc='Reference pressure')
    b.temperature_ref = Param(mutable=True,
                             default=298.15,
                             #units=pyunits.K,
                             doc='Reference temperature')

    dh_vap = {'benzene': 3.387e4,
              'toluene': 3.8262e4,
              'hydrogen': 0,
              'methane': 0}
    b.dh_vap = Param(b.Comp,
                    mutable=False,
                    #units=pyunits.J/pyunits.mol,
                    initialize=extract_data(dh_vap),
                    doc="heat of vaporization")