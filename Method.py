import numpy as np
import pandas as pd
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

def Initialize(b):
    #calculate_variable_from_constraint(b.temperature_bubble, b.eq_T_bubble)
    #b.temperature_bubble.fix()
    #b.del_component(b.eq_T_bubble)
    #b.temperature_bubble.fix(50)

    calculate_variable_from_constraint(b.temperature_dew, b.eq_T_dew)
    b.temperature_dew.fix()
    b.del_component(b.eq_T_dew)

    calculate_variable_from_constraint(b._t1, b._t1_constraint)
    b._t1.fix()
    b.del_component(b._t1_constraint)

    calculate_variable_from_constraint(b._teq, b._teq_constraint)
    b._teq.fix()
    b.del_component(b._teq_constraint)

    if 'R' in b.name:
        for c in b.Comp:
            calculate_variable_from_constraint(b.outlet.flow['Liq',c],b.rule_mass_balance_liq[c])
            b.outlet.flow['Liq',c].fix()
            calculate_variable_from_constraint(b.outlet.flow['Vap',c],b.rule_mass_balance_vap[c])
            b.outlet.flow['Vap',c].fix()
            
        b.del_component(b.rule_mass_balance_liq)
        b.del_component(b.rule_mass_balance_vap)

#def Arc(b1, b2):
#    m.b1.outlet()