import numpy as np
import pandas as pd
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util.initialization import propagate_state


def Initialize(b, Arc=None):

    solver = pyo.SolverFactory('ipopt')
    
    if 'R' not in b.name:
        propagate_state(arc=Arc)

    if 'R' in b.name:
        for c in b.Comp:
            calculate_variable_from_constraint(b.outlet.flow['Liq',c],b.rule_mass_balance_liq[c])
        #    b.outlet.flow['Liq',c].fix()
            calculate_variable_from_constraint(b.outlet.flow['Vap',c],b.rule_mass_balance_vap[c])
        #    b.outlet.flow['Vap',c].fix()
            
        #b.del_component(b.rule_mass_balance_liq)
        #b.del_component(b.rule_mass_balance_vap)
    
    if 'M' not in b.name:
        calculate_variable_from_constraint(b.temperature_dew, b.eq_T_dew)
        #b.temperature_dew.fix()
        #b.del_component(b.eq_T_dew)

        calculate_variable_from_constraint(b._t1, b._t1_constraint)
        #b._t1.fix()
        #b.del_component(b._t1_constraint)
    
        calculate_variable_from_constraint(b._teq, b._teq_constraint)
        #b._teq.fix()
        #b.del_component(b._teq_constraint)

        b.Inlet.fix()
        flag = True
        while flag:
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
            
        b.Inlet.unfix()
        
    else:
        b.Recycle_Feed.fix()
        flag = True
        while flag:
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
            
        b.Recycle_Feed.unfix()
