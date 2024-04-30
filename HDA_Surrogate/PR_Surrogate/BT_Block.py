import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

def State_Block(b):
    
    # Define a set for components
    b.Phase = Set(initialize = ['Liq','Vap'])
    b.Comp = Set(initialize = ['benzene','toluene'])
    
    def state_rule(b):
        # Define a set for components
        b.Phase = Set(initialize = ['Liq','Vap'])
        b.Comp = Set(initialize = ['benzene','toluene'])
        
        b.flow = Var(
            b.Phase,
            b.Comp,
            initialize = 0.1,
            bounds = (1e-12, 100)
            #unit = mol/s
        )

        b.T = Var( # Temperature
            initialize = 298.15,
            bounds = (25, 2000)
            #unit = K
        )

        b.P = Var( # Pressure
            initialize = 101325,
            bounds = (100, 1e+8)
            #unit = Pa
        )
        
        # z:总分率 / x:液相分率 / y:气相分率
        b.total_flow = Var(b.Comp, initialize = 1, within = NonNegativeReals)
        @b.Constraint(b.Comp)
        def _total_flow(b,c):
            return b.total_flow[c] == sum(b.flow[p,c] for p in b.Phase)

        b.z = Var(b.Comp, initialize = 0.5, bounds = (0,1))
        @b.Constraint(b.Comp)
        def _z(b,c): # Fraction_B
            return b.z[c] == b.total_flow[c] / sum(b.total_flow[c] for c in b.Comp)

        b.x = Var(b.Comp, initialize = 0.5, bounds = (0,1))
        @b.Constraint(b.Comp)
        def _x(b,c):
            return b.x[c] == b.flow['Liq', c]/sum(b.flow['Liq', c] for c in b.Comp)
        
        b.y = Var(b.Comp, initialize = 0.5, bounds = (0,1))
        @b.Constraint(b.Comp)
        def _y(b,c):
            return b.y[c] == b.flow['Vap', c]/sum(b.flow['Vap', c] for c in b.Comp)

        def _flow_phase(b,p):
            return sum(b.flow[p,c] for c in b.Comp)
        b.flow_phase = Expression(b.Phase, rule=_flow_phase)

        b.enth_mol_phase_comp = Var(
            b.Phase,
            b.Comp,
            initialize=7e5,
            #units=pyunits.J/pyunits.mol,
            doc='Phase-component molar specific enthalpies')
    
        b.enth_mol_liq = Var(
            initialize=7e5,
            #units=pyunits.J/pyunits.mol,
            doc='Liquid molar specific enthalpies')

        b.enth_mol_vap = Var(
            initialize=7e5,
            #units=pyunits.J/pyunits.mol,
            doc='Vapor molar specific enthalpies')

        b.total_enth = Var(
            initialize=7e5,
            #units=pyunits.J/pyunits.mol,
            doc='total molar specific enthalpies')

        def rule_enth_mol_liq(b):
            return b.enth_mol_liq == sum(\
                    b.enth_mol_phase_comp['Liq', c] * \
                    b.x[c] for c in b.Comp)
        b.eq_enth_mol_liq = Constraint(rule=rule_enth_mol_liq)

        def rule_enth_mol_vap(b):
            return b.enth_mol_vap == sum(\
                    b.enth_mol_phase_comp['Vap', c] * \
                    b.y[c] for c in b.Comp)
        b.eq_enth_mol_vap = Constraint(rule=rule_enth_mol_vap)

        def rule_total_enth(b):
            return b.total_enth == b.enth_mol_liq * b.flow_phase['Liq'] + b.enth_mol_vap * b.flow_phase['Vap']
        b.eq_enth_mol = Constraint(rule = rule_total_enth)
    
    b.inlet = Block(rule = state_rule)
    b.outlet = Block(rule = state_rule)

#--------------------------热量衡算------------------------------------
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
                ('Vap', 'toluene', '5'): 0}
    
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
                'toluene': 3.8262e4}
    b.dh_vap = Param(b.Comp,
                    mutable=False,
                    #units=pyunits.J/pyunits.mol,
                    initialize=extract_data(dh_vap),
                    doc="heat of vaporization")
    
    # 液相中各组分的摩尔焓
    def _inlet_enth_mol_comp_liq(b, c):
        return b.inlet.enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.cp_ig_5['Liq', c] / 5) *
            (b.inlet.T**5 - b.temperature_ref**5)
            + (b.cp_ig_4['Liq', c] / 4) *
            (b.inlet.T**4 - b.temperature_ref**4)
            + (b.cp_ig_3['Liq', c] / 3) *
            (b.inlet.T**3 - b.temperature_ref**3)
            + (b.cp_ig_2['Liq', c] / 2) *
            (b.inlet.T**2 - b.temperature_ref**2)
            + b.cp_ig_1['Liq', c] *
            (b.inlet.T - b.temperature_ref))
    b.rule_inlet_enth_mol_comp_liq = Constraint(\
        b.Comp,rule = _inlet_enth_mol_comp_liq)
    

    def _outlet_enth_mol_comp_liq(b, c):
        return b.outlet.enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.cp_ig_5['Liq', c] / 5) *
            (b.outlet.T**5 - b.temperature_ref**5)
            + (b.cp_ig_4['Liq', c] / 4) *
            (b.outlet.T**4 - b.temperature_ref**4)
            + (b.cp_ig_3['Liq', c] / 3) *
            (b.outlet.T**3 - b.temperature_ref**3)
            + (b.cp_ig_2['Liq', c] / 2) *
            (b.outlet.T**2 - b.temperature_ref**2)
            + b.cp_ig_1['Liq', c] *
            (b.outlet.T - b.temperature_ref))
    b.rule_outlet_enth_mol_comp_liq = Constraint(b.Comp,rule = _outlet_enth_mol_comp_liq)


    # 气相中各组分的摩尔焓
    def _inlet_enth_mol_comp_vap(b, j):
        return b.inlet.enth_mol_phase_comp['Vap', j] == b.dh_vap[j] + \
                ((b.cp_ig_5['Vap', j] / 5) *
                    (b.inlet.T**5 - b.temperature_ref**5)
                    + (b.cp_ig_4['Vap', j] / 4) *
                      (b.inlet.T**4 - b.temperature_ref**4)
                    + (b.cp_ig_3['Vap', j] / 3) *
                      (b.inlet.T**3 - b.temperature_ref**3)
                    + (b.cp_ig_2['Vap', j] / 2) *
                      (b.inlet.T**2 - b.temperature_ref**2)
                    + b.cp_ig_1['Vap', j] *
                      (b.inlet.T - b.temperature_ref))
    b.rule_inlet_enth_mol_comp_vap = Constraint(b.Comp,rule = _inlet_enth_mol_comp_vap)

    def _outlet_enth_mol_comp_vap(b, j):
        return b.outlet.enth_mol_phase_comp['Vap', j] == b.dh_vap[j] + \
                ((b.cp_ig_5['Vap', j] / 5) *
                    (b.outlet.T**5 - b.temperature_ref**5)
                    + (b.cp_ig_4['Vap', j] / 4) *
                      (b.outlet.T**4 - b.temperature_ref**4)
                    + (b.cp_ig_3['Vap', j] / 3) *
                      (b.outlet.T**3 - b.temperature_ref**3)
                    + (b.cp_ig_2['Vap', j] / 2) *
                      (b.outlet.T**2 - b.temperature_ref**2)
                    + b.cp_ig_1['Vap', j] *
                      (b.outlet.T - b.temperature_ref))
    b.rule_outlet_enth_mol_comp_vap = Constraint(b.Comp,rule = _outlet_enth_mol_comp_vap)
