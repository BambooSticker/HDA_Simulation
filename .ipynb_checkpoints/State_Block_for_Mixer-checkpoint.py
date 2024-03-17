import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

def Create_State_Block_for_Mixer(b):

    # Define a set for components
    b.Phase = pyo.Set(initialize=['Liq','Vap'])
    b.Comp = pyo.Set(initialize=['benzene','toluene','methane','hydrogen'])

    def state_rule(b):
        # Define a set for components
        b.Phase = Set(initialize=['Liq','Vap'])
        b.Comp = Set(initialize=['benzene','toluene','methane','hydrogen'])

        b.flow = Var(
            b.Phase,
            b.Comp,
            initialize=0.12,
            bounds=(1e-12, 5)
            #unit=mol/s
        )

        b.T = Var(
            initialize=298.15,
            bounds=(25, 1200)
            #unit=K
        )

        b.P = Var(
            initialize=101325,
            bounds=(100, 1e+10)
            #unit=Pa
        )

        # z:总分率 / x:液相分率 / y:气相分率
        def _total_flow(b,c):
            return sum(b.flow[p,c] for p in b.Phase)
        b.total_flow = Expression(b.Comp, rule = _total_flow)
        
        def _z(b,c):
            return sum(b.flow[p,c] for p in b.Phase) / sum(b.flow[p,c] for p in b.Phase for c in b.Comp)
        b.z = Expression(b.Comp, rule=_z)

        def _x(b,c):
            return b.flow['Liq', c]/sum(b.flow['Liq', c] for c in b.Comp)
        b.x = Expression(b.Comp, rule=_x)

        def _y(b,c):
            return b.flow['Vap', c]/sum(b.flow['Vap', c] for c in b.Comp)
        b.y = Expression(b.Comp, rule=_y)

        def _flow_phase(b,p):
            return sum(b.flow[p,c] for c in b.Comp)
        b.flow_phase = Expression(b.Phase, rule=_flow_phase)

        b.enth_mol_phase_comp = Var(
            b.Phase,
            b.Comp,
            #bounds=(1e-10,1e8),
            initialize=7e5,
            #units=pyunits.J/pyunits.mol,
            doc='Phase-component molar specific enthalpies')

        def rule_enth_mol_liq(b):
            return sum(b.enth_mol_phase_comp['Liq', c] * \
                    b.x[c] for c in b.Comp)
        b.enth_mol_liq = Expression(rule=rule_enth_mol_liq)

        def rule_enth_mol_vap(b):
            return sum(b.enth_mol_phase_comp['Vap', c] * \
                    b.y[c] for c in b.Comp)
        b.enth_mol_vap = Expression(rule=rule_enth_mol_vap)

        def rule_total_enth(b):
            return b.enth_mol_liq * b.flow_phase['Liq'] + b.enth_mol_vap * b.flow_phase['Vap']
        b.total_enth = Expression(rule = rule_total_enth)

    b.inlet_list = Set(initialize=['Toluene','Hydrogen','Recycle'])
    b.inlet_feed = Block(b.inlet_list,rule = state_rule)
    b.outlet = Block(rule = state_rule)


# ----------------------- 构造用于network的Port ------------------------
    b.Toluene_Feed = Port(initialize=[(b.inlet_feed['Toluene'].T, Port.Extensive),
                                         (b.inlet_feed['Toluene'].P, Port.Extensive),
                                          (b.inlet_feed['Toluene'].flow, Port.Extensive)])
    b.Hydrogen_Feed = Port(initialize=[(b.inlet_feed['Hydrogen'].T, Port.Extensive),
                                         (b.inlet_feed['Hydrogen'].P, Port.Extensive),
                                          (b.inlet_feed['Hydrogen'].flow, Port.Extensive)])
    b.Recycle_Feed = Port(initialize=[(b.inlet_feed['Recycle'].T, Port.Extensive),
                                         (b.inlet_feed['Recycle'].P, Port.Extensive),
                                          (b.inlet_feed['Recycle'].flow, Port.Extensive)])
    b.Outlet = Port(initialize=[(b.outlet.T, Port.Extensive),
                                         (b.outlet.P, Port.Extensive),
                                          (b.outlet.flow, Port.Extensive)])

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
    
    # 液相中各组分的摩尔焓
    def _inlet_enth_mol_comp_liq(b, f, c):
        return b.inlet_feed[f].enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.cp_ig_5['Liq', c] / 5) *
            (b.inlet_feed[f].T**5 - b.temperature_ref**5)
            + (b.cp_ig_4['Liq', c] / 4) *
            (b.inlet_feed[f].T**4 - b.temperature_ref**4)
            + (b.cp_ig_3['Liq', c] / 3) *
            (b.inlet_feed[f].T**3 - b.temperature_ref**3)
            + (b.cp_ig_2['Liq', c] / 2) *
            (b.inlet_feed[f].T**2 - b.temperature_ref**2)
            + b.cp_ig_1['Liq', c] *
            (b.inlet_feed[f].T - b.temperature_ref))
    b.rule_inlet_enth_mol_comp_liq = Constraint(b.inlet_list, b.Comp,rule = _inlet_enth_mol_comp_liq)
    

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
    def _inlet_enth_mol_comp_vap(b, f, c):
        return b.inlet_feed[f].enth_mol_phase_comp['Vap', c] == b.dh_vap[c] + \
                ((b.cp_ig_5['Vap', c] / 5) *
                    (b.inlet_feed[f].T**5 - b.temperature_ref**5)
                    + (b.cp_ig_4['Vap', c] / 4) *
                      (b.inlet_feed[f].T**4 - b.temperature_ref**4)
                    + (b.cp_ig_3['Vap', c] / 3) *
                      (b.inlet_feed[f].T**3 - b.temperature_ref**3)
                    + (b.cp_ig_2['Vap', c] / 2) *
                      (b.inlet_feed[f].T**2 - b.temperature_ref**2)
                    + b.cp_ig_1['Vap', c] *
                      (b.inlet_feed[f].T - b.temperature_ref))
    b.rule_inlet_enth_mol_comp_vap = Constraint(b.inlet_list, b.Comp,rule = _inlet_enth_mol_comp_vap)

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
