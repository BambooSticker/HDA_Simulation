import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from Parameter_Ideal import BT_Params

def Create_State_Block_for_Mixer(b):

    b.param = Block(rule = BT_Params)

    # Define a set for components
    b.Phase = pyo.Set(initialize=['Liq','Vap'])
    b.Comp = pyo.Set(initialize=['benzene','toluene','methane','hydrogen'])

    b.Q = Var(initialize = 0)

    def state_rule(b):
        # Define a set for components
        b.Phase = Set(initialize=['Liq','Vap'])
        b.Comp = Set(initialize=['benzene','toluene','methane','hydrogen'])

        b.flow = Var(b.Phase,b.Comp,initialize=0.12,bounds=(1e-12, 5)) # kmol/s
        b.T = Var(initialize=298.15, bounds=(25, 1200)) # K
        b.P = Var(initialize=101325,bounds=(100, 1e+10)) # Pa

        # z:总分率 / x:液相分率 / y:气相分率
        @b.Expression(b.Comp)
        def comp_flow(b,c): # 某组分的总流率
            return sum(b.flow[p,c] for p in b.Phase)
            
        @b.Expression(b.Phase)
        def phase_flow(b,p): # 某相的总流率
            return sum(b.flow[p,c] for c in b.Comp)
            
        b.total_flow = Var(initialize=1, bounds=(1e-12,10))
        @b.Constraint()
        def _total_flow(b): # 总流率
            return b.total_flow == sum(b.flow[p,c] for p in b.Phase for c in b.Comp)

        b.z = Var(b.Comp,initialize=0.5, bounds=(1e-12,1))
        @b.Constraint(b.Comp)
        def _z(b,c):
            return b.z[c] == b.comp_flow[c] / b.total_flow
            
        b.x = Var(b.Comp,initialize=0.5, bounds=(1e-12,1))
        @b.Constraint(b.Comp)
        def _x(b,c):
            return b.x[c] == b.flow['Liq', c ]/ b.phase_flow['Liq']

        b.y = Var(b.Comp,initialize=0.5, bounds=(1e-12,1))
        @b.Constraint(b.Comp)
        def _y(b,c):
            return b.y[c] == b.flow['Vap', c] / b.phase_flow['Vap']

        b.enth_mol_phase_comp = Var(b.Phase,b.Comp,initialize=7e5,
            doc='Phase-component molar specific enthalpies')
            # pyunits.J/pyunits.mol

        @b.Expression()
        def enth_mol_liq(b):
            return sum(b.enth_mol_phase_comp['Liq', c] * \
                    b.x[c] for c in b.Comp)

        @b.Expression()
        def enth_mol_vap(b):
            return sum(b.enth_mol_phase_comp['Vap', c] * \
                    b.y[c] for c in b.Comp)

        @b.Expression()
        def total_enth(b):
            return b.enth_mol_liq * b.phase_flow['Liq'] + b.enth_mol_vap * b.phase_flow['Vap']

    b.inlet_list = Set(initialize=['Toluene','Hydrogen','Recycle'])
    b.inlet_feed = Block(b.inlet_list,rule = state_rule)
    b.outlet = Block(rule = state_rule)


# ----------------------- 构造用于network的Port ------------------------
    b.Toluene_Feed = Port(initialize=[(b.inlet_feed['Toluene'].T, Port.Equality),
                                         (b.inlet_feed['Toluene'].P, Port.Equality),
                                          (b.inlet_feed['Toluene'].flow, Port.Equality)])
    b.Hydrogen_Feed = Port(initialize=[(b.inlet_feed['Hydrogen'].T, Port.Equality),
                                         (b.inlet_feed['Hydrogen'].P, Port.Equality),
                                          (b.inlet_feed['Hydrogen'].flow, Port.Equality)])
    b.Recycle_Feed = Port(initialize=[(b.inlet_feed['Recycle'].T, Port.Equality),
                                         (b.inlet_feed['Recycle'].P, Port.Equality),
                                          (b.inlet_feed['Recycle'].flow, Port.Equality)])
    b.Outlet = Port(initialize=[(b.outlet.T, Port.Equality),
                                         (b.outlet.P, Port.Equality),
                                          (b.outlet.flow, Port.Equality)])

#--------------------------热量衡算------------------------------------
    # 液相中各组分的摩尔焓
    def _inlet_enth_mol_comp_liq(b, f, c):
        return b.inlet_feed[f].enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.param.cp_ig_5['Liq', c] / 5) *
            (b.inlet_feed[f].T**5 - b.param.temperature_ref**5)
            + (b.param.cp_ig_4['Liq', c] / 4) *
            (b.inlet_feed[f].T**4 - b.param.temperature_ref**4)
            + (b.param.cp_ig_3['Liq', c] / 3) *
            (b.inlet_feed[f].T**3 - b.param.temperature_ref**3)
            + (b.param.cp_ig_2['Liq', c] / 2) *
            (b.inlet_feed[f].T**2 - b.param.temperature_ref**2)
            + b.param.cp_ig_1['Liq', c] *
            (b.inlet_feed[f].T - b.param.temperature_ref))
    b.rule_inlet_enth_mol_comp_liq = Constraint(b.inlet_list, b.Comp,rule = _inlet_enth_mol_comp_liq)
    

    def _outlet_enth_mol_comp_liq(b, c):
        return b.outlet.enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.param.cp_ig_5['Liq', c] / 5) *
            (b.outlet.T**5 - b.param.temperature_ref**5)
            + (b.param.cp_ig_4['Liq', c] / 4) *
            (b.outlet.T**4 - b.param.temperature_ref**4)
            + (b.param.cp_ig_3['Liq', c] / 3) *
            (b.outlet.T**3 - b.param.temperature_ref**3)
            + (b.param.cp_ig_2['Liq', c] / 2) *
            (b.outlet.T**2 - b.param.temperature_ref**2)
            + b.param.cp_ig_1['Liq', c] *
            (b.outlet.T - b.param.temperature_ref))
    b.rule_outlet_enth_mol_comp_liq = Constraint(b.Comp,rule = _outlet_enth_mol_comp_liq)


    # 气相中各组分的摩尔焓
    def _inlet_enth_mol_comp_vap(b, f, c):
        return b.inlet_feed[f].enth_mol_phase_comp['Vap', c] == b.param.dh_vap[c] + \
                ((b.param.cp_ig_5['Vap', c] / 5) *
                    (b.inlet_feed[f].T**5 - b.param.temperature_ref**5)
                    + (b.param.cp_ig_4['Vap', c] / 4) *
                      (b.inlet_feed[f].T**4 - b.param.temperature_ref**4)
                    + (b.param.cp_ig_3['Vap', c] / 3) *
                      (b.inlet_feed[f].T**3 - b.param.temperature_ref**3)
                    + (b.param.cp_ig_2['Vap', c] / 2) *
                      (b.inlet_feed[f].T**2 - b.param.temperature_ref**2)
                    + b.param.cp_ig_1['Vap', c] *
                      (b.inlet_feed[f].T - b.param.temperature_ref))
    b.rule_inlet_enth_mol_comp_vap = Constraint(b.inlet_list, b.Comp,rule = _inlet_enth_mol_comp_vap)

    def _outlet_enth_mol_comp_vap(b, j):
        return b.outlet.enth_mol_phase_comp['Vap', j] == b.param.dh_vap[j] + \
                ((b.param.cp_ig_5['Vap', j] / 5) *
                    (b.outlet.T**5 - b.param.temperature_ref**5)
                    + (b.param.cp_ig_4['Vap', j] / 4) *
                      (b.outlet.T**4 - b.param.temperature_ref**4)
                    + (b.param.cp_ig_3['Vap', j] / 3) *
                      (b.outlet.T**3 - b.param.temperature_ref**3)
                    + (b.param.cp_ig_2['Vap', j] / 2) *
                      (b.outlet.T**2 - b.param.temperature_ref**2)
                    + b.param.cp_ig_1['Vap', j] *
                      (b.outlet.T - b.param.temperature_ref))
    b.rule_outlet_enth_mol_comp_vap = Constraint(b.Comp,rule = _outlet_enth_mol_comp_vap)
