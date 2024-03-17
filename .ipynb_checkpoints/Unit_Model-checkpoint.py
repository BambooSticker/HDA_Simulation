import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.network import *
from State_Block import Create_State_Block
from State_Block_for_Mixer import Create_State_Block_for_Mixer

def M101_creater(m):
    m.M101 = Block()
    Create_State_Block_for_Mixer(m.M101)

    # 物料守恒
    def rule_mass_balance(b,p,c):
        return m.M101.outlet.flow[p,c] == sum(m.M101.inlet_feed[f].flow[p,c] for f in m.M101.inlet_list)
    m.M101.rule_mass_balance = Constraint(m.Phase, m.Comp,rule=rule_mass_balance)
    
    # 能量守恒
    def rule_energy_balance(b):
        return sum(m.M101.inlet_feed[f].total_enth for f in m.M101.inlet_list) == m.M101.outlet.total_enth
    m.M101.rule_energy_balance = Constraint(rule = rule_energy_balance)
    
    # 等压
    m.M101.rule_eq_pres = Constraint(expr = m.M101.outlet.P == 350000)

    return m.M101

#--------------------------- Flash 101 -----------------------------------
def F101_creater(m):
    m.F101 = Block()
    Create_State_Block(m.F101)
    m.F101.del_component(m.F101.Outlet)

    # 分出两股物料
    m.F101.virtual_stream_vap = Var(m.Phase,m.Comp,initialize = 1e-8)
    def rule_virtual_stream_vap(b,p,c):
        if p == 'Vap':
            return m.F101.virtual_stream_vap[p,c] == m.F101.outlet.flow[p,c]
        else:
            return m.F101.virtual_stream_vap[p,c] == 1e-8
    m.F101.rule_virtual_stream_vap = Constraint(m.Phase,m.Comp,rule = rule_virtual_stream_vap)
    
    m.F101.virtual_stream_liq = Var(m.Phase,m.Comp,initialize = 1e-8)
    def rule_virtual_stream_liq(b,p,c):
        if p == 'Liq':
            return m.F101.virtual_stream_liq[p,c] == m.F101.outlet.flow[p,c]
        else:
            return m.F101.virtual_stream_liq[p,c] == 1e-8
    m.F101.rule_virtual_stream_liq = Constraint(m.Phase,m.Comp,rule = rule_virtual_stream_liq)
    
    m.F101.Liq_Outlet = Port(initialize=[(m.F101.outlet.T, Port.Equality),
                                         (m.F101.outlet.P, Port.Equality)])
    m.F101.Liq_Outlet.add(m.F101.virtual_stream_liq,name="flow",rule=Port.Equality)
    
    m.F101.Vap_Outlet = Port(initialize=[(m.F101.outlet.T, Port.Equality),
                                         (m.F101.outlet.P, Port.Equality)])
    m.F101.Vap_Outlet.add(m.F101.virtual_stream_vap,name="flow",rule=Port.Equality)

    
    # 物料守恒
    def rule_mass_comp(b, c):
        return m.F101.inlet.total_flow[c] == m.F101.outlet.total_flow[c]
    m.F101.rule_mass_balance = Constraint(m.Comp,rule=rule_mass_comp)

    # 能量守恒
    m.F101.Q = Var(initialize = 0)

    def rule_energy_balance(b):
        return m.F101.inlet.total_enth + m.F101.Q == m.F101.outlet.total_enth
    m.F101.rule_energy_balance = Constraint(rule = rule_energy_balance)
    return m.F101


#-------------------------------- Flash 102 -------------------------------------
def F102_creater(m):
    m.F102 = Block()
    Create_State_Block(m.F102)
    m.F102.del_component(m.F102.Outlet)

    # 分出两股物料
    m.F102.virtual_stream_vap = Var(m.Phase,m.Comp,initialize = 1e-8)
    def rule_virtual_stream_vap(b,p,c):
        if p == 'Vap':
            return m.F102.virtual_stream_vap[p,c] == m.F102.outlet.flow[p,c]
        else:
            return m.F102.virtual_stream_vap[p,c] == 1e-8
    m.F102.rule_virtual_stream_vap = Constraint(m.Phase,m.Comp,rule = rule_virtual_stream_vap)
    
    m.F102.virtual_stream_liq = Var(m.Phase,m.Comp,initialize = 1e-8)
    def rule_virtual_stream_liq(b,p,c):
        if p == 'Liq':
            return m.F102.virtual_stream_liq[p,c] == m.F102.outlet.flow[p,c]
        else:
            return m.F102.virtual_stream_liq[p,c] == 1e-8
    m.F102.rule_virtual_stream_liq = Constraint(m.Phase,m.Comp,rule = rule_virtual_stream_liq)
    
    m.F102.Liq_Outlet = Port(initialize=[(m.F102.outlet.T, Port.Equality),
                                         (m.F102.outlet.P, Port.Equality)])
    m.F102.Liq_Outlet.add(m.F102.virtual_stream_liq,name="flow",rule=Port.Equality)
    
    m.F102.Vap_Outlet = Port(initialize=[(m.F102.outlet.T, Port.Equality),
                                         (m.F102.outlet.P, Port.Equality)])
    m.F102.Vap_Outlet.add(m.F102.virtual_stream_vap,name="flow",rule=Port.Equality)

    # 物料守恒
    def rule_mass_comp(b, c):
        return m.F102.inlet.total_flow[c] == m.F102.outlet.total_flow[c]
    m.F102.rule_mass_balance = Constraint(m.Comp,rule=rule_mass_comp)

    # 能量守恒
    m.F102.Q = Var(initialize = 0)

    def rule_energy_balance(b):
        return m.F102.inlet.total_enth + m.F102.Q == m.F102.outlet.total_enth
    m.F102.rule_energy_balance = Constraint(rule = rule_energy_balance)
    
    return m.F102

#-------------------------------- Heater 101 -------------------------------------
def H101_creater(m):
    m.H101 = Block()
    Create_State_Block(m.H101)
    m.H101.del_component(m.H101.eq_phase_equilibrium)
    
    # 物料守恒(出口气相)
    def rule_mass_balance_vap(b,c):
        return m.H101.outlet.flow['Vap',c] == sum(m.H101.inlet.flow[p,c] for p in m.Phase)
    m.H101.rule_mass_balance_vap = Constraint(m.Comp,rule=rule_mass_balance_vap)

    # 物料守恒（出口液相为0）
    def rule_mass_balance_liq(b,c):
        return m.H101.outlet.flow['Liq',c] == 1e-8
    m.H101.rule_mass_balance_liq = Constraint(m.Comp,rule=rule_mass_balance_liq)

    # 能量守恒
    m.H101.Q = Var(initialize = 0)

    def rule_energy_balance(b):
        return m.H101.inlet.total_enth + m.H101.Q == m.H101.outlet.total_enth
    m.H101.rule_energy_balance = Constraint(rule = rule_energy_balance)

    return m.H101    

#--------------------------------- Reactor 101 ------------------------------------
def R101_creater(m):
    m.R101 = Block()
    Create_State_Block(m.R101)
    m.R101.del_component(m.R101.eq_phase_equilibrium)

    # Reaction Index
    m.R101.rate_reaction_idx = Set(initialize=["R1"])

    # Reaction Stoichiometry
    m.R101.rate_reaction_stoichiometry = {("R1", "Vap", "benzene"): 1,
                                    ("R1", "Vap", "toluene"): -1,
                                    ("R1", "Vap", "hydrogen"): -1,
                                    ("R1", "Vap", "methane"): 1,
                                    ("R1", "Liq", "benzene"): 0,
                                    ("R1", "Liq", "toluene"): 0,
                                    ("R1", "Liq", "hydrogen"): 0,
                                    ("R1", "Liq", "methane"): 0}

    # x_toluene = 0.75
    m.R101.extent_of_reaction = Param(default = 0.75)

    def rule_mass_balance_vap(b,c):
        return b.outlet.flow['Vap',c] == \
            b.inlet.flow['Vap',c] + (-b.rate_reaction_stoichiometry["R1",'Vap',c])\
        *b.rate_reaction_stoichiometry["R1",'Vap','toluene'] *b.inlet.flow['Vap','toluene']* b.extent_of_reaction
    m.R101.rule_mass_balance_vap = Constraint(m.Comp, rule=rule_mass_balance_vap)

    def rule_mass_balance_liq(b,c):
        return b.outlet.flow['Liq',c] == b.inlet.flow['Liq',c]
    m.R101.rule_mass_balance_liq = Constraint(m.Comp, rule=rule_mass_balance_liq)

    # 动量守恒
    def momentum_balance_rule(b):
        return b.outlet.P == b.inlet.P
    m.R101.rule_momentum_balance = Constraint(rule = momentum_balance_rule)

    # 能量守恒
    m.R101.Q = Var(initialize = 0)

    dh_rxn_dict = {"R1": 1.08e5}
    m.R101.dh_rxn = Param(m.R101.rate_reaction_idx,
                     initialize=dh_rxn_dict,
                     #units=pyunits.J/pyunits.mol,
                     doc="Heat of reaction")

    def rule_energy_balance(b):
        return m.R101.inlet.total_enth + m.R101.Q + \
        m.R101.extent_of_reaction*b.inlet.flow['Vap','toluene']*m.R101.dh_rxn['R1']== m.R101.outlet.total_enth
    m.R101.rule_energy_balance = Constraint(rule = rule_energy_balance)

    return m.R101

#----------------------------- Spliter 101 --------------------------------------
def S101_creater(m):
    m.S101 = Block()
    Create_State_Block(m.S101)
    m.S101.del_component(m.S101.eq_phase_equilibrium)
    m.S101.del_component(m.S101.Outlet)

    # 分出两股物料
    m.S101.purge = Var(m.Phase,m.Comp,initialize = 1e-5)
    def rule_purge(b,p,c):
        return m.S101.purge[p,c] == 0.2 * m.S101.inlet.flow[p,c]
    m.S101.rule_purge = Constraint(m.Phase,m.Comp,rule = rule_purge)
    
    m.S101.recycle = Var(m.Phase,m.Comp,initialize = 1e-8)
    def rule_recycle(b,p,c):
        return m.S101.recycle[p,c] == 0.8 * m.S101.inlet.flow[p,c]
    m.S101.rule_recycle = Constraint(m.Phase,m.Comp,rule = rule_recycle)
    
    m.S101.Purge = Port(initialize=[(m.S101.outlet.T, Port.Equality),
                                         (m.S101.outlet.P, Port.Equality)])
    m.S101.Purge.add(m.S101.purge, name="flow", rule=Port.Equality)
    
    m.S101.Recycle = Port(initialize=[(m.S101.outlet.T, Port.Equality),
                                         (m.S101.outlet.P, Port.Equality)])
    m.S101.Recycle.add(m.S101.recycle,  name="flow", rule=Port.Equality)
    
    # 物料守恒(出口气相)
    def rule_mass_balance(b,p,c):
        return m.S101.outlet.flow[p,c] == m.S101.inlet.flow[p,c]
    m.S101.rule_mass_balance = Constraint(m.Phase, m.Comp,rule=rule_mass_balance)

    # 等温等压
    m.S101.rule_eq_temp = Constraint(expr = m.S101.inlet.T == m.S101.outlet.T)
    m.S101.rule_eq_pres = Constraint(expr = m.S101.inlet.P == m.S101.outlet.P)
    
    return m.S101  

#----------------------------- Mixer 101 --------------------------------------


