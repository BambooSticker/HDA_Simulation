import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.network import *
from State_Block import Create_State_Block

#-----------------------------------------------------------------------------
def F101_creater(m):
    m.F101 = Block()
    Create_State_Block(m.F101)

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

#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
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

    
