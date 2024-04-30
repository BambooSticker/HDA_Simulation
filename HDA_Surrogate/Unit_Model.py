import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.network import *
from State_Block_Ideal import Create_State_Block
from State_Block_for_Mixer import Create_State_Block_for_Mixer
from Surrogate_Flash import load_surrogate
from datetime import datetime

Phase = ['Liq','Vap']
Comp = ['benzene','toluene','methane','hydrogen']
Inlet = ['Toluene','Hydrogen']
Outlet = ['Purge','Benzene','Toluene']

def Inlet_Feed_creater(b,i):
    
    b.T = Var(initialize = 298.15, bounds=(25, 1200))
    b.P = Var(initialize = 101325, bounds=(100, 1e+10))
    b.flow = Var(Phase, Comp, initialize=0.12, bounds=(1e-12, 5))
    feed_var_list = [(b.flow, Port.Equality),
                     (b.T, Port.Equality),
                     (b.P, Port.Equality)]
    b.Port = Port(initialize = feed_var_list)

        
def Outlet_Flow_creater(b,o):
    b.T = Var(initialize = 298.15, bounds=(25, 1200))
    b.P = Var(initialize = 101325, bounds=(100, 1e+10))
    b.flow = Var(Phase, Comp, initialize=0.12, bounds=(1e-12, 5))
    outlet_var_list = [(b.flow, Port.Equality),
                       (b.T, Port.Equality),
                       (b.P, Port.Equality)]
    b.Port = Port(initialize = outlet_var_list)

def M101_creater(b):
    Create_State_Block_for_Mixer(b)

    # 物料守恒
    def rule_mass_balance(b,p,c):
        return b.outlet.flow[p,c] == sum(b.inlet_feed[f].flow[p,c] for f in b.inlet_list)
    b.rule_mass_balance = Constraint(Phase, Comp,rule=rule_mass_balance)
    
    # 能量守恒
    def rule_energy_balance(b):
        return sum(b.inlet_feed[f].total_enth for f in b.inlet_list) == b.outlet.total_enth
    b.rule_energy_balance = Constraint(rule = rule_energy_balance)
    
    # 等压
    b.rule_eq_pres = Constraint(expr = b.outlet.P == 350000)

#--------------------------- Flash 101 -----------------------------------
def F101_creater(b):
    Create_State_Block(b, VLE=True)
    b.del_component(b.Outlet)

    # 分出两股物料
    b.virtual_stream_vap = Var(b.Phase,b.Comp,initialize = 1e-8)
    
    @b.Constraint(b.Phase,b.Comp)
    def rule_virtual_stream_vap(b,p,c):
        if p == 'Vap':
            return b.virtual_stream_vap[p,c] == b.outlet.flow[p,c]
        else:
            return b.virtual_stream_vap[p,c] == 1e-8

    b.virtual_stream_liq = Var(b.Phase,b.Comp,initialize = 1e-8)
    
    @b.Constraint(b.Phase,b.Comp)
    def rule_virtual_stream_liq(b,p,c):
        if p == 'Liq':
            return b.virtual_stream_liq[p,c] == b.outlet.flow[p,c]
        else:
            return b.virtual_stream_liq[p,c] == 1e-8
    
    b.Liq_Outlet = Port(
        initialize=[(b.outlet.T, Port.Equality),
                    (b.outlet.P, Port.Equality)])
    
    b.Liq_Outlet.add(b.virtual_stream_liq,
                          name="flow",
                          rule=Port.Equality)
    
    b.Vap_Outlet = Port(
        initialize=[(b.outlet.T, Port.Equality),
                    (b.outlet.P, Port.Equality)])
    
    b.Vap_Outlet.add(b.virtual_stream_vap,
                          name="flow",
                          rule=Port.Equality)

    # 物料守恒
    @b.Constraint(b.Comp)
    def rule_mass_comp(b, c):
        return b.inlet.comp_flow[c] == b.outlet.comp_flow[c]

    # 能量守恒
    @b.Constraint()
    def rule_energy_balance(b):
        return b.inlet.total_enth + b.Q == b.outlet.total_enth
    
#-------------------------------- Flash 102  -------------------------------------
### Surrogate Model
def F102_creater(b):
    Create_State_Block(b)
    load_surrogate(b)
    b.del_component(b.Outlet)

    # 分出气液两股物料
    b.virtual_stream_vap = Var(Phase,Comp,initialize = 1e-8)
    b.virtual_stream_liq = Var(Phase,Comp,initialize = 1e-8)

    @b.Constraint(Phase,Comp)
    def rule_virtual_stream_vap(b,p,c):
        if p == 'Vap':
            return b.virtual_stream_vap[p,c] == b.outlet.flow[p,c]
        else:
            return b.virtual_stream_vap[p,c] == 1e-8
    
    @b.Constraint(Phase,Comp)
    def rule_virtual_stream_liq(b,p,c):
        if p == 'Liq':
            return b.virtual_stream_liq[p,c] == b.outlet.flow[p,c]
        else:
            return b.virtual_stream_liq[p,c] == 1e-8
    
    b.Liq_Outlet = Port(
        initialize=[(b.outlet.T, Port.Equality),
                    (b.outlet.P, Port.Equality)])
    b.Liq_Outlet.add(b.virtual_stream_liq,
                          name="flow",
                          rule=Port.Equality)
    
    b.Vap_Outlet = Port(
        initialize=[(b.outlet.T, Port.Equality),
                    (b.outlet.P, Port.Equality)])
    b.Vap_Outlet.add(b.virtual_stream_vap,
                          name="flow",
                          rule=Port.Equality)

    # 物料守恒
    @b.Constraint(Comp)
    def rule_mass_comp(b, c):
        return b.inlet.comp_flow[c] == b.outlet.comp_flow[c]

    # 能量守恒
    @b.Constraint()
    def rule_energy_balance(b):
        return b.inlet.total_enth + b.Q == b.outlet.total_enth
    
#-------------------------------- Heater 101 -------------------------------------
def H101_creater(b):
    Create_State_Block(b)
    
    # 物料守恒(出口气相为总, 液相为0)
    @b.Constraint(Phase, Comp)
    def rule_mass_balance_vap(b,p,c):
        if p == 'Vap':
            return b.outlet.flow[p,c] == b.inlet.comp_flow[c]
        else:
            return b.outlet.flow[p,c] == 1e-8

    # 能量守恒
    @b.Constraint()
    def rule_energy_balance(b):
        return b.inlet.total_enth + b.Q == b.outlet.total_enth

#--------------------------------- Reactor 101 ------------------------------------
def R101_creater(b):
    Create_State_Block(b)
    
    # x_toluene = 0.75, need to be determined
    b.extent_of_reaction = Var(initialize = 0.75, bounds=(1e-12,1))
    
    # Reaction Index
    b.rate_reaction_idx = Set(initialize=["R1"])

    # Reaction Stoichiometry
    b.rate_reaction_stoichiometry = {("R1", "Vap", "benzene"): 1,
                                    ("R1", "Vap", "toluene"): -1,
                                    ("R1", "Vap", "hydrogen"): -1,
                                    ("R1", "Vap", "methane"): 1,
                                    ("R1", "Liq", "benzene"): 0,
                                    ("R1", "Liq", "toluene"): 0,
                                    ("R1", "Liq", "hydrogen"): 0,
                                    ("R1", "Liq", "methane"): 0}
    @b.Constraint(Comp)
    def rule_mass_balance_vap(b,c):
        return b.outlet.flow['Vap',c] == \
            b.inlet.flow['Vap',c] + (-b.rate_reaction_stoichiometry["R1",'Vap',c])\
        *b.rate_reaction_stoichiometry["R1",'Vap','toluene'] *b.inlet.flow['Vap','toluene']* b.extent_of_reaction

    @b.Constraint(Comp)
    def rule_mass_balance_liq(b,c):
        return b.outlet.flow['Liq',c] == b.inlet.flow['Liq',c]

    # 动量守恒
    @b.Constraint()
    def momentum_balance_rule(b):
        return b.outlet.P == b.inlet.P

    # 能量守恒
    dh_rxn_dict = {"R1": 1.08e5}
    b.dh_rxn = Param(b.rate_reaction_idx,
                     initialize=dh_rxn_dict,
                     #units=pyunits.J/pyunits.mol,
                     doc="Heat of reaction")

    @b.Constraint()
    def rule_energy_balance(b):
        return b.inlet.total_enth + b.Q + \
        b.extent_of_reaction*b.inlet.flow['Vap','toluene']*b.dh_rxn['R1']== b.outlet.total_enth

#----------------------------- Spliter 101 --------------------------------------
def S101_creater(b):
    Create_State_Block(b)
    b.del_component(b.Outlet)

    # 分出两股物料
    b.purge = Var(Phase,Comp,initialize = 1e-5)
    b.recycle_ratio = Var(initialize = 0.5, bounds=(0.5,1))
    def rule_purge(b,p,c):
        return b.purge[p,c] == (1-b.recycle_ratio) * b.inlet.flow[p,c]
    b.rule_purge = Constraint(Phase,Comp,rule = rule_purge)
    
    b.recycle = Var(Phase,Comp,initialize = 1e-8)
    @b.Constraint(Phase, Comp)
    def rule_recycle(b,p,c):
        return b.recycle[p,c] == b.recycle_ratio * b.inlet.flow[p,c]
    
    b.Purge = Port(initialize=[(b.outlet.T, Port.Equality),
                                    (b.outlet.P, Port.Equality)])
    b.Purge.add(b.purge, name="flow", rule=Port.Equality)
    
    b.Recycle = Port(initialize=[(b.outlet.T, Port.Equality),
                                      (b.outlet.P, Port.Equality)])
    b.Recycle.add(b.recycle, name="flow", rule=Port.Equality)
    
    # 物料守恒(出口气相)
    @b.Constraint(Phase, Comp)
    def rule_mass_balance(b,p,c):
        return b.outlet.flow[p,c] == b.inlet.flow[p,c]

    # 等温等压
    b.rule_eq_temp = Constraint(expr = b.inlet.T == b.outlet.T)
    b.rule_eq_pres = Constraint(expr = b.inlet.P == b.outlet.P)
    

def create_unit(m):
    m.units_list = ['Inlet_Feed','R101','H101','F101','F102','M101','S101','Outlet_Flow']

    def time_log():
        return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    for unit in m.units_list:
        if unit == 'Inlet_Feed':
            setattr(m, unit, Block(m.Inlet,rule=globals()[f'{unit}_creater']))
        elif unit == 'Outlet_Flow':
            setattr(m, unit, Block(m.Outlet,rule=globals()[f'{unit}_creater']))
        else:
            setattr(m, unit, Block(rule=globals()[f'{unit}_creater']))

        print(f'[{time_log()}] {unit} Block is created successfully!')
        
    print()
    print('Units list:', m.units_list)
        



    


