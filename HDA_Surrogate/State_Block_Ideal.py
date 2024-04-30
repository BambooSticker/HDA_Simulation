import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from Parameter_Ideal import BT_Params

def Create_State_Block(b, VLE = False):

    b.param = Block(rule = BT_Params)

    # Define a set for components
    b.Phase = Set(initialize=['Liq','Vap'])
    b.Comp = Set(initialize=['benzene','toluene','methane','hydrogen'])

    b.Q = Var(initialize = 0, units = units.kW)

    def state_rule(b):
        # Define a set for components
        b.Phase = Set(initialize=['Liq','Vap'])
        b.Comp = Set(initialize=['benzene','toluene','methane','hydrogen'])

        b.flow = Var(b.Phase,b.Comp,initialize=0.12,bounds=(1e-12, 5)) # kmol/s
        b.T = Var(initialize=298.15, bounds=(25, 1200)) # K
        b.P = Var(initialize=101325,bounds=(100, 1e+10), units=units.Pa) # Pa

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

    b.inlet = Block(rule=state_rule)
    b.outlet = Block(rule=state_rule)

# ----------------------- 构造用于network的Port ------------------------
    inlet_var_list = [(b.inlet.flow, Port.Equality),
                      (b.inlet.T, Port.Equality),
                      (b.inlet.P, Port.Equality)]
    outlet_var_list = [(b.outlet.flow, Port.Equality),
                       (b.outlet.T, Port.Equality),
                       (b.outlet.P, Port.Equality)]
    
    b.Inlet = Port(initialize = inlet_var_list)
    b.Outlet = Port(initialize = outlet_var_list)
# ---------------------------------------------------------------------
    if VLE == True: 
        b.pressure_sat = Var(b.Comp,initialize=101325,bounds=(100, 1e+10)) #unit=pa
    
        # 相平衡约束（理想物系）
        @b.Constraint(b.Comp)
        def rule_Equilibrium(b, c):
            return b.outlet.x[c]*b.pressure_sat[c] == b.outlet.y[c]*b.outlet.P
        #-----------------------------------------------------------------------------    
        # 泡点温度
        b.temperature_bubble = Param(default=50) # K,
        
        # 露点温度
        b.temperature_dew = Var(initialize=370, bounds = (100, 500)) # K

        @b.Expression(b.Comp)
        def _p_sat_dewT(b, j):
            return 1e5*10**(b.param.pressure_sat_coeff_A[j] -
                           b.param.pressure_sat_coeff_B[j] /
                           (b.temperature_dew +
                            b.param.pressure_sat_coeff_C[j]))

        @b.Constraint()
        def rule_dew_T(b):
            return sum(b.inlet.z[c]/b._p_sat_dewT[c] \
                      for c in b.Comp) * b.outlet.P - 1 == 0
    # -----------------------------------------------------------------------------   
        b._teq = Var(initialize=400, # K,
                    doc='Temperature for calculating phase equilibrium')
        b._t1 = Var(initialize=400, # K,
                    doc='Intermediate temperature for calculating Teq')
        
        b.eps_1 = Param(default=0.01, # K
                        mutable=True,
                        doc='Smoothing parameter for T1')
        b.eps_2 = Param(default=0.0005, # K,
                        mutable=True,
                        doc='Smoothing parameter for Teq')
        
        @b.Constraint()
        def rule_t1(b): 
            return b._t1 == 0.5*(
                    b.outlet.T + b.temperature_bubble +
                    sqrt((b.outlet.T-b.temperature_bubble)**2 +
                         b.eps_1**2))
        
        # TODO : Add option for supercritical extension
        @b.Constraint()
        def rule_teq(b):
            return b._teq == 0.5*(
                  b._t1 + b.temperature_dew -
                  sqrt((b._t1-b.temperature_dew)**2 +
                       b.eps_2**2))
        
        # 安托因方程
        @b.Constraint(b.Comp)
        def rule_P_s(b, c):
            return ((log10(b.pressure_sat[c]*1e-5) -
                     b.param.pressure_sat_coeff_A[c]) *
                    (b._teq + b.param.pressure_sat_coeff_C[c])) == \
                   -b.param.pressure_sat_coeff_B[c]
        

#--------------------------热量衡算------------------------------------
    # 液相中各组分的摩尔焓
    @b.Constraint(b.Comp)
    def rule_inlet_enth_mol_comp_liq(b, c):
        return b.inlet.enth_mol_phase_comp['Liq', c] * 1E3 == \
            ((b.param.cp_ig_5['Liq', c] / 5) *
            (b.inlet.T**5 - b.param.temperature_ref**5)
            + (b.param.cp_ig_4['Liq', c] / 4) *
            (b.inlet.T**4 - b.param.temperature_ref**4)
            + (b.param.cp_ig_3['Liq', c] / 3) *
            (b.inlet.T**3 - b.param.temperature_ref**3)
            + (b.param.cp_ig_2['Liq', c] / 2) *
            (b.inlet.T**2 - b.param.temperature_ref**2)
            + b.param.cp_ig_1['Liq', c] *
            (b.inlet.T - b.param.temperature_ref))

    @b.Constraint(b.Comp)
    def rule_outlet_enth_mol_comp_liq(b, c):
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


    # 气相中各组分的摩尔焓
    @b.Constraint(b.Comp)
    def rule_inlet_enth_mol_comp_vap(b, c):
        return b.inlet.enth_mol_phase_comp['Vap', c] == b.param.dh_vap[c] + \
                ((b.param.cp_ig_5['Vap', c] / 5) *
                    (b.inlet.T**5 - b.param.temperature_ref**5)
                    + (b.param.cp_ig_4['Vap', c] / 4) *
                      (b.inlet.T**4 - b.param.temperature_ref**4)
                    + (b.param.cp_ig_3['Vap', c] / 3) *
                      (b.inlet.T**3 - b.param.temperature_ref**3)
                    + (b.param.cp_ig_2['Vap', c] / 2) *
                      (b.inlet.T**2 - b.param.temperature_ref**2)
                    + b.param.cp_ig_1['Vap', c] *
                      (b.inlet.T - b.param.temperature_ref))

    @b.Constraint(b.Comp)
    def rule_outlet_enth_mol_comp_vap(b, c):
        return b.outlet.enth_mol_phase_comp['Vap', c] == b.param.dh_vap[c] + \
                ((b.param.cp_ig_5['Vap', c] / 5) *
                    (b.outlet.T**5 - b.param.temperature_ref**5)
                    + (b.param.cp_ig_4['Vap', c] / 4) *
                      (b.outlet.T**4 - b.param.temperature_ref**4)
                    + (b.param.cp_ig_3['Vap', c] / 3) *
                      (b.outlet.T**3 - b.param.temperature_ref**3)
                    + (b.param.cp_ig_2['Vap', c] / 2) *
                      (b.outlet.T**2 - b.param.temperature_ref**2)
                    + b.param.cp_ig_1['Vap', c] *
                      (b.outlet.T - b.param.temperature_ref))

