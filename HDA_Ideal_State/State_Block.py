import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from idaes.core.util.misc import extract_data
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

def Create_State_Block(b):

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

    b.inlet = pyo.Block(rule=state_rule)
    b.outlet = pyo.Block(rule=state_rule)

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

    b.pressure_sat_coeff_A = pyo.Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'A'},
        doc='Parameter A to compute saturated pressure'
        #units=dimensionless
    )

    b.pressure_sat_coeff_B = pyo.Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'B'},
        doc='Parameter B to compute saturated pressure'
        #units=K
    )

    b.pressure_sat_coeff_C = pyo.Param(
        b.Comp,
        mutable=False,
        initialize={c: v for (c, j), v in pressure_sat_coeff_data.items() if j == 'C'},
        doc='Parameter C to compute saturated pressure'
        #units=K
    )

    b.pressure_sat = pyo.Var(b.Comp,
                                 initialize=101325,
                                 bounds=(100, 1e+10),
                                 #unit=pa
                                 doc='Vapor pressure')

    # 相平衡约束（理想物系）
    def rule_Eq(b, c):
        return b.outlet.x[c]*b.pressure_sat[c] == b.outlet.y[c]*b.outlet.P
    b.eq_phase_equilibrium = Constraint(b.Comp, rule=rule_Eq)
    #-----------------------------------------------------------------------------    
    # 泡点温度
    #b.temperature_bubble = Var(initialize=100,
    #                           bounds = (15, 500),
    #                                #units = K,
    #                                doc="Bubble point temperature")
    b.temperature_bubble = Param(initialize=50,
                                    #units = K,
                                    doc="Bubble point temperature")
    #def rule_psat_bubble(b, j):
    #    return 1e5*10**(b.pressure_sat_coeff_A[j] -
    #                   b.pressure_sat_coeff_B[j] /
    #                   (b.temperature_bubble +
    #                    b.pressure_sat_coeff_C[j]))
    #b._p_sat_bubbleT = Expression(b.Comp,
    #                          rule=rule_psat_bubble)

#    def rule_bubble_T(b):
#        return sum(b.inlet.z[c]/(b.inlet.z['toluene']+b.inlet.z['benzene'])*b._p_sat_bubbleT[c] \
#                  for c in ['benzene','toluene']) - b.outlet.P == 0

#    def rule_bubble_T(b):
#        return sum(b.inlet.z[c]*b._p_sat_bubbleT[c] \
#                  for c in b.Comp) - b.outlet.P == 0
    
    #b.eq_T_bubble = Constraint(rule = rule_bubble_T)
    
    # 露点温度
    b.temperature_dew = Var(initialize=370,
                            bounds = (100, 500),
                            #units = K,
                            doc="Dew point temperature")

    def rule_psat_dew(b, j):
        return 1e5*10**(b.pressure_sat_coeff_A[j] -
                       b.pressure_sat_coeff_B[j] /
                       (b.temperature_dew +
                        b.pressure_sat_coeff_C[j]))
    b._p_sat_dewT = Expression(b.Comp, rule=rule_psat_dew)
    
    def rule_dew_T(b):
        return sum(b.inlet.z[c]/b._p_sat_dewT[c] \
                  for c in b.Comp) * b.outlet.P - 1 == 0
    b.eq_T_dew = Constraint(rule = rule_dew_T)
# -----------------------------------------------------------------------------   
    
    b._teq = Var(initialize=400,
                    #unit = K,
                    doc='Temperature for calculating phase equilibrium')
    b._t1 = Var(initialize=400,
                   #unit = K,
                   doc='Intermediate temperature for calculating Teq')
    
    b.eps_1 = Param(default=0.01,
                       #unit = K,
                       mutable=True,
                       doc='Smoothing parameter for T1')
    b.eps_2 = Param(default=0.0005,
                       #unit = K,
                       mutable=True,
                       doc='Smoothing parameter for Teq')
    
    # PSE paper Eqn 13
    def rule_t1(b): 
        return b._t1 == 0.5*(
                b.outlet.T + b.temperature_bubble +
                sqrt((b.outlet.T-b.temperature_bubble)**2 +
                     b.eps_1**2))
    b._t1_constraint = Constraint(rule=rule_t1)  
    
    # PSE paper Eqn 14
    # TODO : Add option for supercritical extension
    def rule_teq(b):
        return b._teq == 0.5*(
              b._t1 + b.temperature_dew -
              sqrt((b._t1-b.temperature_dew)**2 +
                   b.eps_2**2))
    b._teq_constraint = Constraint(rule=rule_teq)
    
    # 安托因方程
    def rule_P_s(b, j):
        return ((log10(b.pressure_sat[j]*1e-5) -
                 b.pressure_sat_coeff_A[j]) *
                (b._teq + b.pressure_sat_coeff_C[j])) == \
               -b.pressure_sat_coeff_B[j]
    b.eq_pressure_sat = Constraint(b.Comp,
                                   rule=rule_P_s)

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
