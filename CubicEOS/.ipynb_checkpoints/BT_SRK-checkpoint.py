from pyomo.environ import *
import os
import idaes
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from idaes.core.util.math import safe_log
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Parameters
component_list = ["benzene","toluene"]
phase_list = ["Liq", "Vap"]
pressure_crit_data = {"benzene": 48.9e5, "toluene": 41.0e5}
temperature_crit_data = {"benzene": 562.2, "toluene": 591.8}
kappa_data = {
            ("benzene", "benzene"): 0.0000,
            ("benzene", "toluene"): 0.0000,
            ("toluene", "benzene"): 0.0000,
            ("toluene", "toluene"): 0.0000,
        }
omega_data = {"benzene": 0.212, "toluene": 0.263}
R = 8.314462618

def BT_SRK(m):
    # define external function for finding cubic roots
    m.crl = ExternalFunction(
        library=os.path.join(idaes.bin_directory, "cubic_roots.so"),
        function="cubic_root_l_ext",
    )
    m.crh = ExternalFunction(
        library=os.path.join(idaes.bin_directory, "cubic_roots.so"),
        function="cubic_root_h_ext",
    )
    
    m.Comp = Set(initialize=component_list)
    m.Phase = Set(initialize=phase_list)
    
    m.Tc = Param(m.Comp,initialize = temperature_crit_data)
    m.Pc = Param(m.Comp,initialize = pressure_crit_data)
    m.kappa = Param(m.Comp, m.Comp, initialize = kappa_data)
    m.omega = Param(m.Comp, initialize = omega_data)
    m.antoine_coeff_A = Param(m.Comp,mutable=False,
                              initialize={"benzene": 4.202, "toluene": 4.216})
    m.antoine_coeff_B = Param(m.Comp,mutable=False,
                              initialize={"benzene": 1322, "toluene": 1435})
    m.antoine_coeff_C = Param(m.Comp,mutable=False,
                              initialize={"benzene": -38.56, "toluene": -43.33})
    
    # for SRK cEOS
    m.OmegaA = Param(default=0.42748)
    m.OmegaB = Param(default=0.08664)
    m.EoS_u = Param(default=1)
    m.EoS_w = Param(default=0)
    m.EoS_p = sqrt(m.EoS_u**2 - 4 * m.EoS_w)
    def func_fw(m,j):
        return 0.48 + 1.574 * m.omega[j] - 0.176 * m.omega[j] ** 2
    m.fw = Expression(m.Comp, rule = func_fw)
    
    # Define state variables
    m.flow_mol = Var(initialize=1.0,bounds=(1e-8, 120)) # F
    m.mole_frac_comp = Var(m.Comp, bounds=(0, None), initialize=0.5) # z
    m.flow_mol_phase = Var(m.Phase, initialize=0.5, domain=NonNegativeReals) # L, V
    m.mole_frac_phase_comp = Var(m.Phase, m.Comp, initialize=0.5,bounds=(0, None)) # x, y
    
    m.T = Var(initialize = 298.15)
    m.P = Var(initialize = 101325,units=units.Pa)
    
    m.temperature_bubble = Var(initialize=400)
    m._mole_frac_tbub = Var(m.Comp,initialize=1 / len(m.Comp),bounds=(0, None))
    m.temperature_dew = Var(initialize=400)
    m._mole_frac_tdew = Var(m.Comp,initialize=1 / len(m.Comp),bounds=(0, None))
    
    def rule_total_mass_balance(m): # L + V = F
        return m.flow_mol_phase["Liq"] + m.flow_mol_phase["Vap"] == m.flow_mol
    m.total_flow_balance = Constraint(rule=rule_total_mass_balance)
    
    def rule_comp_mass_balance(m, c): # L*x + V*y = F*z
        return (
            m.flow_mol * m.mole_frac_comp[c]
            == m.flow_mol_phase["Liq"] * m.mole_frac_phase_comp["Liq", c]
            + m.flow_mol_phase["Vap"] * m.mole_frac_phase_comp["Vap", c]
        )
    m.component_flow_balances = Constraint(m.Comp, rule=rule_comp_mass_balance)
    
    def rule_mole_frac(m): # sum(x) - sum(y) = 1
        return (
            sum(m.mole_frac_phase_comp["Liq", c] for c in m.Comp)
            - sum(m.mole_frac_phase_comp["Vap", c] for c in m.Comp)
            == 0
        )
    m.sum_mole_frac = Constraint(rule=rule_mole_frac)
    
    # 用逸度估算泡露点时，需要分别计算 泡点/露点温度下 的 气相和液相的逸度
    #-----------------------液相泡点温度---------------------------
    def _bubble_temp_liq(m, j):
        def a(k):
            return (
                m.OmegaA * ((R * m.Tc[k]) ** 2 / m.Pc[k])
                * ((1 + m.fw[k] * (1 - sqrt(m.temperature_bubble / m.Tc[k]))) ** 2)
            )
    
        def b(k):
            return (m.OmegaB * R * m.Tc[k] / m.Pc[k])
            
        am = sum(
            sum(
                m.mole_frac_comp[i]
                * m.mole_frac_comp[j]
                * sqrt(a(i) * a(j))
                * (1 - m.kappa[i, j])
                for j in m.Comp
            )
            for i in m.Comp
        )
        bm = sum(m.mole_frac_comp[i] * b(i) for i in m.Comp)
    
        A = am * m.P / (R * m.temperature_bubble) ** 2
        B = bm * m.P / (R * m.temperature_bubble)
    
        delta = (
            2 * sqrt(a(j)) / am * sum(
                m.mole_frac_comp[i] * sqrt(a(i)) * (1 - m.kappa[j, i])
                for i in m.Comp
            )
        )
    
        Z = m.crl.evaluate(args=((-(1+B-m.EoS_u*B)), 
                                 (A-m.EoS_u*B-(m.EoS_u-m.EoS_w)*B**2), 
                                 (-A*B-m.EoS_w*B**2-m.EoS_w*B**3)))
    
        return exp(
            (
                b(j) / bm * (Z - 1) * (B * m.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * m.EoS_p)
                + A
                * (b(j) / bm - delta)
                * safe_log(
                    (2 * Z + B * (m.EoS_u + m.EoS_p))
                    / (2 * Z + B * (m.EoS_u - m.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * m.EoS_p)
        )
    m.bubble_temp_liq = Expression(m.Comp, rule = _bubble_temp_liq)
    
    #-----------------------液相露点温度---------------------------
    def _dew_temp_liq(m, j):
        def a(k):
            return (
                m.OmegaA * ((R * m.Tc[k]) ** 2 / m.Pc[k])
                * ((1 + m.fw[k] * (1 - sqrt(m.temperature_dew / m.Tc[k]))) ** 2)
            )
        def b(k):
            return (m.OmegaB * R * m.Tc[k] / m.Pc[k])
            
        am = sum(
            sum(
                m._mole_frac_tdew[i]
                * m._mole_frac_tdew[j]
                * sqrt(a(i) * a(j))
                * (1 - m.kappa[i, j])
                for j in m.Comp
            )
            for i in m.Comp
        )
        bm = sum(m._mole_frac_tdew[i] * b(i) for i in m.Comp)
        
        A = am * m.P / (R * m.temperature_dew) ** 2
        B = bm * m.P / (R * m.temperature_dew)
        
        delta = (
            2 * sqrt(a(j)) / am
            * sum(m._mole_frac_tdew[i] * sqrt(a(i)) * (1 - m.kappa[j, i]) for i in m.Comp
            )
        )
        
        Z = m.crl.evaluate(args=((-(1+B-m.EoS_u*B)), 
                                 (A-m.EoS_u*B-(m.EoS_u-m.EoS_w)*B**2), 
                                 (-A*B-m.EoS_w*B**2-m.EoS_w*B**3)))
    
        return exp(
                (
                    b(j) / bm * (Z - 1) * (B * m.EoS_p)
                    - safe_log(Z - B, eps=1e-6) * (B * m.EoS_p)
                    + A
                    * (b(j) / bm - delta)
                    * safe_log(
                        (2 * Z + B * (m.EoS_u + m.EoS_p))
                        / (2 * Z + B * (m.EoS_u - m.EoS_p)),
                        eps=1e-6,
                    )
                )
                / (B * m.EoS_p)
            )
    m.dew_temp_liq = Expression(m.Comp, rule = _dew_temp_liq)
    
    #-----------------------气相泡点温度---------------------------
    def _bubble_temp_vap(m, j):
        def a(k):
            return (
                m.OmegaA * ((R * m.Tc[k]) ** 2 / m.Pc[k])
                * ((1 + m.fw[k] * (1 - sqrt(m.temperature_bubble / m.Tc[k]))) ** 2)
            )
        def b(k):
            return (m.OmegaB * R * m.Tc[k] / m.Pc[k])
            
        am = sum(
            sum(
                m._mole_frac_tbub[i]
                * m._mole_frac_tbub[j]
                * sqrt(a(i) * a(j))
                * (1 - m.kappa[i, j])
                for j in m.Comp
            )
            for i in m.Comp
        )
        bm = sum(m._mole_frac_tbub[i] * b(i) for i in m.Comp)
    
        A = am * m.P / (R * m.temperature_bubble) ** 2
        B = bm * m.P / (R * m.temperature_bubble)
    
        delta = (
            2 * sqrt(a(j)) / am
            * sum(m._mole_frac_tbub[i] * sqrt(a(i)) * (1 - m.kappa[j, i]) for i in m.Comp)
        )
    
        Z = m.crh.evaluate(args=((-(1+B-m.EoS_u*B)), 
                                 (A-m.EoS_u*B-(m.EoS_u-m.EoS_w)*B**2), 
                                 (-A*B-m.EoS_w*B**2-m.EoS_w*B**3)))
    
        return exp(
            (
                b(j) / bm * (Z - 1) * (B * m.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * m.EoS_p)
                + A
                * (b(j) / bm - delta)
                * safe_log(
                    (2 * Z + B * (m.EoS_u + m.EoS_p))
                    / (2 * Z + B * (m.EoS_u - m.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * m.EoS_p)
        )
    m.bubble_temp_vap = Expression(m.Comp, rule = _bubble_temp_vap)
        
    #-----------------------气相露点温度---------------------------
    def _dew_temp_vap(m, j):
        def a(k):
            return (
                m.OmegaA * ((R * m.Tc[k]) ** 2 / m.Pc[k])
                * ((1 + m.fw[k] * (1 - sqrt(m.temperature_dew / m.Tc[k]))) ** 2)
            )
            
        def b(k):
            return (m.OmegaB * R * m.Tc[k] / m.Pc[k])
    
        am = sum(
            sum(
                m.mole_frac_comp[i]
                * m.mole_frac_comp[j]
                * sqrt(a(i) * a(j))
                * (1 - m.kappa[i, j])
                for j in m.Comp
            )
            for i in m.Comp
        )
        bm = sum(m.mole_frac_comp[i] * b(i) for i in m.Comp)
    
        A = am * m.P / (R * m.temperature_dew) ** 2
        B = bm * m.P / (R * m.temperature_dew)
    
        delta = (
            2 * sqrt(a(j)) / am
            * sum(m.mole_frac_comp[i] * sqrt(a(i)) * (1 - m.kappa[j, i]) for i in m.Comp)
        )
    
        Z = m.crh.evaluate(args=((-(1+B-m.EoS_u*B)), 
                                 (A-m.EoS_u*B-(m.EoS_u-m.EoS_w)*B**2), 
                                 (-A*B-m.EoS_w*B**2-m.EoS_w*B**3)))
    
        return exp(
            (
                b(j) / bm * (Z - 1) * (B * m.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * m.EoS_p)
                + A
                * (b(j) / bm - delta)
                * safe_log(
                    (2 * Z + B * (m.EoS_u + m.EoS_p))
                    / (2 * Z + B * (m.EoS_u - m.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * m.EoS_p)
        )
    m.dew_temp_vap = Expression(m.Comp, rule = _dew_temp_vap)
    
    ###################### 分别计算泡露点温度
    m._sum_mole_frac_tbub = Constraint(expr=1e3== \
                                       1e3 * sum(m._mole_frac_tbub[j] for j in m.Comp))
    # ln(xi(T_dew))+ln(f_liq(T_dew)) = ln(xi)+ln(f_vap(T_bub))
    def rule_bubble_temp(m, j):
        return log(m.mole_frac_comp[j]) + log(m.bubble_temp_liq[j]) == \
        log(m._mole_frac_tbub[j]) + log(m.bubble_temp_vap[j])
    m.eq_temperature_bubble = Constraint(m.Comp, rule=rule_bubble_temp)
    
    m._sum_mole_frac_tdew = Constraint(expr=1e3 == \
                                       1e3 * sum(m._mole_frac_tdew[j] for j in m.Comp))
    # ln(xi)+ln(f_liq(T_bub)) = ln(xi(T_bub))+ln(f_vap(T_dew))
    def rule_dew_temp(m, j):
        return log(m._mole_frac_tdew[j]) + log(m.dew_temp_liq[j]) == \
        log(m.mole_frac_comp[j]) + log(m.dew_temp_vap[j])
    m.eq_temperature_dew = Constraint(m.Comp, rule=rule_dew_temp)
    ####################### -------------------------------------------
    
    # 算出泡露点后，可以计算 T_eq 并带入计算了
    m._teq = Var(initialize=m.T.value)
    m._t1 = Var(initialize=m.T.value)
    m.eps_1 = Param(default=0.01,mutable=True)
    m.eps_2 = Param(default=0.0005,mutable=True)
    
    def rule_t1(b):
        return b._t1 == 0.5 * (b.T + b.temperature_bubble
            + sqrt((b.T - b.temperature_bubble) ** 2 + b.eps_1**2))
    m._t1_constraint = Constraint(rule=rule_t1)
    
    def rule_teq(b):
        return b._teq == 0.5 * (b._t1 + b.temperature_dew
            - sqrt((b._t1 - b.temperature_dew) ** 2 + b.eps_2**2))
    m._teq_constraint = Constraint(rule=rule_teq)
    
    def func_a(m, j):
        return (m.OmegaA * ((R * m.Tc[j]) ** 2 / m.Pc[j])
            * ((1 + m.fw[j] * (1 - sqrt(m._teq / m.Tc[j]))) ** 2))
    m.a = Expression(m.Comp,rule=func_a)
    
    def func_b(m, j):
        return (
            m.OmegaB
            * R * m.Tc[j] / m.Pc[j]
        )
    m.b = Expression(m.Comp,rule=func_b)
    
    def rule_am(m, p):
        return sum(
        sum(
            m.mole_frac_phase_comp[p, i]
            * m.mole_frac_phase_comp[p, j]
            * sqrt(m.a[i] * m.a[j])
            * (1 - m.kappa[i, j])
            for j in m.Comp
        )
        for i in m.Comp
        )
    m.am = Expression(m.Phase, rule=rule_am)
    
    def rule_bm(m, p):
        return sum(m.mole_frac_phase_comp[p,i] * m.b[i] for i in m.Comp)
    m.bm = Expression(m.Phase, rule=rule_bm)
    
    def func_A(m,p):
        return m.am[p]*m.P/(R*m._teq)**2
    m.A = Expression(m.Phase, rule=func_A)
    
    def func_B(m,p):
        return m.bm[p]*m.P/(R*m._teq)
    m.B = Expression(m.Phase, rule=func_B)
    
    # 定义三次方程一般形式的系数
    def rule_cubic_coef_b(m,p):
        return (-(1+m.B[p]-m.EoS_u*m.B[p]))
    m.cubic_coef_b = Expression(m.Phase,rule=rule_cubic_coef_b)
    
    def rule_cubic_coef_c(m,p):
        return (m.A[p]-m.EoS_u*m.B[p]-(m.EoS_u-m.EoS_w)*m.B[p]**2)
    m.cubic_coef_c = Expression(m.Phase,rule=rule_cubic_coef_c)
    
    def rule_cubic_coef_d(m,p):
        return (-m.A[p]*m.B[p]-m.EoS_w*m.B[p]**2-m.EoS_w*m.B[p]**3)
    m.cubic_coef_d = Expression(m.Phase,rule=rule_cubic_coef_d)
    
    def evaluate_root(m,p):
        if p == 'Vap':
            return m.crh.evaluate(args=(m.cubic_coef_b[p], m.cubic_coef_c[p], m.cubic_coef_d[p]))
        else:
            return m.crl.evaluate(args=(m.cubic_coef_b[p], m.cubic_coef_c[p], m.cubic_coef_d[p]))
    m.external_expr = Expression(m.Phase, rule=evaluate_root)
    
    m.Z = Var(m.Phase, initialize = 1)
    def rule_cubic_Z(m,p):
        return m.Z[p] == m.external_expr[p]
    m.cubic_z = Constraint(m.Phase, rule = rule_cubic_Z)
    
    # 用压缩因子->逸度系数计算相平衡
    def rule_delta(m, p, i):
        # See pg. 145 in Properties of Gases and Liquids
        return (2 * sqrt(m.a[i]) / m.am[p]
            * sum(
                m.mole_frac_phase_comp[p, j]
                * sqrt(m.a[j])
                * (1 - m.kappa[i, j])
                for j in m.Comp
            )
        )
    m.delta = Expression(m.Phase, m.Comp, rule=rule_delta)
    
    def ln_fug_coeff_cubic(b, p, j):
        # See pg. 145 in Properties of Gases and Liquids
        return (
                b.b[j] / b.bm[p] * (b.Z[p] - 1) * (b.B[p] * b.EoS_p)
                - safe_log(b.Z[p] - b.B[p], eps=1e-6)
                * (m.B[p] * m.EoS_p)
                + b.A[p] * (b.b[j] / b.bm[p] - b.delta[p, j])
                * safe_log(
                    (2 * b.Z[p] + b.B[p] * (b.EoS_u + b.EoS_p))
                    / (2 * b.Z[p] + b.B[p] * (b.EoS_u - b.EoS_p)),
                    eps=1e-6
                )
            ) / (b.B[p] * b.EoS_p) + log(b.mole_frac_phase_comp[p, j])
    m.ln_fug_coeff = Expression(m.Phase, m.Comp, rule = ln_fug_coeff_cubic)
    
    # 相平衡约束
    def rule_VLE(m,c):
        return m.ln_fug_coeff['Vap',c] - m.ln_fug_coeff['Liq',c] ==0
    m.rule_equilibrium = Constraint(m.Comp, rule = rule_VLE)


