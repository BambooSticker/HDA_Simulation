import numpy as np
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util.initialization import propagate_state


def Initialize(b, Arc=None, Order=None):
    '''Initialize units in HDA flowsheet

    Args:
        b: unit as a pyo.Block
        Arc: an pyo.Arc with the b.Inlet as destination 
        Order: initialization order calculated by SD tool
    
    Returns:
        None
    
    '''
    solver = pyo.SolverFactory('ipopt')
    
    if b.name == Order:
        # propagate the stream propeties from the former unit
        propagate_state(arc=Arc)

    if 'R' in b.name:
        for c in b.Comp:
            calculate_variable_from_constraint(b.outlet.flow['Liq',c],b.rule_mass_balance_liq[c])
            calculate_variable_from_constraint(b.outlet.flow['Vap',c],b.rule_mass_balance_vap[c])

    
    if 'M' not in b.name:
        # Flash calculation: T_bub -> T_dew -> T_1 -> T_eq
        calculate_variable_from_constraint(b.temperature_dew, b.eq_T_dew)
        calculate_variable_from_constraint(b._t1, b._t1_constraint)
        calculate_variable_from_constraint(b._teq, b._teq_constraint)

        b.Inlet.fix()
        flag = True
        while flag:
            solver.options['print_user_options'] = 'no'
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
        b.Inlet.unfix()
        
    else:
        b.Recycle_Feed.fix()
        flag = True
        while flag:
            solver.options['print_user_options'] = 'no'
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
            
        b.Recycle_Feed.unfix()

    print("Initialization complete for '%s'! " %(b.name))


def print_outcome(m):
    '''Print the optimal condition of unfixed units

    Args:
        m: flowsheet pyo.model
    
    Returns:
        None
    
    '''
    print("Optimal Values: ")
    print("\t H101 outlet temperature = ", value(m.H101.Outlet.T), "K")
    print("\t R101 outlet temperature = ", value(m.R101.Outlet.T), "K")
    print("\t F101 outlet temperature = ", value(m.F101.Vap_Outlet.T), "K")
    print("\t F102 outlet temperature = ", value(m.F102.Vap_Outlet.T), "K")
    print("\t F102 outlet pressure = ", value(m.F102.Vap_Outlet.P), "Pa")
    print()
    print("Optimal Outcome: ")
    print("\t operating cost  = $", value(m.operating_cost))
    print("\t benzene purity  =", value(m.purity))
    print("\t benzene product =", value(m.product))
    print("\t Overhead loss   =", value(m.loss))


Phase = ['Vap', 'Liq']
Comp = ['benzene','toluene','methane','hydrogen']
mw_comp_data = {'benzene': 78.1136E-3,
                    'toluene': 92.1405E-3,
                    'hydrogen': 2.016e-3,
                    'methane': 16.043e-3}
def stream_report(b):
    '''Get information of a single stream

    Args:
        b: stream as a pyo.Arc
    
    Returns:
        None
    
    '''
    # name of stream
    print("%s: " %(b.name))

    # temperature & pressure
    print("Temperature[˚C]: ", round(value(b.dest.T)-273.15, 2))
    print("Pressure[kPa]: ", round(value(b.dest.P)*1e-3,2))

    # phase type
    if value(sum(b.dest.flow['Vap',c] for c in Comp)) >= 1e-3:
        if value(sum(b.dest.flow['Liq',c] for c in Comp)) >= 1e-3:
            print("Phase type: Both")
        else:
            print("Phase type: Vap")
    else:
        print("Phase type: Liq")
        
    # mole flow rate
    print("flow rate[kmol/hr]: ", )
    for p in Phase:
        for c in Comp:
            print('\t [',p,',',c,']: ',round(value(b.dest.flow[p,c])*3.6,3))
            
    # mass flowrate
    print("Mass flowrate[kg/hr]:")
    mass_flow = {}
    for c in Comp:
        val = value(sum(b.dest.flow[p,c] for p in Phase))* \
                         mw_comp_data[c]*3600
        mass_flow[c] = val
        print("\t [ %s ]:"%(c), round(val,3) )

    # mass fraction
    total_mass_flow = 0
    for mass in mass_flow.values(): 
        total_mass_flow += mass
    print("Mass fraction[-]:")
    for c in Comp:
        print("\t [ %s ]:"%(c), round(mass_flow[c]/total_mass_flow,6)) 
    print("****************************************")

def stream_information(m):
    '''Get information of every stream in HDA flowsheet

    Args:
        m: flowsheet pyo.model
    
    Returns:
        pd.DataFrame presenting the stream information
    
    '''
    data = {f's{i}': [] for i in range(1, 12)}
    stream = [f's{i}' for i in range(1,12)]
    Item = ['温度/˚C', '压力/kPa', '相态','摩尔流量kmol/hr','质量流率kg/hr', \
            '->苯','->甲苯','->甲烷','->氢气','质量分率','->苯','->甲苯','->甲烷','->氢气']

    def PFD_creater(b):
        # name of stream
        s = b.name
    
        # temperature & pressure
        temp = round(value(b.dest.T)-273.15, 2)
        data[s].append(temp)
    
        pres = round(value(b.dest.P)*1e-3,2)
        data[s].append(pres)
    
        # phase type
        if value(sum(b.dest.flow['Vap',c] for c in Comp)) >= 1e-3:
            if value(sum(b.dest.flow['Liq',c] for c in Comp)) >= 1e-3:
                data[s].append('气液共存')
                
            else:
                data[s].append('气相')
        else:
            data[s].append('液相')
    
        # mole flow rate
        mole_flowrate = round(value(sum(b.dest.flow[p,c] for p in Phase for c in Comp))*3.6,3)
        data[s].append(mole_flowrate)
    
        data[s].append('-')
        
        # mass flowrate
        mass_flow = {}
        for c in Comp:
            val = round(value(sum(b.dest.flow[p,c] for p in Phase))* \
                             mw_comp_data[c]*3600,2)
            mass_flow[c] = val
            data[s].append(val)
            
        data[s].append('-')
        # mass fraction
        total_mass_flow = 0
        for mass in mass_flow.values(): 
            total_mass_flow += mass
        for c in Comp:
            mass_frac = round(mass_flow[c]/total_mass_flow,6)
            data[s].append(mass_frac)
    
    for block in m.component_data_objects(Arc, descend_into=False):
        PFD_creater(block)
        
    updated_data = {key: value + [0]*14 if not value else value for key, value in data.items()}
    df = pd.DataFrame(updated_data,index=Item)

    df_s1 = [303.15-273.15,350.0,'气相',f'{0.32*3.60:.3f}']
    total_mf = 0
    mass_flow = {}
    
    df_s1.append('-')
    for c in Comp:
        mf = value(sum(m.M101.Hydrogen_Feed.flow[p,c] for p in Phase))* \
                                 mw_comp_data[c]*3600
        mass_flow[c] = round(mf,2)
        df_s1.append(mass_flow[c])
        
    df_s1.append('-')
    for mass in mass_flow.values(): 
        total_mf += mass
    for c in Comp:
        mass_frac = round(mass_flow[c]/total_mf,6)
        df_s1.append(mass_frac)

    df['s1']=df_s1
    

    df_s2 = [303.15-273.15,350.0,'气相',f'{0.30*3.6:.3f}']
    total_mf = 0
    mass_flow = {}
    
    df_s2.append('-')
    for c in Comp:
        mf = value(sum(m.M101.Toluene_Feed.flow[p,c] for p in Phase))* \
                                 mw_comp_data[c]*3600
        mass_flow[c] = round(mf,2)
        df_s2.append(mass_flow[c])
    
    df_s2.append('-')
    for mass in mass_flow.values(): 
        total_mf += mass
    for c in Comp:
        mass_frac = round(mass_flow[c]/total_mf,6)
        df_s2.append(mass_frac)
            
    df['s2']=df_s2

    df_s8 = [f'{value(m.S101.Purge.T-273.15):.2f}',
             value(m.S101.Purge.P*1e-3),
             '气相',
             f'{value(sum(m.S101.Purge.flow[p,c] for p in Phase for c in Comp)*3.6):.3f}'
            ]
    total_mf = 0
    mass_flow = {}
    
    df_s8.append('-')
    for c in Comp:
        mf = value(sum(m.S101.Purge.flow[p,c] for p in Phase))* \
                                 mw_comp_data[c]*3600
        mass_flow[c] = round(mf,2)
        df_s8.append(mass_flow[c])
    
    df_s8.append('-')
    for mass in mass_flow.values(): 
        total_mf += mass
    for c in Comp:
        mass_frac = round(mass_flow[c]/total_mf,6)
        df_s8.append(mass_frac)
            
    df['s8']=df_s8

    df_s10 = [f'{value(m.F102.Vap_Outlet.T-273.15):.2f}',
             value(m.F102.Vap_Outlet.P*1e-3),
             '气相',
             f'{value(sum(m.F102.Vap_Outlet.flow[p,c] for p in Phase for c in Comp)*3.6):.3f}'
            ]
    total_mf = 0
    mass_flow = {}
    
    df_s10.append('-')
    for c in Comp:
        mf = value(sum(m.F102.Vap_Outlet.flow[p,c] for p in Phase))* \
                                 mw_comp_data[c]*3600
        mass_flow[c] = round(mf,2)
        df_s10.append(mass_flow[c])
    
    df_s10.append('-')
    for mass in mass_flow.values(): 
        total_mf += mass
    for c in Comp:
        mass_frac = round(mass_flow[c]/total_mf,6)
        df_s10.append(mass_frac)
            
    df['s10']=df_s10

    df_s11 = [f'{value(m.F102.Liq_Outlet.T-273.15):.2f}',
             value(m.F102.Liq_Outlet.P*1e-3),
             '液相',
             f'{value(sum(m.F102.Liq_Outlet.flow[p,c] for p in Phase for c in Comp)*3.6):.3f}'
            ]
    total_mf = 0
    mass_flow = {}
    
    df_s11.append('-')
    for c in Comp:
        mf = value(sum(m.F102.Liq_Outlet.flow[p,c] for p in Phase))* \
                                 mw_comp_data[c]*3600
        mass_flow[c] = round(mf,2)
        df_s11.append(mass_flow[c])
    
    df_s11.append('-')
    for mass in mass_flow.values(): 
        total_mf += mass
    for c in Comp:
        mass_frac = round(mass_flow[c]/total_mf,6)
        df_s11.append(mass_frac)
            
    df['s11']=df_s11

    return df