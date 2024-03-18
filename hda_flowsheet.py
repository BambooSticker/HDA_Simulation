import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
from idaes.core.util.model_statistics import degrees_of_freedom
from Method import *
from Unit_Model import *
mw = {'benzene': 78.1136E-3,
    'toluene': 92.1405E-3,
    'hydrogen': 2.016e-3,
    'methane': 16.043e-3}
def hda_flowsheet(TEE = True):
    # build flowsheet
    print('Building flowsheet...')
    print()
    # Create a concrete model
    m = ConcreteModel()
    
    # Define a set for components
    m.Phase = Set(initialize = ['Liq','Vap'])
    m.Comp = Set(initialize = ['benzene','toluene','methane','hydrogen'])
    
    # construct every unit operation
    m.R101 = R101_creater(m)
    m.H101 = H101_creater(m)
    m.F101 = F101_creater(m)
    m.F102 = F102_creater(m)
    m.M101 = M101_creater(m)
    m.S101 = S101_creater(m)

    # Add the Arc to connect
    m.s4 = Arc(source=m.H101.Outlet, destination=m.R101.Inlet)
    m.s5 = Arc(source=m.R101.Outlet, destination=m.F101.Inlet)
    m.s9 = Arc(source=m.F101.Liq_Outlet, destination=m.F102.Inlet)
    m.s6 = Arc(source=m.F101.Vap_Outlet, destination=m.S101.Inlet)
    m.s7 = Arc(source=m.S101.Recycle, destination=m.M101.Recycle_Feed)
    m.s3 = Arc(source=m.M101.Outlet, destination=m.H101.Inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print('Setting inputs...')
    print()
    # Setting for heater
    m.H101.outlet.P.fix(350000)
    m.H101.outlet.T.fix(600)
    
    # Setting for flash vessel
    m.F101.outlet.P.fix(350000)
    m.F101.outlet.T.fix(325)
    
    m.F102.outlet.P.fix(150000)
    m.F102.outlet.T.fix(375)
    
    # Setting for reactor heat duty
    m.R101.Q.fix(0)
    
    # two inlet feeds set
    m.M101.inlet_feed['Toluene'].flow['Liq','benzene'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Liq','toluene'].fix(0.30)
    m.M101.inlet_feed['Toluene'].flow['Liq','methane'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Liq','hydrogen'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Vap','benzene'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Vap','toluene'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Vap','methane'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].flow['Vap','hydrogen'].fix(1e-8)
    m.M101.inlet_feed['Toluene'].T.fix(303.2)
    m.M101.inlet_feed['Toluene'].P.fix(350000)
    
    m.M101.inlet_feed['Hydrogen'].flow['Liq','benzene'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Liq','toluene'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Liq','methane'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Liq','hydrogen'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Vap','benzene'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Vap','toluene'].fix(1e-8)
    m.M101.inlet_feed['Hydrogen'].flow['Vap','methane'].fix(0.02)
    m.M101.inlet_feed['Hydrogen'].flow['Vap','hydrogen'].fix(0.30)
    m.M101.inlet_feed['Hydrogen'].T.fix(303.2)
    m.M101.inlet_feed['Hydrogen'].P.fix(350000)

    print('Initializing flowsheet...')
    print()

    m.cooling_cost = Expression(
        expr=0.212e-7 * (-m.F101.Q) + 0.212e-7 * (-m.R101.Q)
    )
    
    # set heating cost for H101 & F102
    m.heating_cost = Expression(
        expr=2.2e-7 * m.H101.Q + 1.9e-7 * m.F102.Q
    )
    
    m.operating_cost = Expression(
        expr=(3600 * 8000 * (m.heating_cost + m.cooling_cost))
    )
    
    m.purity = m.F102.outlet.flow['Vap','benzene']*mw['benzene']/sum(m.F102.outlet.flow['Vap',c]*mw[c] for c in m.Comp)
    m.loss = m.S101.Purge.flow['Vap','benzene']/m.R101.Outlet.flow['Vap','benzene']
    m.product = m.F102.Vap_Outlet.flow["Vap", "benzene"]

    # Initialization with Sequential Decomposition method
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic" 
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 10
    
    # which stream to tear / what's the calculation order
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    cal_order = seq.calculation_order(G)
    
    # Initialization function for each unit
    def function(unit, order = cal_order[0][0].name):
        if 'R' in unit.name:
            Initialize(unit,Arc=m.s4, Order=order)
        elif 'F101' in unit.name:
            Initialize(unit,Arc=m.s5, Order=order)
        elif 'F102' in unit.name:
            Initialize(unit,Arc=m.s9, Order=order)
        elif 'S' in unit.name:
            Initialize(unit,Arc=m.s6, Order=order)
        elif 'M' in unit.name:
            Initialize(unit,Arc=m.s7, Order=order)
        elif 'H' in unit.name:
            Initialize(unit,Arc=m.s3, Order=order)
    
    # guess for tear stream (s1: H101 -> R101)
    tear_guesses = {
        'flow': {
            ("Vap", "benzene"): 1e-5,
            ("Vap", "toluene"): 0.3,
            ("Vap", "hydrogen"): 0.30,
            ("Vap", "methane"): 0.02,
            ("Liq", "benzene"): 1e-5,
            ("Liq", "toluene"): 1e-5,
            ("Liq", "hydrogen"): 1e-5,
            ("Liq", "methane"): 1e-5,
        },
        'T': 600,
        'P': 350000,
    }
    
    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.R101.Inlet, tear_guesses)
    
    seq.run(m,function)

    print('Solving flowsheet...')
    print()
    # 模型求解
    solver = SolverFactory('ipopt')
    solver.options['tol'] = 1e-4
    results = solver.solve(m,tee=TEE)
    assert results.solver.termination_condition == TerminationCondition.optimal

    return m





