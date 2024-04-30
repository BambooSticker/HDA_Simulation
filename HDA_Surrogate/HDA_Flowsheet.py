from pyomo.environ import (Constraint,
                           Var,
                           Set,
                           Param,
                           Expression,
                           ConcreteModel,
                           SolverFactory,
                           TransformationFactory,
                           units as pyunits,
                           TerminationCondition)

from pyomo.network import Arc, SequentialDecomposition
from Method import (Initialize,
                   cost_vessel, 
                   cost_fired_heater)
from Unit_Model import create_unit

def HDA_Flowsheet():

    # build flowsheet
    print('Building flowsheet...')
    print()
    
    # Create a concrete model
    m = ConcreteModel()
    # Define a set for components
    m.Phase = Set(initialize = ['Liq','Vap'])
    m.Comp = Set(initialize = ['benzene','toluene','methane','hydrogen'])
    m.Inlet = Set(initialize = ['Toluene','Hydrogen'])
    m.Outlet = Set(initialize = ['Purge','Benzene','Toluene'])

    create_unit(m)
    # Add the Arc to connect
    m.s1 = Arc(source=m.Inlet_Feed['Hydrogen'].Port, destination=m.M101.Hydrogen_Feed)
    m.s2 = Arc(source=m.Inlet_Feed['Toluene'].Port, destination=m.M101.Toluene_Feed)
    m.s3 = Arc(source=m.M101.Outlet, destination=m.H101.Inlet)
    m.s4 = Arc(source=m.H101.Outlet, destination=m.R101.Inlet)
    m.s5 = Arc(source=m.R101.Outlet, destination=m.F101.Inlet)
    m.s6 = Arc(source=m.F101.Vap_Outlet, destination=m.S101.Inlet)
    m.s7 = Arc(source=m.S101.Recycle, destination=m.M101.Recycle_Feed)
    m.s8 = Arc(source=m.S101.Purge, destination=m.Outlet_Flow['Purge'].Port)
    m.s9 = Arc(source=m.F101.Liq_Outlet, destination=m.F102.Inlet)
    m.s10 = Arc(source=m.F102.Vap_Outlet, destination=m.Outlet_Flow['Benzene'].Port)
    m.s11 = Arc(source=m.F102.Liq_Outlet, destination=m.Outlet_Flow['Toluene'].Port)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # set inputs
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
    m.R101.extent_of_reaction.fix(0.75)

    m.S101.recycle_ratio.fix(0.8)

    # two inlet feeds set # kmol*s
    m.Inlet_Feed['Toluene'].flow['Liq','benzene'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Liq','toluene'].fix(0.03665)
    m.Inlet_Feed['Toluene'].flow['Liq','methane'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Liq','hydrogen'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Vap','benzene'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Vap','toluene'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Vap','methane'].fix(1e-8)
    m.Inlet_Feed['Toluene'].flow['Vap','hydrogen'].fix(1e-8)
    m.Inlet_Feed['Toluene'].T.fix(303.2)
    m.Inlet_Feed['Toluene'].P.fix(350000)
    
    m.Inlet_Feed['Hydrogen'].flow['Liq','benzene'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Liq','toluene'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Liq','methane'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Liq','hydrogen'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Vap','benzene'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Vap','toluene'].fix(1e-8)
    m.Inlet_Feed['Hydrogen'].flow['Vap','methane'].fix(0.00618)
    m.Inlet_Feed['Hydrogen'].flow['Vap','hydrogen'].fix(0.0599)
    m.Inlet_Feed['Hydrogen'].T.fix(303.2)
    m.Inlet_Feed['Hydrogen'].P.fix(350000)

    # initialize flowsheet
    print('Initializing flowsheet...')
    print()

    # Initialization with Sequential Decomposition method
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic" 
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5

    print('Limiting Wegstein tear to 5 iterations to obtain initial solution,'
          ' if not converged IPOPT will pick up and continue.')
    print()
    
    # which stream to tear / what's the calculation order
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

    tear_guesses = {
            "flow_mol_phase_comp": {
                    (0, "Vap", "benzene"): 1e-5,
                    (0, "Vap", "toluene"): 1e-5,
                    (0, "Vap", "hydrogen"): 0.03665,
                    (0, "Vap", "methane"): 0.00618,
                    (0, "Liq", "benzene"): 1e-5,
                    (0, "Liq", "toluene"): 0.0599,
                    (0, "Liq", "hydrogen"): 1e-5,
                    (0, "Liq", "methane"): 1e-5},
            "temperature": {0: 303},
            "pressure": {0: 350000}}
    
    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.H101.Inlet, tear_guesses)

    # Initialization function for each unit
    def function(unit):
        if 'R' in unit.name:
            Initialize(unit,Arc=m.s4)
        elif 'F101' in unit.name:
            Initialize(unit,Arc=m.s5)
        elif 'F102' in unit.name:
            Initialize(unit,Arc=m.s9)
        elif 'S' in unit.name:
            Initialize(unit,Arc=m.s6)
        elif 'M' in unit.name:
            Initialize(unit,Arc=m.s7)
        elif 'H' in unit.name:
            Initialize(unit,Arc=m.s3)

    seq.run(m, function)

    # define economical items and constraints
    mw_comp_data = {'benzene': 78.1136,
                    'toluene': 92.1405,
                    'hydrogen': 2.016,
                    'methane': 16.043}
    @m.Expression()
    def purity(m):
        return m.F102.outlet.flow['Vap','benzene']*mw_comp_data['benzene']/sum(m.F102.outlet.flow['Vap',c] * mw_comp_data[c] for c in m.Comp)
    
    @m.Expression()
    def loss(m):
        return m.S101.Purge.flow['Vap','benzene']/m.R101.Outlet.flow['Vap','benzene']
    
    @m.Expression()
    def product(m):
        return 3600 * m.F102.outlet.phase_flow["Vap"]  # kmol/hr
    
    # set cooling cost for F101 & R101
    # 单价 元/kJ
    # 0.081 × 0.01/3600 = $0.22/GJ
    @m.Expression()
    def cooling_cost(m):
        return 0.22e-6 * (-m.F101.Q) + 0.22e-6 * (-m.R101.Q)
    
    # set heating cost for H101 & F102
    @m.Expression()
    def heating_cost(m):
        return 4e-6 * m.H101.Q + 1.9e-6 * m.F102.Q
    
    @m.Expression()
    def operating_cost(m):
        return 3600 * 8000 * (m.heating_cost + m.cooling_cost)
    
    # define the capital cost of 4 main units
    cost_vessel(m.R101,
                vertical=False,
                material_type='StainlessSteel304',
                shell_thickness=1.25, # inch
                weight_limit=1,
                aspect_ratio_range=1,
                vessel_diameter=9.53,
                vessel_length=57, #ft
                number_of_units=1,
                R_Ft = 1.6)
    cost_vessel(m.F101,
                vertical=True,
                material_type='CarbonSteel',
                shell_thickness=1.25, # inch
                weight_limit=1,
                aspect_ratio_range=1,
                vessel_volume=40,
                number_of_units=1)
    cost_vessel(m.F102,
                vertical=True,
                material_type='CarbonSteel',
                shell_thickness=1.25, # inch
                weight_limit=1,
                aspect_ratio_range=1,
                vessel_volume=100,
                number_of_units=1)
    cost_fired_heater(m.H101,
                      heat_source='Fuel',
                      material_type='StainlessSteel',
                      integer=True)
    
    # solve model
    print('Solving flowsheet...')
    print()

    # 模型求解
    solver = SolverFactory('ipopt')
    solver.options['tol'] = 1e-3
    results = solver.solve(m,tee=True)

    print('Complete.')
    print()

    return m