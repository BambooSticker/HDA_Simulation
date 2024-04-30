import numpy as np
import pandas as pd
from pyomo.environ import *
from pyomo.network import *
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.util.initialization import propagate_state

pi = 3.1415926535
Phase = ['Liq','Vap']
Comp = ['benzene','toluene','methane','hydrogen']


# -------------------------------- initialization ------------------------------------
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
    solver.options['print_user_options'] = 'no'

    
    #if b.name == Order:
    # propagate the stream propeties from the former unit
    
    #propagate_state(arc=Arc)

    if 'R' in b.name:
        for c in Comp:
            calculate_variable_from_constraint(b.outlet.flow['Liq',c],b.rule_mass_balance_liq[c])
            calculate_variable_from_constraint(b.outlet.flow['Vap',c],b.rule_mass_balance_vap[c])

    elif 'F101' in b.name:
        # Flash calculation: T_bub -> T_dew -> T_1 -> T_eq
        calculate_variable_from_constraint(b.temperature_dew, b.rule_dew_T)
        calculate_variable_from_constraint(b._t1, b.rule_t1)
        calculate_variable_from_constraint(b._teq, b.rule_teq)

        b.Inlet.fix()
        flag = True
        while flag:
            solver.options['print_user_options'] = 'no'
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
        b.Inlet.unfix()

    elif 'Inlet_Feed' in b.name:
        n=0
        
    elif 'M' in b.name:
        b.Recycle_Feed.fix()
        flag = True
        while flag:
            solver.options['print_user_options'] = 'no'
            status = solver.solve(b)
            if (status.solver.termination_condition == TerminationCondition.optimal or 
                status.solver.status == SolverStatus.ok):
                flag = False
    
        b.Recycle_Feed.unfix()
        
    else:
        b.Inlet.fix()
        
        status = solver.solve(b)
            
        b.Inlet.unfix()
    

    #print("Initialization complete for '%s'! " %(b.name))

# -------------------------------- get optimal outcome ------------------------------------
def print_outcome(m):
    '''Print the optimal condition of unfixed units

    Args:
        m: flowsheet pyo.model
    
    Returns:
        None
    
    '''
    print("Decision Variables: ")
    print(f"     H101 outlet temperature = {value(m.H101.Outlet.T):0.2f}K")
    print(f"     R101 outlet temperature = {value(m.R101.Outlet.T):0.2f}K")
    print(f"     F101 outlet temperature = {value(m.F101.Vap_Outlet.T):0.2f}K")
    print(f"     F102 outlet temperature = {value(m.F102.Vap_Outlet.T):0.2f}K")
    print(f"     F102 outlet pressure = {value(m.F102.Vap_Outlet.P)/1e3:0.2f}kPa")
    print(f"With total inlet feed {0.10273*3600:0.2f} kmol/hr, ")
    print(f'Cost ${value(m.Unit_cost):0.2f} per tonne benzene produced.')
    print()
    print("Outcome: ")
    print(f"\t Total Annual Cost = ${value(m.Total_annual_cost):0.2f}")
    print(f"\t Benzene Purity  = {value(m.purity):0.2f}")
    print(f"\t Benzene Product = {value(m.product):0.2f}")
    print(f"\t Overhead Loss   = {value(m.loss):0.2f}")
    print()
    #production = m.F102.Vap_Outlet.flow["Vap", "benzene"] * 8000 * 3600 * 78 / 1e3

# -------------------------------- basic info ------------------------------------
Phase = ['Vap', 'Liq']
Comp = ['benzene','toluene','methane','hydrogen']
mw_comp_data = {'benzene': 78.1136,
                    'toluene': 92.1405,
                    'hydrogen': 2.016,
                    'methane': 16.043} # kmol/kg

# -------------------------------- report a single stream ------------------------------------
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
        
    # mole flowrate
    print("flow rate[kmol/hr]: ", )
    for p in Phase:
        for c in Comp:
            print('\t [',p,',',c,']: ',round(value(b.dest.flow[p,c])*3600,3))
            
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


# -------------------------------- get stream information ------------------------------------
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
        mole_flowrate = value(sum(b.dest.flow[p,c] for p in Phase for c in Comp))*3600
        data[s].append(round(mole_flowrate,3))
    
        data[s].append('-')
        
        # mass flowrate
        mass_flow = {}
        total_mass_flow = sum(b.dest.flow[p,c] * mw_comp_data[c] for p in Phase for c in Comp)*3600
        for c in Comp:
            val = round(value(sum(b.dest.flow[p,c] for p in Phase))* \
                             mw_comp_data[c]*3600,2)
            mass_flow[c] = val
            data[s].append(val)
            
        data[s].append('-')
        # mass fraction
        for c in Comp:
            mass_frac = round(value(sum(b.dest.flow[p,c]*mw_comp_data[c] for p in Phase)*3600/total_mass_flow),5)
            data[s].append(mass_frac)
    
    for block in m.component_data_objects(Arc, descend_into=False):
        PFD_creater(block)
        
    updated_data = {key: value + [0]*14 if not value else value for key, value in data.items()}
    df = pd.DataFrame(updated_data,index=Item)

    return df


# -------------------------------- cost for vessels ------------------------------------
def cost_vessel(
    blk,
    vertical=False,
    material_type='CarbonSteel',
    shell_thickness=1.25, # inch
    weight_limit=1,
    aspect_ratio_range=1,
    vessel_diameter=None,
    vessel_length=None,
    vessel_volume=None,
    number_of_units=1,
    R_Ft = 2.1,
    sensitivity = False):
    """
    Generic vessel costing method.

    Args:
        vertical: alignment of vessel; vertical if True, horizontal if
            False (default=False).
        material_type: VesselMaterial Enum indicating material of
            construction, default = VesselMaterial.CarbonSteel.
        shell_thickness: thickness of vessel shell, including pressure
            allowance. Default = 1.25 inches.
        weight_limit: 1: (default) 1000 to 920,000 lb, 2: 4200 to 1M lb.
            Option 2 is only valid for vertical vessels.
        aspect_ratio_range: vertical vessels only, default = 1;
            1: 3 < D < 21 ft, 12 < L < 40 ft, 2: 3 < D < 24 ft; 27 < L < 170 ft.
        include_platforms_ladders: whether to include platforms and
            ladders in costing , default = True.
        vessel_diameter: Pyomo component representing vessel diameter.
            If not provided, assumed to be named "diameter"
        vessel_length: Pyomo component representing vessel length.
            If not provided, assumed to be named "length".
        number_of_units: Integer or Pyomo component representing the
            number of parallel units to be costed, default = 1.
        number_of_trays: Pyomo component representing the number of
            distillation trays in vessel (default=None)
    """
    # Build generic costing variables
    blk.base_cost_per_unit = Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        # units=pyo.units.USD_CE500,
        doc="Base cost per unit",
    )

    blk.capital_cost = Var(
        initialize=1e4,
        domain=pyo.NonNegativeReals,
        bounds=(0, None),
        # units=pyo.units.USD_CE500,
        doc="Capital cost of all units",
    )
    
    blk.vessel_diameter = Var(initialize = 1, units=units.ft)
    if vessel_diameter != None:
        blk.vessel_diameter.fix(vessel_diameter)

    blk.vessel_length = Var(initialize = 1, units=units.ft)
    if vessel_length != None:
        blk.vessel_length.fix(vessel_length)

    D_in = pyo.units.convert(blk.vessel_diameter, to_units=pyo.units.inch)
    L_in = pyo.units.convert(blk.vessel_length, to_units=pyo.units.inch)
    
    blk.vessel_volume = Var(initialize = 1, units=units.ft**3)
    if vessel_volume != None:
        blk.vessel_volume.fix(vessel_volume)
        blk.L_over_D = Constraint(expr=L_in == 3 * D_in)

    V_in3 = pyo.units.convert(blk.vessel_volume, to_units=pyo.units.inch**3)

    blk.volume_eq = Constraint(expr=V_in3 == L_in * D_in**2 * 0.25 * pi)

        
    # Material densities in lb/cubic inch
    material_factor_dict = {
        'CarbonSteel': {"factor": 1.0, "density": 0.284}, ### 用这个
        'LowAlloySteel': {"factor": 1.2, "density": 0.271},
        'StainlessSteel304': {"factor": 1.7, "density": 0.270},
        'StainlessSteel316': {"factor": 2.1, "density": 0.276},
        'Carpenter20CB3': {"factor": 3.2, "density": 0.29},
        'Nickel200': {"factor": 5.4, "density": 0.3216},
        'Monel400': {"factor": 3.6, "density": 0.319},
        'Inconel600': {"factor": 3.9, "density": 0.3071},
        'Incoloy825': {"factor": 3.7, "density": 0.2903},
        'Titanium': {"factor": 7.7, "density": 0.1628},
    }
    material_factors = material_factor_dict[material_type]

    blk.temperature_factor = Param(default = 1.6, mutable=True)
    if 'R' in blk.name:
        blk.temperature_factor = R_Ft

        ##########################
    if sensitivity == True:
        blk.del_component(blk.temperature_factor)
        blk.BigM = Param(default=800)
    
        blk.y1 = Var(domain=Binary)
        blk.y2 = Var(domain=Binary)
    
        @blk.Constraint()
        def temperature_factor_ref1(blk):
            return blk.outlet.T <= 400 + blk.BigM*(1-blk.y1)
            
        @blk.Constraint()
        def temperature_factor_ref2(blk):
            return blk.outlet.T >= 400 - blk.BigM*(1-blk.y2)
    
        @blk.Constraint()
        def sum_y(blk):
            return blk.y1+blk.y2 == 1
    
        @blk.Expression()
        def temperature_factor(b):
            return 1.6+blk.y2*0.5
        ##########################
    #@blk.Expression()
    #def temperature_factor(b):
    #    if 'R' in b.name:
    #        return 2.1
    #    else:
    #        return 1.6
    
    blk.pressure_factor = Param(default = 1)
    # Users should calculate the pressure design based shell thickness
    # Pressure factor assumed to be included in thickness
    blk.shell_thickness = Param(
        mutable=True, doc="Shell thickness", 
        units=units.inch
    )
    blk.shell_thickness.set_value(shell_thickness)
    blk.material_factor = Param(
        initialize=material_factors["factor"],
        mutable=True,
        doc="Construction material correction factor",
    )
    blk.material_density = Param(
        initialize=material_factors["density"],
        mutable=True,
        doc="Density of the metal",
        units=units.pound / units.inch**3,
    )

    # Calculate weight of vessel
    blk.weight = Var(
        initialize=1000,
        domain=pyo.NonNegativeReals,
        doc="Weight of vessel in lb",
        # units=pyo.units.pound,
    )

    @blk.Constraint()
    def weight_eq(blk):
        return blk.weight == (
            pi
            * (D_in + blk.shell_thickness)
            * (L_in + 0.8 * D_in)
            * blk.shell_thickness
            * blk.material_density
        )

    # Base Vessel cost
    # Alpha factors for correlation
    # 2nd key is weight_limit option
    alpha_dict = {
        "H": {1: {1: 8.9552, 2: -0.2330, 3: 0.04333}, 2: None},
        "V": {
            1: {1: 7.0132, 2: 0.18255, 3: 0.02297},
            2: {1: 7.2756, 2: 0.18255, 3: 0.02297},
        },
    }

    if vertical:
        alpha = alpha_dict["V"][weight_limit]
    else:
        alpha = alpha_dict["H"][weight_limit]

    @blk.Constraint()
    def base_cost_constraint(blk):
        return blk.base_cost_per_unit == (
            exp(
                alpha[1]
                + alpha[2] * (log(blk.weight / units.pound))
                + alpha[3] * (log(blk.weight / units.pound) ** 2)
            )
             # * units.USD_CE500
        )

    # Total capital cost of vessel and ancillary equipment
    @blk.Constraint()
    def capital_cost_constraint(blk):
        cost_expr = blk.material_factor * blk.base_cost_per_unit * blk.temperature_factor * blk.pressure_factor
        return blk.capital_cost == cost_expr * number_of_units

# -------------------------------- creat common cost variables ------------------------------------
def _make_common_vars(blk, integer=True):
    # Build generic costing variables (most costing models need these vars)
    blk.base_cost_per_unit = Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        #units=pyo.units.USD_CE500,
        doc="Base cost per unit",
    )

    blk.capital_cost = Var(
        initialize=1e4,
        domain=pyo.NonNegativeReals,
        bounds=(0, None),
        #units=pyo.units.USD_CE500,
        doc="Capital cost of all units",
    )

    if integer is True:
        domain = Integers
    else:
        domain = NonNegativeReals
    blk.number_of_units = Var(
        initialize=1, domain=domain, bounds=(1, 100), doc="Number of units to install."
    )
    blk.number_of_units.fix(1)

# -------------------------------- cost for fired heater ------------------------------------
def cost_fired_heater(
    blk,
    heat_source='Fuel',
    material_type='CarbonSteel',
    integer=True,
):
    """
    Generic costing method for fired heaters.

    Args:
        heat_source: HeaterSource Enum indicating type of source of heat,
            default = HeaterSource.Fuel.
        material_type: HeaterMaterial Enum indicating material of
            construction, default = HeaterMaterial.CarbonSteel.
        integer: whether the number of units should be constrained to be
            an integer or not (default = True).
    """

    # Build generic costing variables
    _make_common_vars(blk, integer)
    
    # Convert pressure to psi,g
    #psi_to_Pa = 6896.757
    #blk.pressrue = Expression(expr=(blk.inlet.P - 101325)/psi_to_Pa)

    P = units.convert(
        blk.outlet.P,
        to_units=units.psi,
    ) - units.convert(1 * units.atm, to_units=units.psi)

    # Convert heat duty (kW) to BTU/hr
    #BTUhr_to_kW = 0.0855
    #blk.heat_duty = Expression(expr=(blk.Q/BTUhr_to_kW
    #    )
    #    / blk.number_of_units
    #)

    Q = (
        units.convert(
            blk.Q, to_units=units.BTU / units.hr
        )
        / blk.number_of_units
    )

    # Material factor
    material_factor_dict = {
        'CarbonSteel': 1.0,
        'CrMoSteel': 1.4,
        'StainlessSteel': 1.7,
    }
    blk.material_factor = Param(
        initialize=material_factor_dict[material_type],
        domain=pyo.NonNegativeReals,
        doc="Construction material correction factor",
    )

    # Pressure design factor calculation
    blk.pressure_factor = Var(
        initialize=1.1, domain=pyo.NonNegativeReals, doc="Pressure design factor"
    )

    @blk.Constraint()
    def pressure_factor_eq(blk):
        return blk.pressure_factor == (
            0.986
            - 0.0035 * (P / (500.00 * units.psi))
            + 0.0175 * (P / (500.00 * units.psi)) ** 2
        )

    @blk.Constraint()
    def base_cost_per_unit_eq(blk):
        if heat_source == 'Fuel':
            bc_expr = pyo.exp(
                0.32325 + 0.766 * log(Q / units.BTU * units.hr)
            )
        elif heat_source == 'Reformer':
            bc_expr = 0.859 * (Q / units.BTU * units.hr) ** 0.81
        elif heat_source == 'Pyrolysis':
            bc_expr = 0.650 * (Q / units.BTU * units.hr) ** 0.81
        elif heat_source == 'HotWater':
            bc_expr = exp(
                9.593
                - 0.3769 * log((Q / units.BTU * units.hr))
                + 0.03434 * log((Q / units.BTU * units.hr)) ** 2
            )
        elif heat_source == 'Salts':
            bc_expr = 12.32 * (Q / units.BTU * units.hr) ** 0.64
        elif heat_source == 'DowthermA':
            bc_expr = 12.74 * (Q / units.BTU * units.hr) ** 0.65
        elif heat_source == 'steamBoiler':
            bc_expr = 0.367 * (Q / units.BTU * units.hr) ** 0.77

        return blk.base_cost_per_unit == bc_expr

    @blk.Expression(doc="Base cost for all units installed")
    def base_cost(blk):
        return blk.base_cost_per_unit * blk.number_of_units

    # Total capital cost of heater(s)
    @blk.Constraint()
    def capital_cost_constraint(blk):
        return blk.capital_cost == (
            blk.material_factor * blk.pressure_factor * blk.base_cost
        )


# -------------------------------- get info of vessels ------------------------------------
def vessel_info(unit):
    '''
    get the information of a vessel unit (R101, F101 or F102).

    Arg:
        unit: single unit as pyo.Block. (e.g. m.R101)

    return:
        Information as DataFrame.
    '''
    
    info = {}
    # vertical or horizontal
    if 'R' in unit.name:
        # vertical or horizontal
        info['Installation'] = 'horizontal'
        # material_type
        info['Material Type'] = 'Stainless Steel 304'
        
    else:
        info['Installation'] = 'vertical'
        # material_type
        info['Material Type'] = 'CarbonSteel'
    
    # diameter
    info['Diameter'] = f'{unit.vessel_diameter.value:0.2f} ft'
    # length
    info['Length'] = f'{unit.vessel_length.value:0.2f} ft'
    # volume
    info['Volume'] = f'{unit.vessel_volume.value:0.2f} ft3'
    # shell_thickness
    info['Shell Thickness'] = f'{unit.shell_thickness.value:0.2f} inch'
    # cost
    info['Capital Cost'] = f'${unit.capital_cost.value:0.2f}' 
    vessel_data = pd.DataFrame(data=info, index=[0]).transpose()
    vessel_data.columns=[f'{unit.name}']
    return vessel_data

# -------------------------------- get info of heater ------------------------------------
def heater_info(unit):
    '''
    get the information of a vessel unit (H101 only).

    Arg:
        unit: single unit as pyo.Block. (e.g. m.H101)

    return:
        Information as DataFrame.
    '''
    info = {}
    
    # material_type
    info['Heat Source'] = 'Fuel'
    info['Material Type'] = 'Stainless Steel'
    
    # diameter
    info['Heat Duty'] = f'{unit.Q.value:0.2f} kW'
    # length
    info['Pressure'] = f'{unit.outlet.P.value/1e3:0.2f} kPa'
    # cost
    info['Capital Cost'] = f'${unit.capital_cost.value:0.2f}' 
    heater_data = pd.DataFrame(data=info, index=[0]).transpose()
    heater_data.columns=[f'{unit.name}']
    return heater_data