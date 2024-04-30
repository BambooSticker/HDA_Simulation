from pyomo.environ import *
import idaes
import numpy as np
import math

# Initialize bubble and dew point calculations, t1 and teq
def Initialization(m):
    # Antoine equation
    def antoine_P(b, j, T):
        return units.convert_value(
                    value(
                        10
                        ** (
                            b.antoine_coeff_A[j]
                            - b.antoine_coeff_B[j]
                            / (T + b.antoine_coeff_C[j])
                        )
                    ),
                    from_units=units.bar,
                    to_units=units.Pa,
                )
    
    for k in m.values():
        # Bubble temperature initialization
        if hasattr(k, "_mole_frac_tbub"):
            Tbub0 = 0
            for j in k.Comp:
                    Tbub0 += value(
                        k.mole_frac_comp[j]
                        * (
                            k.antoine_coeff_B[j]
                            / (
                                k.antoine_coeff_A[j]
                                - math.log10(
                                    value(
                                        units.convert(
                                            k.P, to_units=units.bar
                                        )
                                    )
                                )
                            )
                            - k.antoine_coeff_C[j]
                        )
                    )
    
            err = 1
            counter = 0
    
            while err > 1e-2 and counter < 100:
                f = value(
                    sum(
                        antoine_P(k, j, Tbub0) * k.mole_frac_comp[j]
                        for j in k.Comp
                    )
                    - k.P
                )
                df = value(
                    sum(
                        k.mole_frac_comp[j]
                        * k.antoine_coeff_B[j]
                        * math.log(10)
                        * antoine_P(k, j, Tbub0)
                        / (Tbub0 + k.antoine_coeff_C[j]) ** 2
                        for j in k.Comp
                    )
                )
    
                if f / df > 20:
                    Tbub1 = Tbub0 - 20
                elif f / df < -20:
                    Tbub1 = Tbub0 + 20
                else:
                    Tbub1 = Tbub0 - f / df
    
                err = abs(Tbub1 - Tbub0)
                Tbub0 = Tbub1
                counter += 1
    
            k.temperature_bubble.value = Tbub0
    
            for j in k.Comp:
                k._mole_frac_tbub[j].value = value(
                    k.mole_frac_comp[j] * antoine_P(k, j, Tbub0) / k.P
                )
    
        # Dew temperature initialization
        if hasattr(k, "_mole_frac_tdew"):
            Tdew0 = 0
            for j in k.Comp:
                Tdew0 += value(
                    k.mole_frac_comp[j]
                    * (
                        k.antoine_coeff_B[j]
                        / (
                            k.antoine_coeff_A[j]
                            - math.log10(
                                value(
                                    units.convert(
                                        k.P, to_units=units.bar
                                    )
                                )
                            )
                        )
                        - k.antoine_coeff_C[j]
                    )
                )
    
            err = 1
            counter = 0
    
            while err > 1e-2 and counter < 100:
                f = value(
                    k.P
                    * sum(
                        k.mole_frac_comp[j] / antoine_P(k, j, Tdew0)
                        for j in k.Comp
                    )
                    - 1
                )
                df = -value(
                    k.P
                    * math.log(10)
                    * sum(
                        k.mole_frac_comp[j]
                        * k.antoine_coeff_B[j]
                        / (
                            (Tdew0 + k.antoine_coeff_C[j]) ** 2
                            * antoine_P(k, j, Tdew0)
                        )
                        for j in k.Comp
                    )
                )
    
                if f / df > 20:
                    Tdew1 = Tdew0 - 20
                elif f / df < -20:
                    Tdew1 = Tdew0 + 20
                else:
                    Tdew1 = Tdew0 - f / df
    
                err = abs(Tdew1 - Tdew0)
                Tdew0 = Tdew1
                counter += 1
    
            k.temperature_dew.value = Tdew0
    
            for j in k.Comp:
                k._mole_frac_tdew[j].value = value(
                    k.mole_frac_comp[j] * k.P / antoine_P(k, j, Tdew0)
                )
                
    m._t1.value = max(m.T.value, m.temperature_bubble.value)
    m._teq.value = min(m._t1.value, m.temperature_dew.value)