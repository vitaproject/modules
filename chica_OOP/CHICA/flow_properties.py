
from math import sqrt, pi
from CoolProp.CoolProp import PropsSI as SI
from sympy import symbols, Eq, solve
from scipy.interpolate import interp1d
from tqdm import tqdm

def mach_number_sympy(Tinput, Ma = None, u = None, R = 8.3145, M = 4.003E-3, gamma = 1.667):

    """
    Mach number equation, may be used to calculate Mach number or velocity
    
    :param float Tinput: Coolant temperature
    :param float Ma: Mach number 
    :param float u: Jet velocity
    :param float R: Universal gas constant
    :param float M: Molar mass
    :param float gamma: Isentropic expansion factor
    """
    
    Ma_calc, u_calc, gamma_calc, Ti_calc, R_calc, M_calc = \
        symbols("Ma_calc u_calc gamma_calc Ti_calc R_calc M_calc")
    
    eq = Eq(u_calc / ((gamma_calc * Ti_calc * (R_calc/M_calc))**0.5), Ma_calc)
    
    if Ma == None:
        Ma_calculated = solve(eq, Ma_calc)
        Ma = float(Ma_calculated[0].subs([(gamma_calc, gamma),\
                                          (u_calc, u), \
                                          (Ti_calc, Tinput), \
                                          (R_calc, R), \
                                          (M_calc, M)]))
    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(gamma_calc, gamma),\
                                        (Ma_calc, Ma), \
                                        (Ti_calc, Tinput), \
                                        (R_calc, R), \
                                        (M_calc, M)]))
    return Ma, u

def mass_flow_sympy(Ti, Pi, D = 0.001, massflow = None, u = None):
    
    """
    Mass flow equation, may be used to calculate mass flow or velocity
    
    :param float Ti: Coolant temperature
    :param float Pi: Coolant pressure
    :param float D: Jet diameter
    :param float massflow: Mass flow of coolant
    :param float u: Jet velocity
    :param float Ajet: Jet cross-sectional area
    """
    
    density = SI("D", "T", Ti, "P", Pi, "helium")
    Ajet = pi * ((D/2)**2) # * 13 # < - group of inlets in an inlet
    
    roe, u_calc, massflow_calc, Ajet_calc = \
        symbols("roe u_calc massflow_calc Ajet_calc")
    
    eq = Eq(roe * u_calc * Ajet_calc, massflow_calc)
    
    if massflow == None:
        massflow_calculated = solve(eq, massflow_calc)
        massflow = float(massflow_calculated[0].subs([(roe, density),\
                                                      (u_calc, u),\
                                                      (Ajet_calc, Ajet)]))
    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(roe, density),\
                                        (massflow_calc, massflow),\
                                        (Ajet_calc, Ajet)]))
    return massflow, u, Ajet, density, D

def mstar_sympy(Tinput, Pinput, A, D = 0.001, t = 0.001, ks = 340, mstar = None, massflow = None):
    
    """
    Non-dimensional mass flow equation, may be used to calculate \
        Non-dimensional mass flow or mass flow. Also returns specific heat \
            capacity
    
    :param float Tinput: Coolant temperature
    :param float Pinput: Coolant pressure
    :param float A: Heat input cross-sectional area
    :param float D: Jet diameter
    :param float t: Armour tile thickness
    :param float ks: Armour tile thermal conduction
    :param float mstar: Non-dimensional mass flow
    :param float massflow: Coolant mass flow
    :param float cp: Specific heat capacity
    """
    
    cp = SI("C", "T", Tinput, "P", Pinput, "helium")
    
    mstar_calc, massflow_calc, cp_calc, ks_calc, A_calc, t_calc = \
        symbols("mstar_calc massflow_calc cp_calc ks_calc A_calc t_calc")
    
    eq = Eq((massflow_calc * cp_calc * t_calc) / (ks_calc * A_calc), mstar_calc)
    
    if mstar == None:
        mstar_calculated = solve(eq, mstar_calc)
        mstar = float(mstar_calculated[0].subs([(massflow_calc, massflow),\
                                                (cp_calc, cp), \
                                                (ks_calc, ks), \
                                                (A_calc, A), \
                                                (t_calc, t)]))
    else:
        massflow_calculated = solve(eq, massflow_calc)
        massflow = float(massflow_calculated[0].subs([(mstar_calc, mstar),\
                                                      (cp_calc, cp), \
                                                      (ks_calc, ks), \
                                                      (A_calc, A), \
                                                      (t_calc, t)]))
    return mstar, massflow, cp

def taus_sympy(cell_type, taus = None, mstar = None):
    
    """
    Non-dimensional temperature coefficient equation, may be used to calculate \
        Non-dimensional temperature or Non-dimensional mass flow
    
    :param float mstar: Non-dimensional mass flow 
    :param float taus: Non-dimensional temperature coefficient
    """
    
    if cell_type == "DLH":
    
        mstar_calc, taus_calc = symbols("mstar_calc taus_calc")
        
        eq = Eq((4.415 * (mstar_calc ** 0.6242)) + 0.7828, taus_calc)
    
        if taus == None:
            taus_calculated = solve(eq, taus_calc)
            taus = float(taus_calculated[0].subs(mstar_calc, mstar))
        else:
            mstar_calculated = solve(eq, mstar_calc)
            mstar = float(mstar_calculated[0].subs(taus_calc, taus))
    
    elif cell_type == "JIVC":
        
        mstar_calc, taus_calc = symbols("mstar_calc taus_calc")
        
        eq = Eq((5.83 * (mstar_calc ** 0.7648)) + 1.743, taus_calc)
    
        if taus == None:
            taus_calculated = solve(eq, taus_calc)
            taus = float(taus_calculated[0].subs(mstar_calc, mstar))
        else:
            mstar_calculated = solve(eq, mstar_calc)
            mstar = float(mstar_calculated[0].subs(taus_calc, taus))
    
    else:
        raise ValueError("Error, no cell type assigned")
    
    return taus, mstar

def metal_temperature_sympy(Tinput, taus, mass_flow, Q, Cp):

    """
    heat sink surface temperature equation, may be used to calculate heat sink\
        surface temperature or developed coolant temperature

    :param float Tinput: Inlet coolant temperature
    :param float taus: Non-dimensional temperature coefficient
    :param float Tnew: Developed coolant temperature
    :param float Tpeak: Peak temperature of the heat sink surface
    """

    Tpeak = ((taus * Q) / (mass_flow * Cp)) + Tinput

    return Tpeak

def pdrop_sympy(density, Eu = None, u = None, delta_p = None):
    
    """
    Pressure drop equation, may be used to calculate pressure drop, Euler \
        number or velocity
    
    :param float density: Coolant density
    :param float Eu: Euler number
    :param float u: Jet velocity
    :param float delta_p: Pressure drop
    """
    
    delta_p_calc, Eu_calc, u_calc = symbols("delta_p_calc Eu_calc u_calc")
    
    eq = Eq(Eu_calc * density * (u_calc ** 2), delta_p_calc)
    
    if Eu == None:
        Eu_calculated = solve(eq, Eu_calc)
        Eu = float(Eu_calculated[0].subs([(u_calc, u), \
                                         (delta_p_calc, delta_p)]))
    elif delta_p == None:
        delta_p_calculated = solve(eq, delta_p_calc)
        delta_p = float(delta_p_calculated[0].subs([(Eu_calc, Eu), \
                                                   (u_calc, u)]))
    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(Eu_calc, Eu), \
                                       (delta_p_calc, delta_p)]))
        
    return delta_p, Eu, u

def Reynolds_sympy(density, mu, D = 0.001, u = None, Re = None):
    
    """
    Reynolds number equation, may be used to calculate Reynolds number or \
        velocity
    
    :param float density: Coolant density
    :param float mu: Coolant viscosity
    :param float D: Jet diameter
    :param float u: Jet velocity
    :param float Re: Reynolds number
    """
    
    Re_calc, u_calc = symbols("Re_calc u_calc")
    
    eq = Eq((density * u_calc * D) / mu, Re_calc)
    
    if Re == None:
        Re_calculated = solve(eq, Re_calc)
        Re = float(Re_calculated[0].subs(u_calc, u))
    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs(Re_calc, Re))
    
    return Re, u

def Euler_sympy(cell_type, Eu = None, Re = None):
    
    """
    Euler number equation, may be used to calculate Euler number or Reynolds \
        number
    
    :param float Eu: Euler number
    :param float Re: Reynolds number
    """
    if cell_type == "DLH":
        Eu_calc, Re_calc = symbols("Eu_calc Re_calc")
        
        eq = Eq((1.486 * (Re_calc ** -0.0721)) - 0.001133, Eu_calc)
        
        if Eu == None:
            Eu_calculated = solve(eq, Eu_calc)
            Eu = float(Eu_calculated[0].subs(Re_calc, Re))
        else:
            Re_calculated = solve(eq, Re_calc)
            Re = float(Re_calculated[0].subs(Eu_calc, Eu))
    
    elif cell_type == "JIVC":
        Eu_calc, Re_calc = symbols("Eu_calc Re_calc")
        
        eq = Eq(11.85 * (Re_calc ** -0.2351), Eu_calc)
        
        if Eu == None:
            Eu_calculated = solve(eq, Eu_calc)
            Eu = float(Eu_calculated[0].subs(Re_calc, Re))
        else:
            Re_calculated = solve(eq, Re_calc)
            Re = float(Re_calculated[0].subs(Eu_calc, Eu))
    
    else:
        raise ValueError("Error, no cell type assigned")
    
    return Eu, Re

def coolant_temperature_sympy(taus, Tinput, metal_temperature):
    
    """
    Developed coolant temperature equation, may be used to calculate the \
        developed coolant temperature or mass flow
    
    :param float 
    """

    Tnew = Tinput + ((metal_temperature - Tinput) / taus)
    
    return Tnew

def massflow_nextrows_DLH(coordinates, massflow, cp, Ahex, Tpeak, Tinput, HFf, \
                          massflow_actual, t = 0.001, ks = 340):
    
    """
    Mass flow per hexagon equation for rows beyond the first row, may be used \
        to calculate the number of hexagons in the next row and the mass flow \
            per hexagon
    
    :param float coordinates: Hexagon coordinates
    :param float massflow: Coolant mass flow
    :param float cp: Coolant specific heat capacity
    :param float Ahex: Area of a hexagon
    :param float Tpeak: Peak temperature of the armour tile
    :param float Tinput: Inlet coolant temperature 
    :param float massflow_actual: Total mass flow inboard/outbaord dependant \
        on the selected split
    :param float t: Armour tile thickness
    :param float ks: Armour tile thermal conduction
    :param scipy.interpolate.interpolate.interp1d HFf: Object containing Heat \
        flux versus radial point 
    """
    
    Ri = sqrt(coordinates[0]**2 + coordinates[2]**2)
    HFi = HFf(Ri) * Ahex
    
    A = massflow
    B = 4.415 * (((cp * t) / (ks * Ahex))**0.6242)
    C = 0.7828
    D = HFi / (cp * (Tpeak - Tinput))
    
    Anew = A - ((A - ((A**0.6242) * D * B) - (D * C))/(1 - (0.6242 * (A **-0.3758) * D * B)))
    
    if D == 0:
        massflow = 0.0
        ni = 0.0
    else:   
        while sqrt((Anew - A)**2) >= 1e-10:
            A = Anew
            Anew = A - ((A - ((A**0.6242) * D * B) - (D * C))/(1 - (0.6242 * (A **-0.3758) * D * B)))
        massflow = Anew
        ni = massflow_actual / massflow
    
    return massflow, ni

def mass_flow_next_row_JIVC(coordinates, cross_sectional_area, HFf, mass_flow, \
                            cp, Tinput, Tpeak, mass_flow_actual, cell_type, t = 0.001, ks = 340):
    
    if cell_type == "DLH":
        phi, theta, omega = [4.415, 0.6242, 0.7828]
    else:
        phi, theta, omega = [5.83, 0.7648, 1.743]
    
    Ri = sqrt(coordinates[0]**2 + coordinates[2]**2)
    HFi = HFf(Ri) * cross_sectional_area
    
    A = mass_flow
    B = phi * (((cp * t) / (ks * cross_sectional_area))**theta)
    C = omega
    D = HFi / (cp * (Tpeak - Tinput))
    
    Anew = lambda A, B, C, D: A - ((A - ((A**theta) * D * B) - (D * C))/(1 - (theta * (A **-(1-theta)) * D * B)))
    A_updated = Anew(A, B, C, D)
    
    if D == 0:
        massflow = 0.0
        ni = 0.0
    else:   
        while sqrt((A_updated - A)**2) >= 1e-10:
            A = A_updated
            A_updated = Anew(A, B, C, D)
        massflow = A_updated
        ni = mass_flow_actual / massflow ## in this case, number of jets total
    
    return massflow, ni