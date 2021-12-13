
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
    Ajet = pi * ((D/2)**2) # < - group of inlets in an inlet
    
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

def Reynolds_sympy(density, mu, D, u = None, Re = None):
    
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

def massflow_nextrows(coordinates, massflow, cp, Ahex, Tpeak, Tinput, HFf, massflow_actual, t = 0.001, ks = 340):
    
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

def pdrop(T, P, D, u, density):
    
    """
    Pressure drop equation, used to calculate the pressure drop, Reynolds \
        number and Euler number
    
    :param float T: Coolant temperature
    :param float P: Coolant pressure
    :param float D: Jet diameter
    :param float u: Jet velocity
    :param float density: Coolant density
    :param float delta_p: Pressure drop
    :param float Re: Reynolds number
    :param float Eu: Euler number
    """
    
    Re = Reynolds(density, P, T, D, u)
    Eu = Euler(Re)
    # u = velocity(massflow, P, T, D)
    delta_p = Eu * density * (u ** 2)
    
    return delta_p, Re, Eu

def coolant_temperature(massFlow, channels, Ti, Pi, inner_radius, outer_radius, \
                      cp, hcp_final, A, D, q, s, t):
    
    """
    Flow properties calculator, returning the developed coolant temperature, \
        pressure, armour tile peak temperature, mass flow per hexagon, non-\
            dimensional mass flow, Reynolds number, coolant density and jet \
                velocity. This is calculated for all hexagons in all channels.
    
    :param float massFlow: Total mass flow for a single plate
    :param float Ti: Coolant temperature
    :param float Pi: Coolant pressure
    :param float inner_radius: Radial location of the plate inboard edge
    :param float outer_radius: Radial location of the plate outboard edge
    :param float cp: Coolant specific heat capacity
    :param float A: Heat input area
    :param float D: Jet diameter
    :param float t: Armour tile thickness
    :param list channels: Coordinate locations for hexagons in given channels, \
        separated by channel
    :param list hcp_final: Coordinate locations for all hexagons
    :param list q: Heat flux density
    :param list s: Segmented distance from R = 0 to outer_radius, defined by \
        the locations of the heat flux density values held in q
    """
    values_per_channel = []
    
    massFlow = massFlow/2 # to account for the fact that the calc is using half the plate
    Ajet = pi*((D/2)**2)
    
    HFf = interp1d(s, q)
    no_hexagons = 0
    total_power = 0
    
    # calc temperature of the coolant
    
    for direction in channels:
        Vjet = 25.0
        count = 0
        massflow_total = 0
        Tinput = Ti
        Pinput = Pi
        for channel in direction:
            if count == 0:
                density = SI('D', 'T', Tinput, 'P', Pinput, 'helium')
                for coordinates in channel:
                    Ri = sqrt(coordinates[0]**2 + coordinates[1]**2)
                    HFi = HFf(Ri) * A
                    Pdrop, massflow, Re, Vjet = pressure_drop(Tinput, Pinput, D, Ajet, density, "V", u = Vjet)
                    Tnew = Tinput + HFi/(massflow * SI('C', 'T', Tinput, 'P', Pinput, 'helium'))
                    mstar = m_star(massflow, t, Ajet, Tinput, Pinput)
                    Tmetal = metal_temperature(mstar, Tinput, Tnew)
                    coordinates.append([Tnew, Pdrop, Tmetal, massflow, mstar, \
                                        Re, density, Vjet])
                    check_cases(Re, mstar)
                    no_hexagons += 1
                    total_power += HFi
                    massflow_total += massflow
                count += 1
            else:
                massflow = massflow_total / 2 / len(channel)
                density = SI('D', 'T', Tinput, 'P', Pinput, 'helium')
                for coordinates in tqdm(channel):
                    Ri = sqrt(coordinates[0]**2 + coordinates[1]**2)
                    HFi = HFf(Ri) * A
                    Pdrop, massflow, Re, Vjet = pressure_drop(Tinput, Pinput, D, Ajet, density, "m", massflow = massflow)
                    Tnew = Tinput + HFi/(massflow * SI('C', 'T', Tinput, 'P', Pinput, 'helium'))
                    mstar = m_star(massflow, t, Ajet, Tinput, Pinput)
                    Tmetal = metal_temperature(mstar, Tinput, Tnew)
                    coordinates.append([Tnew, Pdrop, Tmetal, massflow, mstar, \
                                        Re, density, Vjet])
                    # check_cases(Re, mstar)
                    no_hexagons += 1
                    total_power += HFi
            
            values_per_channel.append([Tnew, Pdrop, massflow, mstar, \
                                        Re, density])
            Tinput = sum([temp[2][0] for temp in channel]) / len(channel)
            Pinput = Pinput - (sum([pressure[2][1] for pressure in channel])  / \
                 len(channel))

    return channels, s, q, values_per_channel

def mass_flow(mstar, t, A, T, P):
    
    """
    Mass flow equation, used to calculate the mass flow per hexagon
    
    :param float T: Coolant temperature
    :param float P: Coolant pressure
    :param float mstar: Non-dimensional mass flow
    :param float t: Armour tile thickness
    :param float A: Jet cross-sectional area
    """
    
    massflow = (mstar * 400 * A) / (SI('C', 'T', T, 'P', P, 'helium') * t)
    
    return massflow

def pressure_drop(T, P, D, Ajet, density, version, u=0, massflow=0):
    
    """
    Pressure drop equation, may be used to calculate the pressure drop, mass \
        flow, Reynolds number or velocity
    
    :param float T: Coolant temperature
    :param float P: Coolant pressure
    :param float D: Jet diameter
    :param float Ajet: Jet cross-sectional area
    :param float density: Coolant density
    :param float u: Jet velocity
    :param float massflow: Coolant mass flow
    """
    
    if version == "V":
        massflow = density * u * Ajet
        Re = Reynolds(density, P, T, D, u)
        Eu = Euler(Re)
        # u = velocity(massflow, P, T, D)
        delta_p = Eu * density * (u ** 2)
    elif version == "m":
        u = massflow / (density * Ajet)
        Re = Reynolds(density, P, T, D, u)
        Eu = Euler(Re)
        # u = velocity(massflow, P, T, D)
        delta_p = Eu * density * (u ** 2)
    else:
        print("please enter a version type")
    
    return delta_p, massflow, Re, u

def velocity(massflow, roe, A):
    
    """
    Velocity equation, used to calculate the velocity
    
    :param float massflow: Coolant mass flow
    :param float roe: Coolant density
    :param float A: Jet cross-sectional area
    """
    
    u = massflow / (roe * A)
    
    return u

def metal_temperature(mstar, Tinput, Tnew):
    
    """
    Armour tile temperature equation, used to calculate the Armour tile \
        temperature
    
    :param float mstar: Non-dimensional mass flow
    :param float Tinput: Coolant temperature
    :param float Tnew: Developed coolant temperature
    :param float taus: Non-dimensional temperature coefficient
    :param float Tpeak: Peak temperature of the armour tile
    """
    
    taus = (4.415 * (mstar ** 0.6242)) + 0.7828
    Tpeak = (taus * (Tnew - Tinput)) + Tinput
    
    return Tpeak

def Reynolds(density, P, T, D, u):
    
    """
    Reynolds number equation, used to calculate the Reynolds number
    
    :param float density: Coolant density
    :param float T: Coolant temperature
    :param float P: Coolant pressure
    :param float D: Jet diameter
    :param float u: Jet velocity
    :param float Re: Reynolds number
    """
    
    Re = (density * u * D) / SI('V', 'T', T, 'P', P, 'helium')
    
    return Re

def Euler(Re):
    
    """
    Euler number equation, used to calculate the Euler number
    
    :param float Eu: Euler number
    :param float Re: Reynolds number
    """
    
    Eu = (1.486 * (Re ** -0.0721)) - 0.001133
    
    return Eu

def m_star(massflow, t, A, T, P):
    
    """
    Non-dimensional mass flow equation, used to calculate the non-dimensional \
        mass flow
    
    :param float massflow: Coolant mass flow
    :param float T: Coolant temperature
    :param float P: Coolant pressure
    :param float t: Armour tile thickness
    :param float A: Heat input cross-sectional area
    """
    
    mstar = (massflow * SI('C', 'T', T, 'P', P, 'helium') * t) / (\
            400 * A)
    
    return mstar

def check_cases(Re, mstar, Ma):
    
    """
    Range of validity check for the calculated Reynolds number and non-\
        dimensional mass flow
    
    :param float mstar: Non-dimensional mass flow
    :param float Re: Reynolds number
    """
    
    if Re >= 40000 and Re <= 120000:
        None
    else:
        print("failed Reynolds check " + str(Re))
    
    if mstar >= 0.45 and mstar <= 1.22:
        None
    else:
        print("failed mass flow check " + str(mstar))
    
    if Ma <= 0.3:
        print
    else:
        print("failed Ma check " + str(Ma))
    
    print([Re, mstar, Ma])
    
    return