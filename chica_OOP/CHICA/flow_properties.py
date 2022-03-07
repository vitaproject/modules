from math import sqrt, pi
from sympy import symbols, Eq, solve

def mach_number_sympy(group_data, trigger = "mach_number", R = 8.3145, M = 4.003E-3, gamma = 1.667):

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

    if trigger == "mach_number":
        Ma_calculated = solve(eq, Ma_calc)
        Ma = float(Ma_calculated[0].subs([(gamma_calc, gamma),\
                                          (u_calc, group_data["jet_velocity"]), \
                                          (Ti_calc, group_data["input_temperature"]), \
                                          (R_calc, R), \
                                          (M_calc, M)]))
        group_data["mach_number"] = Ma

    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(gamma_calc, gamma),\
                                        (Ma_calc, group_data["mach_number"]), \
                                        (Ti_calc, group_data["input_temperature"]), \
                                        (R_calc, R), \
                                        (M_calc, M)]))
        group_data["jet_velocity"] = u

def mass_flow_sympy(group_data, trigger = "mass_flow", D = 0.001):
    
    """
    Mass flow equation, may be used to calculate mass flow or velocity
    
    :param float Ti: Coolant temperature
    :param float Pi: Coolant pressure
    :param float D: Jet diameter
    :param float massflow: Mass flow of coolant
    :param float u: Jet velocity
    :param float Ajet: Jet cross-sectional area
    """
    
    # density = SI("D", "T", group_data["input_pressure"], "P",\
    #               group_data["input_temperature"], "helium")
    density = group_data["fluid_properties"].density
    Ajet = pi * ((D/2)**2) # * 13 # < - group of inlets in an inlet
    
    roe, u_calc, mass_flow_calc, Ajet_calc = \
        symbols("roe u_calc mass_flow_calc Ajet_calc")
    
    eq = Eq(roe * u_calc * Ajet_calc, mass_flow_calc)
    
    if trigger == "mass_flow":
        mass_flow_calculated = solve(eq, mass_flow_calc)
        mass_flow = float(mass_flow_calculated[0].subs([(roe, density),\
                                                        (u_calc, group_data["jet_velocity"]),\
                                                        (Ajet_calc, Ajet)]))
        group_data["mass_flow"] = mass_flow

    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(roe, density),\
                                        (mass_flow_calc, group_data["mass_flow"]),\
                                        (Ajet_calc, Ajet)]))
        group_data["jet_velocity"] = u
    
    group_data["jet_cross_sectional_area"] = Ajet
    group_data["density"] = density
    group_data["jet_diameter"] = D
        
def m_star_sympy(group_data, trigger = "non_dimensional_mass_flow", D = 0.001, t = 0.001, ks = 340):
    
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

    cp = group_data["fluid_properties"].specific_heat_capacity

    m_star_calc, mass_flow_calc, cp_calc, ks_calc, A_calc, t_calc = \
        symbols("m_star_calc mass_flow_calc cp_calc ks_calc A_calc t_calc")

    eq = Eq((mass_flow_calc * cp_calc * t_calc) / (ks_calc * A_calc), m_star_calc)

    if trigger == "non_dimensional_mass_flow":

        m_star_calculated = solve(eq, m_star_calc)
        m_star = float(m_star_calculated[0].subs([(mass_flow_calc, group_data["mass_flow"]),
            (cp_calc, cp), (ks_calc, ks),
            (A_calc, group_data["cross_sectional_area"]), (t_calc, t)]))

        group_data["non_dimensional_mass_flow"] = m_star

    else:

        mass_flow_calculated = solve(eq, mass_flow_calc)
        mass_flow = float(mass_flow_calculated[0].subs([(m_star_calc, group_data["non_dimensional_mass_flow"]),
            (cp_calc, cp), (ks_calc, ks),
            (A_calc, group_data["cross_sectional_area"]), (t_calc, t)]))

        group_data["mass_flow"] = mass_flow

def taus_sympy(group_data, trigger = "non_dimensional_temperature"):

    """
    Non-dimensional temperature coefficient equation, may be used to calculate \
        Non-dimensional temperature or Non-dimensional mass flow

    :param float mstar: Non-dimensional mass flow 
    :param float taus: Non-dimensional temperature coefficient
    """

    m_star_calc, taus_calc = symbols("m_star_calc taus_calc")

    if group_data["cell_type"] == "DLH":

        eq = Eq((4.415 * (m_star_calc ** 0.6242)) + 0.7828, taus_calc)

        if trigger == "non_dimensional_temperature":

            taus_calculated = solve(eq, taus_calc)
            taus = float(taus_calculated[0].subs(m_star_calc,
                group_data["non_dimensional_mass_flow"]))
            group_data["non_dimensional_temperature"] = taus

        else:

            m_star_calculated = solve(eq, m_star_calc)
            m_star = float(m_star_calculated[0].subs(taus_calc,
                group_data["non_dimensional_temperature"]))
            group_data["non_dimensional_mass_flow"] = m_star

    elif group_data["cell_type"] == "JIVC":

        eq = Eq((5.83 * (m_star_calc ** 0.7648)) + 1.743, taus_calc)

        if taus == None:

            taus_calculated = solve(eq, taus_calc)
            taus = float(taus_calculated[0].subs(m_star_calc,
                group_data["non_dimensional_mass_flow"]))
            group_data["non_dimensional_temperature"] = taus

        else:

            m_star_calculated = solve(eq, m_star_calc)
            m_star = float(m_star_calculated[0].subs(taus_calc,
                group_data["non_dimensional_temperature"]))
            group_data["non_dimensional_mass_flow"] = m_star

    else:
        raise ValueError("Error, no cell type assigned")

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

def pdrop_sympy(group_data, trigger = "pressure_drop"):
    
    """
    Pressure drop equation, may be used to calculate pressure drop, Euler \
        number or velocity
    
    :param float density: Coolant density
    :param float Eu: Euler number
    :param float u: Jet velocity
    :param float delta_p: Pressure drop
    """
    
    density = group_data["fluid_properties"].density
    
    delta_p_calc, Eu_calc, u_calc = symbols("delta_p_calc Eu_calc u_calc")
    
    eq = Eq(Eu_calc * density * (u_calc ** 2), delta_p_calc)
    
    if trigger == "euler_number":

        Eu_calculated = solve(eq, Eu_calc)
        Eu = float(Eu_calculated[0].subs([(u_calc, group_data["jet_velocity"]),\
                                         (delta_p_calc, group_data["pressure_drop"])]))
        group_data["euler_number"] = Eu

    elif trigger == "pressure_drop":

        delta_p_calculated = solve(eq, delta_p_calc)
        delta_p = float(delta_p_calculated[0].subs([(Eu_calc, group_data["euler_number"]), \
                                                   (u_calc, group_data["jet_velocity"])]))
        group_data["pressure_drop"] = delta_p

    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs([(Eu_calc, group_data["euler_number"]), \
                                       (delta_p_calc, group_data["pressure_drop"])]))
        group_data["jet_velocity"] = u

def reynolds_sympy(group_data, trigger = "reynolds_number", D = 0.001):

    """
    Reynolds number equation, may be used to calculate Reynolds number or \
        velocity

    :param float density: Coolant density
    :param float mu: Coolant viscosity
    :param float D: Jet diameter
    :param float u: Jet velocity
    :param float Re: Reynolds number
    """

    density = group_data["fluid_properties"].density
    viscostiy = group_data["fluid_properties"].viscosity

    Re_calc, u_calc = symbols("Re_calc u_calc")

    eq = Eq((density * u_calc * D) / viscostiy, Re_calc)

    if trigger == "reynolds_number":
        Re_calculated = solve(eq, Re_calc)
        Re = float(Re_calculated[0].subs(u_calc, group_data["jet_velocity"]))

        group_data["reynolds_number"] = Re

    else:
        u_calculated = solve(eq, u_calc)
        u = float(u_calculated[0].subs(Re_calc, group_data["reynolds_number"]))

        group_data["jet_velocity"] = u

def euler_sympy(group_data, trigger = "euler_number"):

    """
    Euler number equation, may be used to calculate Euler number or Reynolds \
        number

    :param float Eu: Euler number
    :param float Re: Reynolds number
    """

    if group_data["cell_type"] == "DLH":

        Eu_calc, Re_calc = symbols("Eu_calc Re_calc")
        eq = Eq((1.486 * (Re_calc ** -0.0721)) - 0.001133, Eu_calc)

        if trigger == "euler_number":

            Eu_calculated = solve(eq, Eu_calc)
            Eu = float(Eu_calculated[0].subs(Re_calc,
                                             group_data["reynolds_number"]))
            group_data["euler_number"] = Eu

        else:

            Re_calculated = solve(eq, Re_calc)
            Re = float(Re_calculated[0].subs(Eu_calc,
                                             group_data["euler_number"]))
            group_data["reynolds_number"] = Re

    elif group_data["cell_type"] == "JIVC":

        Eu_calc, Re_calc = symbols("Eu_calc Re_calc")
        eq = Eq(11.85 * (Re_calc ** -0.2351), Eu_calc)

        if trigger == "euler_number":

            Eu_calculated = solve(eq, Eu_calc)
            Eu = float(Eu_calculated[0].subs(Re_calc,
                                             group_data["reynolds_number"]))
            group_data["euler_number"] = Eu

        else:

            Re_calculated = solve(eq, Re_calc)
            Re = float(Re_calculated[0].subs(Eu_calc, Eu))
            group_data["reynolds_number"] = Re

    else:
        raise ValueError("Error, no cell type assigned")


def coolant_temperature_sympy(taus, Tinput, metal_temperature):
    
    """
    Developed coolant temperature equation, may be used to calculate the \
        developed coolant temperature or mass flow
    
    :param float 
    """

    Tnew = Tinput + ((metal_temperature - Tinput) / taus)
    
    return Tnew

def mass_flow_next_row(group_data, cells, t = 0.001, ks = 340):

    if group_data["cell_type"] == "DLH":
        phi, theta, omega = [4.415, 0.6242, 0.7828]
    else:
        phi, theta, omega = [5.83, 0.7648, 1.743]

    # Ri = sqrt(coordinates[0]**2 + coordinates[2]**2)    
    Ri = cells.heat_cells[group_data["direction"]][0].R3
    HFi = group_data["HFf"](Ri) * group_data["cross_sectional_area"]

    A = group_data["mass_flow"]
    B = phi * (((group_data["fluid_properties"].specific_heat_capacity * t) / (ks * group_data["cross_sectional_area"]))**theta)
    C = omega
    D = HFi / (group_data["fluid_properties"].specific_heat_capacity * (group_data["max_temperature"] - group_data["input_temperature"]))

    Anew = lambda A, B, C, D: A - ((A - ((A**theta) * D * B) - (D * C))/(1 - (theta * (A **-(1-theta)) * D * B)))
    A_updated = Anew(A, B, C, D)

    if D == 0:
        mass_flow = 0.0
        ni = 0.0
    else:   
        while sqrt((A_updated - A)**2) >= 1e-10:
            A = A_updated

            A_updated = Anew(A, B, C, D)
        mass_flow = A_updated
        ni = group_data["mass_flow_total"] / mass_flow ## in this case, number of jets total
    
    group_data["mass_flow"] = mass_flow
    group_data["specified_number_of_hexagons"] = ni
    
    # return massflow, ni