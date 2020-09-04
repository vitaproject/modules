
from chica.coolant_geometry import toroidal_flow_rectangle, \
                                   poloidal_flow_rectangle, poloidal_flow_circle, \
                                   Qdot_toroidal_flow
from chica.non_dimensional import nusselt, film_coeff
from math import pi

def initial_setup(section_0, input_temperature, input_pressure, VelInput, \
                  input_power, MassFlow, input_rho, epsi, section_1, n, m, \
                  channel_type, m_min, AR, orientation):
    
    if orientation == "poloidal" and channel_type == "rectangle":
        A1, A2, dh, input_power, a, b, rows, \
            m_row, m, freeCAD_input_angle = poloidal_flow_rectangle(\
            MassFlow, input_rho, VelInput, section_0, section_1, input_power, \
            n, m, channel_type, m_min, AR)
    
    if orientation == "poloidal" and channel_type == "circle":
        A1, A2, deltaz, dh, input_power, phi, a, b, rows, hmin, \
            m_row, m, thetap = poloidal_flow_circle(\
            MassFlow, input_rho, VelInput, section_0, section_1, input_power, \
            n, m, channel_type, m_min, AR)
    
    if orientation == "toroidal" and channel_type == "rectangle":
        A1, A2, dh, input_power, a, b, rows, m_row, m \
            = toroidal_flow_rectangle(\
            MassFlow, input_rho, VelInput, section_0, section_1, input_power, \
            n, m, channel_type, m_min, AR)
    
    # print(rows, m, m_row, freeCAD_input_angle*180/pi)
    print(rows, m, m_row)
    
    htc_0, h_f, v_secc, T_ref, T_metal, P_secc, hf_tot, Re, Pr, Nu, Ma = \
        flow_properties(section_0, input_temperature, input_pressure, VelInput)
    
    # input definitions
    # non-dimensional terms at node point, don't use v_s from this calc, just shows that v_secc[0] = y
    nu, re, pr, h_1, v_s = nusselt(MassFlow/(n*m), T_ref[0], A1, input_pressure, \
                            dh, epsi, deltaz) # A1i refers to discrete pipe
      
    # The below should be assigned in non_dimensional.nusselt() to clean it up
    h_f.append(h_1)
    Re.append(re)
    Pr.append(pr)
    Nu.append(nu)
    # compute the first film coefficient
    htc_0.append(film_coeff(dh, input_pressure, T_ref[0], nu) \
                  * A2)

    return htc_0, Re, Pr, Nu, h_f, dh, v_secc, T_ref, T_metal, \
        P_secc, hf_tot, deltaz, input_power, Ma, A1, A2, phi, a, b, rows, hmin, m_row, m, thetap

def looper_setup(section_0, input_temperature, P_secc, y, input_power, h, MassFlow, \
                 input_rho, section_1, epsi, deltaz, Ag, dh):

    # setup parameters for first run

    # Lists of data in files and initial parameters:
    htc_0 = [] # Film transfer coefficient in W/m2/K
    h_f = [] # head loss in pipe due to friction, m
    v_secc = [] # velocity at each node
    T_ref = [0.0 for i in range(len(section_0))] # temperature of free stream
    T_metal = [0.0 for i in range(len(section_0))] # temperature of metal (in this case copper)
    hf_tot = [0.0 for i in range(len(section_0))] # total head loss at each node
    T_ref[0] = input_temperature # Initial temperature, K
    v_secc.append(y)
    moodyf = []
    Re = []
    Pr = []
    Nu = []
    Ma = []
    # input definitions
    # non-dimensional terms at node point, don't use v_s from this calc, just shows that v_secc[0] = y
    nu, re, pr, h_1, v_s = nusselt(MassFlow, T_ref[0], Ag[0], P_secc[0], \
                            dh, epsi, deltaz[0])

    h_f.append(h_1)
    Re.append(re)
    Pr.append(pr)
    Nu.append(nu)

    # compute the first film coefficient
    htc_0.append(film_coeff(dh, P_secc[0], T_ref[0], nu))

    return htc_0, Re, Pr, Nu, h_f, v_secc, T_ref, T_metal, hf_tot, moodyf, Ma

def flow_properties(section_0, input_temperature, input_pressure, VelInput):
    """
    List of flow properties and how they develop along a single pipe

    :param list htc_0: Film transfer coefficient, W/m^2/K
    :param list h_f: Pressure head loss in pipe, m
    :param list v_secc: Velocity along pipe length, m/s
    :param list T_ref: Temperature along pipe length, K
    :param list T_metal: Metal temperature along pipe length, K
    :param list P_secc: Total pressure along pipe, Pa
    :param list hf_tot: Accumulated head loss along pipe, m
    :param list Re: Reynolds number
    :param list Pr: Prandtl number
    :param list Nu: Nusselt number
    :param list Ma: Mach number
    :rtype float:
    """
    # Lists of data in files and initial parameters:
    htc_0 = [] # Film transfer coefficient in W/m2/K
    h_f = [] # Friction coefficient
    v_secc = [] # velocity at each node
    T_ref = [0.0 for i in range(len(section_0))] # temperature of free stream
    T_metal = [0.0 for i in range(len(section_0))] # temperature of metal (in this case copper)
    P_secc = [0.0 for i in range(len(section_0))] # pressure at each node
    hf_tot = [0.0 for i in range(len(section_0))] # total head loss at each node
    T_ref[0] = input_temperature # Initial temperature, K
    P_secc[0] = input_pressure # Initial pressure value, Pa
    v_secc.append(VelInput)
    Re = []
    Pr = []
    Nu = []
    Ma = []
    
    return htc_0, h_f, v_secc, T_ref, T_metal, P_secc, hf_tot, Re, Pr, Nu, Ma