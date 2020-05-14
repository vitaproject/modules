from coolant_geometry import coolant_geometry
from non_dimensional import non_dimensional

class setup:
    
    def __init__(self):
        
        return
    
    def initial_setup(section_0, input_temperature, input_pressure, VelInput, \
                     input_power, h, MassFlow, input_rho, \
                     epsi, section_1, n, m, channel_type, m_min, AR):
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
        # dimension terms
        A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, \
            m_row, m, thetap = coolant_geometry.discrete_pipes_poloidal_flow(\
            MassFlow, input_rho, VelInput, section_0, section_1, input_power, \
            h, n, m, channel_type, m_min, AR)
        # input definitions
        # non-dimensional terms at node point, don't use v_s from this calc, just shows that v_secc[0] = y
        nu, re, pr, h_1, v_s = non_dimensional.nusselt(MassFlow/(n*m), T_ref[0], A1, input_pressure, \
                                dh, epsi, deltaz) # A1i refers to discrete pipe
        h_f.append(h_1)
        Re.append(re)
        Pr.append(pr)
        Nu.append(nu)
        # compute the first film coefficient
        htc_0.append(non_dimensional.film_coeff(dh, input_pressure, T_ref[0], nu) \
                      * A2)
        
        return htc_0, Re, Pr, Nu, h_f, dh, v_secc, T_ref, T_metal, \
            P_secc, hf_tot, deltaz, input_power, Ma, A1, A2, phi, a, b, FC_input, rows, hmin, m_row, m, thetap

    def looper_setup(section_0, input_temperature, P_secc, y, input_power, h, MassFlow, \
                     input_rho, section_1, epsi, deltaz, Ag, dh):
        
        # setup parameters for first run
        
        # Lists of data in files and initial parameters:
        htc_0 = [] # Film transfer coefficient in W/m2/K
        h_f = [] # Friction coefficient
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
        nu, re, pr, h_1, v_s = non_dimensional.nusselt(MassFlow, T_ref[0], Ag[0], P_secc[0], \
                                dh, epsi, deltaz[0])
        
        h_f.append(h_1)
        Re.append(re)
        Pr.append(pr)
        Nu.append(nu)
        
        # compute the first film coefficient
        htc_0.append(non_dimensional.film_coeff(dh, P_secc[0], T_ref[0], nu))
        
        return htc_0, Re, Pr, Nu, h_f, v_secc, T_ref, T_metal, hf_tot, moodyf, Ma
