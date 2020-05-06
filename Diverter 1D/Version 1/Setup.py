from Coolant_Geometry import coolant_geometry
from Non_Dimensional import non_dimensional

class setup:
    
    def __init__(self):
        
        return
    
    def initial_setup(section_0, input_temperature, input_pressure, y, \
                     input_power, h, MassFlow, input_rho, \
                     epsi, section_1):
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
        v_secc.append(y)
        Re = []
        Pr = []
        Nu = []
        Ma = []
        # dimension terms
        Ag, tg, dh, deltaz, Alist, input_power = coolant_geometry.a_r_poloidal_flow(MassFlow, \
                                                input_rho, y, section_0, section_1, input_power, h)
        # input definitions
        # non-dimensional terms at node point, don't use v_s from this calc, just shows that v_secc[0] = y
        nu, re, pr, h_1, v_s = non_dimensional.nusselt(MassFlow, T_ref[0], Ag[0], input_pressure, \
                                dh, epsi, deltaz[0])
        h_f.append(h_1)
        Re.append(re)
        Pr.append(pr)
        Nu.append(nu)
        # compute the first film coefficient
        htc_0.append(non_dimensional.film_coeff(dh, input_pressure, T_ref[0], nu) \
                     * Alist[0])
        return htc_0, Re, Pr, Nu, h_f, Ag, tg, dh, v_secc, T_ref, T_metal, \
            P_secc, hf_tot, Alist, deltaz, input_power, Ma

    def looper_setup(section_0, input_temperature, P_secc, y, input_power, h, MassFlow, \
                     input_rho, section_1, epsi, deltaz):
        
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
        # dimension terms
        # Ag, tg, dh, deltaz, Alist, input_power = coolant_geometry.a_r_poloidal_flow(MassFlow, \
        #                                         input_rho, y, section_0, section_1, input_power, h)
        
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
        
        return htc_0, Re, Pr, Nu, h_f, Ag, tg, dh, v_secc, T_ref, T_metal, hf_tot, moodyf, input_power
