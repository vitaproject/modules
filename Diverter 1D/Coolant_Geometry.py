# computes the tangential area at each node point between two adjacent nodes
# area of the gap Ag [m**2]
# thickness of the gap tg [m]

from math import *

class coolant_geometry:
    
    def __init__(self):
        
        return 

    def a_r_poloidal_flow(MassFlow, input_rho, y, section_0, section_1, input_power, h):
        
        # dimension terms
        deltaz = []
        Ag = [] # tangential area of the gap
        Alist = [] # parallel area of the gap
        
        for i in range(len(section_1)-1):
            A = pi * (section_0[i] + section_0[i+1]) * sqrt(((section_1[i]-section_1[i+1])**2) \
                + ((section_0[i]-section_0[i+1])**2)) # for a frustum
            # A = 2 * pi() * section_0[i] * (section_1[i+1] - section_1[i]) # for a cylinder
            Alist.append(A)
            input_power[i] = input_power[i] * A * \
                (section_0[i+1]/(section_0[i+1] + h)) # scaling for copper thickness
        
        # thickness of gap
        step1 = ((MassFlow)/(input_rho*y*pi)) + (section_0[0]**2)
        a = 2 * section_0[0]
        tg = (-a + sqrt((a**2)-(4*((section_0[0]**2)-step1)))) / 2 
        
        # tg = 0.01 # tg fixed
        
        # deltaz needs to be calculated for each section
        for i in range( len(section_1)-1):
            deltaz.append(sqrt((section_1[i+1] - section_1[i])**2 + (section_0[i+1] - section_0[i])**2)) #Length in [m], lenght of pipe section being considered
        
        # compute for each node past the first node
        for i, x in enumerate(section_1):
            Ag.append(pi*((section_0[i] + tg)**2) - pi*(section_0[i]**2))
        
        # hydraulic diameter
        dh = 2 * tg # for annular cylinder
        # dh = diameter of the tubes for perfectly circular tubes
        
        return Ag, tg, dh, deltaz, Alist, input_power
    
    def a_r_toroidal_flow(MassFlow, input_rho, y, internal_radius):
        
        # Need to update tg and Ag calcs
        
        # dimension terms
        step1 = ((MassFlow)/(input_rho*y*pi)) + (internal_radius**2)
        a = 2 * internal_radius
        tg = (-a + sqrt((a**2)-(4*((internal_radius**2)-step1)))) / 2 #thickness of gap
        # tg = 0.01 # tg fixed
        Ag = pi*((internal_radius + tg)**2) - pi*(internal_radius**2)
        
        return Ag, tg