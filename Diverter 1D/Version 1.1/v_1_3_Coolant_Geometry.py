# computes the tangential area at each node point between two adjacent nodes
# area of the gap Ag [m**2]
# thickness of the gap tg [m]

from math import *
import sympy

class coolant_geometry:
    
    def __init__(self):
        
        return 

    def annulus_poloidal_flow(MassFlow, input_rho, y, section_0, section_1, input_power, h):
        
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
    
    def discrete_pipes_poloidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n, m):
        
        # A2, C2, AS2, xi, xii, deltazi, a, b, ASG2, A2i, AS2i, ASG2i, ai, bi, C2i, ASGP2, phi, phii = \
        #     sympy.symbols("A2, C2, AS2, xi, xii, deltazi, a, b, ASG2, A2i, AS2i, ASG2i, ai, bi, C2i, ASGP2i, phi, phii")
        
        deltaz = sqrt((section_1[1] - section_1[0])**2 + (section_0[1] - section_0[0])**2) #Length in [m], lenght of pipe section being considered
        theta = 2*pi/(n*m) # phi
        a = 2 * section_0[0] * sin(theta/2) # a
        b = (MassFlow/(n*m)) / (input_rho * VelInput * a) # b
        A1 = a * b # A1, area of the pipe with fluid flowing
        A2 = a * deltaz # A2, area of the pipe facing the power input
        dh = 2 * a * b / (a + b)
        phi = asin((a/2)/section_0[0])
        
        FC_input = []
        z = -2.82733
        
        for i in range(n*m):
            
            thetai = i * 2*theta
                
            x1 = (section_0[0] * cos(thetai + phi)) + h * cos(thetai) # x1
            y1 = (section_0[0] * sin(thetai + phi)) + h * sin(thetai) # y1
            x2 = (section_0[0] * cos(thetai + phi) + b * cos(thetai)) + h * cos(thetai) # x2
            y2 = (section_0[0] * sin(thetai + phi) + b * sin(thetai)) + h * sin(thetai) # y2
            x3 = (section_0[0] * cos(thetai - phi) + b * cos(thetai)) + h * cos(thetai) # x4
            y3 = (section_0[0] * sin(thetai - phi) + b * sin(thetai)) + h * sin(thetai) # y4
            x4 = (section_0[0] * cos(thetai - phi)) + h * cos(thetai) # x3
            y4 = (section_0[0] * sin(thetai - phi)) + h * sin(thetai) # y3
            
            FC_input.append([x1, y1, z])
            FC_input.append([x2, y2, z])
            FC_input.append([x3, y3, z])
            FC_input.append([x4, y4, z])
            FC_input.append([x1, y1, z])
        
        # first loop to calculate a and b, second to calulate copper area as pipes progress
        
        for i in FC_input:
            for j, y in enumerate(i):
                i[j] = y * 1000
        
        for i in range(len(section_1)-1):
            input_power[i] = input_power[i] * A2 * \
                (section_0[i]/(section_0[i] + h)) # scaling for copper thickness
        
        return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input
            