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
    
    def discrete_pipes_poloidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n):
        
        A2, C2, AS2, xi, xii, deltazi, a, b, ASG2, A2i, AS2i, ASG2i, ai, bi = \
            sympy.symbols("A2, C2, AS2, xi, xii, deltazi, a, b, ASG2, A2i, AS2i, ASG2i, ai, bi")
        
        # n is the number of segments, i.e. number of pipes required
        
        equation1 = pi * deltazi * (xi + xii) # A2
        equation2 = A2 / n # AS2
        equation3 = AS2 - (C2 * AS2) # ASG2
        equation4 = ASG2/deltazi # a
        equation5 = MassFlow/(input_rho * VelInput * a * n) # b
        
        eq1 = sympy.Eq(equation1, A2)
        eq2 = sympy.Eq(equation2, AS2)
        eq3 = sympy.Eq(equation3, ASG2)
        eq4 = sympy.Eq(equation4, a)
        eq5 = sympy.Eq(equation5, b)
        
        output_A2 = []
        output_AS2 = []
        output_ASG2 = []
        output_a = []
        output_b = []
        output_dh = []
        output_A1i = []
        output_A1 = []
        deltaz = []
        
        for i in range( len(section_1)-1):
            deltaz.append(sqrt((section_1[i+1] - section_1[i])**2 + (section_0[i+1] - section_0[i])**2)) #Length in [m], lenght of pipe section being considered
        
        for i in range(len(section_1)-1):
            A2i = sympy.solve(eq1.subs({deltazi:deltaz[i], xi:section_0[i], xii:section_0[i+1]}), A2)
            AS2i = sympy.solve(eq2.subs(A2, float(A2i[0])), AS2)
            ASG2i = sympy.solve(eq3.subs({AS2:float(AS2i[0]), C2:0.0}), ASG2) # C2 this is a user input
            ai = sympy.solve(eq4.subs({ASG2:float(ASG2i[0]), deltazi:deltaz[i]}), a)
            bi = sympy.solve(eq5.subs({a:float(ai[0])}))
            dhi = 2*float(ai[0])*float(bi[0]) / (float(ai[0]) + float(bi[0]))
            A1i = float(ai[0]) * float(bi[0])
            A1 = n * A1i
            
            input_power[i] = input_power[i] * float(ASG2i[0]) * \
                (section_0[i+1]/(section_0[i+1] + h)) # scaling for copper thickness
            
            output_A2.append(float(A2i[0]))
            output_AS2.append(float(AS2i[0]))
            output_ASG2.append(float(ASG2i[0]))
            output_a.append(float(ai[0]))
            output_b.append(float(bi[0]))
            output_dh.append(dhi)
            output_A1i.append(A1i)
            output_A1.append(A1)
        
        return output_A2, output_AS2, output_ASG2, output_a, output_b, output_dh, \
            output_A1i, output_A1, input_power, deltaz
        
        
    
    