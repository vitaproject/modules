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
    
    def discrete_pipes_poloidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n, m, channel_type, m_min, AR):
        
        deltaz = sqrt((section_1[1] - section_1[0])**2 + (section_0[1] - section_0[0])**2) #Length in [m], lenght of pipe section being considered
        thetat = 2*pi/(n*m) # phi
        thetap = 2*pi/n - 1 # angle of the panel
        FC_input = []
        z = -2.82733
        
        if channel_type == "rectangle":
            
            r2 = section_0[0] + m_min
            x = r2
            y = x * tan(thetat/2)
            R = sqrt(x ** 2 + y ** 2)
            thetam = 2 * asin((m_min/2)/R)
            rows = 1
            
            theta = thetat - thetam
            
            # Minimum material thickness between channels check
            while theta < 0.0:
                
                thetat = 2*pi/(n*m//rows)
                theta = thetat - thetam
                channels = m/(rows + 1)
                rows += 1
            
            a = 2 * R * sin(theta/2) # a
            b = (MassFlow/(n*m)) / (input_rho * VelInput * a) # b
            dh = 2 * a * b / (a + b)
            A2 = a * deltaz # A2, area of the pipe facing the power input
            phi = asin((a/2)/R)
            ARi = b/a
            rowsi = rows # set new variable rowsi to update in next loop
            m_row = m/rows # number of channels per row
            m_new = m
            b_new = b
            
            # Channel aspect ratio check
            while ARi > AR:
                rowsi += 1
                b_new = b * (rows / rowsi)
                ARi = b_new/a
                m_new += m_row
            
            rows = rowsi # updated number of rows
            b = b_new # updated channel width
            m = m_new # updated number of channels
            A1 = a * b # A1, area of the pipe with fluid flowing in a single row
            
            thetai = 0.0 # for one output, symmetric about the x axis
            #thetai = i * 2*theta  #for multiple outputs
            
            # for i in range(0, rows):
            
                # x1 = (R * cos(thetai + phi) + (i * ((b + m_min) * cos(thetai)))) # x1
                # y1 = (R * sin(thetai + phi) + (i * ((b + m_min) * sin(thetai)))) # y1
                # x2 = (R * cos(thetai + phi) + b * cos(thetai)) + (i * ((b + m_min) * cos(thetai))) # x2
                # y2 = (R * sin(thetai + phi) + b * sin(thetai)) + (i * ((b + m_min) * sin(thetai))) # y2
                # x3 = (R * cos(thetai - phi) + b * cos(thetai)) + (i * ((b + m_min) * cos(thetai))) # x4
                # y3 = (R * sin(thetai - phi) + b * sin(thetai)) + (i * ((b + m_min) * sin(thetai))) # y4
                # x4 = (R * cos(thetai - phi)) + (i * ((b + m_min) * cos(thetai))) # x3
                # y4 = (R * sin(thetai - phi)) + (i * ((b + m_min) * sin(thetai))) # y3
            
            x1 = (R * cos(thetai + phi)) # x1
            y1 = (R * sin(thetai + phi)) # y1
            x2 = (R * cos(thetai + phi) + b * cos(thetai)) # x2
            y2 = (R * sin(thetai + phi) + b * sin(thetai)) # y2
            x3 = (R * cos(thetai - phi) + b * cos(thetai)) # x4
            y3 = (R * sin(thetai - phi) + b * sin(thetai)) # y4
            x4 = (R * cos(thetai - phi))# x3
            y4 = (R * sin(thetai - phi))# y3
                
            FC_input.append([x1, y1, z])
            FC_input.append([x2, y2, z])
            FC_input.append([x3, y3, z])
            FC_input.append([x4, y4, z])
            FC_input.append([x1, y1, z])
            
            hmin = ((x2 - x1) * rows) + (m_min * (rows - 1))
            
        else: # circular geometry
            
            dh = 2 * sqrt((MassFlow/(n*m))/(pi*input_rho*VelInput))
            R = section_0[0] + m_min + dh/2
            thetam = 2 * asin(((dh/2) + (0.5 * m_min))/((dh/2) + R))
            rows = 1
            
            theta = thetat - thetam
            
            while theta < 0:
                
                thetat = 2*pi/(n*m//rows)
                theta = thetat - thetam
                channels = m/(rows + 1)
                rows += 1
            
            m_row = m/rows
            A1 = pi * ((dh/2) ** 2)
            A2 = dh * deltaz
            
            FC_input.append([R + (dh/2), 0, 0]) # centre point
            FC_input.append([R + (dh), 0, 0]) # radius point
            
            # values need to be set for the output, needs to be cleaned up
            phi = 0 # not required for circular output
            a = 0 # not required for circular output
            b = 0 # not required for circular output
            hmin = (dh * rows) + (m_min * (rows - 1))
        
        for i in FC_input:
            for j, y in enumerate(i):
                i[j] = y * 1000
        
        for i in range(len(section_1)-1):
            input_power[i] = input_power[i] * A2 # power scaled for incident area or pipe
        
        f = open("Example_divertor_complex//coolant_geometry.asc", "w")

        for line in FC_input:
            for i in line:
                f.write(str(i) + "\t")
            f.write("\n")
        
        f.close()
        
        return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, m
            