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
        FC_input = []
        z = -2.82733
        
        if channel_type == "rectangle":
            
            thetap = 2*pi/n - (1*pi/180) # angle of the panel
            thetat = thetap / m # phi
            r2 = section_0[0] + m_min
            x = r2
            y = x * tan(thetat/2)
            R = sqrt(x ** 2 + y ** 2)
            thetam = 2 * asin((m_min/2)/R)
            
            thetap_old = 0
            theta = 0
            
            while sqrt((thetap_old - thetap) ** 2) > 1E-30:
                
                thetap_old = thetap
                # thetap = 2*pi/n - (1*pi/180) - (thetat + thetam) # angle of the panel
                thetat = thetap / m # phi
                r2 = section_0[0] + m_min
                x = r2
                y = x * tan(thetat/2)
                R = sqrt(x ** 2 + y ** 2)
                thetam = 2 * asin((m_min/2)/R)
                thetap = 2*pi/n - (1*pi/180) - (thetat + thetam) # angle of the panel
            
            theta = thetat - thetam
            rows = 1
            m_new = m
            
            # Minimum material thickness between channels check
            while theta < 0.0:
                
                rows += 1
                m_new = m / rows
                thetat = thetap/(m_new)
                theta = thetat - thetam
                
            a = 2 * R * sin(theta/2) # a
            b = (MassFlow/(n*m)) / (input_rho * VelInput * a) # b
            dh = 2 * a * b / (a + b)
            A2 = a * deltaz # A2, area of the pipe facing the power input
            phi = asin((a/2)/R)
            ARi = b/a
            rowsi = rows # set new variable rowsi to update in next loop
            m_row = m_new
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
            
            #     x1 = (R * cos(thetai + phi) + (i * ((b + m_min) * cos(thetai)))) # x1
            #     y1 = (R * sin(thetai + phi) + (i * ((b + m_min) * sin(thetai)))) # y1
            #     x2 = (R * cos(thetai + phi) + b * cos(thetai)) + (i * ((b + m_min) * cos(thetai))) # x2
            #     y2 = (R * sin(thetai + phi) + b * sin(thetai)) + (i * ((b + m_min) * sin(thetai))) # y2
            #     x3 = (R * cos(thetai - phi) + b * cos(thetai)) + (i * ((b + m_min) * cos(thetai))) # x4
            #     y3 = (R * sin(thetai - phi) + b * sin(thetai)) + (i * ((b + m_min) * sin(thetai))) # y4
            #     x4 = (R * cos(thetai - phi)) + (i * ((b + m_min) * cos(thetai))) # x3
            #     y4 = (R * sin(thetai - phi)) + (i * ((b + m_min) * sin(thetai))) # y3
            
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
            
            guide_wire = [[section_0[0], section_1[0], 0], [section_0[-1], section_1[-1], 0]]
            
            hmin = ((x2 - x1) * rows) + (m_min * (rows - 1))
            
        else: # circular geometry
            
            thetap = 2*pi/n - (1*pi/180) # angle of the panel
            thetat = thetap / m # phi
            dh = 2 * sqrt((MassFlow/(n*m))/(pi*input_rho*VelInput))
            R = section_0[0] + m_min + dh/2
            thetam = 2 * asin(((dh/2) + (0.5 * m_min))/((dh/2) + R))
            
            thetap_old = 0
            theta = 0
            
            while sqrt((thetap_old - thetap) ** 2) > 1E-30:
                
                thetap_old = thetap
                thetat = thetap / m # phi
                dh = 2 * sqrt((MassFlow/(n*m))/(pi*input_rho*VelInput))
                R = section_0[0] + m_min + dh/2
                thetam = 2 * asin(((dh/2) + (0.5 * m_min))/((dh/2) + R))
                thetap = 2*pi/n - (1*pi/180) - (thetat + thetam) # angle of the panel
            
            print(thetap*180/pi)
            
            rows = 1
            theta = thetat - thetam
            m_new = m
            
            while theta < 0:
                
                rows += 1
                m_new = m / rows
                thetat = thetap/(m_new)
                theta = thetat - thetam
            
            m = m_new * rows
            A1 = pi * ((dh/2) ** 2)
            A2 = dh * deltaz
            m_row = m/rows
            
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
        
        for i in guide_wire:
            for j, y in enumerate(i):
                i[j] = y * 1000
        
        # scale power input
        for i in range(len(input_power)):
            
            A0 = pi * deltaz * (section_0[i] + section_0[i+1]) # area of a frustum between nodes
            P0 = A0 * input_power[i] # power in [W] for the defined frustum
            Pp = P0 / n # divide by number of plates
            Ppp = Pp / m # divide by number of pipes per plate
            input_power[i] = Ppp # updates power input list
            
        f = open("Example_divertor_complex//coolant_geometry.asc", "w")

        for line in FC_input:
            for i in line:
                f.write(str(i) + "\t")
            f.write("\n")
        
        f.close()
        
        f = open("Example_divertor_complex//coolant_guide.asc", "w")

        for line in guide_wire:
            for i in line:
                f.write(str(i) + "\t")
            f.write("\n")
        
        f.close()
        
        return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, m, thetap
            
    
    
    def discrete_pipes_toroidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n, m, channel_type, m_min, AR):
        
        print(m)
        
        A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, thetap = [0 for i in range(13)]
        FC_input = []
        
        h = section_1[-1] - section_1[0]
        a = (h - ((m+1)*m_min)) / m
        b = (MassFlow / (m * n)) / (input_rho * a * VelInput)
        theta = atan((h)/(section_0[-1] - section_0[0]))
        z = 0.0
        A1 = a * b
        
        x1 = (section_0[0] + m_min + (a + m_min)/tan(theta))
        y1 = - 0.5*a
        x2 = x1
        y2 = - y1
        x3 = x1 + b
        y3 = y2
        x4 = x3
        y4 = y1
        
        FC_input.append([x1, y1, z])
        FC_input.append([x2, y2, z])
        FC_input.append([x3, y3, z])
        FC_input.append([x4, y4, z])
        FC_input.append([x1, y1, z])
        
        for i in FC_input:
            for j, y in enumerate(i):
                i[j] = y * 1000
        
        f = open("Example_divertor_complex//coolant_geometry.asc", "w")

        for line in FC_input:
            for i in line:
                f.write(str(i) + "\t")
            f.write("\n")
        
        return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, m, thetap
