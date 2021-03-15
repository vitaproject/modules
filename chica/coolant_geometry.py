from math import sin, cos, tan, asin, atan, sqrt, pi
from numpy import array

def poloidal_flow_rectangle(MassFlow, input_rho, VelInput, section_0, section_1, \
                            input_power, n, m, channel_type, m_min, AR):
    """
    Creates pipe profile geometry file for importing in to freeCAD.
    Creates pipe geometry information for computing flow properties.
    Assumes rectanglar pipe geometry and poloidal flow.
     
    :rtype float:
    """
    A0 = MassFlow / (input_rho * VelInput)
    sectionCoord = array([[section_0[0], section_1[0]], [section_0[-1], section_1[-1]]])
    d = sqrt((section_1[-1] - section_1[0]) ** 2 + (section_0[-1] - section_0[0]) ** 2)
    # ^ Length in [m], lenght of pipe section being considered
    
    thetap = 2*pi/n - (1*pi/180) # angle of the panel
    thetat = thetap / m          # incorrect
    x = section_0[0] + m_min
    y = x * tan(thetat/2)
    R = sqrt(x ** 2 + y ** 2)
    thetam = 2 * atan((m_min/2)/R) # swept angle of minimum material thickness
    theta = (thetap - thetam * (m+1))/m # swept angle of a single pipe
    a = 2 * R * sin(theta/2) # single pipe width
    b = (MassFlow/(n*m)) / (input_rho * VelInput * a) # single pipe depth
    
    freeCAD_input_angle = thetap - 2 * thetam - theta
    
    rows = 1
    b_orig = b
    m_orig = m

    a, b, rows, m = AR_check(a, b, m, b_orig, m_orig, rows, AR)
    freeCAD_input_angle = freeCAD_angle_back_calc(a, m, R, thetap, thetat)
    FC_input, guide_wire = coordinates("poloidal", m, rows, None, a, sectionCoord, m_min, b, R, theta)
    A1, A2, dh, m_row = pipe_geometry(a, b, d, m, rows)
    write_files(FC_input, guide_wire, rows, m)
    input_power = scale_power(input_power, A2)
    
    return A1, A2, dh, input_power, a, b, rows, m_row, m, freeCAD_input_angle

def toroidal_flow_rectangle(MassFlow, input_rho, VelInput, section_0, section_1, \
                            input_power, n, m, channel_type, m_min, AR):
    """
    Creates pipe profile geometry file for importing in to freeCAD.
    Creates pipe geometry information for computing flow properties.
    Assumes rectanglar pipe geometry and toroidal flow.
     
    :rtype float:
    """
    A0 = MassFlow / (input_rho * VelInput)
    sectionCoord = array([[section_0[0], section_1[0]], [section_0[-1], section_1[-1]]])
    d = sqrt((section_1[-1] - section_1[0]) ** 2 + (section_0[-1] - section_0[0]) ** 2)
    a = (d - (m_min * (m + 1))) / m
    b = A0 / (a * m * n)

    rows = 1
    b_orig = b
    m_orig = m

    a, b, rows, m = AR_check(a, b, m, b_orig, m_orig, rows, AR)
    FC_input, guide_wire = coordinates("toroidal", m, rows, d, a, sectionCoord, m_min, b, None, None)
    A1, A2, dh, m_row = pipe_geometry(a, b, d, m, rows)
    write_files(FC_input, guide_wire, rows, m)
    input_power = scale_power(input_power, A2)
    
    return A1, A2, dh, input_power, a, b, rows, m_row, m

def add_rows(rows, b, b_orig, m, m_orig, a, AR):
    """
    First check case for AR, this ensures a > b when checking user AR against \
    a/b by adding rows, thus reducing b. This is done as, if b > a and so \
    a/b < AR, to adhere to AR and keep the A0, a would need to \
    increase, therefore failing the m_min requirement. Multiple rows must be \
    added and b reduced for each, to keep A0 constant whilst adhereing to the \
    m_min requirement.
     
    :rtype float:
    """
    while a / b < AR:
        rows += 1
        b = b_orig * (1 / rows)
        m += m_orig  # total number of pipes
    return rows, b, m


def adjust_AR(a, b, AR):
    """
    Second check case for aspect ratio, ensures final aspect ratio matches \
    user input
    
    :param float A1: Cross sectional area of a single pipe 
    :rtype float:
    """
    A1 = a * b
    b = sqrt(A1 / AR)
    a = AR * b
    return a, b

def coordinates(orientation, m, rows, d, a, sectionCoord, m_min, b, R, theta):
    """
    Assigns pipe profile geometry points to freeCAD input file

    :param list FC_input: list of points to be used in freeCAD for creating \
                          pipe profile
    :param list guide_wire: list of points to be used in freeCAD for creating \
                            guide wire to extrude pipe profile
    :rtype float:
    """
    FC_input = []

    if orientation == "poloidal":
        z = -2.82733
        x1 = R * cos(theta/2) # x1
        y1 = R * sin(theta/2) # y1
        x2 = x1 + b
        y2 = y1
        x3 = x1
        y3 = -y1
        x4 = x2
        y4 = -y1
        guide_wire = [[sectionCoord[0][0], sectionCoord[0][1], 0], [sectionCoord[1][0], sectionCoord[1][1], 0]]
        
    if orientation == "toroidal":
        z = 0
        m_row = m / rows
        m_act = (d - ((m_row) * a)) / (m_row + 1)
        theta = atan((sectionCoord[0][1] - sectionCoord[1][1]) / (sectionCoord[0][0] - sectionCoord[1][0]))
        x1 = sectionCoord[0][0] + m_act * cos(theta) + m_min * cos(theta - pi / 2)
        y1 = sectionCoord[0][1] + m_act * sin(theta) - m_min * cos(theta - pi / 2)
        x2 = x1 + b * cos(pi / 2 - theta)
        y2 = y1 - b * sin(pi / 2 - theta)
        x3 = x1 + a * cos(theta)
        y3 = y1 + a * sin(theta)
        x4 = x2 + a * cos(theta)
        y4 = y2 + a * sin(theta)
        guide_wire = [[sectionCoord[0][0], sectionCoord[0][1], 0], [sectionCoord[1][0], sectionCoord[1][1], 0]]

    FC_input.append([x1, y1, z]) # vertex 1
    FC_input.append([x2, y2, z]) # vertex 2
    FC_input.append([x4, y4, z]) # vertex 3
    FC_input.append([x3, y3, z]) # vertex 4
    FC_input.append([x1, y1, z]) # vertex 1, closes shape in freeCAD
    
    return FC_input, guide_wire

def write_files(FC_input, guide_wire, rows, m):
    
    """
    Assigns pipe profile geometry points to freeCAD input file

    :param list FC_input: list of points to be used in freeCAD for creating \
                          pipe profile
    :param list guide_wire: list of points to be used in freeCAD for creating \
                            guide wire to extrude pipe profile
    :rtype float:
    """
    
    for i in FC_input:
        for j, y in enumerate(i):
            i[j] = y * 1000 # convert to mm

    for i in guide_wire:
        for j, y in enumerate(i):
            i[j] = y * 1000 # convert to mm

    f = open("coolant_geometry.asc", "w")

    for line in FC_input:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")

    f.close()

    f = open("coolant_guide.asc", "w")

    for line in guide_wire:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")

    f.close()
    
    f = open("important_info.asc", "w")

    f.write(str(m/rows) + "\n" + str(rows))

    f.close()

def AR_check(a, b, m, b_orig, m_orig, rows, AR):
    
    """
    Conducts an aspect ratio check on the current pipe geometry, will return 
    pipe geometry that conforms to the user input aspect ratio 

    :rtype float:
    """
    
    while a / b != AR:

        if AR > 1 and a > b and a / b > AR:
            print("Gate 1")
            a, b = adjust_AR(a, b, AR)

        if AR > 1 and a > b and a / b < AR:
            print("Gate 2")
            rows, b, m = add_rows(rows, b, b_orig, m, m_orig, a, AR)
            a, b = adjust_AR(a, b, AR)

        if AR > 1 and a < b and a / b < AR:
            print("Gate 3")
            rows, b, m = add_rows(rows, b, b_orig, m, m_orig, a, AR)
            a, b = adjust_AR(a, b, AR)

        if AR < 1 and a > b and a / b > AR:
            print("Gate 4")
            a, b = adjust_AR(a, b, AR)

        if AR < 1 and a < b and a / b > AR:
            print("Gate 5")
            a, b = adjust_AR(a, b, AR)

        if AR < 1 and a < b and a / b < AR:
            print("Gate 6")
            rows, b, m = add_rows(rows, b, b_orig, m, m_orig, a, AR)
            a, b = adjust_AR(a, b, AR)

        if AR == 1 and a > b:
            print("Gate 7")
            a, b = adjust_AR(a, b, AR)

        if AR == 1 and a < b:
            print("Gate 8")
            rows, b, m = add_rows(rows, b, b_orig, m, m_orig, a, AR)
            a, b = adjust_AR(a, b, AR)

        Error = a / b
        if sqrt((AR - Error) ** 2.0) < 1E-5:
            break

        else:
            print("SOMETHING HAS GONE TERRIBLY WRONG")
            break
    
    return a, b, rows, m

def freeCAD_angle_back_calc(a, m, R, thetap, thetat):
    
    """
    Identifies the swept angle between the centre point of pipes at opposite
    ends of single plate. Required input for freeCAD geometry production
    for poloidal pipes.

    :rtype float:
    """
    
    theta = 2 * asin((a/2)/R)
    thetam = (thetat - theta) * (m/(m+1))
    freeCAD_input_angle = thetap - 2 * thetam - theta
        
    return freeCAD_input_angle

def pipe_geometry(a, b, d, m, rows):
    
    """
    Assigns pipe profile geometry points to freeCAD input file

    :rtype float:
    """
    
    A1 = a * b
    A2 = a * d
    dh = 2 * a * b / (a + b)
    m_row = m/rows
    
    return A1, A2, dh, m_row

def scale_power(input_power, A2):
    
    """
    Scales the power input to an individual pipe based on the pipe area 
    incident to the plasma
    
    :rtype float:
    """
    
    # # scale power input
    for i in range(len(input_power)):
        P0 = A2 * input_power[i]  # assigns power in [W]
        input_power[i] = P0  # updates power input list
    
    return input_power

#     # A0 = pi * deltaz * (section_0[i] + section_0[i+1]) # area of a frustum between nodes
#     # P0 = A0 * input_power[i] # power in [W] for the defined frustum
#     # Pp = P0 / n # divide by number of plates
#     # Ppp = Pp / m # divide by number of pipes per plate
#     # input_power[i] = Ppp # updates power input list

""" below is work in progress"""

def poloidal_flow_circle(MassFlow, input_rho, VelInput, section_0, section_1, \
                            input_power, n, m, channel_type, m_min, AR):
    
    deltaz = sqrt((section_1[1] - section_1[0])**2 + (section_0[1] - section_0[0])**2) 
    # ^ Length in [m], lenght of pipe section being considered
    FC_input = []
    
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

    # scale power input
    for i in range(len(input_power)):

        # A0 = pi * deltaz * (section_0[i] + section_0[i+1]) # area of a frustum between nodes
        # P0 = A0 * input_power[i] # power in [W] for the defined frustum
        # Pp = P0 / n # divide by number of plates
        # Ppp = Pp / m # divide by number of pipes per plate
        # input_power[i] = Ppp # updates power input list

        P0 = A2 * input_power[i]  # assigns power in [W]
        input_power[i] = P0  # updates power input list

    f = open("example_data/divertor_complex/coolant_geometry.asc", "w")

    for line in FC_input:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")

    f.close()

    return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, m, thetap
    
def Qdot_toroidal_flow(SxW, SyH, m, section_0, section_1, AR, rows, MassFlow, input_rho, VelInput, fixed_parameter,
                       m_min):

    # geometry parameters are driving, so flow parameters calculated after
    # either mass flow or velicty fixed
    # AR = W/H, will likely be provided as H/W

    d = sqrt((section_1[-1] - section_1[0]) ** 2 + (section_0[-1] - section_0[0]) ** 2)
    Sx = d / m
    a = Sx / SxW
    b = a * (1 / AR)

    sectionCoord = array([[section_0[0], section_1[0]], [section_0[-1], section_1[-1]]])

    coordinates("toroidal", m, rows, d, a, sectionCoord, m_min, b)

    A = a * b
    A0 = A * m * rows

    if fixed_parameter is "MassFlow":
        VelInput = MassFlow / (A0 * input_rho)

    elif fixed_parameter is "Velocity":
        MassFlow = VelInput * input_rho * A0
    else:
        print("SOMETHING HAS GONE TERRIBLY WRONG")

    return