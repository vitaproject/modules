# computes the tangential area at each node point between two adjacent nodes
# area of the gap Ag [m**2]
# thickness of the gap tg [m]

from math import *
import numpy


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


def poloidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n, m, channel_type, m_min, AR):

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

    sectionCoord = numpy.array([[section_0[0], section_1[0]], [section_0[-1], section_1[-1]]])

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


def toroidal_flow(MassFlow, input_rho, VelInput, section_0, section_1, input_power, h, n, m, channel_type, m_min,
                  AR):

    A0 = MassFlow / (input_rho * VelInput)
    sectionCoord = numpy.array([[section_0[0], section_1[0]], [section_0[-1], section_1[-1]]])
    d = sqrt((section_1[-1] - section_1[0]) ** 2 + (section_0[-1] - section_0[0]) ** 2)
    a = (d - (m_min * (m + 1))) / m
    b = A0 / (a * m)

    rows = 1
    b_orig = b
    m_orig = m

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

    m_act, theta = coordinates("toroidal", m, rows, d, a, sectionCoord, m_min, b)

    return A1, A2, deltaz, dh, input_power, phi, a, b, FC_input, rows, hmin, m_row, m, thetap


def add_rows(rows, b, b_orig, m, m_orig, a, AR):
    while a / b < AR:
        rows += 1
        b = b_orig * (1 / rows)
        m += m_orig  # total number of pipes
    return rows, b, m


def adjust_AR(a, b, AR):
    A1 = a * b
    b = sqrt(A1 / AR)
    a = AR * b
    return a, b


def coordinates(orientation, m, rows, d, a, sectionCoord, m_min, b):
    FC_input = []

    if orientation == "poloidal":
        return
    else:
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

    FC_input.append([x1, y1, z])
    FC_input.append([x2, y2, z])
    FC_input.append([x4, y4, z])
    FC_input.append([x3, y3, z])
    FC_input.append([x1, y1, z])

    for i in FC_input:
        for j, y in enumerate(i):
            i[j] = y * 1000

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

    return m_act, theta
