
# import matplotlib.pyplot as plt
from os import path
from CoolProp.CoolProp import PropsSI as SI
import numpy
# import csv
from chica.utility import get_example_data_path
from chica.runner import initial
from chica.geometry_setup import initial_setup


#################################
# Input definitions:
isection_path = get_example_data_path("divertor_complex/Divertor_lower_LFS_sink.asc")

# Sensitivity study parameters:
q_r = numpy.linspace(5, 50, 1)  # Mass flow rate, kg/s
v_r_input = numpy.linspace(5, 50, 1)  # Flow velocity, m/s

# Overall parameter definitions:
# Tokamak:
n_input = 12  # total number of plates
m_input = 50  # number of pipes per plate
input_power = [10E6 for i in range(9)]  # W/m^2
# Cooling channel
epsi = 0.00000015  # Surface Roughness, m
input_pressure = 8E6  # input pressure, Pa
input_temperature = 373.15  # input temperature, K
# Channel geometry
# channel_type = "rectangle"
channel_type = "circle" 
orientation = "poloidal"
# orientation = "toroidal"
m_min = 1E-4  # Minimum material between channels, m
AR = 1  # Aspect ratio, width / depth

################################
# Initial calculated parameters
input_rho = SI('D', 'P', input_pressure, 'T', input_temperature, "helium")  # kg/m3

pressure_output = []
htc_0_output = []
temp_output = []
vel_output = []
mass_flow_input = []
velocity_input = []
reynolds_output = []
prandtl_output = []
nusselt_output = []
T_metal_output = []
moodyf_output = []
Ma_output = []
a_output = []
b_output = []
A1_output = []
A2_output = []
hmin_output = []
rows_output = []
thetap_output = []

section_0i = []  # r in m
section_1i = []  # z in m

# txt and asc files
with open(isection_path) as isection:
    count = 0
    for line in isection:
        if count < 2:
            data = line.split()
            section_0i.append(float(data[0]) / 1000)
            section_1i.append(float(data[2]) / 1000)
            count += 1

section_0 = numpy.flip(numpy.linspace(section_0i[0], section_0i[1], 10))
section_1 = numpy.flip(numpy.linspace(section_1i[0], section_1i[1], 10))

count_global = 0

for massflow in q_r:

    for velocityinput in v_r_input:
        mass_flow_input.append(massflow)  # add first mass flow term to input, needs to be in loop as changes with loop
        velocity_input.append(
            velocityinput)  # add first velocity term to input, needs to be in loop as changes with loop

        htc_0, Re, Pr, Nu, h_f, dh, v_secc, T_ref, T_metal, P_secc, hf_tot, \
        deltaz, input_power, Ma, A1, A2, phi, a, b, rows, hmin, m_row, m, thetap = \
            initial_setup(section_0, input_temperature, input_pressure, \
                          velocityinput, input_power, massflow, input_rho, epsi, section_1, \
                          n_input, m_input, channel_type, m_min, AR, orientation)

        # -------------------------------Start of first procedure------------------------------- #

        P_secc, v_secc, hf_tot, T_metal, T_ref = initial(T_ref, section_0, \
                                                         section_1, massflow / (n_input * m), input_power,
                                                         input_pressure, A1, dh, epsi, \
                                                         deltaz, T_metal, Re, Pr, h_f, hf_tot, P_secc, v_secc,
                                                         Nu, htc_0, A2, Ma)

        # -------------------------------setup definitions to update pressure------------------- #

        # reset parameters
        #     htc_0, Re, Pr, Nu, h_f, v_secc, T_ref, T_metal, \
        #     hf_tot, moodyf, Ma = setup.looper_setup(section_0, input_temperature, \
        #     P_secc, VelocityInput, input_power, h, MassFlow, input_rho, section_1, epsi, deltaz, Ag, dh)

        #     P_secc, v_secc, hf_tot, T_metal, T_ref = solver.looper(T_ref, section_0, \
        #     section_1, MassFlow, input_power, Ag, dh, epsi, deltaz, T_metal, Re, \
        #     Pr, h_f, hf_tot, P_secc, v_secc, Nu, htc_0, Alist, Ma)

        pressure_output.append(P_secc[len(P_secc) - 1])
        htc_0_output.append(htc_0[len(htc_0) - 1])
        temp_output.append(T_ref[len(T_ref) - 1])
        vel_output.append(v_secc[len(v_secc) - 1])
        Ma_output.append(Ma)
        reynolds_output.append(Re[len(Re) - 1])
        prandtl_output.append(Pr[len(Pr) - 1])
        nusselt_output.append(Nu[len(Nu) - 1])
        T_metal_output.append(T_metal[len(T_metal) - 1])
        a_output.append(a)
        b_output.append(b)
        A1_output.append(A1)
        A2_output.append(A2)
        hmin_output.append(hmin)
        rows_output.append(rows)
        thetap_output.append(thetap)
        # moodyf_output.append(moodyf[len(moodyf)-1])
        # count_global += 1
        # print(count_global)
    # if count_global == 30000:
    #     break