#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
CHICA - Cooling and Heating Interaction and Corrosion Analysis
        Calculates the temperature evolution for a thermofluid flowing along an axis-symmetric channel
"""
__author__ = "Marcos Parro, Jack Taylor and Daniel Iglesias"
__copyright__ = "Copyright 2012, Marcos Parro"
__credits__ = ["Marcos Parro", "Jack Taylor", "Daniel Iglesias"]
__license__ = "LGPL"
__version__ = "1.0.1"
__maintainer__ = "Jack Taylor"
__email__ = "Jack.Taylor@tokamakenergy.co.uk"
__status__ = "Experimental"

# import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as SI
import numpy
#import csv
from runner import solver
from setup import setup

# File definitions:
iPot = open("input_power.txt")
iSection = open("Divertor_lower_LFS_sink.asc")

q_r = numpy.linspace(5, 50, 1)# kg/s, mass flow water, needs to vary between 5 kg/s and 50 kg/s
v_r_input = numpy.linspace(5, 50, 1) # m/s, Velocity helium, remember 50 m/s being mentione this is not fixed)

pressure_output = []
htc_0_output = []
Temp_output = []
Vel_output = []
Mass_Flow_input = []
Velocity_input= []
reynolds_output = []
prandtl_output =[]
nusselt_output = []
T_metal_output = []
moodyf_output = []
Ma_output = []
a_output = []
b_output = []
A1_output = []
A2_output = []

section_0i = [] # r in m
section_1i = [] # z in m

# csv files
# with open("input_cross_section_vertical_actual.csv") as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     count = 0
#     for row in csv_reader:
#         if count == 0:
#             section_0.append(float(row[0][3:]))
#             section_1.append(float(row[1]))
#             count += 1
#         else:
#             section_0.append(float(row[0]))
#             section_1.append(float(row[1]))

#txt and asc files
count = 0
for line in iSection :
    if count < 2:
        data = line.split()
        section_0i.append(float(data[0])/1000)
        section_1i.append(float(data[2])/1000)
        count += 1

section_0 = numpy.linspace(section_0i[0], section_0i[1], 10)
section_1 = numpy.linspace(section_1i[0], section_1i[1], 10)

# Parameter definitions:
h = 200E-3 # m, thickness of the copper, this is not fixeed
epsi = 0.00000015 # Value in m, surface roughness, only considered in Nu calc
input_pressure = 8E6 # input pressure in Pa
input_temperature = 373.15 # input temperature in K
input_rho =  SI('D', 'P', input_pressure, 'T', input_temperature, 'Helium') # Helium density in kg/m3
n = 12 # total number of plates
m = 2000 # number of pipes per plate
channel_type = "rectangle" # options: rectangle, circle
m_min = 0.5E-3 # minimum material between channels, [m]
AR = 4

count_global = 0

for MassFlow in q_r:
    
    for VelocityInput in v_r_input:
        
        input_power = [10E6 for i in range(len(section_0)-1)] # W/m2
        Mass_Flow_input.append(MassFlow) # add first mass flow term to input, needs to be in loop as changes with loop
        Velocity_input.append(VelocityInput) # add first velocity term to input, needs to be in loop as changes with loop
        
        htc_0, Re, Pr, Nu, h_f, dh, v_secc, T_ref, T_metal, P_secc, hf_tot, \
        deltaz, input_power, Ma, A1, A2, phi, a, b, FC_input, rows, hmin = \
        setup.initial_setup(section_0, input_temperature, input_pressure, \
        VelocityInput, input_power, h, MassFlow, input_rho, epsi, section_1, \
        n, m, channel_type, m_min, AR)
        
        # -------------------------------Start of first procedure------------------------------- #
        
        # P_secc, v_secc, hf_tot, T_metal, T_ref = solver.initial(T_ref, section_0, \
        # section_1, MassFlow/(n*m), input_power, input_pressure, A1, dh, epsi, deltaz, T_metal, Re, \
        # Pr, h_f, hf_tot, P_secc, v_secc, Nu, htc_0, A2, Ma)
        
        # -------------------------------setup definitions to update pressure------------------- #
        
        #reset parameters
    #     htc_0, Re, Pr, Nu, h_f, v_secc, T_ref, T_metal, \
    #     hf_tot, moodyf, Ma = setup.looper_setup(section_0, input_temperature, \
    #     P_secc, VelocityInput, input_power, h, MassFlow, input_rho, section_1, epsi, deltaz, Ag, dh)
        
    #     P_secc, v_secc, hf_tot, T_metal, T_ref = solver.looper(T_ref, section_0, \
    #     section_1, MassFlow, input_power, Ag, dh, epsi, deltaz, T_metal, Re, \
    #     Pr, h_f, hf_tot, P_secc, v_secc, Nu, htc_0, Alist, Ma)
        
        # pressure_output.append(P_secc[len(P_secc)-1])
        # htc_0_output.append(htc_0[len(htc_0)-1])
        # Temp_output.append(T_ref[len(T_ref)-1])
        # Vel_output.append(v_secc[len(v_secc)-1])
        # Ma_output.append(Ma)
        # reynolds_output.append(Re[len(Re)-1])
        # prandtl_output.append(Pr[len(Pr)-1])
        # nusselt_output.append(Nu[len(Nu)-1])
        # T_metal_output.append(T_metal[len(T_metal)-1])
        # a_output.append(a)
        # b_output.append(b)
        # A1_output.append(A1)
        # A2_output.append(A2)
        # moodyf_output.append(moodyf[len(moodyf)-1])
        # count_global += 1
        # print(count_global)
    # if count_global == 30000:
    #     break

# f = open("coolant_geometry.txt", "w")

# for line in FC_input:
#     f.write(line + "\n")
    
# f.close()
