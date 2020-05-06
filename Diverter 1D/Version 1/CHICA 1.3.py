#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
CHICA - Cooling and Heating Interaction and Corrosion Analysis
        Calculates the temperature evolution for a thermofluid flowing along an axis-symmetric channel
"""
__author__ = "Marcos Parro and Daniel Iglesias"
__copyright__ = "Copyright 2012, Marcos Parro"
__credits__ = ["Marcos Parro", "Daniel Iglesias"]
__license__ = "LGPL"
__version__ = "1.0.1"
__maintainer__ = "Daniel Iglesias"
__email__ = "Daniel.Iglesias@tokamakenergy.co.uk"
__status__ = "Experimental"

# import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as SI
import numpy
from Runner_v_1_1 import solver
from Setup import setup

# File definitions:
iPot = open("input_power.txt")
iSection = open("input_cross_section_vertical_actual.txt")

q_r = numpy.linspace(5, 50, 10)# kg/s, mass flow water, needs to vary between 5 kg/s and 50 kg/s
v_r_input = numpy.linspace(5, 50, 10) # m/s, Velocity helium, remember 50 m/s being mentione this is not fixed)

pressure_output = []
htc_0_output = []
Temp_output = []
Vel_output = []
tg_output = []
Ag_output = []
Mass_Flow_input = []
Velocity_input= []
Ag_output = []
noMicroChannels_output = []
reynolds_output = []
prandtl_output =[]
nusselt_output = []
T_metal_output = []
moodyf_output = []
Ma_output = []

# pot_0 = [] # z in W/m^2
# pot_1 = [] # DP in W/m^2
section_0 = [] # r in m
section_1 = [] # z in m

# for line in iPot :
#   data = line.split()
#   pot_0.append( float(data[0]))
#   pot_1.append( float(data[1])*1E6)
for line in iSection :
  data = line.split()
  section_0.append(float(data[0]))
  section_1.append(float(data[1]))

# Parameter definitions:
h = 200E-3 # m, thickness of the copper, this is not fixeed
 # Cooling type
#ref_type = "constant_film"
ref_type = "axial" #i.e. staight along the nodes, no out of plane
#ref_type = "spiral" #out of plane, no node to node transport by convection
constant_film_value = 2. #constant film not used
spiral_tube = "rectangular" # A=a*b #For spiral tube, shape defines hydraulic diameter
rectangular_base = 10E-2 # m, =a
g = 9.81 #Gravity
epsi = 0.00000015 # Value in m, surface roughness, only considered in Nu calc
input_pressure = 8E6 # input pressure in Pa
input_temperature = 373.15 # input temperature in K
input_rho =  SI('D', 'P', input_pressure, 'T', input_temperature, 'Helium') # Helium density in kg/m3
# internal_radius = 1.156 #internal radius of the plates
# input_power = 0

count_global = 0

for MassFlow in q_r:
    
    for y in v_r_input:
        
        input_power = [10E6 for i in range(len(section_0)-1)] # W/m2
        Mass_Flow_input.append(MassFlow) # add first mass flow term to input, needs to be in loop as changes with loop
        Velocity_input.append(y) # add first velocity term to input, needs to be in loop as changes with loop
        
        htc_0, Re, Pr, Nu, h_f, Ag, tg, dh, v_secc, T_ref, T_metal, P_secc, hf_tot, \
        Alist, deltaz, input_power, Ma = setup.initial_setup(section_0, input_temperature, input_pressure, \
        y, input_power, h, MassFlow, input_rho, epsi, section_1)
        
        # -------------------------------Start of first procedure------------------------------- #
        
        P_secc, v_secc, hf_tot, T_metal, T_ref = solver.initial(T_ref, section_0, \
        section_1, MassFlow, input_power, input_pressure, Ag, dh, epsi, deltaz, T_metal, Re, \
        Pr, h_f, hf_tot, P_secc, v_secc, Nu, htc_0, Alist, Ma)
        
        # # -------------------------------setup definitions to update pressure------------------------------- #
        
        #reset parameters
        # htc_0, Re, Pr, Nu, h_f, Ag, tg, dh, v_secc, T_ref, T_metal, \
        # hf_tot, moodyf, input_power = setup.looper_setup(section_0, input_temperature, \
        # P_secc, y, input_power, h, MassFlow, input_rho, section_1, epsi, deltaz)
        
        # P_secc, v_secc, hf_tot, T_metal, T_ref = solver.looper(T_ref, section_0, \
        # section_1, MassFlow, input_power, Ag, dh, epsi, deltaz, T_metal, Re, \
        # Pr, h_f, hf_tot, P_secc, v_secc, Nu, htc_0, Alist)
        
        pressure_output.append(P_secc[len(P_secc)-1])
        htc_0_output.append(htc_0[len(htc_0)-1])
        Temp_output.append(T_ref[len(T_ref)-1])
        Vel_output.append(v_secc[len(v_secc)-1])
        tg_output.append(tg)
        Ag_output.append(Ag[len(Ag)-1])
        Ma_output.append(Ma)
        # noMicroChannels_output.append(Ag/1E-6)
        reynolds_output.append(Re[len(Re)-1])
        prandtl_output.append(Pr[len(Pr)-1])
        nusselt_output.append(Nu[len(Nu)-1])
        T_metal_output.append(T_metal[len(T_metal)-1])
        # moodyf_output.append(moodyf[len(moodyf)-1])
        # count_global += 1
        # print(count_global)
    # if count_global == 30000:
    #     break

