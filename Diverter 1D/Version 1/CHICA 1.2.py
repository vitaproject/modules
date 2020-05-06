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

from math import *
import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import PropsSI as SI
import numpy

# File definitions:
iPot = open("input_power.txt")
iSection = open("input_cross_section_vertical_ideal.txt")

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
input_rho_list = []

pot_0 = [] # z in W/m^2
pot_1 = [] # DP in W/m^2
section_0 = [] # r in m
section_1 = [] # z in m

for line in iPot :
  data = line.split()
  pot_0.append( float(data[0]))
  pot_1.append( float(data[1])*1E6)
for line in iSection :
  data = line.split()
  section_0.append( float(data[0]))
  section_1.append( float(data[1]))

count_global = 0

for MassFlow in q_r:
    
    for j, y in enumerate(v_r_input):
        
        Mass_Flow_input.append(MassFlow)
        Velocity_input.append(y)
        
        
        # reset parameter
        # pot_1[0] = 6413394.48 # W/m2
        pot_1[0] = 8E6
        
        # Parameter definitions:
        h = 200E-3 # m, thickness of the copper, this is not fixeed
         # Cooling type
        #ref_type = "constant_film"
        ref_type = "axial" #i.e. staight along the nodes, no out of plane
        #ref_type = "spiral" #out of plane, no node to node transport by convection
        constant_film_value = 2. #constant film not used
        spiral_tube = "rectangular" # A=a*b #For spiral tube, shape defines hydraulic diameter
        rectangular_base = 10E-2 # m, =a
        deltaz = section_1[1] #Length in [m], lenght of pipe section being considered
        g = 9.81 #Gravity
        epsi = 0.00000015 # Value in m, surface roughness, only considered in Nu calc
        hf_tot = 0 #Total pressure head loss from friction and linear losses
        input_pressure = 8E6 # input pressure in Pa
        input_temperature = 373.15 # input temperature in K
        input_rho =  SI('D', 'P', input_pressure, 'T', input_temperature, 'Helium') # Helium density in kg/m3
        internal_radius = 1.156 #internal radius of the plates
        
        input_rho_list.append(input_rho)
        
        # dimension terms
        step1 = ((MassFlow)/(input_rho*y*pi)) + (internal_radius**2)
        a = 2 * internal_radius
        tg = (-a + sqrt((a**2)-(4*((internal_radius**2)-step1)))) / 2 #thickness of gap
        # tg = 0.01 # tg fixed
        Ag = pi*((internal_radius + tg)**2) - pi*(internal_radius**2)
        
        # Lists of data in files:
        htc_0 = [] # Film transfer coefficient in W/K
        h_f = [] # Friction coefficient
        v_secc = [] # Speed array
        P_secc = [] # Pressure array
        hf_tot = [] # Charge loss array
        
        # read files and store data
        
        # Lists of temperatures:
        T_ref = [0.0 for i in range( len(section_0) ) ]
        T_metal = [0.0 for i in range( len(section_0) ) ]
        P_secc = [0.0 for i in range( len(section_0) ) ]
        hf_tot = [0.0 for i in range( len(section_0) ) ]
        
        # print(len(section_0))
        
        T_ref[0] = input_temperature # Initial temperature, K
        P_secc[0] = input_pressure # Initial pressure value, Pa
        v_secc.append(y)
        
        ####### first guess definitions
        
        # Coolant section area (m**2)
        def a_r () :
          if ref_type == "spiral" :
            return a*b # this will be the toroidal version, need to think
          elif ref_type == "axial" :
            return Ag # axial version i.e. poloidal channel
        
        def dh() :
          # Hydraulic diameter
            # Dh = 4 Area / wet perimeter
          # Case of spiral rectangular tube:
            # Dh = (2ab)/(a+b)
          if ref_type == "spiral" :
            return ((2*a*b)/(a+b))
          elif ref_type == "axial" :
            return 2*tg
          else :
            print ("WARNING: NOT IMPLEMENTED!!!!\n")
            return 1.
        
        input_re = input_rho*y*dh()/SI('V', 'T', input_temperature, 'P', input_pressure, 'helium')
        input_pr = SI('V', 'T', input_temperature, 'P', input_pressure, 'helium') \
              * SI('C', 'T', input_temperature, 'P', input_pressure, 'helium') \
              / SI('L', 'T', input_temperature, 'P', input_pressure, 'helium')
        
          # Correlation:
        def film_coef(index) :
          if ref_type == "constant_film" :
            return constant_film_value
          else :
            # Nu*K/Dh
            film = nusselt(index)*SI('L', 'T', T_ref[index], 'P', input_pressure, 'helium') \
                    / dh()
            
            htc_0.append(film * (2.0*pi*internal_radius*section_1[1]))
            return film
        
        moodyf = []
        Nu = []
        Re = []
        Re.append(input_re)
        Pr = []
        Pr.append(input_pr)
        # adjust power density of first point to outer surface:
        
        pot_1[0] *=  section_0[0]/(section_0[0]+h)
        
        def nusselt(index) :
          # Nu=f(Re,Pr) ,, Depends on correlation
            # Reynold's number, Re = rho*v_s*Dh/mu
            # v_s = q_r/(rho*A) -> Bulk velocity (away from boundary layer)
          v_s = MassFlow/(SI('D', 'T', T_ref[index], 'P', input_pressure, 'helium')*a_r())
          re = SI('D', 'T', T_ref[index], 'P', input_pressure, 'helium') \
                * v_s * dh() / SI('V', 'T', T_ref[index], 'P', input_pressure, 'helium')
          Re.append(re)
          # fricc = 0.25 / (log10 ((epsi /(3.71*dh())) + (5.74 / (re)**0.9)))**2
          a = 2/(log(10))
          d = log(10)*re/5.2
          b = epsi/(dh()*3.7)
          s = b*d + log(d)
          r = s**(s/(s+1))
          m = b*d + log(d/r)
          p = log(r/m)
          DLA = p*(m/(m+1))
          DCFA = DLA*(1+((p/2)/(((m+1)**2)+((p/3)*(2*m-1)))))
          RHS = a*((log(d/r)) + DCFA)
          fricc = (1/RHS)**2
          
          h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(g) * dh())
          h_f.append(h_1)
          
          
            # smooth tube:
          # friction=(1.82*log10(re)-1.64)**(-2.)
          prandtl = SI('V', 'T', T_ref[index], 'P', input_pressure, 'helium') \
              * SI('C', 'T', T_ref[index], 'P', input_pressure, 'helium') \
              / SI('L', 'T', T_ref[index], 'P', input_pressure, 'helium')
          Pr.append(prandtl)
          # oVars.write(str(v_s)+"\t"+str(re)+"\t")
          # oVarsSI.write(str(v_s)+"\t"+str(re)+"\t")
          v_secc.append(v_s)
          
            # Petukhov:
          return ( (fricc/8.)*re*prandtl/(1.07+12.7*sqrt(fricc/8.) \
                  * (prandtl**(2./3.)-1)) )
        
        input_htc_0 = nusselt(0)*SI('L', 'T', input_temperature, 'P', input_pressure, 'helium') \
                    / dh()
        
        htc_0.append(input_htc_0)
        
        # -------------------------------Start of first procedure------------------------------- #
          # Computing the temperature of coolant:
        
        #1st loop
        for i in range( len(section_1)-1 ) :
          
          T_ref[i+1] = T_ref[i]+(2.0*pi*internal_radius*section_1[1]*section_0[1]) \
                      / (MassFlow*SI('C', 'T', T_ref[i], 'P', input_pressure, 'helium')) \
                        *pot_1[0]  
          
          # Computing the temperature of interface:
          T_metal[i+1] = T_ref[i+1] + pot_1[0]/film_coef(i+1)
        
        hf_tot[0] = h_f[0]
        
        ##2nd loop
        for x in range(len(h_f)-1) :
          hf_tot[x+1] = float(h_f[x]) + hf_tot[x]
        hf_bar = hf_tot[len(h_f)-1] * float(g) * SI('D', 'T', T_ref[x], 'P', input_pressure, 'helium')\
          / 100000
        
        h_tot = hf_bar
        
        ##3rd loop
        for i in range( len(section_0)-1) :
            
            P_secc[i+1] = P_secc[i] + (0.5*SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium') \
            * (v_secc[i]**2 - v_secc[i+1]**2))\
            + (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g))\
            * float(section_1[i]-section_1[i+1])\
            - (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)*float(h_f[i]))
        
        # -------------------------------setup definitions to take updated pressure------------------------------- #
        ####### iterative definitions
        
        def film_coef_iterable(index) :
          if ref_type == "constant_film" :
            return constant_film_value
          else :
            # Nu*K/Dh
            film = nusselt_iterable(index)*SI('L', 'T', T_ref[index], 'P', P_secc[index], 'helium') \
                    / dh()
            
            htc_0[index] = film * (2.0*pi*internal_radius*section_1[1])
            return film
        
        def nusselt_iterable(index) :
          # Nu=f(Re,Pr) ,, Depends on correlation
            # Reynold's number, Re = rho*v_s*Dh/mu
            # v_s = q_r/(rho*A) -> Bulk velocity (away from boundary layer)
            
          v_s = MassFlow/(SI('D', 'T', T_ref[index], 'P', P_secc[index], 'helium')*a_r())
          re = SI('D', 'T', T_ref[index], 'P', P_secc[index], 'helium') \
                * v_s * dh() / SI('V', 'T', T_ref[index], 'P', P_secc[index], 'helium')
          Re[index] = re
          
          # fricc = 0.25 / (log10 ((epsi /(3.71*dh())) + (5.74 / (re)**0.9)))**2
          a = 2/(log(10))
          d = log(10)*re/5.2
          b = epsi/(dh()*3.7)
          s = b*d + log(d)
          r = s**(s/(s+1))
          m = b*d + log(d/r)
          p = log(r/m)
          DLA = p*(m/(m+1))
          DCFA = DLA*(1+((p/2)/(((m+1)**2)+((p/3)*(2*m-1)))))
          RHS = a*((log(d/r)) + DCFA)
          fricc = (1/RHS)**2
          
          h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(g) * dh())
          h_f[index] = h_1
          
          moodyf.append(fricc)
          
            # smooth tube:
          # friction=(1.82*log10(re)-1.64)**(-2.)
          prandtl = SI('V', 'T', T_ref[index], 'P', P_secc[index], 'helium') \
              * SI('C', 'T', T_ref[index], 'P', P_secc[index], 'helium') \
              / SI('L', 'T', T_ref[index], 'P', P_secc[index], 'helium')
          Pr[index] = prandtl
          v_secc[index] = v_s
            # Petukhov: (page 53 in PHD paper)
          
          Nu.append((fricc/8.)*re*prandtl/(1.07+12.7*sqrt(fricc/8.) \
                  * (prandtl**(2./3.)-1)) )
          
          return ( (fricc/8.)*re*prandtl/(1.07+12.7*sqrt(fricc/8.) \
                  * (prandtl**(2./3.)-1)) )
        
        # -------------------------------Start of second procedure------------------------------- #
        ## repeat above steps with pressure estimates until pressure and temperature converge
        
        P_secc_new = [8E6 for i in range(len(section_0))]
        
        
        
        count = 0
        # while sqrt((1-(P_secc_new[len(P_secc_new)-1]/P_secc[len(P_secc)-1]))**2) > 1E-30:
        while count <= 100:
            
            #1st loop
            for i in range( len(section_1)-1 ) :
                
              T_ref[i+1] = T_ref[i]+(2.0*pi*internal_radius*section_1[1]*section_0[1]) \
                          / (MassFlow*SI('C', 'T', T_ref[i], 'P', P_secc[i], 'helium')) \
                            *pot_1[0]
              
              # Computing the temperature of interface:
              T_metal[i+1] = T_ref[i+1] + pot_1[0]/film_coef_iterable(i+1)
            
            hf_tot[0] = h_f[0]
            
            ##2nd loop
            for x in range(len(h_f)-1) :
              hf_tot[x+1] = float(h_f[x]) + hf_tot[x]
            hf_bar = hf_tot[len(h_f)-1] * float(g) * SI('D', 'T', T_ref[x], 'P', P_secc[x], 'helium')\
              / 100000
            
            h_tot = hf_bar
            
            ##3rd loop
            for i in range( len(section_0)-1) :
                
                P_secc_new[i+1] = P_secc[i] + (0.5*SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium') \
                * (v_secc[i]**2 - v_secc[i+1]**2))\
                + (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g))\
                * float(section_1[i]-section_1[i+1])\
                - (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)*float(h_f[+1]))
        
            count += 1
            
            
            if sqrt((1-(P_secc_new[len(P_secc_new)-2]/P_secc[len(P_secc)-2]))**2) < 1E-8:
                break
            else:
                P_secc = P_secc_new
        
            # print(pot_1)
        
        pressure_output.append(P_secc[len(P_secc)-1])
        htc_0_output.append(htc_0[len(htc_0)-1])
        Temp_output.append(T_ref[len(T_ref)-1])
        Vel_output.append(v_secc[len(v_secc)-1])
        tg_output.append(tg)
        Ag_output.append(Ag)
        noMicroChannels_output.append(Ag/1E-6)
        reynolds_output.append(Re[len(Re)-1])
        prandtl_output.append(Pr[len(Pr)-1])
        nusselt_output.append(Nu[len(Nu)-1])
        T_metal_output.append(T_metal[len(T_metal)-1])
        moodyf_output.append(moodyf[len(moodyf)-1])
        count_global += 1
        # print(count_global)
    if count_global == 30000:
        break

