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

# File definitions:
iPot = open("input_power.txt")
iKMat = open("input_conductivity_Cu.txt")
iCpRef = open("input_specific_heat_water.txt")
iKRef = open("input_conductivity_water.txt")
iDensRef = open("input_density_water.txt")
iPrandtl = open("input_Prandtl_number.txt")
iVisco = open("input_dynamic_viscosity_water.txt")
iSection = open("input_cross_section.txt")

oTemp = open("output_temperatures.txt", "w")
oVars = open("output_variables.txt", "w")
oTempSI = open("output_temperaturesSI.txt", "w")
oVarsSI = open("output_variablesSI.txt", "w")
oPresSI = open("output_pressureSI.txt", "w")
oVelcrit = open("output_VMiller.txt", "w")
oRet = open("output_Vret.txt", "w")
oCrit = open("output_CriticalFlux.txt", "w")

oVars.write("# z(cm)\t Velocity(m/s)\t Re\t h(W/(cm^2*C)) \t thickness (cm) \n")
oVarsSI.write("# z(m)\t Velocity(m/s)\t Re\t h(W/(m^2*C)) \t thickness (m) \n")
oPresSI.write("# z(m)\t Pressure (Pa)\n")
oVelcrit.write("# R(m)\t Critical velocity (m/s)\t Curved critical velocity curvada \n")
oRet.write("# Req(m)\t Velocity (m/s)\t h(W/m^2*ºC)\t Temperatura(ºC) \n")
oCrit.write("z(m) \t q* (n/u)\t Critical flux (W/m^2) \n")

# Parameter definitions:
h = 3.5E-3 # m, thickness
q_r = 30 # kg/s, mass flow water
v_r = 4 # m/s, Velocity water
thickness_min = 0.4 # cm
pipe_init_radius = 3. # cm
 # Cooling type
#ref_type = "constant_film"
ref_type = "axial"
#ref_type = "spiral"
constant_film_value = 2.
thickness_type = "fernando_1"
spiral_tube = "rectangular" # A=a*b
rectangular_base = 10E-2 # m, =a
deltaz = 0.05 #Lenght of each section in m
g = 9.81 #Gravity
epsi = 0.00000015 # Value in m
hf_tot = 0
h_2 = 0 # Cumulative friction parameter for pressure calculation
rho = 1000 # Water density in kg/m3

#prints a and b
if spiral_tube == "rectangular" :
  a = rectangular_base
  b = q_r*1E-3/(v_r*a)
  print ("a="+str(a)+"m")
  print ("b="+str(b)+"m")
else : # if is squared
  a = b = math.sqrt(q_r*1E-3/v_r)
#a_r = a*b # only spiral!

# Lists of data in files:
kMat_0 = [] # temp in C
kMat_1 = [] # K in W(m C)^-1
cpRef_0 = [] # temp in C
cpRef_1 = [] # cp in J(kg C)^-1
kRef_0 = [] # temp in C 
kRef_1 = [] # K in W(m C)^-1
pot_0 = [] # z in mm
pot_1 = [] # DP in W/cm^2
densRef_0 = [] # temp in C
densRef_1 = [] # density in kg(m)^-3
prandtl_0 = [] # temp in C
prandtl_1 = [] # Prandtl number (adimensional)
visco_0 = [] # temp in C
visco_1 = [] # dynamic viscosity in kg(m s)^-1
section_0 = [] # z in cm
section_1 = [] # r in cm
htc_0 = [] # Film transfer coefficient in W/m2/K
h_f = [] # Friction coefficient
v_secc = [] # Speed array
P_secc = [] # Pressure array
hf_tot = [] # Charge loss array
secc_ret = [] # Section of the return water channel

# read files and store data
for line in iKMat :
  data = line.split()
  kMat_0.append( float(data[0]) )
  kMat_1.append( float(data[1]) )
for line in iCpRef :
  data = line.split()
  cpRef_0.append( float(data[0]) )
  cpRef_1.append( float(data[1]) )
for line in iKRef :
  data = line.split()
  kRef_0.append( float(data[0]) )
  kRef_1.append( float(data[1]) )
for line in iPot :
  data = line.split()
  pot_0.append( float(data[0])*1E-2 )
  pot_1.append( float(data[1])*1E4 )
for line in iDensRef :
  data = line.split()
  densRef_0.append( float(data[0]) )
  densRef_1.append( float(data[1]) )
for line in iPrandtl :
  data = line.split()
  prandtl_0.append( float(data[0]) )
  prandtl_1.append( float(data[1]) )
for line in iVisco :
  data = line.split()
  visco_0.append( float(data[0]) )
  visco_1.append( float(data[1]) )
for line in iSection :
  data = line.split()
  section_0.append( float(data[0])*1E-2 )
  section_1.append( float(data[1])*1E-2 )

# Lists of temperatures:
T_ref = [0.0 for i in range( len(section_0) ) ]
T_metal = [0.0 for i in range( len(section_0) ) ]
P_secc = [0.0 for i in range( len(section_0) ) ]
hf_tot = [0.0 for i in range( len(section_0) ) ]

T_ref[0] = 31 # Initial temperature, C
P_secc[0] = 3*1E5 # Initial pressure value, Pa
#V = q_r/(interpolation(densRef_0,densRef_1,T_ref[49])* pi *(0.18**2 - (section_1[49] + thickness(49))**2))

oTemp.write(str(section_0[0]))
oTemp.write("\t")
oTemp.write(str(T_ref[0]))
oTemp.write("\n")
oPresSI.write(str(section_1[0]))
oPresSI.write("\t")
oPresSI.write(str(P_secc[0]))
oPresSI.write("\n")

# Definition of functions:
 # Interpolation:
def interpolation(list_0, list_1, value_0) :
  for j in range( len(list_0)-1 ) :
    if list_0[j+1] > value_0 :
#      print list_0[j],list_0[j+1],value_0,list_1[j],list_1[j+1],\
#            list_1[j] + (list_1[j+1]-list_1[j])/ \
#               (list_0[j+1]-list_0[j])*(value_0-list_0[j])
      # y2 = y0 + (y1-y0)/(x1-x0)*x2
      return float( list_1[j] + (list_1[j+1]-list_1[j])/ \
               (list_0[j+1]-list_0[j])*(value_0-list_0[j]) )
  # Over maximum value...
  return list_1[j+1]

 # Gap thickness (m)
def thickness(index_in) :
  if thickness_type == "fernando_1" :
    if section_0[index_in] <= 10.8333E-2 : return 2.266*1E-2
    elif section_0[index_in] <= 68.815E-2 : 
      return ( 2.266 + (section_0[index_in] - 10.8333E-2)*(0.9734-2.266)/(68.815E-2-10.8333E-2) )*1E-2
    elif section_0[index_in] <= 123.462E-2 : 
      return ( 0.9734+ (section_0[index_in] - 68.815E-2)*(.6-0.9734)/(123.462E-2-68.815E-2) )*1E-2
    else : return 0.6*1E-2
  else :
    return max( pipe_init_radius-section_1[index_in], thickness_min )

# Coolant section area (m**2)
def a_r (index_in) :
  if ref_type == "spiral" :
    return a*b
  elif ref_type == "axial" :
    return pi*( (section_1[index_in]+thickness(index_in) )**2 -\
		(section_1[index_in] )**2 )

 # Correlation:
def film_coef(index) :
  if ref_type == "constant_film" :
    return constant_film_value
  else :
    oVars.write(str(section_0[index]*1E2)+"\t")
    oVarsSI.write(str(section_0[index])+"\t")
    # Nu*K/Dh
    film = nusselt(index)*interpolation(kRef_0,kRef_1,T_ref[index]) \
           / dh(index)
    htc_0.append(film)
    oVars.write(str(film*1E-4)+"\t"+str(thickness(index)*1E2)+"\n")
    oVarsSI.write(str(film)+"\t"+str(thickness(index))+"\n")
    # print (film)
    return film

def dh(index_in) :
  # Hydraulic diameter
    # Dh = 4 Area / wet perimeter
  # Case of spiral rectangular tube:
    # Dh = (2ab)/(a+b)
  if ref_type == "spiral" :
    return ((2*a*b)/(a+b))
  elif ref_type == "axial" :
    return (2*thickness(index_in))
  else :
    print ("WARNING: NOT IMPLEMENTED!!!!\n")
    return 1.

Re = []
Pr = []

def nusselt(index) :
  # Nu=f(Re,Pr) ,, Depends on correlation
   # Reynold's number, Re = rho*v_s*Dh/mu
    # v_s = q_r/(rho*A) -> Bulk velocity (away from boundary layer)
  
  v_s = q_r/(interpolation(densRef_0,densRef_1,T_ref[index])*a_r(index))
  re = interpolation(densRef_0,densRef_1,T_ref[index]) \
        * v_s * dh(index) / interpolation(visco_0,visco_1,T_ref[index])
  Re.append(re)
  fricc = 0.25 / (log10 ((epsi /(3.71*dh(index))) + (5.74 / (re)**0.9)))**2
  h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(g) * thickness(index))
  h_f.append(h_1)
  
  print(re)
  
  #print q_r, v_s, dh(index), interpolation(visco_0,visco_1,T_ref[index])
  #print re

   # smooth tube:
  friction=(1.82*log10(re)-1.64)**(-2.)
  prandtl = interpolation(prandtl_0,prandtl_1,T_ref[index])
  Pr.append(prandtl)
  oVars.write(str(v_s)+"\t"+str(re)+"\t")
  oVarsSI.write(str(v_s)+"\t"+str(re)+"\t")
  v_secc.append(v_s)
  #print h_f
   # Petukhov:
  return ( (friction/8.)*re*prandtl/(1.07+12.7*sqrt(friction/8.) \
         * (prandtl**(2./3.)-1)) )


  # adjust power density of first point to outer surface:
  pot_1[0] *=  section_1[0]/(section_1[0]+h)

# -------------------------------Start of iterative procedure:# ------------------------------- #
 # Computing the temperature of coolant:
#print(pot_1)

##1st loop
for i in range( len(section_0)-1 ) :
  # adjust power density to outer surface:
  pot_1[i+1] *=  section_1[i+1]/(section_1[i+1]+h)
  
  # print(pot_1)
  
  T_ref[i+1] = T_ref[i]+2.0*3.1416*(section_0[i+1]-section_0[i]) \
             / (q_r*interpolation(cpRef_0,cpRef_1,section_0[i])) \
               *pot_1[i]*section_1[i]
  oTemp.write(str(section_0[i+1]*1E2) +"\t"+str(T_ref[i+1]))
  oTempSI.write(str(section_0[i+1])+"\t"+str(T_ref[i+1]))
    
  
 # Computing the temperature of interface:
  T_metal[i+1] = T_ref[i+1] + pot_1[i+1]/film_coef(i+1)

  oTemp.write("\t")
  oTemp.write(str(T_metal[i+1]))
  oTemp.write("\n")
  oTempSI.write("\t")
  oTempSI.write(str(T_metal[i+1]))
  oTempSI.write("\n")

hf_tot[0] = h_f[0]

##2nd loop
for x in range(len(h_f)-1) :
  hf_tot[x+1] = float(h_f[x]) + hf_tot[x]
hf_bar = hf_tot[len(section_0)-2 ] * float(g) * rho / 100000
#print hf_tot,len(h_f),hf_tot2 

# Charge loss due to the section change

k_1 = 1 - (float(dh(4))**4) /(float(dh(3))**4) #Friction parameter
hs_1 = k_1 * (v_secc[4]**2) * 0.5 / float(g)
k_2 = 1 - (float(dh(15))**4) /(float(dh(14))**4) #Friction parameter
hs_2 = k_2 * (v_secc[15]**2) * 0.5 / float(g)
k_3 = 1 - (float(dh(19))**4) /(float(dh(18))**4) #Friction parameter
hs_3 = k_3 * (v_secc[19]**2) * 0.5 / float(g)
hs_bar = (hs_1 + hs_2 + hs_3) *rho * float(g) / 100000
h_tot = hf_bar + hs_bar

##3rd loop
for i in range( len(section_0)-2 ) :

# Computing the pressure evolution throughout the whole section:
  if section_0[i] <= 10.8333E-2 :
     P_secc[i+1] = P_secc[i] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i])
  elif section_0[i] <= 68.815E-2 :
      P_secc[i+1] = P_secc[i] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i])
      P_secc[14] = P_secc[13] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i] + hs_2)
  elif section_0[i] <= 123.462E-2 :
      P_secc[i+1] = P_secc[i] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i])
      P_secc[19] = P_secc[18] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i] + hs_3)
  else:
      P_secc[i+1] = P_secc[i] + rho*0.5*(v_secc[i]**2 - v_secc[i+1]**2)\
                + rho*float(g)*(float(section_1[i]) - float(section_1[i+1]))*1E-2 - rho*float(g)*float(h_f[i])
      
  P_secc[3] = P_secc[2] + rho*0.5*(v_secc[2]**2 - v_secc[3]**2)\
                + rho*float(g)*(float(section_1[2]) - float(section_1[3]))*1E-2 - rho*float(g)*float(h_f[2] + hs_1)
  
  oPresSI.write(str(section_0[i+1]))
  oPresSI.write("\t")
  oPresSI.write(str(P_secc[i+1]))
  oPresSI.write("\n")
  
#print (P_secc[3],h_tot)

v = []
for i,x in enumerate(section_1):
    v.append(q_r/(interpolation(densRef_0,densRef_1,T_ref[i])*pi))
             #*thickness(i))

# print(v)
print("SUCCESS!!!")

def plot_xy(x_in, y_in, xtitle, ytitle):
  plt.plot(x_in, y_in)
  plt.xlabel(xtitle)
  plt.ylabel(ytitle)
  plt.show(block=True)

fig,axs = plt.subplots(3)
fig.set_size_inches(10.5, 10.5, forward=True)
plt.xlabel("Lenght (m)")
axs[0].grid()
axs[0].set_ylabel("Power density, $\\frac{MW}{m^2}$")
axs[0].plot(section_0, np.array(pot_1)*1E-6)
axs[1].grid()
axs[1].set_ylabel("HTC, $\\frac{W}{m^2 K}$")
axs[1].plot(section_0[1:], htc_0)
axs[2].grid()
axs[2].set_ylabel("Temperature, $C$")
axs[2].plot(section_0[1:], T_metal[1:])
axs[2].plot(section_0[0:], T_ref[0:])
plt.show(block=True)

print(htc_0)