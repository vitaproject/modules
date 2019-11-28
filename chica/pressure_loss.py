#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
CHICA - Cooling and Heating Interaction and Corrosion Analysis
        Calculates the temperature evolution for a thermofluid flowing along an axis-symmetric channel
"""
__author__ = "Marcos Parro"
__copyright__ = "Copyright 2012, Marcos Parro"
__credits__ = ["Marcos Parro", "Daniel Iglesias"]
__license__ = "LGPL"
__version__ = "1.0.1"
__maintainer__ = "Daniel Iglesias"
__email__ = "Daniel.Iglesias@tokamakenergy.co.uk"
__status__ = "Experimental"

from pandas import DataFrame
import pandas as pd
import numpy as np
from math import*
import os

os.chdir("/home/markus/Documentos/Beam Dump/Refrigeración/programas_refrigeracion/Nominal")

q_r = 30
epsi = 0.0000065
g = 9.81
duty_cycle = 1
tb = []
tb = [0.0 for i in range(2499) ]
tb [0] = 31
rho_aux = []
rho_aux = [0.0 for i in range(2499) ]
cp_aux = []
cp_aux = [0.0 for i in range(2499) ]
pr_aux = []
pr_aux = [0.0 for i in range(2499) ]
mu_aux = []
mu_aux = [0.0 for i in range(2499) ]

pot = pd.read_table("curva de potencia1mm.txt", names = ["Secc", "Pot"], sep = "\t")
sect = pd.read_table("curva de seccion1mm.txt", names = ["Secc", "Radius"], sep = "\t")
sect_1 = sect.multiply(1E-3, axis = 1)
basic = pd.concat([sect_1, pot["Pot"]], axis = 1)
rho = pd.read_table("densidad refrigerante.txt", names = ["Temp", "Dens"], sep = "\s+")
cp = pd.read_table("calor especifico refrigerante.txt", names = ["Temp", "Cp"], sep = "\s+")
pr = pd.read_table("numero de Prandtl.txt", names = ["Temp", "Pr"], sep = "\s+")
mu = pd.read_table("viscosidad dinamica refrigerante.txt", names = ["Temp", "Mu"], sep = "\s+")
geom = DataFrame(np.arange(9996).reshape(2499,4), columns = ["Secc", "Thk", "Width","Area"])
temp = DataFrame(np.arange(12495).reshape(2499,5), columns = ["Secc", "Pot", "Tbulk", "Tsb", "Tonb"])
var = DataFrame(np.arange(12495).reshape(2499,5), columns = ["Secc", "Vel", "Re", "Tsb", "Tonb"])

def Eint(row):
    
    if row["Secc"] <= 20*1E-2 : return (1.7E-2 - row["Radius"])
    elif row["Secc"] <= 2000*1E-3 :
     return 5E-3
    elif row["Secc"] <= 2015*1E-3 :
     return 0.1*1E-3 + Eint(row - 1)
    elif row["Secc"] > 2015*1E-3 :
     return 6.5E-3

def interpolation(list_0, list_1, value_0) :
  for j in range( len(list_0)-1 ) :
    if list_0[j+1] > value_0 :

      return float( list_1[j] + (list_1[j+1]-list_1[j])/ \
               (list_0[j+1]-list_0[j])*(value_0-list_0[j]) )
  # Over maximum value...
  return list_1[j+1]
  
def espesor(row) :

    if row["Secc"] <= 20E-2 : return 2.3614*1E-2
    elif row["Secc"] <= 68.8E-2 : 
      return ( 2.3614 + (row["Secc"] - 20E-2)*(1.212-2.3614)/(68.8E-2-20E-2) )*1E-2
    elif row["Secc"] <= 125.13E-2 : 
      return ( 1.212 + (row["Secc"] - 68.8E-2)*(.7-1.212)/(125.13E-2-68.8E-2) )*1E-2
    elif row["Secc"] <= 200E-2 :
      return 0.7E-2
    elif row["Secc"] <= 201.5E-2 :
      return (0.7E-2 - Eint(row) + 5E-3)
    elif row["Secc"] <= 250E-2 :
      return (0.55 + (row["Secc"] - 201.5E-2)*(.7-.55)/(250*1E-2-201.5*1E-2) )*1E-2
      
def a_r(row) :
    
    return pi*( (row["Radius"]+espesor(row) + Eint(row) )**2 -\
          (row["Radius"] + Eint(row))**2 )
          
# Hydraulic diameter (m):

def dh(row) :
  # Hydraulic diameter
    # Dh = 4 Area / wet perimeter
  
    return (2*espesor(row))
    
def vel(row) :
  
   return  q_r/(row["Rho"]*a_r(row))
       
def rey(row) :
  
  re = row["Rho"] * vel(row) * dh(row) / row["Mu"]
  return re

def power(row) :
    
    pot_1 =  duty_cycle*row["Radius"]/(row["Radius"]+Eint(row))
    return pot_1
    
def fricc(row) :

  a = 2/log(10)
  b = epsi/(row["Diameter"] * 3.7) 
  d = (log(10) * row["Re"])/5.2
  s = b * d + log(d)
  q = s**(s/(s+1))
  g = b * d + log(d/q)
  z = log(q/g)
  Dla = (g * z)/(g + 1)
  Dcfa = Dla*(1 + z/(2*(g + 1)**2 + (z/3)*(2*g - 1)))
  friction = (a * (log(d/q) + Dcfa))**(-2)
  
  return friction

# Darcy-Weishbach linear loss height (m):

def carga(row)  :
  
  h_1 = (fricc(row) * (1E-3) * (row["Vel"])**2) / (2 * float(g) * row["Diameter"])
  
  return h_1
    
geom["Secc"] = sect_1["Secc"]
geom["Thk"] = sect_1.apply(Eint, axis = 1)
geom["Width"] = sect_1.apply(espesor, axis = 1)
geom["Area"] = sect_1.apply(a_r, axis = 1)
geom["Diameter"] = sect_1.apply(dh, axis = 1)
geom["InR"] = geom["Thk"] + sect_1["Radius"]
temp["Secc"] = sect["Secc"]
temp["Pot"] = sect.apply(power, axis = 1)*pot["Pot"]*1E4

for i in range(2498) :
    
    tb[i+1] = tb[i] + 2.0*3.1416*1E-3 * (geom.ix[i]["InR"]) * temp.ix[i]["Pot"]\
    / (q_r*4184)
    
    temp["Tbulk"] = tb
    rho_aux[i+1] = interpolation(rho["Temp"], rho["Dens"], temp.ix[i]["Tbulk"])
    cp_aux[i+1] = interpolation(cp["Temp"], cp["Cp"], temp.ix[i]["Tbulk"])
    pr_aux[i+1] = interpolation(pr["Temp"], pr["Pr"], temp.ix[i]["Tbulk"])
    mu_aux[i+1] = interpolation(mu["Temp"], mu["Mu"], temp.ix[i]["Tbulk"])
    basic["Rho"] = rho_aux
    basic["Cp"] = cp_aux
    basic["Pr"] = pr_aux
    basic["Mu"] = mu_aux

basic.ix[0]["Rho"] = basic.ix[1]["Rho"]
basic.ix[0]["Cp"] = basic.ix[1]["Cp"]
basic.ix[0]["Pr"] = basic.ix[1]["Pr"]
basic.ix[0]["Mu"] = basic.ix[1]["Mu"]

basic["Diameter"] = geom["Diameter"]
basic["Vel"] = basic.apply(vel, axis = 1)
basic["Re"] = basic.apply(rey, axis = 1)
basic["Delta"] = basic.apply(carga, axis = 1)

