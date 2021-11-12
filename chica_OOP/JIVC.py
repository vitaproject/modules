# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 08:21:27 2021

@author: jack.taylor

Finds the diamater and mach number of a circular inlet pipe to a JIVC-like
system. Returns the values for a single segment of a single cassette, there 
may be multiple cassettes on the carrier, and multiple adjacent rows and 
columns of JIVC-packs on the cassette.

"""

from pandas import DataFrame, concat
from numpy import pi, sqrt, linspace, array
from CoolProp.CoolProp import PropsSI as SI
from main import Panel
import matplotlib.pyplot as plt

class JIVC_manifold():

    """
    Mach number equation, may be used to calculate Mach number or velocity
    
    :param float density: Coolant density
    :param float Ajet: Cross-sectional area of the JIVC jet
    :param float R: Universal gas constant
    :param float M: Molar mass
    :param float gamma: Isentropic expansion factor
    """

    density = SI("D", "T", 373.15, "P", 100E5, "helium")
    Ajet = pi * (0.0005 ** 2)
    R = 8.3145
    M = 4.003E-3
    gamma = 1.667
    range_of_multipliers = array(linspace(1, 5, 9))
    speed_of_sound = sqrt(gamma * 373.15 * R / M)

    def __init__(self, Ma, n, m):

        self.mach_number = Ma
        self.n = n
        self.m = m
        pass

    def create_panel(self):
        self.panel = Panel(float(self.mach_number), "inputs.xml", "q_adjusted.asc", \
                           n=self.n, m=self.m)
        self.panel("JIVC")

    def repeated_operations(self):

        # - mass flow per jet set by the mach number
        mass_flow_per_jet = self.density * self.Ajet * self.speed_of_sound * \
            self.mach_number

        # - total mass flow set by the number of jets, slabs and mach number
        total_mass_flow = self.panel.n_slabs * self.panel.n_jets_per_slab \
            * mass_flow_per_jet

        # - volumetric flow rate of inlet pipe
        volumetric_flow_rate = total_mass_flow / self.density

        # - total area of the jets
        Aout = self.Ajet * self.panel.n_slabs * self.panel.n_jets_per_slab

        # - area of the inlet pipe
        Ain_range = Aout * self.range_of_multipliers

        # - range of diameters corresponding to range of areas
        Din_range = sqrt(4 * Ain_range / pi)

        # - velocity range for range of areas
        velocity_range = volumetric_flow_rate / Ain_range

        # - mach range for range of areas
        mach_range = velocity_range / self.speed_of_sound

        self.diameters = DataFrame(Din_range)
        self.mach = DataFrame(mach_range)

    def manager(self):
        self.create_panel()
        self.repeated_operations()

if __name__ == "__main__":

    # - vary n, m and Ma

    Diameter_master = []
    Mach_master = []
    n_list = []
    m_list = []

    m = linspace(5, 10, 6, dtype = int)
    n = linspace(5, 10, 6, dtype = int)
    Ma = linspace(0.3, 0.3, 1)
    labels = ["n = %s, m = %s" % (ni, mi) for ni in n for mi in m]
    n_list = [array(linspace(nx, nx, mx, dtype = int)) for ni, nx in enumerate(n) \
              for mi, mx in enumerate(m)]
    m_list = [len(n) for n in n_list]

    # - initialise storage

    for mi, mx in enumerate(m_list):
    
        diameter_subset = DataFrame()
        mach_subset = DataFrame()
    
        for Mai in Ma:

            test = JIVC_manifold(Mai, n_list[mi], mx)
            test.manager()
            diameter_subset = concat([test.diameters, diameter_subset], ignore_index = True)
            mach_subset = concat([test.mach, mach_subset], ignore_index = True)

        Diameter_master.append(diameter_subset)
        Mach_master.append(mach_subset)

    for index, diameter_set in enumerate(Diameter_master):
        plt.scatter(Mach_master[index][0], diameter_set[0], cmap = index, label = labels[index])
    plt.legend()
