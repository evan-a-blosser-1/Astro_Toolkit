# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 23:59:44 2023

@author: Evan A. Blosser 
 Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com
 
 Program:  Bi-elliptical Transfer Delta-V Optimizer
         - For calcualting the requird Delta-V for a bi-elliptical transfer, 
           while optimizing the eccentricity to minimize the total required
           Delta-V.
"""

'''
Packages
'''
#
# Basic
import numpy as np
#
#
'''
Variable Declaration
'''
#
# Gravitational Parameter for Earth in km^3/s^2
# [If working in Canonical Units replace with 1 (DU^3/TU^2)]
mu_e = 3.986004418e5
#
# Equitorial Radius of the Earth in km
# [If working in Canonical Units replace with 1 (DU)]
r_e = 6378.136
#
#
"""
Main
"""
############################
# Initial Orbital Elements #
############################
# These radii are in km
r_i = 200.0 + r_e
r_f = 72559.50 + r_e
# Semi-major axis of transfer
a_t = (r_i + r_f)/2.0
# Energy of transfer
E_t = -mu_e/(2.0*a_t)
######################
# Initial Velocities #
######################
# Initial Velocity of 1st circular orbit
v_i = np.sqrt(mu_e/r_i)
# Final velocity of 2nd circular orbit
v_f = np.sqrt(mu_e/r_f)
###########################
# Optimizing eccentricity #
###########################
for inc in np.arange( 0.0,  1.0, 0.0001):
    e_orb_2 = inc
#####################################
# Transfer Velocities for 2nd Orbit #
#####################################
#
r_orb2_per = r_i 
#
r_orb2_apo = r_orb2_per*(1.0 + e_orb_2)/(1.0 - e_orb_2)
#
a_orb2 = (r_orb2_apo + r_orb2_per)/2.0
#
E_orb2 = -mu_e/(2.0*a_orb2)
#
v_orb2_i = np.sqrt(2.0*(mu_e/r_orb2_per + E_orb2))
#
v_orb2_f = np.sqrt(2.0*(mu_e/r_orb2_apo + E_orb2))
#
P_orb2 = 2.0*np.pi*np.sqrt(a_orb2**3/mu_e)
TOF_orb2 = P_orb2/2.0
#############################
# Eccentricity of 3rd Orbit #
#############################
#
r_orb3_per = r_f
r_orb3_apo = r_orb2_apo
e_orb_3 = (r_orb3_apo - r_orb3_per)/(r_orb3_apo + r_orb3_per)
#####################################
# Transfer Velocities for 3rd Orbit #
#####################################
#
a_orb3 = (r_orb3_apo + r_orb3_per)/2.0
#
E_orb3 = -mu_e/(2.0*a_orb3)
#
v_orb3_f = np.sqrt(2.0*(mu_e/r_orb3_per + E_orb3))
#
v_orb3_i = np.sqrt(2.0*(mu_e/r_orb3_apo + E_orb3))
#
P_orb3 = 2.0*np.pi*np.sqrt(a_orb3**3/mu_e)
TOF_orb3 = P_orb3/2.0
########################
# Delta-V Calculations #
########################
# Delta-V required to initilize transfer orbit
delta_v_1 = v_orb2_i - v_i
# Delt-V required to go from the 2nd to 3rd orbit
delta_v_2 = v_orb3_i - v_orb2_f
# Delta-V required to end transfer and achieve desired orbit
delta_v_3 = v_f - v_orb3_f 
# Total required Delta-V 
delta_v_total = abs(delta_v_3)+abs(delta_v_2)+abs(delta_v_1)



print(f"delta_v_total = {delta_v_total:7.5} km/s")
print(f"e_orb_2 = {e_orb_2:7.5}")