# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 23:59:44 2023

@author: Evan A. Blosser 
 Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com
 
 Program: Hohmann Transfer Delta-V Calculator
         - For calcualting the requird Delta-V for a Hohmann stransfer
           between 2 orbits. Both must be Circular in order to use the 
           Hohmann transfer.
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
r_i = 3000.0 + r_e
r_f = 6000.0 + r_e
# Semi-major axis of transfer
a_t = (r_i+r_f)/2.0
# Energy of transfer
E_t = -mu_e/(2.0*a_t)
###############################
# Initial/Transfer Velocities #
###############################
# Initial Velocity of 1st circular orbit
v_i = np.sqrt(mu_e/r_i)
# Final velocity of 2nd circular orbit
v_f = np.sqrt(mu_e/r_f)
# Transfer velocity at 1st circlular orbit
v_tc1 = np.sqrt(2.0*(mu_e/r_i+E_t))
# Transfer velocity at 2nd circlular orbit
v_tc2 = np.sqrt(2.0*(mu_e/r_f+E_t))
########################
# Delta-V Calculations #
########################
# Delta-V required to initilize transfer orbit
delta_v_1 = v_tc1 - v_i
# Delta-V required to end transfer and achieve desired orbit
delta_v_2 = v_f - v_tc2
# Total required Delta-V 
delta_v_total = abs(delta_v_2)+abs(delta_v_1)
print(f"delta_v_total = {delta_v_total:7.5} km/s")