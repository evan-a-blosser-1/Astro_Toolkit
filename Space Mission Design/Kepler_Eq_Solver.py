# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 21:14:16 2023

@author: Evan A. Blosser 
 Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com
 
 Program:  Kepler's Equation Calculator
         - 
"""
'''
Packages
'''
#
# Basic Math
import numpy as np

'''
Variable Declaration
'''
# Eccentricity calculated from problem 2
e = 0.012756420733831
"""
Main
"""
############
# Below PI #
############
# Mean Anomaly initial guess below pi ~ 3.14 
M_b = 0.1
E_b = M_b + e/2.0

    
f = lambda E: E - e*np.sin(E)-M_b
f_prime = lambda E: 1.0- e*np.cos(E)
newton_root_below_PI = E_b - (f(E_b))/(f_prime(E_b))

print("Newton root below PI =", newton_root_below_PI)


############
# Above PI #
############
# Mean Anomaly initial guess above pi ~ 3.14 
M_a = 4.1
E_a = M_a + e/2.0

f_2 = lambda E: E - e*np.sin(E)-M_a
f_prime2 = lambda E: 1.0- e*np.cos(E)
newton_root_above_PI = E_a - (f(E_a))/(f_prime(E_a))

print("Newton root above PI =", newton_root_above_PI)