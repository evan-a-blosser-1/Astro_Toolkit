##############################################
# Author: Evan A. Blosser                    #
# Email: evan.a.blosser-1@ou.edu             #
#                                            #
# Program: 
#
#############################################
#
############
# Packages #
############
# Basic
import numpy as np
#
import matplotlib as matp
########################
# Variable Declaration #
########################
# Gravitational Parameter for Earth in km^3/s^2
# [If working in Canonical Units replace with 1 (DU^3/TU^2)]
mu_e = 3.986004418e5
# Equitorial Radius of the Earth in km
# [If working in Canonical Units replace with 1 (DU)]
r_e = 6378.136
###############
# User Inputs #
###############
#
p  = float(input('Enter the Periapsis:')) 
e  = float(input('Enter the Eccentricity:'))
nu = float(input('Enter the True Anomaly:'))
i  = float(input('Enter the Inclination:'))
Om = float(input('Enter the Longitude of Ascending Node:'))
w = float(input('Enter the Argumaent of Periapsis:'))
################
# Calculations #
################
# Position magnitue for elements
r = p / (1 + e*np.cos(nu))
# Calculate perifocal vector
r_pf_p = r*np.cos(nu)
r_pf_q = r*np.sin(nu)
# Set perifocal position vector in array
r_pf = np.array([r_pf_p, r_pf_q, 0]) 
# Calculate perifocal velocity
v_pf_p = np.sqrt(mu_e/p)*(-np.sin(nu))
v_pf_q = np.sqrt(mu_e/p)*(e + np.cos(nu))
# Set perifocal velocity vector arry
v_pf = np.array([v_pf_p, v_pf_q, 0])
#
#
##################
# Rotaton matrix #
##################
R_11 =  np.cos(Om)*np.cos(w) - np.sin(Om)*np.sin(w)*np.cos(i)
R_12 = -np.cos(Om)*np.sin(w) - np.sin(Om)*np.cos(w)*np.cos(i)
R_13 =  np.sin(Om)*np.sin(i)
R_21 =  np.sin(Om)*np.cos(w) + np.cos(Om)*np.sin(w)*np.cos(i)
R_22 = -np.sin(Om)*np.sin(w) + np.cos(Om)*np.cos(w)*np.cos(i)
R_23 = -np.cos(Om)*np.sin(i)
R_31 =  np.sin(w)*np.sin(i)
R_32 =  np.cos(w)*np.sin(i)
R_33 =  np.cos(i)

R_matrix = np.array([[R_11, R_12, R_13], [R_21, R_22, R_23], 
                     [R_31, R_32, R_33]])

r_geo = np.matmul(R_matrix,r_pf)

v_geo = np.matmul(R_matrix,v_pf)


print (f'{r_geo}& {v_geo}')


