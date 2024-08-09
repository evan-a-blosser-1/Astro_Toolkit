##############################################
# Author: Evan A. Blosser                    #
# Email: evan.a.blosser-1@ou.edu             #
#                                            #
# Program: Calculating perigee given 2       #
#          altitudes, change in true anomaly #
#          and time passed                   #   
##############################################
#
############
# Packages #
############
# Basic
import numpy as np
########################
# Variable Declaration #
########################
# Gravitational Parameter for Earth in km^3/s^2
# [If working in Canonical Units replace with 1 (DU^3/TU^2)]
mu_e = 3.986004418e5
# Equitorial Radius of the Earth in km
# [If working in Canonical Units replace with 1 (DU)]
r_e = 6378.136
############
### Main ###
############
#
"""
###############
# User Inputs #
###################################################################
# Input & Assignment of Initial Position Vector                   #
r_1_i = float(input("|Enter the initial position Vector, (i): ")) #
r_1_j = float(input("|Enter the initial position Vector, (j): ")) #
r_1_k = float(input("|Enter the initial position Vector, (k): ")) #
r_1  = np.array([r_1_i, r_1_j, r_1_k])                            #
# Input & Assignment of Final Position Vector                     #
r_2_i = float(input("|Enter the final position Vector, (i): "))   #
r_2_j = float(input("|Enter the final position Vector, (j): "))   #
r_2_k = float(input("|Enter the final position Vector, (k): "))   #
r_2  = np.array([r_2_i, r_2_j, r_2_k])                            #
# Input for time passed between vectors                           #
t_m      = float(input('Enter how many minutes passed:'))         #
# Is this a Prograde Orbit                                        #
pro_retro = input("|Is this a Prograde Orbit (y/n):")             #
###################################################################
#
"""
# Debug inputs
r_1  = np.array([-1285.58, 8784.04, 2590.81]) 
r_2  = np.array([-9076.47, 5525.34, 3756.74])
t_m  = 20
pro_retro = 'y'
########
# Time #
########
t = t_m*60
#
################
# Calculations #
################
#
################################
# Position Vector Magnitudes   #
r_i = np.linalg.norm(r_1)      #
r_f = np.linalg.norm(r_2)      #
# Cross product for nu calc.   #
C_12 = np.cross(r_1,r_2)       #
# Dot Product for nu clac.     #
D_12 = np.dot(r_1,r_2)         #
################################
#
################
# True Anomaly #
################
if pro_retro == 'y':
    if C_12[2] >= 0:
        nu = np.arccos(D_12/(r_i*r_f))
    else:
        nu = 360 - np.arccos(D_12/(r_i*r_f))
if pro_retro == 'n':
    if C_12[2] >= 0:
        nu = 360 - np.arccos(D_12/(r_i*r_f))
    else:
        nu = np.arccos(D_12/(r_i*r_f))
###################
#
##############
# Lambert's  #
###################################################
A =  np.sin(nu)*np.sqrt((r_i*r_f)/(1-np.cos(nu))) #
# Set z and z_o                                   #
z   = 0.000000000000001                           #
z_o = 0.0                                         #
###################################################


############################################################
# Approximated as shown in 20230213.pdf slde 31            #
if z == 0:                                                 #
    S = 1/6                                                #
if z < 0:                                                  #
    S = (np.sinh(np.sqrt(z)) - np.sqrt(z))/(np.sqrt(z))**3 #
else:                                                      #
    S = (np.sqrt(z)-np.sin(np.sqrt(z)))/(np.sqrt(z))**3    #
# C approximation                                          #
if z == 0:                                                 #
    C = 1/2                                                #
if z < 0:                                                  #
    C = (np.cosh(np.sqrt(z)) - 1)/(-z)                     #
else:                                                      #
    C = (1 - np.cos(np.sqrt(z)))/z                         #
############################################################
# Substitution functions #
##########################
# Equation 5.38 
y = r_i + r_f + A*((z*S - 1)/np.sqrt(C))
# Equation 5.40
f_z_t= lambda z: ((y/C)**1.5)*S + A*np.sqrt(y) - np.sqrt(mu_e)*t
# Equation 5.43

#####################################################################
# F and F' #
############
if z == 0:
    f_zt_prime = lambda y: (np.sqrt(2)/40)*y**1.5 + \
        (A/8)*(np.sqrt(0) + A*np.sqrt(1/(2*y))) 
else:
    f_zt_prime = lambda z:  (y/C)**1.5*(1/(2*z)*(C -\
                                               (3*S**2)/(4*C) + \
                                               (3*S**2)/(4*C))) + \
        (A/8)*(3*(S/C)*np.sqrt(y) + A*np.sqrt(C/y))
######################################################################
# Iterate that Z! #
###################
#
for j in range(0,200000):
    function_ratio = f_z_t(z)/f_zt_prime(z)
    print(f"Iteration [{j}], Yields Function Ratio = {function_ratio}")
    func_check = np.abs(function_ratio)
    if func_check <= 1.0e-6:
        break
    z   = z_o - function_ratio
    z_o = z
#########################
# Lagrange Coefficiants #
#########################
# f lagrange coeff
f     = 1 - (y/r_i)
# g lagrange coeff
g     = A*np.sqrt(y/mu_e)
# g dot lagrange coeff
g_dot = 1 - (y/r_f)
#################
#
#########################
# Velocity calculations #
###############################
# First Velocity Vector       #
v_1 = (1/g)*(r_2 - f*r_1)     #
# Second Velocity Vector      #
v_2 = (1/g)*(g_dot*r_2 - r_1) #
#                            #
print (f'{v_1},{v_2}')      #
############################
# Final Calculations #
####################################
# Angular momentum                 #
h = np.cross(r_1,v_1)              #
# Magnitude of Ang. Mom.           #
h_mag = np.sqrt((h*h).sum())       #
# Semi-Latus Rectum                #
p = h_mag**2/mu_e                  #
# Eccentricity vector              #
e = np.cross(v_1,h)/mu_e - r_1/r_i #
# Eccentricity Magnitude           # 
e_mag = np.sqrt((e*e).sum())       # 
# Perigee Radius                   #
r_p = p/(1 + e_mag)                #
# Perigee Altitude                 #
p_alt = r_p - r_e                  #
# Specific Energy                  #
E = -(mu_e*(1-e_mag**2))/(2*p)     #
# Inclination                      #
i = np.arccos(h[2]/h_mag)          #
i_deg = np.rad2deg(i)
####################################
Answers = f"""
|The specific energy is: 
|{E} (km^2/s^2)
|
|The perigee altitude is: 
|{p_alt} (km)
|
|The inclination of th eorbit is: 
|{i_deg} (deg)
"""
print(Answers)