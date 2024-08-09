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
# Plot and Animation
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

###############
# User Inputs #
########################################################################
r_i_a    = float(input('Enter the Initial Altitude:'))                 #
r_f_a    = float(input('Enter the Final Altitude:'))                   #
nu       = float(input('Enter the True Anomaly Increase:'))            #
t_m      = float(input('Enter how many minutes passed:'))              #
########################################################################
'''
# Debug inputs
r_i_a    = 400
r_f_a    = 1000
nu       = 120
t_m      = 30
'''
########
# Time #
########
t = t_m*60
#
################
# Calculations #
##################################
# Calculating radius             #
r_i = r_i_a + 6378.136           #
r_f = r_f_a + 6378.136           #
# Assign first altitude as vec   #
r_1 = np.array([r_i,0,0])        #
# Assign second altitude as vec  #
r_f_i = r_f*np.cos(nu)           #
r_f_j = r_f*np.sin(nu)           #
r_2 = np.array([r_f_i ,r_f_j,0]) #
r_1mag = np.sqrt((r_1*r_1).sum())
r_2mag = np.sqrt((r_2*r_2).sum())
##################################
#
##############
# Lambert's  #
###################################################
A =  np.sin(nu)*np.sqrt((r_i*r_f)/(1-np.cos(nu))) #
# Set z and z_o                                   #
z   = 0.0000001                                   #
z_o = 0                                           #
###################################################
# Stumpff Functions #
############################################################
# Approximated as shown in 20230213.pdf slde 31            #
if z == 0:                                                 #
    S = 1/6                                                #
elif z < 0:                                                #
    S = (np.sinh(np.sqrt(z)) - np.sqrt(z))/(np.sqrt(z))**3 #
else:                                                      #
    S = (np.sqrt(z)-np.sin(np.sqrt(z)))/(np.sqrt(z))**3    #
# C approximation                                          #
if z == 0:                                                 #
    C = 1/2                                                #
elif z < 0:                                                #
    C = (np.cosh(np.sqrt(z)) - 1)/(-z)                     #
else:                                                      #
    C = (1 - np.cos(np.sqrt(z)))/z                         #
############################################################
# Substitution functions #
##########################
# Equation 5.38 
y = r_1mag + r_2mag + A*((z*S - 1)/np.sqrt(C))
# Equation 5.40
f_z_t= lambda z: ((y/C)**1.5)*S + A*np.sqrt(y) - np.sqrt(mu_e)*t
# Equation 5.43
#####################################################################
# F and F' #
############
if z == 0:
    f_zt_prime = lambda y: (np.sqrt(2)/40)*y**1.5 + (A/8)*(np.sqrt(0) +\
                                                           A*np.sqrt(1/(2*y))) 
else:
    f_zt_prime = lambda z:  (y/C)**1.5*(1/(2*z)*(C -(3*S**2)/(4*C) +\
                                                 (3*S**2)/(4*C))) +\
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
# Perigee Calculations #
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
####################################
print(f'Perigee altitude is: {p_alt}')


