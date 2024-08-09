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
#####################################
# Mathmatical!!                     #
import numpy as np                  #
import math as mt                   #
# System Time                       #
import time                         #
# System Error/Exit                 #
from sys import exit                #
#####################################

###########################
## Constants Declaration ##
###########################
# Gravitational Parameter for Earth in km^3/s^2
# [If working in Canonical Units replace with 1 (DU^3/TU^2)]
mu_e = 3.986004418e5
# Equitorial Radius of the Earth in km
# [If working in Canonical Units replace with 1 (DU)]
r_e = 6378.136
# Maximum Iterations
nmax = 1000

############
### Main ###
############
#
'''
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
'''
# Debug inputs
r_1  = np.array([-1285.58, 8784.04, 2590.81]) 
r_2  = np.array([-9076.47, 5525.34, 3756.74])
t_m  = 20
pro_retro = 'y'

########
# Time #
########
Delta_t = t_m*60
#%% Calculations 
########################
# Initial Calculations #
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
## Lambert's  ##
################
#
##########################################
# True Anomaly/ Prograde or Reotrograde? #
##########################################
Error_Message_pro_ret = f""" 
{'-'*42}
| Error: Please Enter (Y/N) for  
|         Prograde or Retrograde
|         [Not Case Sensitive]
{'-'*42}
"""
# For Prograde Orbit                                          #
if pro_retro in('y','Y'):                                     #
#### Checking The Z component of Position Vec. Cross Product  #    
    if C_12[2] >= 0:                                          #
######## Assinging Change in True Anomaly Equation            #
        Delta_nu = np.arccos(D_12/(r_i*r_f))                  #
    else:                                                     #
        Delta_nu = 2*np.pi - np.arccos(D_12/(r_i*r_f))        #
###############################################################
# For Retrograde Orbit                                        #
elif pro_retro in('n','N'):                                   #
#### Checking The Z component of Position Vec. Cross Product  # 
    if C_12[2] >= 0:                                          #
######## Assinging Change in True Anomaly Equation            #
        Delta_nu = 2*np.pi - np.arccos(D_12/(r_i*r_f))        #
    else:                                                     #
        Delta_nu = np.arccos(D_12/(r_i*r_f))                  #
else:                                                         #
    print(Error_Message_pro_ret)                              #
    input('| Press Enter To Exit...')                         #
    exit()                                                    #
###############################################################
###################################################
# Constant (A) within Universal Kepler's Equation #############
A =  np.sin(Delta_nu)*np.sqrt((r_i*r_f)/(1-np.cos(Delta_nu))) #
###############################################################
##############
# Define Z_0 #
z_0 = 0      #
##############
##########################################################
# Define Lambert's Problem (Utilizing Kepler's Equation) #
##########################################################
#   These Include the Stumpff Functions C & S, and 
#       the Kepler's Universal Equation as a Function
#       of Z as well as the functions derivitive for 
#       the iteration of Z; (Y is used as a substitution).
#
#   tried implimenting a summation for C & S:
#       (-1**k*z**k)/mt.factorial(2*k + 3)
#   revisit but approx. as shown in 20230213.pdf slde 31
##########################################################
# S Stumpff function approximated to the third term      #
def S(z):
    S_Stumpff = 1/6 - z/120 + z**2/5040 - z**3/362880
    return S_Stumpff
# C Stumpff function approximated to the third term      #
def C(z):
    C_Stumpff = 1/2 - z/24 + z**2/720 - z**3/40320
    return C_Stumpff
# Y substitution function to simplify Lambert's Problem  #
def y(z):
    y_sub = r_i + r_f + A*((z*S(z) - 1)/np.sqrt(C(z)))
    return y_sub
# Kepler's Eq. as a function of z                        #
def F(z):
    F_kepler = ((y(z)/C(z))**1.5)*S(z) + A*np.sqrt(y(z))\
        - np.sqrt(mu_e)*Delta_t
    return F_kepler
# The Derivitive of Kepler's Eq. as a function of z       #
def F_Prime(z):
    if z == 0:
        F_z_Prime = (np.sqrt(2)/40)*y(0)**1.5 + \
            (A/8)*(np.sqrt(0) + A*np.sqrt(1/(2*y(0)))) 
    else:
        F_z_Prime = (y(z)/C(z))**1.5*(1/(2*z)*(C(z) -\
        (3*S(z)**2)/(4*C(z)) + (3*S(z)**2)/(4*C(z)))) + \
            (A/8)*(3*(S(z)/C(z))*np.sqrt(y(z)) + A*np.sqrt(C(z)/y(z)))
    return F_z_Prime
###############################################################################
# Iterate that Z! #
###################
#   Using Newtonian Interpolation
#       -Finding the Roots (ie solutions) to
#        an Implicit Function
############################################
def Newt_Int(x_0):
#### Iteration for loop 
    for i in range(0, nmax):
        ratio = F(x_0)/F_Prime(x_0)
        Newt_Int_Message = f"""
        {'-'*42}
        | Itration {i} 
        | F = {F} 
        | F' = {F_Prime}
        | |F/F'| = {np.abs(ratio)}
        {'-'*42}
        """
        print(Newt_Int_Message)
# stop when achieve the desired precision
        if np.abs(ratio) <= 1.0e-6:
            break
        x = x_0 - ratio
        x_0 = x #update the value of x_0
        print(f"Initial values:\nS(z) = {S(x_0)}, C(z) = {C(x_0)}, y(z) = {y(x_0)}, z = {x_0}\n")
        if i == nmax - 1:
            print(f"Newton: no convergence. Try another initial guess.")
    return x

z = Newt_Int(z_0)
###############################################################################
#########################
# Velocity calculations #
#########################
# Calculating the Velocity
#   vectors using the
# Lagrange Coefficiants 
###############################
# f lagrange coeff            #
f     = 1 - (y(z)/r_i)        #
# g lagrange coeff            #
g     = A*np.sqrt(y(z)/mu_e)  #
# g dot lagrange coeff        #
g_dot = 1 - (y(z)/r_f)        #
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
if z < 0:
    print(f"z = {z}; the orbit is hyperbolic.")
elif z == 0:
    print(f"z = {z}; the orbit is parabolic.")
else:
    print(f"z = {z}; the orbit is elliptic.")

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