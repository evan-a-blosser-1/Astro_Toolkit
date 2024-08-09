# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 00:51:47 2023

@author: Evan A. Blosser 
 Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com
 
 Program: Orbital Calculator  
         - Used with three position vectors measured
           in the Geocentric frame.
"""
##############
## Packages ##
##############
# Basic Math
import numpy as np
#
import sys
from sys import exit
##########################
## Variable Declaration ##
##########################
# Gravitational Parameter for Earth in km^3/s^2
# [If working in Canonical Units replace with 1 (DU^3/TU^2)]
mu_e = 3.986004418e5
# Equitorial Radius of the Earth in km
# [If working in Canonical Units replace with 1 (DU)]
r_e = 6378.136
##############
###  Main  ###
##############
#
#######################################
welcome =f"""
{'-'*42}
            Welcome to the 
Geocentric Equtorial Position Calcualtor
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
 Utilizing the Gibbs Method, to calculate the velocity of
  an orbiting body given 3 position vectors measured 
  in the Geocentric frame. 
{'-'*42}
  Warning: The position vectors should be input as
  each i, j, & k direction seperatly when prompted. 
{'-'*42}
"""                              
print(welcome)         

################################
# Input Position Vectors Below #
##################################################################
# Input & Assignment of First Position Vector                    #
r_1_i = float(input("|Enter the first position Vector, (i): "))  #
r_1_j = float(input("|Enter the first position Vector, (j): "))  #
r_1_k = float(input("|Enter the first position Vector, (k): "))  #
r_1  = np.array([r_1_i, r_1_j, r_1_k])                           #
# Input & Assignment of Second Position Vector                   #
r_2_i = float(input("|Enter the second position Vector, (i): ")) #
r_2_j = float(input("|Enter the second position Vector, (j): ")) #
r_2_k = float(input("|Enter the second position Vector, (k): ")) #
r_2  = np.array([r_2_i, r_2_j, r_2_k])                           #
# Input & Assignment of Third Position Vector                    #
r_3_i = float(input("|Enter the third position Vector, (i): "))  #
r_3_j = float(input("|Enter the third position Vector, (j): "))  #
r_3_k = float(input("|Enter the third position Vector, (k): "))  #
r_3  = np.array([r_3_i, r_3_j, r_3_k])                           #
##################################################################

# Test Values  #################################  DEBUG Values
#r_1  = np.array([5887.0, -3520.0, -1204.0])   # ITs WORKING !!! 
#r_2  = np.array([5572.0, -3457.0, -2376.0])   # ITS WORKING !!!!!
#r_3  = np.array([5088.0, -3289.0, -3480.0])   #
################################################
# Calculate Magnitude of Position Vectors #
r_1_mag = np.sqrt((r_1*r_1).sum())       #
r_2_mag = np.sqrt((r_2*r_2).sum())      #
r_3_mag = np.sqrt((r_3*r_3).sum())     #
#######################################
#######################################
Magnitude =f"""
{'-'*42}
|The megnitude of the first position vector
| is: {r_1_mag} (km)
|
|The megnitude of the second position vector 
| is: {r_2_mag} (km)
|
|The megnitude of the third position vector 
| is: {r_3_mag} (km)
|
{'-'*42}
"""                              
print(Magnitude)                
###################################################
# Cross Products                                  #
#  - For ease of calcualting the vectors N,D, & S #
#    used in the Gibb's Method                    #
###################################################
C_1_2 = np.cross(r_1,r_2) #
C_2_3 = np.cross(r_2,r_3) #
C_3_1 = np.cross(r_3,r_1) #
######################################################
#Conservation Check                                  #
# - determine if the angular momentum is conserved   #
######################################################
C_check = C_2_3/np.sqrt((C_2_3*C_2_3).sum()) #
U_check = r_1/r_1_mag                       #
Conservation = np.dot(U_check,C_check)     #
###########################################
Conserve_message = f"""
{'-'*42}
|The calculation of the dot product of the  
| unit vector in the direction of the first 
| position vector with the unit vector 
| normal to the second and third positon
| vector was found to be: {Conservation}
|
{'-'*42}
"""
Conserve_success = f"""
{'-'*42}
| Calcualtions for the conservation 
| of angular momentum was found to be
| relatively Zero, thus the Gibb's 
| Method can be used with these 3
| given positon vectors.
{'-'*42}
"""
Error= """
| The angular momentum was not 
| conserved between these 3 
| vectors, and Gibb's Method
| can NOT continue.
|
|Exiting.....
"""
print(Conserve_message)
if Conservation < 10**-4:
    print(Conserve_success)
else:
    print(Error)
    input("| Press any key to Exit...")
    exit()
##########################
# N, D, & S Calculations #
##########################
N     = r_1_mag*C_2_3 + r_2_mag*C_3_1 + r_3_mag*C_1_2
N_mag = np.sqrt((N*N).sum())
#
D     = C_1_2 + C_2_3 + C_3_1
D_mag = np.sqrt((D*D).sum())
#
Delta_r_23 = r_2_mag - r_3_mag
Delta_r_31 = r_3_mag - r_1_mag
Delta_r_12 = r_1_mag - r_2_mag
S     = r_1*Delta_r_23 + r_2*Delta_r_31 + r_3*Delta_r_12
#################################
# Velocity Vecctor Calculations #
#################################
D_x_r1 = np.cross(D,r_1)
# Velocity of first positon
v_1 = np.sqrt(mu_e/(N_mag*D_mag))*((D_x_r1/r_1_mag) + S) 
#
D_x_r2 = np.cross(D,r_2)
# Velocity of first positon
v_2 = np.sqrt(mu_e/(N_mag*D_mag))*((D_x_r2/r_2_mag) + S) 
#
D_x_r3 = np.cross(D,r_3)
# Velocity of first positon
v_3 = np.sqrt(mu_e/(N_mag*D_mag))*((D_x_r3/r_3_mag) + S) 
#######################
# Velocity Magnitudes #
####################################
v_1_mag = np.sqrt((v_1*v_1).sum()) #
v_2_mag = np.sqrt((v_2*v_2).sum()) #
v_3_mag = np.sqrt((v_3*v_3).sum()) #
####################################
Vector_message = f"""
{'-'*42}
| The velocity of the first position 
|  vector is: 
|  {v_1[0]} (i)  
|  {v_1[1]} (j) 
|  {v_1[2]} (k)
|
| The magnitude of this velocity 
|  vector is: {v_1_mag} (km/s)
{'-'*42}
| The velocity of the second position 
|  vector is: 
|  {v_2[0]} (i)  
|  {v_2[1]} (j) 
|  {v_2[2]} (k)
|
| The magnitude of this velocity 
|  vector is: {v_2_mag} (km/s)
{'-'*42}
| The velocity of the third position 
|  vector is: 
|  {v_3[0]} (i)  
|  {v_3[1]} (j)
|  {v_3[2]} (k)
|
| The magnitude of this velocity 
|  vector is: {v_3_mag} (km/s)
{'-'*42}
"""
print(Vector_message)
#%%
#################################################
## Elemetns From Positon and Velocity Vectors  ##
##  - This uses the previously found velocity  ##
##     vectors along with the position vectors ##
##     to calculate the orbital elements       ##
##                                             ##
#################################################
#
####################
# Angular Momentum #
####################################
# Angular momentum of 1st position #
h_1 = np.cross(r_1,v_1)            #
h_1_mag = np.sqrt((h_1*h_1).sum()) #
# Angular momentum of 2nd position #
h_2 = np.cross(r_2,v_2)            #
h_2_mag = np.sqrt((h_2*h_2).sum()) #
# Angular momentum of 3rd position #
h_3 = np.cross(r_3,v_3)            #
h_3_mag = np.sqrt((h_3*h_3).sum()) #
# Averaging the angular momentum   #
#                        ??        #
####################################
#
#####################
# Semi-Latus Rectum # 
#######################
# 1st semi-latus rec  #
p_1 = h_1_mag**2/mu_e #
# 2nd semi-latus rec  #
p_2 = h_2_mag**2/mu_e #
# 3rd semi-latus rec  #
p_3 = h_3_mag**2/mu_e #
#######################
#
################
# Eccentricity #
############################################
# 1st eccentricity vector                  #
e_1 = np.cross(v_1,h_1)/mu_e - r_1/r_1_mag #
# 1st eccentricity vector magnitude        #
e_1_mag = np.sqrt((e_1*e_1).sum())         #
# 2nd eccentricity vector                  #
e_2 = np.cross(v_2,h_2)/mu_e - r_2/r_2_mag #
# 1st eccentricity vector magnitude        #
e_2_mag = np.sqrt((e_2*e_2).sum())         #
# 3rd eccentricity vector                  #
e_3 = np.cross(v_3,h_3)/mu_e - r_3/r_3_mag #
# 3rd eccentricity vector magnitude        #
e_3_mag = np.sqrt((e_3*e_3).sum())         #
# Averaging the eccentricity               #
#                                          #
############################################
#
###################
# Semi-Major Axis #
############################
# 1st semi-major axis      #
a_1 = p_1/(1 - e_1_mag**2) #
# 2nd semi-major axis      #
a_2 = p_2/(1 - e_2_mag**2) #
# 3rd semi-major axis      #
a_3 = p_3/(1 - e_3_mag**2) #
############################
#
#################
# Normal Vector #
####################################
# Sets a normalizing vector        #
n_calc = np.array([0, 0, 1])       #
# 1st normal vector                #
n_1 = np.cross(n_calc,h_1)         #
# 1st normal vector magnitude      #
n_1_mag = np.sqrt((n_1*n_1).sum()) #
# 2nd normal vector                #
n_2 = np.cross(n_calc,h_2)         #
# 2nd normal vector magnitude      #
n_2_mag = np.sqrt((n_2*n_2).sum()) #
# 3rd normal vector                #
n_3 = np.cross(n_calc,h_3)         #
# 3rd normal vector magnitude      #
n_3_mag = np.sqrt((n_3*n_3).sum()) #
####################################
#
###############
# Inclination #
################################
# 1st inclination              #
i_1_calc = h_1[2]/h_1_mag      #
i_1      = np.arccos(i_1_calc) #
i_1_deg  = np.rad2deg(i_1)     #
# 2nd inclination              #
i_2_calc = h_2[2]/h_2_mag      #
i_2      = np.arccos(i_2_calc) #
i_2_deg  = np.rad2deg(i_2)     #
# 3rd inclination              #
i_3_calc = h_3[2]/h_3_mag      #
i_3      = np.arccos(i_3_calc) #
i_3_deg  = np.rad2deg(i_3)     #
################################
#
###############################
# Longitude of Ascending Node #
############################################
# 1st Longitude (Omega)                    #
Ome_1_calc  = n_1[0]/n_1_mag               #
Omega_1     = np.arccos(Ome_1_calc)        #
# Quadrant check                           #
if n_1[1] >= 0:                            #
    Omega_1_deg = np.rad2deg(Omega_1)      #
else:                                      #
    Omega_1_quad = 2*np.pi - Omega_1       #
    Omega_1_deg = np.rad2deg(Omega_1_quad) #
# 2nd Longitude (Omega) ####################
Ome_2_calc  = n_2[0]/n_2_mag               #
Omega_2     = np.arccos(Ome_2_calc)        #
# Quadrant check                           #
if n_2[1] >= 0:                            #
    Omega_2_deg = np.rad2deg(Omega_2)      #
else:                                      #
    Omega_2_quad = 2*np.pi - Omega_2       #
    Omega_2_deg = np.rad2deg(Omega_2_quad) #
# 3rd Longitude (Omega)#####################
Ome_3_calc  = n_3[0]/n_3_mag               #
Omega_3     = np.arccos(Ome_3_calc)        #
# Quadrant check                           #
if n_3[1] >= 0:                            #
    Omega_3_deg = np.rad2deg(Omega_3)      #
else:                                      #
    Omega_3_quad = 2*np.pi - Omega_3       #
    Omega_3_deg = np.rad2deg(Omega_3_quad) #
############################################
#
##################################
# Argumaent of Periapsis (omega) #
#################################################
# 1st Argumaent of Periapsis (omega)            #
ome_1_calc  = np.dot(n_1,e_1)/(n_1_mag*e_1_mag) #
omega_1     = np.arccos(ome_1_calc)             #
# Quadrant check                                #
if e_1[2] >= 0:                                 #
    omega_1_deg = np.rad2deg(omega_1)           #
else:                                           #
    omega_1_quad = 2*np.pi - omega_1            #
    omega_1_deg = np.rad2deg(omega_1_quad)      #
# 2nd Argumaent of Periapsis (omega)#############
ome_2_calc  = np.dot(n_2,e_2)/(n_2_mag*e_2_mag) #
omega_2     = np.arccos(ome_2_calc)             #
omega_2_deg = np.rad2deg(omega_2)               #
# Quadrant check                                #
if e_2[2] >= 0:                                 #
    omega_2_deg = np.rad2deg(omega_2)           #
else:                                           #
    omega_2_quad = 2*np.pi - omega_2            #
    omega_2_deg = np.rad2deg(omega_2_quad)      #
# 3rd Argumaent of Periapsis (omega)#############
ome_3_calc  = np.dot(n_3,e_3)/(n_3_mag*e_3_mag) #
omega_3     = np.arccos(ome_3_calc)             #
omega_3_deg = np.rad2deg(omega_3)               #
# Quadrant check                                #
if e_3[2] >= 0:                                 #
    omega_3_deg = np.rad2deg(omega_3)           #
else:                                           #
    omega_3_quad = 2*np.pi - omega_3            #
    omega_3_deg = np.rad2deg(omega_3_quad)      #
#################################################
#
#####################
# True Anomaly (nu) #
###############################################
# 1st true anomaly                            #
nu_1_calc = np.dot(e_1,r_1)/(e_1_mag*r_1_mag) #
nu_1      = np.arccos(nu_1_calc)              #
# Quadrant check                              #
if np.dot(r_1,v_1) >= 0:                      #             
    nu_1_deg  = np.rad2deg(nu_1)              #
else:                                         #
    nu_1_quad = 2*np.pi - nu_1                #
    nu_1_deg = np.rad2deg(nu_1_quad)          #
# 2nd true anomaly ############################
nu_2_calc = np.dot(e_2,r_2)/(e_2_mag*r_2_mag) #
nu_2      = np.arccos(nu_2_calc)              #
nu_2_deg  = np.rad2deg(nu_2)                  #
# Quadrant check                              #
if np.dot(r_2,v_2) >= 0:                      #
   nu_2_deg  = np.rad2deg(nu_2)               #
else:                                         #
    nu_2_quad = 2*np.pi - nu_2                #
    nu_2_deg = np.rad2deg(nu_2_quad)          #
# 3rd true anomaly ############################
nu_3_calc = np.dot(e_3,r_3)/(e_3_mag*r_3_mag) #
nu_3      = np.arccos(nu_3_calc)              #
nu_3_deg  = np.rad2deg(nu_3)                  #
# Quadrant check                              #
if np.dot(r_3,v_3)  >= 0:                     #
   nu_3_deg  = np.rad2deg(nu_3)               #
else:                                         #
    nu_3_quad = 2*np.pi - nu_3                #
    nu_3_deg = np.rad2deg(nu_3_quad)          #
###############################################
#
############################
# Eccentricity Anomaly (E) #
##################################################################
# 1st eccentricity anomaly                                       #
E_1_calc = (e_1_mag + np.cos(nu_1))/(1.0 + e_1_mag*np.cos(nu_1)) #
E_1      = np.arccos(E_1_calc)                                   #   
E_1_deg  = np.rad2deg(E_1)                                       #
# 2nd eccentricity anomaly                                       #
E_2_calc = (e_2_mag + np.cos(nu_2))/(1.0 + e_2_mag*np.cos(nu_2)) #
E_2      = np.arccos(E_2_calc)                                   #   
E_2_deg  = np.rad2deg(E_2)                                       #
# 3rd eccentricity anomaly                                       #
E_3_calc = (e_3_mag + np.cos(nu_3))/(1.0 + e_3_mag*np.cos(nu_3)) #
E_3      = np.arccos(E_3_calc)                                   #   
E_3_deg  = np.rad2deg(E_3)                                       #
##################################################################
#
####################
# Mean Anomaly (M) #
#####################################
# 1st mean anomaly                  #
M_1     = E_1 - e_1_mag*np.sin(E_1) #
M_1_deg = np.rad2deg(M_1)           #
# 2nd mean anomaly                  #
M_2     = E_2 - e_2_mag*np.sin(E_2) #
M_2_deg = np.rad2deg(M_2)           #
# 3rd mean anomaly                  #
M_3     = E_3 - e_3_mag*np.sin(E_3) #
M_3_deg = np.rad2deg(M_3)           #
#####################################
Orbital_Elements_report = f"""
{'-'*42}
| The orbital elements were calcualted with
| each pair of position and velocity vectors
{'-'*42}
| The Angular Momentums are:
|  h_1 = {h_1_mag} (km^2/s)
|  h_2 = {h_2_mag} (km^2/s)
|  h_3 = {h_3_mag} (km^2/s)
{'-'*42}
| The Semi-Latus Rectums are:
|  p_1 = {p_1} (km)
|  p_2 = {p_2} (km)
|  p_3 = {p_3} (km)
{'-'*42}
| The Eccentricities are:
|  e_1 = {e_1_mag}
|  e_2 = {e_2_mag}
|  e_3 = {e_3_mag}
{'-'*42}
| The Semi-Major Axii are:
|  a_1 = {a_1} (km)
|  a_2 = {a_2} (km)
|  a_3 = {a_3} (km)
{'-'*42}
| The Inclinations are:
|  i_1 = {i_1_deg} (degrees)
|  i_2 = {i_2_deg} (degrees)
|  i_3 = {i_3_deg} (degrees)
{'-'*42}
| The Longitude of Ascending Nodes are:
|  Omega_1 = {Omega_1_deg} (degrees)
|  Omega_2 = {Omega_2_deg} (degrees)
|  Omega_3 = {Omega_3_deg} (degrees)
|{'-'*42}
| The Argumaent of Periapsis are:
|  omega_1 = {omega_1_deg} (degrees)
|  omega_2 = {omega_2_deg} (degrees)
|  omega_3 = {omega_3_deg} (degrees)
|{'-'*42}
| The True Anomalies are:
|  nu_1 = {nu_1_deg} (degrees)
|  nu_2 = {nu_2_deg} (degrees)
|  nu_3 = {nu_3_deg} (degrees)
|{'-'*42}
| The Eccentricity Anomalies are:
|  E_1 = {E_1_deg} (degrees)
|  E_2 = {E_2_deg} (degrees)
|  E_3 = {E_3_deg} (degrees)
|{'-'*42}
| The Mean Anomalies are:
|  M_1 = {M_1_deg} (degrees)
|  M_2 = {M_2_deg} (degrees)
|  M_3 = {M_3_deg} (degrees)
|{'-'*42}
"""
print(Orbital_Elements_report)































