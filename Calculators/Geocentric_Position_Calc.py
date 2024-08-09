##############################################
# Author: Evan A. Blosser                    #
# Email: evan.a.blosser-1@ou.edu             #
#                                            #
# Program: Using orbital elements to plot    #
#           an orbit around earth            #
#                                            # 
##############################################
#%% Packages & Constants
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
#####################################
#########################
# Constants Declaration #
##############################################################
# Gravitational Parameter for Earth in km^3/s^2              #
# [If working in Canonical Units replace with 1 (DU^3/TU^2)] #
mu_e = 3.986004418e5                                         #
# Equitorial Radius of the Earth in km                       #
# [If working in Canonical Units replace with 1 (DU)]        #
r_e = 6378.136                                               #
##############################################################
############
### Main ###
############
welcome =f"""
{'-'*42}
            Welcome to the 
Geocentric Equtorial Position Calcualtor
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
| Calculating orbital elements from a given position
|    and velocity vector in the Geocentric Frame 
{'-'*42}
|  Warning: The position vectors should be input as
|  each i, j, & k direction seperatly when prompted. 
{'-'*42}
"""                              
print(welcome)   
###############
# User Inputs #
################################
# Input Position Vectors Below #
#########################################################
# Input & Assignment of First Position Vector           #
r_i = float(input("|Enter the position Vector, (i): ")) #
r_j = float(input("|Enter the position Vector, (j): ")) #
r_k = float(input("|Enter the position Vector, (k): ")) #
r_vec  = np.array([r_i, r_j, r_k])                      #
# Input & Assignment of Second Position Vector          #
v_i = float(input("|Enter the velocity Vector, (i): ")) #
v_j = float(input("|Enter the velocity Vector, (j): ")) #
v_k = float(input("|Enter the velocity Vector, (k): ")) #
v_vec  = np.array([v_i, v_j, v_k])                      #
#########################################################
'''
### Test Values 
r_vec  = np.array([5887.0, -3520.0, -1204.0])
v_vec  = np.array([-1.420789549655011, 0.06108240648049763, -7.469331602357601]) 
'''
#####################
# Start Calculating #
#####################
r_mag = np.linalg.norm(r_vec)
v_mag = np.linalg.norm(v_vec)
####################
# Angular Momentum #
####################################
# Angular momentum of 1st position #
h = np.cross(r_vec,v_vec)          #
h_mag = np.linalg.norm(h)          #
####################################
#
#####################
# Semi-Latus Rectum # 
#######################
# semi-latus rec      #
p = h_mag**2/mu_e     #
#######################
#
################
# Eccentricity #
############################################
# eccentricity vector                      #
e = np.cross(v_vec,h)/mu_e - r_vec/r_mag   #
# eccentricity vector magnitude            #
e_mag = np.sqrt((e*e).sum())               #
############################################
#
###################
# Semi-Major Axis #
##################################
a = h_mag**2/mu_e/(1 - e_mag**2) #
##################################
#
#################
# Normal Vector #
######################################
# Sets a normalizing vector          #
n_calc = np.array([0, 0, 1])         #
# Normal vector                      #
n_vec = np.cross(n_calc,h)           #
# Normal vector magnitude            #
n_mag = np.sqrt((n_vec*n_vec).sum()) #
######################################
#
###############
# Inclination #
############################
# 1st inclination          #
i_calc = h[2]/h_mag        #
i      = np.arccos(i_calc) #
i_deg  = np.rad2deg(i)     #
############################
#
###############################
# Longitude of Ascending Node #
########################################
# 1st Longitude (Omega)                #
Ome_calc  = n_vec[0]/n_mag             #
Omega     = np.arccos(Ome_calc)        #
# Quadrant check                       #
if n_vec[1] > 0:                       #
    Omega_deg = np.rad2deg(Omega)      #
else:                                  #
    Omega_quad = 2*np.pi - Omega       #
    Omega_deg = np.rad2deg(Omega_quad) #
########################################
#
##################################
# Argumaent of Periapsis (omega) #
############################################
# 1st Argumaent of Periapsis (omega)       #
ome_calc  = np.dot(n_vec,e)/(n_mag*e_mag)  #
omega     = np.arccos(ome_calc)            #
# Quadrant check                           #
if e[2] > 0:                               #
    omega_deg = np.rad2deg(omega)          #
else:                                      #
    omega_quad = 2*np.pi - omega           #
    omega_deg = np.rad2deg(omega_quad)     #
############################################
#
#####################
# True Anomaly (nu) #
###########################################
# 1st true anomaly                        #
nu_calc = np.dot(e,r_vec)/(e_mag*r_mag)   #
nu      = np.arccos(nu_calc)              #
# Quadrant check                          #
if np.dot(r_vec,v_vec) > 0:               #             
    nu_deg  = np.rad2deg(nu)              #
else:                                     #
    nu_quad = 2*np.pi - nu                #
    nu_deg = np.rad2deg(nu_quad)          #
###########################################
#
Orbital_Elements_report = f"""
{'-'*42}
| The orbital elements were calcualted with
| each pair of position and velocity vectors
{'-'*42}
| The Angular Momentum is:
|  h = {h_mag} (km^2/s)
{'-'*42}
| The Semi-Latus Rectum is:
|  p = {p} (km
{'-'*42}
| The Eccentricity is:
|  e = {e_mag}
{'-'*42}
| The Semi-Major Axis is:
|  a = {a} (km)
{'-'*42}
| The Inclination is:
|  i = {i_deg} (degrees)
{'-'*42}
| The Longitude of Ascending Node is:
|  Omega = {Omega_deg} (degrees)
|{'-'*42}
| The Argumaent of Periapsis is:
|  omega = {omega_deg} (degrees)
|{'-'*42}
| The True Anomaly is:
|  nu = {nu_deg} (degrees)
|{'-'*42}

"""
print(Orbital_Elements_report)
input("| Press Enter To Exit...")