#########################################
# Author: Evan A. Blosser               #
# Email: evan.a.blosser-1@ou.edu        #
#                                       #
# Program: Mission to Titan             #
#            - X-fer calculations and   #
#              charting of a patched    #
#               conic trajectory        # 
#########################################

############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
from sys import exit                #
# Plotting & Animation              #
#####################################
#
#############
# Constants #
##############################
# Grav. Parameter for Earth  #  
mu_e = 3.986004418e5          #
# Equitorial Radius of the Earth in km #
r_e = 6378.136               #
##############################
# Earth Smi-major axis       #
R1 = 149.6e6                 #
# Jupiter Smi-major axis     #
R2 = 778.6e6                 #
# GRav. Perameter of the Sun #
mu_s = 1333                  #
# Accel. of Gravity (km/s^2) #
g_0 = 9.81e-3                #
##############################
# Replace with menu
mu_1 = mu_e
#
################
## User Input ##
################################################################
I_sp    = float(input("|Enter to Specific Impulse in Sec.: ")) #
Alt_loc = float(input('|Enter your current Altitde: '))        #
################################################################
#
######################
# Altitude to Radius #
r_p = Alt_loc + r_e  #
######################
#
###########################
# Hyperbolic Excess Speed #
#########################################################################
if R1 < R2:                                                             #
#  Inner to Outer Planet                                                #
    v_hyper_depart = np.sqrt((mu_s/R1)*(np.sqrt((2*R2)/(R1 + R2)) - 1)) #
#########################################################################
else:                                                                   #
#  Outer to Inner Planet                                                #
    v_hyper_depart = np.sqrt((mu_s/R1)*(1 - np.sqrt((2*R2)/(R1 + R2)))) #
#########################################################################

#######################
# Velocity of perigge #
##########################################################
v_per_depart = np.sqrt(v_hyper_depart**2 + (2*mu_1)/r_p) #
##########################################################
#
#####################################
# Current Veloctiy before departure #
#####################################
v_cur_depart = np.sqrt(mu_1/r_p) #
##################################
#
##################
# Delta-V Total # for inner to outer orbit
############################################################################
Delta_V = v_cur_depart*(np.sqrt(2 + (v_hyper_depart/v_cur_depart)**2) - 1) #
############################################################################
#
######################################
# Fuel Req. to Spacecraft mass ratio ##############
Fuel_Mass_Ratio = 1 - np.exp(-Delta_V/(I_sp*g_0)) #
###################################################
#


