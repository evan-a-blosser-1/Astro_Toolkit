#########################################
# Author: Evan A. Blosser               #
# Email: evan.a.blosser-1@ou.edu        #
#                                       #
# Program: Sphere of infulance Calc.    # 
#########################################
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
# System Error/Exit                 #
from sys import exit                #
#####################################
#########################
# Constants Declaration #
###################################################
# Gravitational Parameter for the Sun in km^3/s^2 #
mu_sun = 132712440018                             #
###################################################
#########################
# Semi-Major axii       #
#(circular orbit approx)#
#     a = r_c (km)      #
#########################
# Mercury               #
r_mercury = 57.91e6     #
# Venus                 #
r_venus = 108.2e6       #
# Earth                 #
r_earth = 149.6e6       #
# Mars                  #
r_mars  = 227.9e6       #
# Jupiter/Titan?        #
r_jupiter = 778.6e6     #
# Saturn                #
r_saturn = 1.433e9      #
# Uranus                #
r_uranus = 2.872e9      #
# Neptune               #
r_neptune = 4.495e9     # 
#########################
#######################
# Planet Mass (kg)    #
#######################
# The Sun             #
M_Sun     = 1.989e30  #
# Mercury             #
M_mercury = 330.2e21  #
# Venus               #
M_venus   = 4.869e24  #
# Earth               #
M_earth   = 5.974e24  #
# Mars                #
M_mars    =  641.9e21 #
# Jupiter             #
M_jupiter = 1.899e27  #
# Saturn              #
M_saturn  = 568.5e24  #
# Uranus              #
M_uranus  = 86.83e24  #
# Neptune             #
M_neptune = 102.4e24  # 
#######################
####################
# Eccentricity     #
####################
# Mercury          #
e_mercury = 0.2056 #
# Venus            #
e_venus   = 0.0067 #
# Earth            #
e_earth   = 0.0167 #
# Mars             #
e_mars    = 0.0935 #
# Jupiter          #
e_jupiter = 0.0489 #
# Saturn           #
e_saturn  = 0.0565 #
# Uranus           #
e_uranus  = 0.0457 #
# Neptune          #
e_neptune = 0.0113 # 
####################
############
### Main ###
############
#
##################
# Print Messages #
##################
welcome =f"""
{'-'*42}
            Welcome to the Sphere of
             Influence Calculator
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
| Calculates the Sphere of Influence and 
| the Hill's Raduis for Stability Zone 
{'-'*42}
| User Inputs: enter the letter of the Planet
| when prompted. Not case sensitive.
| 
| Mercury is Me or me
"""                              
Error =f"""
{'-'*42}
| User Input Error!!
| Please Input: Me for Mercury
| V for Venus, E for Earth
| M for Mars, J for Jupiter,
| S for Saturn, U for Uranus,
| or N for Neptune
{'-'*42}
| Please Restart Program
{'-'*42}
"""
Error_silly =f"""
{'-'*42}
| Your already there silly! :)
{'-'*42}
| Please Restart Program
{'-'*42}
"""
print(welcome)      
###############
# User Inputs #
###########################################
Planet  = input('| Where are you at?: ')  #
###########################################
# Location Selection #
######################################
# Set Mercury to be the target mass  #
if Planet == 'Me' or Planet =='me':  #
    M_p = M_mercury                  #
    R   = r_mercury                  #
    e   = e_mercury                  #
    sel = 'Mercury'                  #
# Set Venus to be the target mass    #
elif Planet == 'V' or Planet =='v':  #
    M_p = M_venus                    #
    R   = r_venus                    #
    e   = e_venus                    #
    sel = 'Venus'                    #
# Set Earth to be the target mass    #
elif Planet =='E' or Planet =='e':   #
    M_p = M_earth                    #
    R   = r_earth                    #
    e   = e_earth                    #
    sel = 'Earth'                    #
# Set Mars to be the target mass     #
elif Planet =='M' or Planet =='m':   #
    M_p = M_mars                     #
    R   = r_mars                     #
    e   = e_mars                     #
    sel = 'Mars'                     #
# Set Jupiter to be the target mass  #
elif Planet =='J' or Planet =='j':   #
    M_p = M_jupiter                  #
    R   = r_jupiter                  #
    e   = e_jupiter                  #
    sel = 'Jupiter'                  #
# Set Saturn to be the target mass   #
elif Planet =='S' or Planet =='s':   #
    M_p = M_saturn                   #
    R   = r_saturn                   #
    e   = e_saturn                   #
    sel = 'Saturn'                   #
# Set Uranus to be the target mass   #
elif Planet =='U' or Planet =='u':   #
    M_p = M_uranus                   #
    R   = r_uranus                   #
    e   = e_uranus                   #
    sel = 'Uranus'                   #
# Set Neptune to be the target mass #
elif Planet =='N' or Planet =='n': #
    M_p = M_neptune               #
    R   = r_neptune              #
    e   = e_neptune             #
    sel = 'Neptune'            #
else:                         #
    exit(Error)              #
#############################
####################################
# Sphere of Influence Calculations #
r_sol = R*(M_p/M_Sun)**(2/5) #######
####################################
########################
# Hill's Radius/Sphere #################
R_h = R*(1 - e)*(M_p/(3*M_Sun))**(1/3) #
########################################
##################
# Stability Zone ##########
Stability_Pro = 0.529*R_h #
Stability_Ret = 0.693*R_h #
###########################
############
## Output ##
############
Data_Output= f"""
{'-'*42}
| The Sphere of influence 
|   of {sel} is:
|   {r_sol} km
{'-'*42}
| The Hills' Radius/Sphere
|  of {sel} is:
|   {R_h} km
{'-'*42}
| The Stability Zones
|   of {sel} are:
|
| Progade Orbit:
|   {Stability_Pro} km
|
| Retrograde Orbit:
|   {Stability_Ret} km
{'-'*42}
"""
print(Data_Output)
