#########################################
# Author: Evan A. Blosser               #
# Email: evan.a.blosser-1@ou.edu        #
#                                       #
# Program: Interplanetary Hohmann       #
#            X-fer Calculator           # 
#########################################
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
from sys import exit                #
#####################################
#########################
# Constants Declaration #
####################################################
# Gravitational Parameter for Earth in km^3/s^2   #
mu_e = 3.986004418e5                              #
# Gravitational Parameter for the Sun in km^3/s^2 #
mu_sun = 132712440018                             #
###################################################
#########################
# Semi-Major axii       #
#(circular orbit approx)#
#     a = r_c (km)      #
#########################
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
#########################
# Orbit sidreal period  #
#        (days)         #
#########################
# Venus                 #
T_venus = 224.7         #
# Earth                 #
T_earth = 365.256       #
# Mars                  #
T_mars  = 687.0456      #
# Jupiter/Titan?        #
T_jupiter = 4331.93616  #
# Saturn                #
T_saturn = 10760.44176  #
# Uranus                #
T_uranus = 30685.15656  #
# Neptune               #
T_neptune = 60194.1888  # 
#########################
############
### Main ###
############
#
##################
# Print Messages #
##################
welcome =f"""
{'-'*42}
            Welcome to the Interplanetary
             Hohmann Transfer Calculator
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
| Used as a good approximation for Delta-V 
| req needed for interplanetary travel: 
{'-'*42}
| Note: Mercury and Pluto are excluded as their 
| eccentricity is too large to be approximated 
| as zero (Hohmann X-fer will not work!)
{'-'*42}
| User Inputs: enter the letter of the Planet
| when prompted. Not case sensitive.
"""                              
Error =f"""
{'-'*42}
| User Input Error!!
| Please Input:
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
################################################
Initial_loc  = input('| Where are you at?: ')  #
Destination  = input('| Where to?: ')          #
################################################
# Location Selection #
#############################################
# Set Venus to be the initial radius        #
if Initial_loc == 'V' or Initial_loc =='v': #
    r_i = r_venus                           #
    T_i = T_venus                           #
    i_sel = 'Venus'                         #
# Set Earth to be the initial radius        #
elif Initial_loc =='E' or Initial_loc =='e':#
    r_i = r_earth                           #
    T_i = T_earth                           #
    i_sel = 'Earth'                         #
# Set Mars to be the initial radius         #
elif Initial_loc =='M' or Initial_loc =='m':#
    r_i = r_mars                            #
    T_i = T_mars                            #
    i_sel = 'Mars'                          #
# Set Jupiter to be the initial radius      #
elif Initial_loc =='J' or Initial_loc =='j':#
    r_i = r_jupiter                         #
    T_i = T_jupiter                         #
    i_sel = 'Jupiter'                       #
# Set Saturn to be the initial radius       #
elif Initial_loc =='S' or Initial_loc =='s':#
    r_i = r_saturn                          #
    T_i = T_saturn                          #
    i_sel = 'Saturn'                        #
# Set Uranus to be the initial radius       #
elif Initial_loc =='U' or Initial_loc =='u':#
    r_i = r_uranus                          #
    T_i = T_uranus                          #
    i_sel = 'Uranus'                        #
# Set Neptune to be the initial radius      #
elif Initial_loc =='N' or Initial_loc =='n':#
    r_i = r_neptune                         #
    T_i = T_neptune                        #
    i_sel = 'Neptune'                     #
else:                                    #
    exit(Error)                         #
########################################
#########################
# Destination Selection #
##############################################
# Set Venus to be the final radius           #
if Destination == 'V' or Destination =='v':  # 
    r_f = r_venus                            # 
    T_f = T_venus                            #
    f_sel = 'Venus'                          #
# Set Earth to be the final radius           #
elif Destination =='E' or Destination =='e': #
    r_f = r_earth                            #
    T_f = T_earth                            #
    f_sel = 'Earth'                          #
# Set Mars to be the final radius            #
elif Destination == 'M' or Destination =='m':#
    r_f = r_mars                             #
    T_f = T_mars                             #
    f_sel = 'Mars'                           #
# Set Jupiter to be the final radius         #
elif Destination == 'J' or Destination =='j':#
    r_f = r_jupiter                          #
    T_f = T_jupiter                          #
    f_sel = 'Jupiter'                        #
# Set Saturn to be the final radius          #
elif Destination == 'S' or Destination =='s':#
    r_f = r_saturn                           #
    T_f = T_saturn                           #
    f_sel = 'Saturn'                         #
# Set Uranus  to be the final radius         #
elif Destination == 'U' or Destination =='u':#
    r_f = r_uranus                           #
    T_f = T_uranus                           #
    f_sel = 'Uranus'                         #
# Set Neptune to be the final radius         #
elif Destination == 'N' or Destination =='n':#
    r_f = r_neptune                          #
    T_f = T_neptune                         #
    f_sel = 'Neptune'                      #
else:                                     #
    exit(Error)                          #
#########################################
##################
# Idiot Proofing #
########################
if r_i == r_f:         #
    print(Error_silly) #
########################
#
################
# Calculations #
################
####################
# Hohmann Transfer #
####################################################################
# Initial Detla-V                                                  #
Delta_V_i = np.sqrt(mu_sun/r_i)*(np.sqrt((2*r_f)/(r_i + r_f)) - 1) #
# Final Delta-V                                                    #
Delta_V_f = np.sqrt(mu_sun/r_f)*(1 - np.sqrt((2*r_i)/(r_i + r_f))) #
# Toatl Delta-V fro Hohmann Transfer                               #
Delta_V_total = abs(Delta_V_i) + abs(Delta_V_f)                    #
####################################################################
#
##################
# Synodic Period #
########################
# Synodic Period of    #
# location/destination ###########
T_syn = (T_i*T_f)/abs(T_i - T_f) #
##################################
############
## Output ##
############
##################
# Delta-V Output #
##################
Data_Output= f"""
{'-'*42}
| The Delta V requird to travel from
|  {i_sel} to {f_sel} is:
|   {Delta_V_total} km/s
{'-'*42}
| The synodic period of {f_sel}
|    with respect to {i_sel} is:
|    {T_syn} days
{'-'*42}
"""
print(Data_Output)
'''
for .exe 
input('Press any key to exit...')
'''


































