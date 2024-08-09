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
##############################
# Accel. of Gravity (km/s^2) #
g_0 = 9.81e-3                #
##############################
#
##############################
# Grav. Parameter (km^3/s^2) #
##############################
# Sun                        #
mu_sun = 132712440018        #
# Mercury                    #
mu_me  = 22032               #
# Venus                      #
mu_v   = 324859              #
# Earth                      #
mu_e   = 3.986004418e5       #
# Mars                       #
mu_m   = 42828               #
# Jupiter/Titan              #
mu_j   = 126686534           #
# Saturn                     #
mu_sat = 37931187            #
# Uranus                     #
mu_u   = 5793939             #
# Neptune                    #
mu_n   = 6836529             #
##############################
######################
# Planet Radius (km) #
######################
# Mercury            #
r_me  = 2440         #
# Venus              #
r_v   = 6052         #
# Earth              #
r_e   = 6378.136     #
# Mars               #
r_m   = 3396         #
# Jupiter            #
r_j   = 71490        #
# Saturn             #
r_sat = 60270        #
# Uranus             #
r_u   = 25560        #
# Neptune            #
r_n   = 24764        # 
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
#########################
# Orbit sidreal period  #
#        (days)         #
#########################
# Mercury               #
T_mercury = 87.97       #
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
|
| Note: Mercury was included as a test
|        with known values and results
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
##############################
# Mass of Space craft in kg  #
Craft_mass = 2000            #
# Specific Impulse in Sec.   #
I_sp = 300                   #
##############################
r_i_park_alt = 150
r_f_park_alt = 150
Initial_loc  = input('| Where are you at?: ')  #
Destination  = input('| Where to?: ')          #
################################################
# Location Selection #
###############################################
# Set Mercury  to be the initial radius       #
if Initial_loc == 'ME' or Initial_loc =='me': #
    r_i   = r_mercury                         #
    Rad_i = r_me                              #
    mu_i  = mu_me                             #
    T_i   = T_mercury                         #
    i_sel = 'Mercury'                         #
# Set Venus to be the initial radius          #
elif Initial_loc == 'V' or Initial_loc =='v': #
    r_i   = r_venus                           #
    Rad_i = r_v                               #
    mu_i  = mu_v                              #
    T_i   = T_venus                           #
    i_sel = 'Venus'                           #
# Set Earth to be the initial radius          #
elif Initial_loc =='E' or Initial_loc =='e':  #
    r_i   = r_earth                           #
    Rad_i = r_e                               #
    mu_i  = mu_e                              #
    T_i   = T_earth                           #
    i_sel = 'Earth'                           #
# Set Mars to be the initial radius           #
elif Initial_loc =='M' or Initial_loc =='m':  #
    r_i   = r_mars                            #
    Rad_i = r_m                               #
    mu_i  = mu_m                              #
    T_i   = T_mars                            #
    i_sel = 'Mars'                            #
# Set Jupiter to be the initial radius        #
elif Initial_loc =='J' or Initial_loc =='j':  #
    r_i   = r_jupiter                         #
    Rad_i = r_j                               #
    mu_i  = mu_j                              #
    T_i   = T_jupiter                         #
    i_sel = 'Jupiter'                         #
# Set Saturn to be the initial radius         #
elif Initial_loc =='S' or Initial_loc =='s':  #
    r_i   = r_saturn                          #
    Rad_i = r_sat                             #
    mu_i  = mu_sat                            #
    T_i   = T_saturn                          #
    i_sel = 'Saturn'                          #
# Set Uranus to be the initial radius         #
elif Initial_loc =='U' or Initial_loc =='u':  #
    r_i   = r_uranus                          #
    Rad_i = r_u                               #
    mu_i  = mu_u                              #
    T_i   = T_uranus                          #
    i_sel = 'Uranus'                          #
# Set Neptune to be the initial radius        #
elif Initial_loc =='N' or Initial_loc =='n':  #
    r_i   = r_neptune                         #
    Rad_i = r_n                               #
    mu_i  = mu_n                              #
    T_i   = T_neptune                         #
    i_sel = 'Neptune'                         #
else:                                         #
    print(Error)                              #
    print(f"| User input was: {Initial_loc}") #
    input("| Please Press Enter To Exit...")  #
    exit()                                    #
###############################################
#########################
# Destination Selection #
###############################################
# Set Mercury to be the final radius          #
if Destination == 'ME' or Destination =='me': # 
    r_f   = r_mercury                         # 
    Rad_f = r_me                              #
    mu_f  = mu_me                             #    
    T_f   = T_mercury                         #
    f_sel = 'Mercury'                         #
# Set Venus to be the final radius            #
elif Destination == 'V' or Destination =='v': # 
    r_f   = r_venus                           # 
    Rad_f = r_v                               #
    mu_f  = mu_v                              #  
    T_f   = T_venus                           #
    f_sel = 'Venus'                           #
# Set Earth to be the final radius            #
elif Destination =='E' or Destination =='e':  #
    r_f   = r_earth                           #
    Rad_f = r_e                               #
    mu_f  = mu_e                              #  
    T_f   = T_earth                           #
    f_sel = 'Earth'                           #
# Set Mars to be the final radius             #
elif Destination == 'M' or Destination =='m': #
    r_f   = r_mars                            #
    Rad_f = r_m                               #
    mu_f  = mu_m                              #  
    T_f   = T_mars                            #
    f_sel = 'Mars'                            #
# Set Jupiter to be the final radius          #
elif Destination == 'J' or Destination =='j': #
    r_f   = r_jupiter                         #
    Rad_f = r_j                               #
    mu_f  = mu_j                              #  
    T_f   = T_jupiter                         #
    f_sel = 'Jupiter'                         #
# Set Saturn to be the final radius           #
elif Destination == 'S' or Destination =='s': #
    r_f   = r_saturn                          #
    Rad_f = r_sat                             #
    mu_f  = mu_sat                            #  
    T_f   = T_saturn                          #
    f_sel = 'Saturn'                          #
# Set Uranus  to be the final radius          #
elif Destination == 'U' or Destination =='u': #
    r_f   = r_uranus                          #
    Rad_f = r_u                               #
    mu_f  = mu_u                              #  
    T_f   = T_uranus                          #
    f_sel = 'Uranus'                          #
# Set Neptune to be the final radius          #
elif Destination == 'N' or Destination =='n': #
    r_f   = r_neptune                         #
    Rad_f = r_n                               #
    mu_f  = mu_n                              #  
    T_f   = T_neptune                         #
    f_sel = 'Neptune'                         #
else:                                         #
    print(Error)                              #
    print(f"| User input was: {Destination}") #
    input("| Please Press Enter To Exit...")  #
    exit()                                    #
###############################################
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
############
# Transfer #
############
####################################################################
# Heliocentric velocity of craft on departure from earth           #
v_helio_craft = np.sqrt(mu_sun*((2/r_i)-(2/(r_i+r_f))))            #
# Earth's Heliocentric Velocity                                    #
v_helio_planet = np.sqrt(mu_sun/r_i)                               #
# The velocity of the departure hyperbola                          #
v_inf_dep =  v_helio_craft - v_helio_planet                        #
# Geocentric Spacecraft velocity                                   #
v_geo = np.sqrt(mu_i/(r_i_park_alt+Rad_i))                         #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.          #
v_geo_dep = np.sqrt(v_inf_dep**2 + (2*mu_i)/(r_i_park_alt+Rad_i))  #
# Departure Delta-V                                                #
Delta_V_dep = v_geo_dep - v_geo                                    #
####################################################################
###########
# Arrival #
###########
#######################################################################
# Heliocentric velocity of craft on departure from earth              #
v_helio_craft_arrive = np.sqrt(mu_sun*((2/r_f)-(2/(r_i+r_f))))        #
# Earth's Heliocentric Velocity                                       #
v_helio_arrive_planet = np.sqrt(mu_sun/r_f)                           #
# The velocity of the departure hyperbola                             #
v_inf_ari = v_helio_craft_arrive - v_helio_arrive_planet              #
# Geocentric Spacecraft velocity                                      #
v_geo_Rende = np.sqrt(mu_f/(r_f_park_alt+Rad_f))                      #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.             #
v_geo_ari = np.sqrt(v_inf_ari**2 + (2*mu_f)/(r_f_park_alt+Rad_f))     #
# Departure Delta-V                                                   #
Delta_V_ari = v_geo_ari - v_geo_Rende                                 #
#######################################################################
###################################################################
# Time of Filght from Earth to Saturn for Hohman                  #
TOF_Hohm_Sec = (np.pi)/(np.sqrt(mu_sun))*((r_i + r_f)/2)**(3/2)   #
# Time of Flight in years                                         #
TOF_Hohm = TOF_Hohm_Sec/(60*60*24*365)                            #
###################################################################
######################
# Total Delta-V Req. ################################
Delta_V_total = abs(Delta_V_dep) + abs(Delta_V_ari) #
#####################################################
# Propellent Req. for the departure ##############################
Fuel_Mass_Ratio = 1 - np.exp(-Delta_V_total/(I_sp*g_0))          #
# Spacecraft Initial mass Equation                               #
Mass_Prop = (Fuel_Mass_Ratio*Craft_mass)/(1 - Fuel_Mass_Ratio)   #
##################################################################
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
| The Time of Flight is:
|     {TOF_Hohm} years
{'-'*42}
| The required fuel mass ratio is
|   {Fuel_Mass_Ratio} x100%
{'-'*42}
| The mass of propellant required for
|   the spacecraft weighing 
|   {Craft_mass} kg is:
|    {Mass_Prop} kg
{'-'*42}
"""
print(Data_Output)
input("| Please Press Enter To Exit...")  



































