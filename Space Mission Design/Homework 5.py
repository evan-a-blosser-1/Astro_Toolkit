#########################################
# Author: Evan A. Blosser               #
# Email: evan.a.blosser-1@ou.edu        #
#                                       #
# Program: Homework 5                   #
#           utilizing my interplanetary # 
#           x-fer calculator            # 
#########################################
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
from sys import exit                #
#####################################
#############
# Constants #
##############################
# Grav. Parameter for Earth  #  
mu_e = 3.986004418e5         #
# Equi.Radius of Earth in km #
r_e = 6378.136               #
##############################
# Earth Smi-major axis       #
RE = 149.6e6                 #
# Grav. P of Sun in km^3/s^2 #
mu_sun = 132712440018        #
##############################
##############################
# Earth Smi-major axis       #
R1 = 149.6e6                 #
# Saturn Smi-major axis      #
R2 = 1.433e9                 #
# Accel. of Gravity (km/s^2) #
g_0 = 9.81e-3                #
##############################
#%% Q1
##############
# Question 1 #
################################
# Earth's position on a day    #
r_1_Earth = 147.4*10**6        #
# The parking orbit given      #
#   alt 200km                  #
r_1_park = 200 + r_e           #
# Assign variables for the new #  
#   orbit aroudn the sun       #
r_1_per = 120*10**6            #
r_1_aph   = r_1_Earth          #
################################
#
#####################################
# Departure Hyperbolic Excess Speed #
###############################################################################
# Hyperbolic Excess Speed Planet Departure from parking orbit #################
v_1_hyper_dep = np.sqrt((mu_sun/r_1_per))*\
    (np.sqrt((2*r_1_aph)/(r_1_per + r_1_aph)) - 1)                            #
# Current Veloctiy at Perigee of Parking orbit ################################
v_1_current = np.sqrt(mu_e/r_1_park)                                          #
# Delta-V Total ###############################################################
Delta_V_1 = v_1_current*(np.sqrt(2 + (v_1_hyper_dep/v_1_current)**2) - 1)     #
Answer_Q_1 = f'''
{'-'*42}
|  Answer to Question 1:
{'-'*42}
|  The total Delta-V Required is:
|       {Delta_V_1} (km/s)
{'-'*42}
|   The Velocity of the departure 
|    hyperbola is: {v_1_hyper_dep} (km/s)
{'-'*42}
'''
print(Answer_Q_1)
#%% Q2
##############
# Question 2 #
##############################
# Mass of Space craft in kg  #
Craft_2_mass = 2000          #
# Parking Radius (alt=180 km)#
r_2_circ_park = 180 + r_e    #
# Specific Impulse in Sec.   #
I_2_sp = 300                 #
##############################
#############
# Departure #
#############
# Going from the parking orbit around earth into a hyperbolic 
#  departure trajectory. 
#########################################################################
# Hyperbolic Excess Speed Planet Departure from parking orbit           #
v_hyper_depart_2 = np.sqrt((mu_sun/R1))*(np.sqrt((2*R2)/(R1 + R2)) - 1) #
# Current Veloctiy at Perigee of Parking orbit ##########################
v_current_2 = np.sqrt(mu_e/r_2_circ_park)                               #
# Delta-V Total ##############################################################
Delta_V_2 = v_current_2*(np.sqrt(2 + (v_hyper_depart_2/v_current_2)**2) - 1) #
# Propellent Req. for the departure ##########################################
Fuel_Mass_Ratio_Depart_2 = 1 - np.exp(-Delta_V_2/(I_2_sp*g_0))               #
# Spacecraft Initial mass Equation                                           #
Mass_Prop_2 = (Fuel_Mass_Ratio_Depart_2*Craft_2_mass)/\
    (1 - Fuel_Mass_Ratio_Depart_2)                                           #
##############################################################################
####################
# Hohmann Transfer # Parts From Project
####################    over mission to titan
###################################################################
# Time of Filght from Earth to Saturn for Hohman                  #
TOF_Hohm_2_Sec = (np.pi)/(np.sqrt(mu_sun))*((R1 + R2)/2)**(3/2)   #
# Time of Flight in years                                         #
TOF_Hohm_2 = TOF_Hohm_2_Sec/(60*60*24*365)                        #
###################################################################
Answer_Q_2 =f'''
{'-'*42}
|  Answer to Question 2:
{'-'*42}
| The Mass of Requird Propellent is:
|    {Mass_Prop_2} (kg) 
| For the given Isp:
|    {I_2_sp} (seconds)
{'-'*42}
| The time requird for the mission is:
|    {TOF_Hohm_2} Years
{'-'*42}
'''
print(Answer_Q_2)
#%% Q3
##############
# Question 3 #
############################
# Perigee Radius km        #
r_3_p = 7000.0             #
# Velocity at perigee km/s #
v_3_p = 9.0                #
# Part A)                  #
Del_3 = 1                  #
# and B)                   #
Del_3_vel = .001
############################
#############################
# Calculate original Apogee #
############################################
# Angular momentum after change in Radius  #
h_3 = r_3_p*v_3_p                          #
# Semi-Latus Rectum after change in Radius #
semi_lat_3 = (h_3**2)/ mu_e                #
# Eccentricity after change in Radius      #
Ecc_3 = (h_3**2)/(mu_e*r_3_p)-1            #
# Apogee Radius after change in Radius     #
r_3_ap_origin = semi_lat_3/(1 - Ecc_3)     #
############################################
##########
# Part A #
############################################
# Change Radius by 1 km                    #
new_r_3_p = r_3_p + Del_3                  #
# Angular momentum after change in Radius  #
h_3_A = new_r_3_p*v_3_p                    #
# Semi-Latus Rectum after change in Radius #
semi_lat_3_A = (h_3_A**2)/ mu_e            #
# Eccentricity after change in Radius      #
Ecc_3_A = (h_3_A**2)/(mu_e*new_r_3_p)-1    #
# Apogee Radius after change in Radius     #
r_3_ap_A = semi_lat_3_A/(1 - Ecc_3_A)      #
# Change in Apogee                         #
Change_3_A = r_3_ap_A - r_3_ap_origin      #
############################################
##########
# Part B #
############################################
# Change Velocity by 1 km                  #
new_v_3_p = v_3_p + Del_3_vel              #
# Angular momen. after change in Velocity  #
h_3_B = r_3_p*new_v_3_p                    #
# Semi-Latus Rec. after change in Velocity #
semi_lat_3_B = (h_3_B**2)/ mu_e            #
# Eccentricity after change in Velocity    #
Ecc_3_B = (h_3_B**2)/(mu_e*r_3_p)-1        #
# Apogee Radius after change in Velocity   #
r_3_ap_B = semi_lat_3_B/(1 - Ecc_3_B)      #
# Change in Apogee                         #
Change_3_B = r_3_ap_B - r_3_ap_origin      #
############################################
Answer_Q_3 = f'''
{'-'*42}
|  Answer to Question 3:
{'-'*42}
| The Apogee Radius after a 1 km change 
|    in the Perigee Radius is:
|    {r_3_ap_A} (km)
{'-'*42}
| The change in Apogee Radius for Part A is:
|    {Change_3_A}
{'-'*42}
| The Apogee Radius after a 1 km/s change 
|    in the Perigee Velocity is:
|    {r_3_ap_B} (km)
{'-'*42}
| The change in Apogee Radius for Part B is:
|    {Change_3_B}
{'-'*42}
'''
print(Answer_Q_3)
#%% Q4
##############
# Question 4 #
#############################
# Earth's Semi_Major Axis   #
r_E_4  = 149.6e6            #
# Mercury's Semi_Major Axis #
r_ME_4 = 57.91e6            #
# Mercury's Radius          #
r_ME   = 2440               #
# Mercury's Grav. Per.      #
mu_ME  = 22032              #
# Circular parking altitude #
#  for Earth and Mercury    #      
r_4_park_alt = 150          #
#############################
############
# Transfer #
############
####################################################################
# Heliocentric velocity of craft on departure from earth           #
v_4_helio_craft = np.sqrt(mu_sun*((2/r_E_4)-(2/(r_E_4+r_ME_4))))   #
# Earth's Heliocentric Velocity                                    #
v_4_helio_earth = np.sqrt(mu_sun/r_E_4)                            #
# The velocity of the departure hyperbola                          #
v_inf_4_dep = v_4_helio_earth - v_4_helio_craft                    #
# Geocentric Spacecraft velocity                                   #
v_4_geo = np.sqrt(mu_e/(r_4_park_alt+r_e))                         #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.          #
v_geo_4_dep = np.sqrt(v_inf_4_dep**2 + (2*mu_e)/(r_4_park_alt+r_e))#
# Departure Delta-V                                                #
Delta_V_4_dep = v_geo_4_dep - v_4_geo                              #
####################################################################
###########
# Arrival #
###########
############################################################################
# Heliocentric velocity of craft on departure from earth                   #
v_4_helio_craft_arrive = np.sqrt(mu_sun*((2/r_ME_4)-(2/(r_E_4+r_ME_4))))   #
# Earth's Heliocentric Velocity                                            #
v_4_helio_merc = np.sqrt(mu_sun/r_ME_4)                               ######
# The velocity of the departure hyperbola                             #
v_inf_4_ari = v_4_helio_craft_arrive - v_4_helio_merc                 #
# Geocentric Spacecraft velocity                                      #
v_4_geo_Rende = np.sqrt(mu_ME/(r_4_park_alt+r_ME))                    #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.             #
v_geo_4_ari = np.sqrt(v_inf_4_ari**2 + (2*mu_ME)/(r_4_park_alt+r_ME)) #
# Departure Delta-V                                                   #
Delta_V_4_ari = v_geo_4_ari - v_4_geo_Rende                           #
#######################################################################
######################
# Total Delta-V Req. ######################################
Delta_V_total_4 = abs(Delta_V_4_dep) + abs(Delta_V_4_ari) #
###########################################################
Answer_Q_4 = f'''
{'-'*42}
|  Answer to Question 4:
{'-'*42}
|    {Delta_V_total_4} (km/s)
{'-'*42}
'''
print(Answer_Q_4)
#%% Q5
##############
# Question 5 #
##############
##################################
# Fly by Altitude around Jupiter #
r_flyby_5 = 200000               #
##################################
# Earth's Semi_Major Axis     #
R_E_5 = 149.6e6               #
# Jupiter's Semi_Major Axis   #
R_J_5 = 778.6e6               #
# Jupiter's Grav. Per.        #
mu_j = 126686534              #
# Jupiter's Radius            #
r_j_5  = 71490                #
###########################################
# Assumed Hohmann transfer from Earth     #
#   - in order to define angular momentum #
###########################################
# Semi-Major Axis of Transfer       #
semi_major_5_hohm = (R_E_5+R_J_5)/2 # 
# Jupiter's Speed                   #
V_Jup_5 = np.sqrt(mu_sun/R_J_5)     #
#  Heliocentric Spacecraft Velocity ###########################
V_Craft_5 = np.sqrt(mu_sun*((2/R_J_5)-(1/semi_major_5_hohm))) #
# Excess Speed upon Jupiter Arrival                           #
V_inf_5 = V_Jup_5 - V_Craft_5                                 #
# Eccentricity of Fly By Entry                                #
e_5 = 1 + ((r_j_5 + r_flyby_5)*(V_inf_5)**2)/mu_j             #
# Turning Angle                                               #
turn_ang_rad = 2*np.arcsin(1/e_5)                             #
turn_ang = np.rad2deg(turn_ang_rad)                           #
# Delta-V Transfered to the spacecraft (!Use Rad!)            #
Delta_V_5 = 2*V_inf_5*np.sin(turn_ang_rad/2)                  #
###############################################################
###################################
# Calculate Exit Orbital Elements #
############################################
# Turning angle for the outbound direction ################## 
Outbound_ang = 180 + turn_ang              
Outbound_ang_rad = np.deg2rad(Outbound_ang)                 
# Outbound Velocity Vector
V_outbound_Uv = V_Jup_5 + V_inf_5*np.cos(Outbound_ang_rad)
V_outbound_Us = V_inf_5*np.sin(Outbound_ang_rad)
# Magnitude/Normal of the outbound Velocity
V_5_outbound = np.sqrt(V_outbound_Uv**2 + V_outbound_Us**2)
#########
# Angular Momentum at exit
#    -assuming we have made a new apogee/perigee
h_5   = R_J_5*V_outbound_Uv
######################
# Eccentricity is Calculated as shown in Lec. 20230410
#  were we use relations between the eccentricity and
#   true anomaly at exit
#####################################################
# e*cos relation                                    #
Ecos = h_5**2/(mu_sun*R_J_5) - 1                    #
# e*sin relation                                    #
Esin = (V_outbound_Us*h_5)/mu_sun                   #
# divide both relations to cancel out eccentricity  #
Tan_theta_5 = Esin/Ecos                             #
# Solve for True anomaly at exit                    #
Theta_5_2   =np.arctan(Tan_theta_5)                 #
# New Eccenmtricity calculated with Theta_2         #
Ecc_5 =  Ecos/np.cos(Theta_5_2)                     #
# New Perihelion of Spacecraft                      #
R_perihelion_5 = (h_5**2/mu_sun)*(1/(1 + Ecc_5))    #
# New Semi-Major Axis Calculations using new        #
#   eccentricity and angular momentumn              #
a_5 =  R_perihelion_5/(1 - Ecc_5)                   #
#####################################################
Answer_Q_5 = f'''
{'-'*42}
|  Answer to Question 5:
{'-'*42}
|    The Semi-Major axis of the post
|        flyby trajectory is:
|            {a_5} (km) 
{'-'*42}
|    The Eccentricity of the post 
|       flyby trajectory is:
|           {Ecc_5}   
{'-'*42}
|    The Delta-V imparted on the
|       spacecraft is:
|         {Delta_V_5} (km/s)
{'-'*42}
'''
print(Answer_Q_5)


