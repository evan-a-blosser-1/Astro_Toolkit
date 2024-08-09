########################################################
# Author: Evan A. Blosser                              #
# Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com #
#                                                      #
# Program: Lambert's_Velocity_Calculator               #
#           -Uses Lambert's Method to compute the      #
#            velocity of an orbiting body.             #
#           -Uses the velocities to calculate the      #
#            Classic Orbital Elements.                 #
########################################################
#
############
# Packages #
############
########################################
# Mathmatical!!                        #
import numpy as np                     #
# System Error/Exit                    #
from sys import exit                   #
# .exe output text color               #
from colorama import init, Fore, Style #
########################################
# Colorama Settings #
########################################
# Forces .exe to accept colorama       #
init()                                 #
# Increase Brightness of Text          #
print(Style.BRIGHT)                    #
########################################
#
###########################
## Constants Declaration ##
####################################
# Grav. Parameter Earth (km^3/s^2) #
mu_e = 3.986004418e5               #
# Grav. Parameter Sun in km^3/s^2  #
mu_sun = 132712440018              #
####################################
#
#######################
# Iteration Tolerance #
Tolerance = 1.0e-6    #
#######################
#
############
### Main ###
############
#
##########################################################
welcome =f"""
{'-'*42}
            Welcome to the 
Velocity & Orbital Element's Calcualtor
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
| Utilizing Lambert's Method, to calculate the velocity of
|  an orbiting body given Position Vectors 
|  (or Magnitudes and change in true anomaly) 
|      and passage of time .
|
{'-'*42}
|  Warning: The position vectors should be input as
|  each i, j, & k direction seperatly when prompted. 
{'-'*42}
"""                              
print(Fore.CYAN + welcome)       
##########################################################
#%% User Inputs 
#
#############################
# Define Mu_Earth or Mu_Sun #
####################################################################
# User Prompt to Input Selection for Earth or Sun                  #
Orbiting_Body_Earth_Helio = input(f"{'-'*42}\n"                    #
        "| Is this a Heliocentric Orbiting Body or Geocentric?: ") #
##############                                                     #
# Selecting Gravitational Parameter by User Input                  #
if Orbiting_Body_Earth_Helio in('h','helio','heliocentric',        #
                                'H','Helio','Heliocentric',        #
                                'HELIO','HELIOCENTRIC'):           #
# Set Gravitational Parameter to The Sun                           #
    mu = mu_sun                                                    #
################                                                   #
elif Orbiting_Body_Earth_Helio in('g','geo','geocentric',          #
                                  'G','Geo','Geocentric',          #
                                  'GEO','GEOCENTRIC'):             # 
# Set Gravitational Parameter to Earth                             #
    mu = mu_e                                                      #
##############                                                     #
else:                                                              #
### Error Message for Wrong Input ##################################
    Error_Message_Helio = f"""
{'-'*42}
| Error: The Gravitational parameter was not set!
|       '{Orbiting_Body_Earth_Helio}' is not an acceptable input
{'-'*42}
| IF: This is a Planet then Enter Heliocentric,
|     Helio, or Just H
|            
| IF: This is an Earth Satallite then select
|     Geocentric
{'-'*42}
    """
    print(Fore.RED + Error_Message_Helio)    ####################
    input(Fore.RED +'| Press Enter To Exit...')                 #
    exit()                                                      #
#################################################################
#
##############################################
# Selection for Position Vector Celculations #
####################################################################
#User Prompt to Input Select Position Vector Calculatoins          #
radius_vec_check = input(f"{'-'*42}\n"                             #  
    "| Do the position vector need to be calculated?(y/n): ")      #
if radius_vec_check in ('y','Y','yes','Yes','YES'):                #
    r_i_a      = float(input("| Enter the Initial Altitude: "))    #
    r_f_a      = float(input("| Enter the Final Altitude: "))      #
    nu_deg_inp = float(input("| Enter the True Anomaly Increase"   #
                           "(in degrees): "))                      #
    nu = np.deg2rad(nu_deg_inp)                                    #
    ################################################################
    #######################
    # Vector Calculations #
    #######################
    # Assuming that the first poistion
    #  is at zero true anomaly in the
    #  Perifocal System 
    #
    # Note: These Coefficiants are not dependant on the 
    #       reference system!
    ##################################
    # Calculating radius             #
    r_i = r_i_a + 6378.136           #
    r_f = r_f_a + 6378.136           #
    # Assign first altitude as vec   #
    r_1 = np.array([r_i,0,0])        #
    # Assign second altitude as vec  #
    #  as a change of the true anom. #
    r_f_i = r_f*np.cos(nu)           #
    r_f_j = r_f*np.sin(nu)           #
    r_2 = np.array([r_f_i ,r_f_j,0]) #
    r_1_mag = np.linalg.norm(r_1)    #
    r_2_mag = np.linalg.norm(r_2)    #
    ################################## 
#########################
# Position Vector Input ################################################
elif radius_vec_check in ('n','N','no','No','NO'):                     #
    # Input & Assignment of Initial Position Vector                    #
    r_1_i = float(input("| Enter the initial position Vector, (i): ")) #
    r_1_j = float(input("| Enter the initial position Vector, (j): ")) #
    r_1_k = float(input("| Enter the initial position Vector, (k): ")) #
    r_1  = np.array([r_1_i, r_1_j, r_1_k])                             #
    # Input & Assignment of Final Position Vector                      #
    r_2_i = float(input("| Enter the final position Vector, (i): "))   #
    r_2_j = float(input("| Enter the final position Vector, (j): "))   #
    r_2_k = float(input("| Enter the final position Vector, (k): "))   #
    r_2  = np.array([r_2_i, r_2_j, r_2_k])                             #
########################################################################
#################################    
# Error Message for Wrong Input #################################
else:                                                           #
    Error_Message_Helio = f"""
{'-'*42}
| Error: The Position Vector/Magnitude was not selected!
|       '{radius_vec_check}' is not an acceptable input
{'-'*42}
    """
    print(Fore.RED + Error_Message_Helio)    ####################
    input(Fore.RED +"| Press Enter To Exit...")                 #
    exit()                                                      #
#################################################################
#####################################
# Time and Prograde/Retrograd Orbit # 
###################################################################
# Input for Time Unit passed between Positions                    #
Time_Scale_Input = input("| Are you entering Time in Minutes, "   #
                         "Hours, Days, or Years?: ")              #
# Select Time Scale up to Years (60*60*24*365.256)                #
if Time_Scale_Input in ('y','Y','Years','YEARS','years'):         #
    Time_Scale = 60*60*24*365.256                                 #
elif Time_Scale_Input in ('d','D','Days','DAYS','days'):          #
    Time_Scale = 60*60*24                                         #
elif Time_Scale_Input in ('h','H','Hours','HOURS','hours'):       #
    Time_Scale = 60*60                                            #
elif Time_Scale_Input in ('m','M','Minutes','MINUTES','minutes'): #
    Time_Scale = 60                                               #
else:                                                             #
    Error_Message_Time_Scale = f"""
{'-'*42}
| Error: The Time Scale was not selected!
|       '{Time_Scale_Input}' is not an acceptable input
{'-'*42}
    """
    print(Fore.RED + Error_Message_Time_Scale)    ###############
    input(Fore.RED +"| Press Enter To Exit...")                 #
    exit()                                                      #
#################################################################
###########################################
# Input for time passed between Positions #########################
t_m      = float(input("| Enter the amount of time passed: "))    #
###################################################################
##################
# Change in Time #########
Delta_t = t_m*Time_Scale #
##########################
###################################
# Prograde/Retrograde Orbit Input ######################
pro_retro = input("| Is this a Prograde Orbit (y/n):") #
########################################################
# Maximum Iterations  #
######################################################################
nmax = int(input(f"{'-'*42}\n"                                       #
    "| Please enter a maximum number of iterations\n"                #
             " for the calculations (if Not Sure start with 100):")) #
######################################################################
'''
# Debug inputs - old ones...remeber the days?
r_1  = np.array([-1285.58, 8784.04, 2590.81]) 
r_2  = np.array([-9076.47, 5525.34, 3756.74])
t_m  = 20
pro_retro = 'y'
'''
#%% Calculations 
########################
# Initial Calculations #
################################
# Position Vector Magnitudes   #
r_1_mag = np.linalg.norm(r_1)  #
r_2_mag = np.linalg.norm(r_2)  #
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
# For Prograde Orbit                                           #
if pro_retro in('y','Y'):                                      #
#### Checking The Z component of Position Vec. Cross Product   #    
    if C_12[2] >= 0:                                           #
######## Assinging Change in True Anomaly Equation             #
        Delta_nu = np.arccos(D_12/(r_1_mag*r_2_mag))           #
    else:                                                      #
        Delta_nu = 2*np.pi - np.arccos(D_12/(r_1_mag*r_2_mag)) #
################################################################
# For Retrograde Orbit                                         #
elif pro_retro in('n','N'):                                    #
#### Checking The Z component of Position Vec. Cross Product   # 
    if C_12[2] >= 0:                                           #
######## Assinging Change in True Anomaly Equation             #
        Delta_nu = 2*np.pi - np.arccos(D_12/(r_1_mag*r_2_mag)) #
    else:                                                      #
        Delta_nu = np.arccos(D_12/(r_1_mag*r_2_mag))           #
else:                                                          #
    print(Fore.RED + Error_Message_pro_ret)                    #
    input(Fore.RED +'| Press Enter To Exit...')                #
    exit()                                                     #
################################################################
###################################################
# Constant (A) within Universal Kepler's Equation #####################
A =  np.sin(Delta_nu)*np.sqrt((r_1_mag*r_2_mag)/(1-np.cos(Delta_nu))) #
#######################################################################
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
    y_sub = r_1_mag + r_2_mag + A*((z*S(z) - 1)/np.sqrt(C(z)))
    return y_sub
# Kepler's Eq. as a function of z                        #
def F(z):
    F_kepler = ((y(z)/C(z))**1.5)*S(z) + A*np.sqrt(y(z))\
        - np.sqrt(mu)*Delta_t
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
#   Using Newtonian Polynomial Interpolation
#       -Finding the Roots (ie solutions) to
#        an Implicit Function
############################################
def Newt_Int(x_0):
#### Iteration for loop 
    for i in range(0, nmax):
######## Define Variables for Each Iteration Value (To Print To User)
        f       = F(x_0)
        f_prime = F_Prime(x_0)
######## Define the Ratio of the functions (F/F')
        Ratio = f/f_prime
######## Break For Loop When Desired Tolerance Achieved
        if np.abs(Ratio) <= Tolerance:
            break
####### Newton's Polynomial Interpolation x_(i) = x_(i-1) - (F/F') 
        x = x_0 - Ratio
####### Set New calculated value as the next to be used for calculations      
        x_0 = x 
####### Print Iteration Data to User
        Newt_Int_Message = f"""
{'-'*42}
| Itration {i} 
|   F = {f} 
|   F' = {f_prime}
|   |F/F'| = {np.abs(Ratio)}
{'-'*42}
| Using the Values:
|  S(z) = {S(x_0)}
|  C(z) = {C(x_0)}
|  y(z) = {y(x_0)}
|  z    = {x_0}
{'-'*42}
        """
        print(Fore.MAGENTA +Newt_Int_Message)
####### Exit If Maximum Iterations are Reached
        if i == nmax - 1:
####### Exit Error Message
            Error_Message_Iteration = f""""
{'-'*42}          
| Error: Newtonian Polynomial Interpolation Did not Converge!!
|          -Please ensure all enterd values are correct.
|           Or Increase Iteration Maximum.
{'-'*42}
            """
            print(Fore.RED + Error_Message_Iteration)                              
            input(Fore.RED + '| Press Enter To Exit...')                         
            exit() 
### Return Value
    return x
##########################
# Run Newt_Int to Find Z #
z = Newt_Int(z_0)        #
##########################
# Determine the Type of Orbit and set String for 
#   Output to User
if z < 0:
    Orbit_Determination = 'HYPERBOLIC'
elif z == 0:
    Orbit_Determination = 'PARABOLIC'
else:
    Orbit_Determination = 'ELLIPTIC'
###############################################################################
#
#########################################
# Velocity calculations Using Lambert's #
#########################################
# Calculating the Velocity
#   vectors using the
# Lagrange Coefficiants 
###############################
# f lagrange coeff            #
f     = 1 - (y(z)/r_1_mag)    #
# g lagrange coeff            #
g     = A*np.sqrt(y(z)/mu)    #
# g dot lagrange coeff        #
g_dot = 1 - (y(z)/r_2_mag)    #
# First Velocity Vector       #
v_1 = (1/g)*(r_2 - f*r_1)     #
# Second Velocity Vector      #
v_2 = (1/g)*(g_dot*r_2 - r_1) #
# First Velocity Magnitude    #
v_1_mag = np.linalg.norm(v_1) #
# Second Velocity Magnitude   #
v_2_mag = np.linalg.norm(v_2) #
###############################
#
###################################
## Orbital Elements Calculations ##
###################################
#   From the Calculated Velocities
#       from Lambert's
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
####################################
#
#####################
# Semi-Latus Rectum # 
#######################
# 1st semi-latus rec  #
p_1 = h_1_mag**2/mu   #
# 2nd semi-latus rec  #
p_2 = h_2_mag**2/mu   #
#######################
#
################
# Eccentricity #
############################################
# 1st eccentricity vector                  #
e_1 = np.cross(v_1,h_1)/mu - r_1/r_1_mag   #
# 1st eccentricity vector magnitude        #
e_1_mag = np.sqrt((e_1*e_1).sum())         #
# 2nd eccentricity vector                  #
e_2 = np.cross(v_2,h_2)/mu - r_2/r_2_mag   #
# 1st eccentricity vector magnitude        #
e_2_mag = np.sqrt((e_2*e_2).sum())         #
############################################
#
###################
# Semi-Major Axis #
############################
# 1st semi-major axis      #
a_1 = p_1/(1 - e_1_mag**2) #
# 2nd semi-major axis      #
a_2 = p_2/(1 - e_2_mag**2) #
############################
#
####################################
# Specific Energy of Each Position #
###########################################
# 1st Specific Energy                     #
Energy_1 = -(mu*(1-e_1_mag**2))/(2*p_1)   #
# 2nd Specific Energy                     #
Energy_2 = -(mu*(1-e_2_mag**2))/(2*p_2)   #
###########################################
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
#####################################
#
###############
# Data Report #
###############
Final_Data_Report = f"""
{'-'*42}
| The Velocities are: (km/s)
| 1st Position: {v_1[0]} i\u0302  {v_1[1]} j\u0302  {v_1[2]} k\u0302
|    Magnitude: {v_1_mag}
| 2nd Position: {v_2[0]} i\u0302  {v_2[1]} j\u0302  {v_2[2]} k\u0302
|    Magnitude: {v_2_mag}
|
| For the Given Positions: (km)
| 1st Position:  {r_1[0]} i\u0302  {r_1[1]} j\u0302  {r_1[2]} k\u0302
|    Magnitude:  {r_1_mag}
| 2nd Position:  {r_2[0]} i\u0302  {r_2[1]} j\u0302  {r_2[2]} k\u0302  
|    Magnitude:  {r_2_mag} 
|
{'-'*42}
| z = {z}
| This is a/an {Orbit_Determination} Orbit
|
{'-'*42}
| The orbital elements were calcualted with
|  both position vectors and the velocity 
|  vectors calculated by Lambert's
|  
|  -The initial and final position are
|   subscripted as 1 & 2 respectivly
{'-'*42}
| The Angular Momentums are:
|  h_1 = {h_1_mag} (km^2/s)
|  h_2 = {h_2_mag} (km^2/s)
{'-'*42}
| The Semi-Latus Rectums are:
|  p_1 = {p_1} (km)
|  p_2 = {p_2} (km)
{'-'*42}
| The Eccentricities are:
|  e_1 = {e_1_mag}
|  e_2 = {e_2_mag}
{'-'*42}
| The Semi-Major Axii are:
|  a_1 = {a_1} (km)
|  a_2 = {a_2} (km)
{'-'*42}
| The Specific Energy of Each Position is:
|  Energy_1 = {Energy_1} (km^2/s^2)
|  Energy_2 = {Energy_2} (km^2/s^2)
{'-'*42}
| The Inclinations are:
|  i_1 = {i_1_deg} (degrees)
|  i_2 = {i_2_deg} (degrees)
{'-'*42}
| The Longitude of Ascending Nodes are:
|  \u03A9_1 = {Omega_1_deg} (degrees)
|  \u03A9_2 = {Omega_2_deg} (degrees)
|{'-'*42}
| The Argumaent of Periapsis are:
|  \u03C9_1 = {omega_1_deg} (degrees)
|  \u03C9_2 = {omega_2_deg} (degrees)
|{'-'*42}
| The True Anomalies are:
|  \u03BD_1 = {nu_1_deg} (degrees)
|  \u03BD_2 = {nu_2_deg} (degrees)
|{'-'*42}
| The Eccentricity Anomalies are:
|  E_1 = {E_1_deg} (degrees)
|  E_2 = {E_2_deg} (degrees)
|{'-'*42}
| The Mean Anomalies are:
|  M_1 = {M_1_deg} (degrees)
|  M_2 = {M_2_deg} (degrees)
|{'-'*42}
"""
print(Fore.CYAN + Final_Data_Report)
###################################
# Keep .exe Open Until User Exits ############
input(Fore.RED + "| Press Enter to Exit...") #
exit()                                       #
##############################################