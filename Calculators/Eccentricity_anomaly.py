##############################################
# Author: Evan A. Blosser                    #
# Email: evan.a.blosser-1@ou.edu             #
#                                            #
# Program: Using eccentricity & mean         #
#            anomaly to calculate the        #
#            eccentricity anomaly            # 
##############################################
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
import math                         #
#####################################
welcome = """
{'-'*42}
| Welcome to the Eccentricity Anomaly Calculator
{'-'*42}
"""
print(welcome)
###############
# User Inputs #
#########################################################################
e = float(input("|Please enter the eccentricity: " ))                   #
M_deg = float(input("|Please enter the mean anomaly (in degrtess): " )) #
nmax = int(input("|Please enter the desired number of iterations: " ))  #
### Convert input degrees to rad:                                       #
M_e = np.deg2rad(M_deg)                                                 #
#########################################################################
''' Debugg Variables
e = 0.8 # initial eccentricity
M_deg =90    #60 degrees
nmax = 1000 # maximum number of iterations
'''
########
# Main #
###########################################################################
# Initial guess for the eccentric anomaly
if M_e < math.pi:
    E_0 = M_e + e/2
elif M_e > math.pi:
    E_0 = M_e - e/2
else:
    E_0 = M_e
# For loop solving for eccentricity anomly using Newtons method
for j in range(0, nmax):
# Newton Root Finding 
    ratio = (E_0 - e*math.sin(E_0) - M_e)/(1 - e*math.cos(E_0))
    print(f"iteration no. {j}, |ratio| = {np.abs(ratio)}")
    if np.abs(ratio) <= 1.0e-6:
        break
    E = E_0 - ratio                                                      
    E_0 = E #update the value of E_0                                     
    E = np.rad2deg(E_0)                                                  
###########################################################################
output = f"""
{'-'*42}
Eccentric anomaly E = {E}
after {j} iterations
{'-'*42}
"""
print(output)