########################################################
# Author: Evan A. Blosser                              #
# Email: evan.a.blosser-1@ou.edu/galactikhan@gmail.com #
#                                                      #
# Program: Julian_Sidereal
#           -For loading ephemris data and calculating
#               the Julian/Sidereal Time 
########################################################
#
import numpy as np
###################
# User Time Input #
###################
# Have the input be split with  dashes '-'
#
#   Calculation Errors:
#       Keeps second digit is always off by 1
#       ie. 2464819.0 days, was calculated as 
#           2464809.0 Days???
###############################################################################
User_Input_Prompt = f"""
{'-'*42}
| Please input desired time in UT 
|   starting with
{'-'*42}
"""
print(User_Input_Prompt)
'''
##############################################################################
# Time Inputs                                                                #
year  = float(input('| Year:  '))                                            #
month = float(input('| Month: '))                                            #
day   = float(input('| Day:   '))                                            #
# Universla Time Input as dash seperated string                              # 
UT_hi, UT_mi, UT_si = input('| UT (dash seperated hr-min-sec): ').split('-') #
'''
year  = 2036                                                                  #
month = 5                                                                     #
day   = 5                               
'''                                                                           # 
# Convert Universla Time Input to float                                       #
UT_H = float(UT_hi)                                                           #
UT_M = float(UT_mi)                                                           #
UT_S = float(UT_si)                  
'''                                        
UT_H = 12                                                                   #
UT_M = 0                                                                    #
UT_S = 0
# Universla Time Calculations                                                #
UT = UT_H + (UT_M/60) + (UT_S/3600)                                          #
##############################################################################
########################
# Julian/Sidereal Time #
##############################################################################
# Simplest Formula fo obtaining J0, Boulet(1991)  ############################
Julian_0 = 367*year - np.fix((7*(year + np.fix((month + 9)/12)))/4) +\
    np.fix((257*month)/9) + day + 1721013.5
# Julian Date given Year, Month, and Day #####################################
Julian_D = Julian_0 + UT/24                                                  #
# T0 in Julian Cenbturies                                                    #
G_T0 = (Julian_0 - 2451545)/36525                                            #
# Greenwich Sidereal Time at 0h UT                                           #
Green_T0 = 100.4606184 + 36000.77004*G_T0 + 0.000387933*G_T0**2 -\
    2.583*10**(-8)*G_T0**3
# Greenwich Sidreal Time #####################################################
Green_Sidreal = Green_T0 + 360.98564724*(UT/24)                              #
##############################################################################