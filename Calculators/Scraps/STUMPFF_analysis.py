# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 00:13:03 2023

@author: Galactikhan
"""
'''
def S(z):
    for k in range (0,stumff_lim):
        S_Stumpff = ((-1)**k*z**k)/mt.factorial(2*k + 3)
    return S_Stumpff
    
def C(z):
    for k in range (0,stumff_lim):
        C_Stumpff = ((-1)**k*z**k)/mt.factorial(2*k + 2)
    return C_Stumpff
'''
import math as mt
z = 1

def S(z):
    S_Stumpff = 0
    for k in range (0,4):
        S_Stumpff += ((-1)**k)*((z)**k)/mt.factorial(2*k + 3)
    return S_Stumpff


S_i = 1/6 - z/120 + z**2/5040 - z**3/362880


out = f'''
| {S(z)}
|
|  {S_i}
'''
print(out)