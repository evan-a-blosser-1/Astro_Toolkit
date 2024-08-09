# Dr. Diogo Shanchez

import math
import numpy as np
mu  = 3.986004418e5
R_e = 6378.136 

# r_1_vec = np.array([5000, 10000, 2100]) # Example 5.2 (Curtis)
# r_1_vec = np.array([3600, 4600, 3600]) # Problem 5.3 (Curtis)
# r_1_vec = np.array([5644, 2830, 4170]) # Problem 5.6 (Curtis)
# r_1_vec = np.array([1, 0, 0])*R_e # Problem 5.10 (Bate)
r_1_vec = np.array([-1285.58, 8784.04, 2590.81])
r_1 = np.linalg.norm(r_1_vec)
print(f"r_1_vec = {r_1_vec} (km), r_1 = {r_1} km")
# r_2_vec = np.array([-14600, 2500, 7000]) # Example 5.2 (Curtis)
# r_2_vec = np.array([-5500, 6240, -5200]) # Problem 5.3 (Curtis -- corrected)
# r_2_vec = np.array([-2240, 7320, 4980]) # Problem 5.6 (Curtis)
# r_2_vec = np.array([1, 1/8, 1/8])*R_e # Problem 5.10 (Bate)
r_2_vec = np.array([-9076.47, 5525.34, 3756.74])
r_2 = np.linalg.norm(r_2_vec)
print(f"r_2_vec = {r_2_vec} (km), r_1 = {r_2} km")
# delta_t = 3600 #sec # Example 5.2 (Curtis)
# delta_t = 30*60 #sec # Problem 5.3 (Curtis)
# delta_t = 20*60 #sec # Problem 5.6 (Curtis)
# delta_t = 1/8*806.8118744 # Problem 5.10 (Bate)
delta_t = 20*60 #sec
print(f"delta_t = {delta_t} sec")
# I will use numpy from now on to perform the vectorial calculations.
# However, you can explicitly calculate everything.
prograde = 1 # I am choosing a prograde trajectory; otherwise, prograde = 0
r1_dot_r2 = np.dot(r_1_vec,r_2_vec) # dot product between r_1 and r_2
print(f"r1_dot_r2 = {r1_dot_r2}")
r1_cross_r2 = np.cross(r_1_vec,r_2_vec) # cross product between r_1 and r_2
print(f"r_1_vec x r_2_vec = {r1_cross_r2}")
argument = r1_dot_r2/(r_1*r_2)
print(f"r1_dot_r2/(r_1*r_2) = {r1_dot_r2/(r_1*r_2)}")
if prograde == 1:
    if r1_cross_r2[2] >= 0:
        delta_theta = math.acos(argument)
    else:
        delta_theta = 2*math.pi - math.acos(argument)
else:
        if r1_cross_r2[2] < 0:
            delta_theta = math.acos(argument)
        else:
            delta_theta = 2*math.pi - math.acos(argument)
print(f"delta_theta = {delta_theta} rad or {np.rad2deg(delta_theta)}")

A = math.sin(delta_theta)*math.sqrt((r_1*r_2)/(1 - math.cos(delta_theta)))
print(f"A = {A}")

# In this problem, I will define all the funcions as python functions to optmize the calculations.
# The student can explicitly define everything.
z_0 = 0 # initial value for z
def S(z):
    S_value = 1/6 - z/120 + z**2/5040 - z**3/362880
    return S_value

def C(z):
    C_value = 1/2 - z/24 + z**2/720 - z**3/40320
    return C_value

def y(z):
    y_value = r_1 + r_2 + A*(z*S(z) - 1)/(math.sqrt(C(z)))
    return y_value 

def F(z):
    F_value = (y(z)/C(z))**1.5*S(z) + A*math.sqrt(y(z)) - math.sqrt(mu)*delta_t
    return F_value
def F_prime(z):
    if z == 0:
        term_1 = math.sqrt(2)/40*y(0)**1.5
        term_2 = math.sqrt(y(0)) + A*math.sqrt(1/(2*y(0)))
    else:
        term_1 = (y(z)/C(z))**1.5*(1/(2*z)*(C(z) - 1.5*S(z)/C(z)) + 0.75*S(z)**2/C(z))
        term_2 = 3*S(z)/C(z)*math.sqrt(y(z)) + A*math.sqrt(C(z)/y(z))
    
    F_prime_value = term_1 + A/8*term_2
    return F_prime_value

print(f"Initial values:\nS(z) = {S(z_0)}, C(z) = {C(z_0)}, y(z) = {y(z_0)}\n")
# Newton Method
def Newton(x_0):
    nmax = 1000
    
    for j in range(0, nmax):
        f = F(x_0)
        f_prime = F_prime(x_0)
        ratio = f/f_prime
        print(f"Newton: iteration {j}")
        print(f"Newton: F = {f}, F' = {f_prime}, |ratio| = {np.abs(ratio)}")
# stop when achieve the desired precision
        if np.abs(ratio) <= 1.0e-6:
            break
        x = x_0 - ratio
        x_0 = x #update the value of x_0
        print(f"Initial values:\nS(z) = {S(x_0)}, C(z) = {C(x_0)}, y(z) = {y(x_0)}, z = {x_0}\n")
        if j == nmax - 1:
            print(f"Newton: no convergence. Try another initial guess.")
    return x

z = Newton(z_0)
if z < 0:
    print(f"z = {z}; the orbit is hyperbolic.")
elif z == 0:
    print(f"z = {z}; the orbit is parabolic.")
else:
    print(f"z = {z}; the orbit is elliptic.")
    
    
f = 1 - y(z)/r_1
g = A*math.sqrt(y(z)/mu)
gdot = 1 - y(z)/r_2
# Calculating the velocities using the Lagrange coefficients:
v_1_vec = 1/g*(r_2_vec - f*r_1_vec)
v_1 = np.linalg.norm(v_1_vec) # magnitude of the vector v_1
v_2_vec = 1/g*(gdot*r_2_vec - r_1_vec)
v_2 = np.linalg.norm(v_2_vec) # magnitude of the vector v_2
print(f"v_1_vec = {v_1_vec} (km/s), v_1 = {v_1} km/s")
print(f"v_2_vec = {v_2_vec} (km/s), v_2 = {v_2} km/s")



# I will create a funtion to calculate the orbital elements given position and velocity:
def orb_elem(r_vec, v_vec):
    r = np.linalg.norm(r_vec) # magnitude of the position vector
    v = np.linalg.norm(v_vec) # magnitude of the velocity vector
    h_vec = np.cross(r_vec, v_vec) # vector angular momentum
    h = np.linalg.norm(h_vec) # magnitude of the angular momentum
# print(f"h = {h} km^2/s")
    p = h**2/mu # semi latus rectum
# print(f"p = {p} km")
    e_vec = 1/mu*((v**2 - mu/r)*r_vec - np.dot(r_vec,v_vec)*v_vec) # eccentricity vector
    e = np.linalg.norm(e_vec) # eccentricity
# print(f"e = {e}")
    a = h**2/mu/(1 - e**2) # semi-major axis
# print(f"a = {a} km")
    i = math.acos(h_vec[2]/h) # inclination
# print(f"i = {np.rad2deg(i)} deg")
# Defining the direction k and the Nodal vector
    k_dir = np.array([0, 0, 1])
    Nodal_vec = np.cross(k_dir, h_vec)
    Nodal = np.linalg.norm(Nodal_vec)
# Longitude of the ascending node
    if Nodal_vec[1] > 0:
        omega = math.acos(Nodal_vec[0]/Nodal)
    else:
        omega = 2*math.pi - math.acos(Nodal_vec[0]/Nodal)
# print(f"\u03A9 = {np.rad2deg(omega)} deg")
# Argument of the perigee
    if e_vec[2] > 0:
        argp = math.acos(np.dot(Nodal_vec,e_vec)/(Nodal*e))
    else:
        argp = 2*math.pi - math.acos(np.dot(Nodal_vec,e_vec)/(Nodal*e))
# print(f"\u03C9 = {np.rad2deg(argp)} deg")
# True anomaly
    if np.dot(r_vec,v_vec) > 0:
        theta = math.acos(np.dot(e_vec,r_vec)/(e*r))
    else:
        theta = 2*math.pi - math.acos(np.dot(e_vec,r_vec)/(e*r))
# print(f"\u03BD = {np.rad2deg(theta)} deg")
# I will return the all the evalues in a list:
    oe = [h, a, e, i, omega, argp, theta]
    return oe
# Calculating the values using the first pair of vectors (r_1_vec and v_1_vec):
print(f"First point in the orbit:")
oe = orb_elem(r_1_vec, v_1_vec)
print(f"h = {oe[0]} km^2/s")
print(f"a = {oe[1]} km")
print(f"e = {oe[2]}")
print(f"i = {np.rad2deg(oe[3])} deg")
print(f"\u03A9 = {np.rad2deg(oe[4])} deg")
print(f"\u03C9 = {np.rad2deg(oe[5])} deg")
print(f"\u03BD = {np.rad2deg(oe[6])} deg")
# Calculate the energy:
E = - mu/(2*oe[1])
print(f"Energy: E = {E} (km/s)^2")
# Perigee altitude
rp = oe[1]*(1 - oe[2]) # perigee
zp = rp - R_e
print(f"Perigee altitude = {zp} km.\n")
# Calculating the values using the second pair of vectors (r_2_vec and v_2_vec):
print(f"Second point in the orbit:")
oe = orb_elem(r_2_vec, v_2_vec)
print(f"h = {oe[0]} km^2/s")
print(f"a = {oe[1]} km")
print(f"e = {oe[2]}")
print(f"i = {np.rad2deg(oe[3])} deg")
print(f"\u03A9 = {np.rad2deg(oe[4])} deg")
print(f"\u03C9 = {np.rad2deg(oe[5])} deg")
print(f"\u03BD = {np.rad2deg(oe[6])} deg")
# Calculate the energy:
E = - mu/(2*oe[1])
print(f"Energy: E = {E} (km/s)^2")
# Perigee altitude
rp = oe[1]*(1 - oe[2]) # perigee
zp = rp - R_e
print(f"Perigee altitude = {zp} km.")