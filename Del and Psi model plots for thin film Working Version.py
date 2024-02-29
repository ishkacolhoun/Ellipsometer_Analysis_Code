# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Coded based on the Ellipsometry and polarised light pg 283 - 288
#Ellipsometric function of a film-substrate system: Applications to the design of reflection-type optical devices and to ellipsometry
# assuming isotropic materials



import numpy as np
import matplotlib
import matplotlib.pyplot as plt


wavelength = 650e-9 # wavelength of light


phi_0_DEG = np.linspace(1,91,100) # incident angle in medium 1 (air), currently sweeping over angles from 0 to 90 deg
#phi_0_DEG = (10,20,30,40,50,60)
phi_0 = np.deg2rad(phi_0_DEG)

N_0 = 1.000 # refractive index of ambient medium (air)
N_1 = 0.15557 + 3.6024j  # refractive index of thin film
N_2 = 1.5   #refractive index of Substrate
d = 10e-9 # Film thickness

# Use Snells law to find phi_1 and phi_2

phi_1 = np.arcsin(N_0*np.sin(phi_0)/N_1)
phi_2 = np.arcsin(N_1*np.sin(phi_1)/N_2)

# Fresnel Equations for s and p polarised light at the different interfaces
r_01p = (N_1*np.cos(phi_0) - N_0*np.cos(phi_1)) / (N_1*np.cos(phi_0) + N_0*np.cos(phi_1))
r_12p = (N_2*np.cos(phi_1) - N_1*np.cos(phi_2)) / (N_2*np.cos(phi_1) + N_1*np.cos(phi_2))
r_01s = (N_0*np.cos(phi_0) - N_1*np.cos(phi_1)) / (N_0*np.cos(phi_0) + N_1*np.cos(phi_1))
r_12s = (N_1*np.cos(phi_1) - N_2*np.cos(phi_2)) / (N_1*np.cos(phi_1) + N_2*np.cos(phi_2))


Beta = 2 * np.pi * (d / wavelength) * N_1 * np.cos(phi_1)

# Reflection Coefficients for p and s waves
R_p =( r_01p + r_12p * np.exp(-2j*Beta) )/( 1 + r_01p * r_12p * np.exp(-2j*Beta))
R_s =( r_01s + r_12s * np.exp(-2j*Beta) )/( 1 + r_01s * r_12s * np.exp(-2j*Beta))

rho = R_p / R_s
#print(rho)
Psi = np.degrees(np.arctan(np.absolute(rho)))
#Delta = np.unwrap(np.angle(rho))
Delta = np.angle(rho,deg=True)
#Delta = np.rad2deg(Delta)  #+24
#Delta = Delta + 24
#Delta = np.unwrap(Delta)
plt.figure(dpi=1200)
#plt.figure(0)
plt.plot(phi_0_DEG,Psi,label='ψ')
plt.xlabel('Angle of Incidence')
plt.ylabel('Psi')

#plt.plot(phi_0,Delta)
#plt.show

#plt.figure(0)
plt.xlabel('Angle of Incidence (deg)')
plt.ylabel('Ellipsometry Angle (deg)')
plt.plot(phi_0_DEG,Delta,label='Δ')
plt.legend()
plt.title('Thin Gold Film')
plt.savefig('Linear Graph Final1.png')
plt.show()

