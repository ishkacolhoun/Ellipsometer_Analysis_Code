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
phi_0 = np.deg2rad(phi_0_DEG)

N_0 = 1.000 # refractive index of ambient medium (air)
N_1 = 0.15557 + 3.6024j  # refractive index of thin film
N_2 = 1.5   #refractive index of Substrate
d = 10e-9 # Film thickness
analyser_angle = np.radians(68.57)
pol_angle = np.radians(45)

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
R_p = r_01p + r_12p * np.exp(-2j*Beta) / 1 + r_01p * r_12p * np.exp(-2j*Beta)
R_s = r_01s + r_12s * np.exp(-2j*Beta) / 1 + r_01s * r_12s * np.exp(-2j*Beta)

rho = R_p / R_s

Psi = np.arctan(np.absolute(rho))
Delta = np.angle(rho)

Psi_prime = np.arctan(np.tan(Psi)/np.tan(pol_angle))

S_0 = 1 

S_1 = - np.cos(2*Psi_prime) 

S_2 = np.sin(2*Psi_prime)*np.cos(Delta)

I = .5*(S_0 + S_1*np.cos(2*analyser_angle)+S_2*np.sin(2*analyser_angle)) # Transmitted Intensity / Incident Intensity




plt.figure(0)
plt.plot(np.degrees(phi_0),I)
plt.xlabel('Incident Angle (deg)')
plt.ylabel('Intensity(%)')

plt.show()