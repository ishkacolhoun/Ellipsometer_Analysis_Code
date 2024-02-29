# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Coded based on the following article
#Ellipsometric function of a film-substrate system: Applications to the design of reflection-type optical devices and to ellipsometry
# assuming isotropic materials



import numpy as np
import matplotlib
import matplotlib.pyplot as plt


wavelength = 546.1e-9 # wavelength of light


phi_0 = np.linspace(1,90,100)

phi_0_radians = np.radians(phi_0) # incident angle in medium 0
N_0 = 1  # refractive index of ambient medium (air)
N_1 = 3.8515 - 0.016460j  # refractive index of substrate


# Use Snells law to find phi_1 and phi_2

phi_1 = np.arcsin(N_0*np.sin(phi_0_radians)/N_1)



# Fresnel Equations for s and p polarised light at the different interfaces
r_01p = (N_1*np.cos(phi_0_radians) - N_0*np.cos(phi_1)) / (N_1*np.cos(phi_0_radians) + N_0*np.cos(phi_1))
r_01s = (N_0*np.cos(phi_0_radians) - N_1*np.cos(phi_1)) / (N_0*np.cos(phi_0_radians) + N_1*np.cos(phi_1))



#D_theta = wavelength/2 *(N_1**2 - np.sin(phi_0)**2)**0.5 # Eq. 13

#X = 0 + np.exp(-1j*2*np.pi*(d/D_theta))  # Eq. 12 


rho = r_01p / r_01s

#Psi = np.rad2deg(np.arctan(rho))
Psi = np.rad2deg(np.arctan(np.absolute(rho)))
Delta = np.angle(rho)
Delta = np.rad2deg(Delta)


plt.figure(dpi=1200)
#plt.figure(0)
plt.plot(phi_0,Psi,label='ψ')
plt.xlabel('Angle of Incidence')
plt.ylabel('Psi')

#plt.plot(phi_0,Delta)
#plt.show

#plt.figure(0)
plt.xlabel('Angle of Incidence (deg)')
plt.ylabel('Ellipsometry Angle (deg)')
plt.plot(phi_0,Delta,label='Δ')
plt.legend()
plt.title('Thick Gold Film')

plt.savefig('Linear Graph Final1.png')

plt.show