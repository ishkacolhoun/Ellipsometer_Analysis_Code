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
import csv
import scipy
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
from symfit import parameters, variables, log, Fit, Model

#x = []
#X_error = 0
Psi_exp = []
Delta_exp = []
#a = []
#b = []
#Y_error = []
header = []
Series_Label = ''
Graph_Title = 'Fixed AoI for Thick Au'
CSV_file_location = 'Exp data.csv'
Dot_Size = 3
N_2_given = 3.8515 - 0.016460j
thickness = 70e-9
print('d:',thickness)
print(N_2_given)

with open(CSV_file_location ,'r') as csvfile:
    header = next(csvfile).split(',')
    plots = csv.reader(csvfile, delimiter=',')
    X_title = ''#header[0]
    Y_title = ''#header[1]
    #Y_error_title = header[9]
    for row in plots:
        Psi_exp.append(float(row[1]))
     #   x.append(float(row[0]))
        Delta_exp.append(float(row[2]))

def Delta(phi_0_DEG, N_1_real, N_1_imag):
    wavelength = 650e-9
    d = thickness
    N_0 = 1
    N_2 = N_2_given
    N_1 = N_1_real + N_1_imag*1j
    #phi_0_DEG = np.linspace(1,91,100) # incident angle in medium 1 (air), currently sweeping over angles from 0 to 90 deg
    phi_0 = np.deg2rad(phi_0_DEG)
    phi_1 = np.arcsin(N_0*np.sin(phi_0)/N_1)
    phi_2 = np.arcsin(N_1*np.sin(phi_1)/N_2)
    
    r_01p = (N_1*np.cos(phi_0) - N_0*np.cos(phi_1)) / (N_1*np.cos(phi_0) + N_0*np.cos(phi_1))
    r_12p = (N_2*np.cos(phi_1) - N_1*np.cos(phi_2)) / (N_2*np.cos(phi_1) + N_1*np.cos(phi_2))
    r_01s = (N_0*np.cos(phi_0) - N_1*np.cos(phi_1)) / (N_0*np.cos(phi_0) + N_1*np.cos(phi_1))
    r_12s = (N_1*np.cos(phi_1) - N_2*np.cos(phi_2)) / (N_1*np.cos(phi_1) + N_2*np.cos(phi_2))
    
    Beta = 2 * np.pi * (d / wavelength) * N_1 * np.cos(phi_1)
    
    # Reflection Coefficients for p and s waves
    R_p =( r_01p + r_12p * np.exp(-2j*Beta) )/( 1 + r_01p * r_12p * np.exp(-2j*Beta))
    R_s =( r_01s + r_12s * np.exp(-2j*Beta) )/( 1 + r_01s * r_12s * np.exp(-2j*Beta))

    rho = R_p / R_s
    #Delta = np.unwrap(np.angle(rho))
    Delta = np.angle(rho,deg=True)
    return(Delta)


def Psi(phi_0_DEG, N_1_real, N_1_imag):
    wavelength = 650e-9
    d = thickness
    N_0 = 1
    N_2 = N_2_given
    #print(N_2)
    N_1 = N_1_real + N_1_imag*1j
    
    #phi_0_DEG = np.linspace(1,91,100) # incident angle in medium 1 (air), currently sweeping over angles from 0 to 90 deg
    phi_0 = np.deg2rad(phi_0_DEG)
    phi_1 = np.arcsin((N_0*np.sin(phi_0))/N_1)
    phi_2 = np.arcsin((N_1*np.sin(phi_1))/N_2)
    
    #print('Phi_1',phi_1)
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
    return(Psi)

phi_0_DEG = (20,25,30,35,40,45,50,55,60)
#print(Psi(phi_0_DEG,.15+3.6j,10e-9))

def residual(x):
    N_1_real = x[0]
    N_1_imag = x[1]
    Psi_pred = Psi(phi_0_DEG, N_1_real,N_1_imag)
    Delta_pred = Delta(phi_0_DEG, N_1_real, N_1_imag)
    #print('Psi_pred:',Psi_pred ,'/n','Psi_exp',Psi_exp)
   # print('Delta_pred:',Delta_pred ,'/n','Delta_exp',Delta_exp)
    err1 = np.mean((Psi_exp - Psi_pred)**2)
    err2 = np.mean((Delta_exp - Delta_pred)**2)
    error = (err1 + err2) / 2
    return error

x0 = [ 3.8515, 0.016460 ] # initial guess for N_1 and d

#print(x0[2])

res = minimize(residual, x0,bounds=[(0,None),(0,None)], method="Nelder-Mead")  
print(' n:',res.x[0],'\n k:',res.x[1])#,'\n d:',res.x[2])




AoI_range = np.linspace(1, 91,100)
Delta_plot = Delta(AoI_range, res.x[0], res.x[1])
Psi_plot = Psi(AoI_range, res.x[0], res.x[1])

plt.figure(dpi=1200)
#plt.figure(0)
plt.plot(AoI_range,Psi_plot,label='ψ Theory')
#plt.xlabel('Angle of Incidence')
#plt.ylabel('Psi')

# #plt.plot(phi_0,Delta)
# #plt.show

#plt.figure(0)
plt.xlabel('Angle of Incidence (deg)')
plt.ylabel('Ellipsometry Angle (deg)')
plt.plot(AoI_range,Delta_plot,label='Δ Theory')
plt.scatter(phi_0_DEG,Psi_exp,label='ψ Measured')
plt.scatter(phi_0_DEG,Delta_exp,label='Δ Measured')
plt.legend()
plt.title('Thick Gold Film on Glass')
plt.savefig('Linear Graph Final1.png')
plt.show()


# #fit = scipy.optimise.curve_fit(Delta,x,y)

