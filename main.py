# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 16:26:03 2016

@author: Sarkar Anwar
"""

import numpy as np
from scipy import interpolate
import constants as const
import matplotlib.pyplot as plt
from fermi_half import fermi_half_numeric, fermi_half_approx
from time import sleep
#import os


#os.system('cls')


T = 300			# Temperature in Kelvin
Ndop = 1e17		# Doping density in /cm3
tox = 2			# EOT in nm
sc_type = 'p'	# Semiconductor type

V_FB = 0.0

fermi_option = "numeric"
verbose = 0

#print("Setup:")
#print("Temperature:", T, "K")
#print("Doping Density:", Ndop, "/cm3")
#print("EOT", tox, "nm")

phi = np.linspace(-2,2,201)
#psi_s = np.linspace(-0.2,1.1,101)
Vg_exp = np.linspace(-2,2,101)

Eg = 1.12         	# Bandgap in eV
Nc = 3.2e19       	# Conduction band density of states /cm3
Nv = 2e19         	# Valence band density of states /cm3
ni = 1e10		    # Intrinsic carrier concentration in /cm3
epsr_s = 11.7		# Silicon dielectric constant
epsr_ox = 3.9

Ndop = 1e6 * Ndop	# Doping density in /m3
Nc = 1e6 * Nc       # Conduction band density of states /m3
Nv = 1e6 * Nv       # Valence band density of states /m3
tox = 1e-9 * tox	# EOT in m
ni = 1e6 * ni       # Intrinsic carrier concentration in /m3

eps_s = epsr_s * const.eps0
#print(eps_s)

Vt = const.kB * T / const.q

if fermi_option == "numeric":
    fermi = np.vectorize(fermi_half_numeric)
else :
    fermi = np.vectorize(fermi_half_approx)

#print("Thermal voltage- ",Vt)
n = 2*Nc/np.sqrt(np.pi) * fermi((phi - Eg/2)/Vt)
p = 2*Nv/np.sqrt(np.pi) * fermi((-phi - Eg/2)/Vt)

if sc_type == 'n':
    Nd = Ndop
    Na = 0
else:
    Nd = 0
    Na = Ndop

total_carrier = p-n+Nd-Na
#plt.semilogy(phi,np.abs(total_carrier))

carrier_2_phi = interpolate.interp1d(total_carrier,phi,assume_sorted = False)
phi_2_carrier = interpolate.interp1d(phi,total_carrier,assume_sorted = False)
phi_b = carrier_2_phi(0)
print('Bulk potential is ',phi_b, 'eV')

#plt.semilogy(phi,n1)
#plt.semilogy(phi,n2)
#plt.semilogy(phi,abs(total_carrier))
#plt.xlabel('Potential')
#plt.ylabel('Carrier density (/m3)')
#plt.show()


# Define Grid

a = 0.21e-9         # Atomic radius in m
xmin = 0
xmax = 0.5e-7
nx = 1000
dx = (xmax-xmin)/nx

x = np.arange(xmin,xmax,dx)
#print(x)

print("Number of grid points is ",np.size(x))
#"""
phi_x = np.zeros(np.size(x))
n_x = np.zeros(np.size(x))
p_x = np.zeros(np.size(x))
carrier_x = np.zeros(np.size(x))
Nd_x = Nd * np.ones(np.size(x))
Na_x = Na * np.ones(np.size(x))

phi_s = np.linspace(-1,1,2)

max_iter = 500


#plt.ion()
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax2 = fig.add_subplot(212)
#ax1.set_ylim([-2,2])
#ax2.set_ylim([-1e23,1e23])
#line1, = ax1.plot(x, phi_x, 'r-')
#line2, = ax2.plot(x,carrier_x, 'b-')


for phi in phi_s:
    print("Solving for phi_s = ",phi)
    phi_x0 = np.linspace(phi,phi_b,nx)
    carrier0 = phi_2_carrier(phi_x0)
    phi_x0[1:nx-1] = phi_b    
    phi_x0[0] = phi
    i = 0

    while (i<max_iter):
        #Solve poisson equation
        phi_x[0] = phi
        for idx, xi in enumerate(x):        
#            phi_x[1:nx-2] = (phi_x0[0:nx-3] + phi_x0[2:nx-1] + dx * dx * const.q*carrier0[1:nx-2]/eps_s)/2 #does not work
            if 0<idx<nx-1:
                phi_x[idx] = (phi_x0[idx-1] + phi_x0[idx+1] + dx**2 * const.q * carrier0[idx]/eps_s)/2
        phi_x[nx-1] = phi_b
        carrier_x = phi_2_carrier(phi_x)
        
        #Calculatue error        
        error_carrier = np.max((np.abs(carrier0 - carrier_x))/carrier0)
        error_phi = np.max(np.abs(phi_x0 - phi_x))
        phi_x0 = phi_x
        carrier0 = carrier_x
        if ((error_phi<1e-19) and (error_carrier<1e-9)):
            break
        i = i+1
    plt.subplot(211)        
    plt.plot(x,phi_x)
    plt.subplot(212)
    plt.plot(x,carrier_x)
plt.show()
