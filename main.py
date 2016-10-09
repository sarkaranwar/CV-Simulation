# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 16:26:03 2016

@author: Sarkar Anwar
@e-mail: sarkar.rm.anwar@gmail.com
"""

import numpy as np
from scipy import interpolate
import constants as const
import matplotlib.pyplot as plt
from fermi_half import fermi_half_numeric, fermi_half_approx
from semiconductor import semiconductor as sc




x = sc("Silicon", 300)


""" Define simulation parameters """
max_iter = 500
error = 1e-9
verbose = True


"""Define device parameters """
T = 300			# Temperature in Kelvin
Ndop = 1e17		# Doping density in /cm3
tox = 2			# EOT in nm
sc_type = 'p'	# Semiconductor type
V_FB = 0.0
fermi_option = "numeric"

Vg_exp = np.linspace(-2,2,101)

"""Define semiconductor parameters """
Eg = 1.12               # Bandgap in eV
Nc = 3.2e19             # Conduction band density of states /cm3
Nv = 2e19               # Valence band density of states /cm3
ni = 1e10               # Intrinsic carrier concentration in /cm3
epsr_s = 11.7           # Silicon dielectric constant
epsr_ox = 3.9

"""Convert semiconductor parameters to SI unit"""
Ndop = 1e6 * Ndop	       # Doping density in /m3
Nc = 1e6 * Nc           # Conduction band density of states /m3
Nv = 1e6 * Nv           # Valence band density of states /m3
tox = 1e-9 * tox        # EOT in m
ni = 1e6 * ni           # Intrinsic carrier concentration in /m3

""" Calculate some regularly used constants """
eps_s = epsr_s * const.eps0
Vt = const.kB * T / const.q

if fermi_option == "numeric":
    fermi = np.vectorize(fermi_half_numeric)
else :
    fermi = np.vectorize(fermi_half_approx)

""" Calculate Fermi Dirac distribution of carriers for a range of potential for easy interpolation """
phi_temp = np.linspace(-2,2,201)
n = 2*Nc/np.sqrt(np.pi) * fermi((phi_temp - Eg/2)/Vt)
p = 2*Nv/np.sqrt(np.pi) * fermi((-phi_temp - Eg/2)/Vt)

if sc_type == 'n':
    Nd = Ndop
    Na = 0
else:
    Nd = 0
    Na = Ndop

total_carrier = p-n+Nd-Na

""" Define interpolation functions for easy interpolation """
carrier_2_phi = interpolate.interp1d(total_carrier,phi_temp,assume_sorted = False)
phi_2_carrier = interpolate.interp1d(phi_temp,total_carrier,assume_sorted = False)

""" Calculate bulk potential """
phi_b = carrier_2_phi(0)


""" Define Grid """
a = 0.21e-9         # Atomic radius in m
xmin = 0
xmax = 0.5e-7
nx = 1000
dx = (xmax-xmin)/nx
x = np.arange(xmin,xmax,dx)
print("Number of grid points is ",np.size(x))

""" Define Surface potential """
phi_start = -1
phi_end = 1
phi_n = 2
phi_s = np.linspace(phi_start,phi_end,phi_n)

phi_x = np.zeros(np.size(x))
carrier_x = np.zeros(np.size(x))

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
            print("Error")
            break
        i = i+1
    if verbose:
        plt.subplot(211)        
        plt.plot(x,phi_x)
        plt.subplot(212)
        plt.plot(x,carrier_x)
plt.show()

if verbose:
    print('Setup:')
    print('Temperature: ',T, 'K')
    print('Doping density: ',Ndop/1e6, '/cm3')
    print('EOT: ',tox/1e-9, 'nm')
    print('Type: ',sc_type)
    print('Flat band voltage: ',V_FB, 'V')
    print('Semiconductor parameter: ')
    print('Bandgap: ',Eg, 'eV')
    print('CB DOS: ',Nc/1e6, '/cm3')    
    print('VB DOS: ',Nv/1e6, '/cm3')    
    print('Bulk potential: ',phi_b, 'eV')