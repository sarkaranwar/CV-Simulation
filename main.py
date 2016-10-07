# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 16:26:03 2016

@author: Sarkar Anwar
"""

import numpy as np
from scipy import interpolate
import constants as const
import matplotlib.pyplot as plt
from fermi_half import fermi_half

#import os


#os.system('cls')


T = 300			# Temperature in Kelvin
Ndop = 1e16		# Doping density in /cm3
tox = 2			# EOT in nm
sc_type = 'p'	# Semiconductor type

V_FB = 0.0


#print("Setup:")
#print("Temperature:", T, "K")
#print("Doping Density:", Ndop, "/cm3")
#print("EOT", tox, "nm")

phi = np.linspace(-2,2,401)
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


Vt = const.kB * T / const.q

#print("Thermal voltage- ",Vt)
n = 2*Nc/np.sqrt(np.pi) * fermi_half((phi - Eg/2)/Vt)
p = 2*Nv/np.sqrt(np.pi) * fermi_half((-phi - Eg/2)/Vt)

if sc_type == 'n':
    Nd = Ndop
    Na = 0
else:
    Nd = 0
    Na = Ndop

total_carrier = p-n+Nd-Na
#plt.semilogy(phi,np.abs(total_carrier))

#phi_b = np.interp(1e24,total_carrier,phi)
f_carrier_phi = interpolate.interp1d(total_carrier,phi,assume_sorted = False)
f_phi_carrier = interpolate.interp1d(phi,total_carrier,assume_sorted = False)
phi_b = f_carrier_phi(0)
print('Bulk potential is ',phi_b, 'eV')

#plt.semilogy(phi_s,n)
#plt.semilogy(phi_s,p)
#plt.xlabel('Potential')
#plt.ylabel('Carrier density (/m3)')
#plt.show()


# Define Grid

a = 0.21e-9         # Atomic radius in m
xmin = 0
xmax = 1e-6
x = np.arange(xmin,xmax,a)
#print(np.size(x))

phi_x = np.zeros(np.size(x))
#print(np.size(phi_x))

phi_s = np.linspace(-2,2,10)

for phi in phi_s:
    phi_x0 = np.linspace(phi,phi_b,np.size(x))
    carrier0 = f_phi_carrier(phi_x0)
#    print(np.size(carrier0))
#    plt.plot(x,carrier0)
    








