# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 06:22:46 2016

@author: Sarkar Anwar
@e-mail: sarkar.rm.anwar@gmail.com

Returns the electron or hole density for a given fermi level

"""
from numpy import exp, sqrt, trapz, arange, pi

def fermi_half_approx(eta):
#    https://nanohub.org/resources/5475
    x = eta
    mu = x**4 + 50 + 33.6*x*( 1 - 0.68*exp(-0.17*(x+1)**2))
    xi = 3*sqrt(pi)/(4*mu**(3/8))
    y = (exp(-x)+xi)**-1
    return y
    
def fermi_half_numeric(eta):
    xi = arange(0,300,0.01)
    y = sqrt(xi)/(1+exp(xi - eta))
    f = trapz(y,xi)
    return f
