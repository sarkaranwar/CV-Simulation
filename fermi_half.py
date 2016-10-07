# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 06:22:46 2016

@author: sarkaranwar
"""
import numpy as np

def fermi_half(eta):
    x = eta
    mu = x**4 + 50 + 33.6*x*( 1 - 0.68*np.exp(-0.17*(x+1)**2))
    xi = 3*np.sqrt(np.pi)/(4*mu**(3/8))
    y = (np.exp(-x)+xi)**-1
    return y