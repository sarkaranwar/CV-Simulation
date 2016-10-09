# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 14:08:11 2016

@author: Sarkar Anwar
@e-mail: sarkar.rm.anwar@gmail.com

Semiconductor class, calculates and returns based on only name of the semiconductor
"""

import constants as const

class semiconductor(object):
    
    def readparams(self):
        sc_filename = self.name + '.txt'
        sc_file = open(sc_filename,"r")
        self.epsr = float(sc_file.readline().split("\t")[1])
        self.affinity = float(sc_file.readline().split("\t")[1])
        self.eg_gamma = float(sc_file.readline().split("\t")[1])
        self.eg_X = float(sc_file.readline().split("\t")[1])
        self.eg_L = float(sc_file.readline().split("\t")[1])
        self.alpha_T = float(sc_file.readline().split("\t")[1])
        self.beta_T = float(sc_file.readline().split("\t")[1])
        self.me_gamma = float(sc_file.readline().split("\t")[1])
        self.me_X = float(sc_file.readline().split("\t")[1])
        self.me_L = float(sc_file.readline().split("\t")[1])
        self.mh = float(sc_file.readline().split("\t")[1])
        self.vt_e = float(sc_file.readline().split("\t")[1])
        self.vt_h = float(sc_file.readline().split("\t")[1])
        sc_file.close()
        
    def calculate_bandgap(self):
        self.eg_gamma = self.eg_gamma - self.alpha_T * self.T**2/(self.T + self.beta_T)
        self.eg_X = self.eg_X - self.alpha_T * self.T**2/(self.T + self.beta_T)
        self.eg_L = self.eg_L - self.alpha_T * self.T**2/(self.T + self.beta_T)

#    def calculate_DOS(self):
#        self.Nc = 
    
    def __init__(self,name,T = 300):
        self.name = name
        self.T = T
        self.readparams()
        self.calculate_bandgap()
    