#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:01:46 2018

@author: Pierrem
"""
import numpy as npy
from scipy import interpolate

#RadialRollerBearing_MaxAxialLoad
punctual_load={'data':[[1.2745097706120987,171.35678277400137],[12.058823356808148,71.85929385708647],
                       [29.80392214397284,42.21105446491762]],'x':'Linear','y':'Linear'}
mixed_load={'data':[[0.9803920599624654,114.07034832200415],[12.05882335680814,47.738690627211724],
                    [29.80392214397284,27.135674477838023]],'x':'Linear','y':'Linear'}
constant_load={'data':[[0.5882354212992595,58.79396770125257],[12.156863211430622,23.61808739733698],
                       [29.90196107198642,14.070348322004122]],'x':'Linear','y':'Linear'}

class RadialRollerBearingSNR:
    def __init__(self, d, D, B, Z, alpha, Dpw=None):
        self.d = d
        self.D = D
        self.B = B
        if Dpw == None:
            Dpw = (d + D)/2
        self.Dpw = Dpw
        self.alpha = alpha
        self.Z = Z
        
    def RuleAxialLoad(self, fr, fa, n, level_axial_load='constant_load'):
        k = 0.04 #variation of this coefficient with serial number
        val_x = ((self.Dpw*1e3)*(n/(2*npy.pi)*60))/(1e4) #unit in mm*tr/min
        if level_axial_load == 'constant_load':
            Pz = self.FunCoeff(val_x,constant_load)
        elif level_axial_load == 'mixed_load':
            Pz = self.FunCoeff(val_x,mixed_load)
        elif level_axial_load == 'punctual_load':
            Pz = self.FunCoeff(val_x,punctual_load)
        Pt = k*(self.d*1e3)**2*max(Pz, 0)
        check_rule_axial_load = True
        if 0.4*fr < Pt:
            check_rule_axial_load = False
        if fa > fr:
            check_rule_axial_load = False
        rule_axial_load = min(0.4*fr-Pt, fr-fa)
        return check_rule_axial_load, rule_axial_load
        
    def FunCoeff(self,x,data_struct,type_x='Linear',type_y='Linear'):
        data = npy.array(data_struct['data'])
        if type_x == 'Log': 
            x = npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol = float(f(x))
        if type_y == 'Log':
            sol = 10**sol
        return sol
    
#R1 = RadialRollerBearingSNR(1,1,1,1,1)
#print(R1.RuleAxialLoad(1,1,1))
    