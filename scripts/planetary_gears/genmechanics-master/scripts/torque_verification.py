# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 17:31:34 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import numpy as npy
import scipy.linalg as linalg

F=100
L=0.1

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')

bearing=genmechanics.linkages.FrictionlessRevoluteLinkage(ground,shaft1,[0,0,0],[0,0,0],'bearing')


load1=loads.KnownLoad(shaft1,[0,L,0],[0,0,0],[0,0,F],[0,0,0],'load')
load2=loads.SimpleUnknownLoad(shaft1,[0,0,0],[0,0,0],[],[0],'output_torque')

mech=genmechanics.Mechanism([bearing], ground,[(shaft1,0,10)],[load1],[load2])


# r=mech.StaticAnalysis()

for l,lv in mech.static_results.items():
    
    for d,v in lv.items():
        print(l.name,d,v)