#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:46:45 2019

@author: ringhausen
"""


import mechanical_components.shafts as shafts
import matplotlib.pyplot as plt
import math
import random


vise = shafts.Accessory([(0.05, 0.02), (0.02, 0.025), (0, 0.03)], (0, 1))
roulement1 = shafts.Accessory([(0, 0.035), (0.03, 0.035), (0.03, 0.04)], (0, -1), [(0, 0.065), (0.03, 0.065)], (0, 1))
pignon1 = shafts.Accessory([(0.03, 0.04), (0.08, 0.04),], (0, -1))
circlip1 = shafts.Accessory([(0.08, 0.04), (0.08, 0.035), (0.085, 0.035), (0.085, 0.04)], (0, -1))
synchro1 = shafts.Accessory([(0.085, 0.04), (0.115, 0.04)], (0, -1))
pignon2 = shafts.Accessory([(0.115, 0.04), (0.165, 0.04)], (0, -1))
circlip2 = shafts.Accessory([(0.165, 0.04), (0.165, 0.035), (0.17, 0.035), (0.17, 0.04)], (0, -1))
pignon3 = shafts.Accessory([(0.17, 0.04), (0.22, 0.04)], (0, -1))
circlip3 = shafts.Accessory([(0.22, 0.04), (0.22, 0.035), (0.225, 0.035), (0.225, 0.04)], (0, -1))
synchro2 = shafts.Accessory([(0.225, 0.04), (0.25, 0.04)], (0, -1))
pignon4 = shafts.Accessory([(0.25, 0.04), (0.30, 0.04), (0.30, 0.05)], (0, -1))

pignon5 = shafts.Accessory([(0.30, 0.05) ,(0.30, 0.06), (0.33, 0.06), (0.33, 0.04)], (0, -1), montage='inbuilt')
pignon6 = shafts.Accessory([(0.41, 0.04), (0.41, 0.05), (0.45, 0.05), (0.45, 0.045)], (0, -1), montage='inbuilt')

roulement2 = shafts.Accessory([(0.45, 0.045), (0.45, 0.035), (0.48, 0.035)], (0, -1))
joint1 = shafts.Accessory([(0.485, 0.03), (0.505, 0.03)], (0, -1))

accessoires = [vise, roulement1, pignon1, circlip1, synchro1, pignon2,
               circlip2, pignon3, circlip3, synchro2, pignon4, pignon5,
               pignon6, roulement2, joint1]

#accessoires.reverse()
#accessoires.sort(key=lambda a: a.functional_points[0][1], reverse=True)
#accessoires = sorted(accessoires, key=lambda a: a.functional_points[0][0])

functional_accessories = shafts.FunctionalAccessories(accessoires)
shaft = shafts.Shaft(functional_accessories)

print(shaft.ShaftCheck())
shaft.Plot2()

















































