#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:38:58 2019

@author: ringhausen
"""

import mechanical_components.shafts as shafts
import mechanical_components.shafts_assembly as shaftsass
import matplotlib.pyplot as plt
import math
import random

#
#borne_sup_accessoire = 0.1
#longueur_accessoire = 0.02
#distance_accessoire = 0.1
#nb_accessoires = 10
#
#y = []
#for i in range(nb_accessoires):
#    y.append(random.random()*borne_sup_accessoire)
#
#external_vector = []
#for i in range(nb_accessoires):
#    if random.random() < 0.5:
#        external_vector.append((0,1))
#    else:
#        external_vector.append((0,-1))
#
#accessory = []
#for i in range(nb_accessoires):
#    accessory.append(shafts.Accessory([((i-1)*longueur_accessoire+(i-1)*distance_accessoire, y[i]),(i*longueur_accessoire+(i-1)*distance_accessoire, y[i])], external_vector[i]))
#    
#functional_accessories = shafts.FunctionalAccessories(accessory)
#shaft = shafts.Shaft(functional_accessories)
#
#
#while shaft.ShaftCheck() != True:
#    
#    y = []
#    for i in range(nb_accessoires):
#        y.append(random.random()*borne_sup_accessoire)
#    
#    external_vector = []
#    for i in range(nb_accessoires):
#        if random.random() < 0.5:
#            external_vector.append((0,1))
#        else:
#            external_vector.append((0,-1))
#    
#    accessory = []
#    for i in range(nb_accessoires):
#        accessory.append(shafts.Accessory([((i-1)*longueur_accessoire+(i-1)*distance_accessoire, y[i]),(i*longueur_accessoire+(i-1)*distance_accessoire, y[i])], external_vector[i]))
#        
#    functional_accessories = shafts.FunctionalAccessories(accessory)
#    shaft = shafts.Shaft(functional_accessories)
#
##    print(shaft.ShaftCheck())
##    print('--------------------------------')
##print(shaft.ShaftDrawingPoints())
#shaft.Plot2()


###############################################################################




