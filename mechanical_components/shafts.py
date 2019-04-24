#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:36:08 2019

@author: ringhausen
"""

import volmdlr as vm
import volmdlr.primitives2D
import volmdlr.primitives3D as vm3d
import math
import numpy as npy
import matplotlib.pyplot as plt
from scipy import interpolate
import json
import pkg_resources
from itertools import chain


class ShaftMaterial:
    def __init__(self, Re, t_yield, C, name='',):
        self.Re = Re            # Ultimate tensile strength (Pa)
        self.t_yield = t_yield  # Yield stength (Pa)
        self.C = C              # Hardness (Brinell Scale) (Bhn)
        self.name = name
        
        

class Accessory:
    def __init__(self, functional_points, external_vector, shaft, linkage_points=None, linkage_external_vector=None, link=None, montage='axial'):
        """
        We suppose that the coordinate origin of the accessories 
            match the origin of the shaft
        The functional_points should be sorted according to increasing x
        The external_vector should be either (0,1) or (0,-1)
        The link can either be with an other accessory (ie. gears), 
            or with a shaft (ie. bearings)
        The montage should be either radial, axial or inbuilt
        """
        self.functional_points = functional_points
        self.external_vector = external_vector
        self.shaft = shaft
        self.linkage_points = linkage_points
        self.linkage_external_vector = linkage_external_vector
        self.link = link        
        self.montage = montage
        
#    def AccessoryBoundary(self):
#        min_point = min(self.functional_points, key=lambda point: point[1])
#        max_point = max(self.functional_points, key=lambda point: point[1])
#        return min_point[1], max_point[1]
        
    
        
class FunctionalAccessories:
    def __init__(self, accessories):
        self.accessories = accessories # list of Accessory
        
    def SortAccessories(self):
        int_acc = [acc for acc in self.accessories if acc.external_vector[1] > 0]
        ext_acc = [acc for acc in self.accessories if acc.external_vector[1] < 0]
        sorted_int_acc = sorted(int_acc, key=lambda a: a.functional_points[0][0])
        sorted_ext_acc = sorted(ext_acc, key=lambda a: a.functional_points[0][0], reverse=True)
        sorted_ext_acc.reverse()
        sorted_accessories = sorted_int_acc + sorted_ext_acc
        return sorted_accessories
        
    
    def FunctionalContour(self): 
        functional_contours = []
        functional_contours_inbuilt = []
        linkage_contours = []
        for accessory in self.accessories:
            if accessory.montage == 'axial' or accessory.montage == 'radial':
                functional_contours.append(vm.Contour2D([vm.primitives2D.RoundedLineSegments2D(accessory.functional_points, {})]))   
                print(accessory.linkage_points)
                if accessory.linkage_points != None:
                    linkage_contours.append(vm.Contour2D([vm.primitives2D.RoundedLineSegments2D(accessory.linkage_points, {})]))
            if accessory.montage == 'inbuilt':
                functional_contours_inbuilt.append(vm.Contour2D([vm.primitives2D.RoundedLineSegments2D(accessory.functional_points, {})]))
        return functional_contours, functional_contours_inbuilt, linkage_contours
    
        

class Shaft:
    def __init__(self, functional_accessories, origin=(0,0), Dmax=None, shaft_material=None):
        self.functional_accessories = functional_accessories 
        self.origin = origin
        self.Dmax = Dmax
        self.shaft_material = shaft_material
        
        self.sorted_functional_accessories = self.functional_accessories.SortAccessories()
        
        self.ext_index = []
        self.int_index = []
        for i in range(len(self.sorted_functional_accessories)):
            if self.sorted_functional_accessories[i].external_vector[1] > 0:
                self.int_index.append(i)
            else:
                self.ext_index.append(i)
        
        self.offset = 0.0
            
        

    def ShaftDrawingPoints(self):
        
        points_inf = []
        for accessoire in self.sorted_functional_accessories:
            if accessoire.external_vector[1] > 0:
                points_inf.extend([tuple(vm.Point2D(p) for p in accessoire.functional_points)])
            else:
                if self.offset != 0:
                    points_offset = []
                    for point in accessoire.functional_points:
                        pt = vm.Point2D(point)
                        points_offset.append(pt.Translation(self.offset * vm.Vector2D(accessoire.external_vector)))
                    points_inf.append(tuple(points_offset))
                    
        points_sup = []
        for accessoire in self.sorted_functional_accessories:
            if accessoire.external_vector[1] < 0:
                points_sup.extend([tuple(vm.Point2D(p) for p in accessoire.functional_points)])
            else:
                if self.offset != 0:
                    points_offset = []
                    for point in accessoire.functional_points:
                        pt = vm.Point2D(point)
                        points_offset.append(pt.Translation(self.offset * vm.Vector2D(accessoire.external_vector)))
                    points_sup.append(tuple(points_offset))

        return points_inf, points_sup 
    
    
    
    def ShaftDrawingContour(self):
        points_inf, points_sup = self.ShaftDrawingPoints()
        print('points inf: ', points_inf, '\n')
        print('points sup: ', points_sup)
        
        points_inf = list(chain.from_iterable(points_inf))
        points_sup = list(chain.from_iterable(points_sup))
        
#        points_sup.reverse()
        points = points_inf + points_sup + [points_inf[0]]
        
        line = vm.primitives2D.RoundedLineSegments2D(points, {})
        contour = vm.Contour2D([line])
        return contour

            
    
    def Plot2(self):
        contour = self.ShaftDrawingContour()
        f, a = contour.MPLPlot()
        contour_accessories, contour_inbuilts, linkage_contours = self.functional_accessories.FunctionalContour()
        for contour_accessory in contour_accessories:
            contour_accessory.MPLPlot(ax=a, style='r')
        for contour_inbuilt in contour_inbuilts:
            contour_inbuilt.MPLPlot(ax=a, style='b')
        for contour_linkage in linkage_contours:
            contour_linkage.MPLPlot(ax=a, style='g')
            
        
    def ShaftCheck(self):
        points_inf, points_sup = self.ShaftDrawingPoints()
        
        for i, accessoire in enumerate(self.sorted_functional_accessories[1:-1], 1):
            if accessoire.montage == 'axial':

                if accessoire.external_vector[1] > 0:
                    points_inf_left = points_inf[: i]
                    points_inf_right = points_inf[i+1 :]

                    if points_inf_left != [] and points_inf_right != []:
                        min_inf_left = min([point[1] for point in list(chain.from_iterable(points_inf_left))])
                        min_inf_right = min([point[1] for point in list(chain.from_iterable(points_inf_right))])
                        r_accessoire = min([point[1] for point in accessoire.functional_points])
                        
                        if r_accessoire > min_inf_left and r_accessoire > min_inf_right:
                            return False
                    
                else:
                    points_sup_left = points_sup[: i]
                    points_sup_right = points_sup[i+1 :]

                    if points_sup_left != [] and points_sup_right != []:
                        max_sup_left = max([point[1] for point in list(chain.from_iterable(points_sup_left))])
                        max_sup_right = max([point[1] for point in list(chain.from_iterable(points_sup_right))])
                        R_accessoire = max([point[1] for point in accessoire.functional_points])
                        
                        if R_accessoire < max_sup_left and R_accessoire < max_sup_right:
                            return False
        return True
            


class ShaftsAssembly:
    def __init__(self, shafts):
        self.shafts = shafts # A list of Shaft objects
    

        

    def TwoShaftsAssembly(self, shaft1, shaft2):
        return
#        if 
        
#        matching_accessories = []
#        for i, accessory1 in enumerate(shaft1.accessories):
#            for j, accessory2 in enumerate(shaft2.accessories):
#                if self.AccessoryComparisonFunctionalToLinkage(accessory1, accessory2):
#                    matching_accessories.append((i, j))
#                elif self.AccessoryComparisonLinkagelToLinkage(accessory1, accessory2):
#                    matching_accessories.append((i, j))
#                
                    
        
    def AccessoryComparisonFunctionalToLinkage(self, accessory1, accessory2):
        f1  = accessory1.functional_points
        v1  = accessory1.external_vector
        l1  = accessory1.linkage_points
        lv1 = accessory1.linkage_external_vector
        f2  = accessory2.functional_points
        v2  = accessory2.external_vector
        l2  = accessory2.linkage_points
        lv2 = accessory2.linkage_external_vector
        
        if f1 == l2 and l1 == f2 and v1 == lv2 and lv1 == v2:
            return True
        else:
            return False
            
            
    def AccessoryComparisonLinkagelToLinkage(self, accessory1, accessory2):
        l1  = accessory1.linkage_points
        lv1 = accessory1.linkage_external_vector
        l2  = accessory2.linkage_points
        lv2 = accessory2.linkage_external_vector
        
        if l1 == l2 and lv1[1] == -lv2[1]:
            return True
        else:
            return False
            
            





























    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    