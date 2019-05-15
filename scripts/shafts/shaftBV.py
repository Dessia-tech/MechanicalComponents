#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:46:45 2019

@author: ringhausen
"""

import mechanical_components.shafts_assembly as ass
import matplotlib.pyplot as plt
import math
import random
import volmdlr as vm


#vis = shafts.TheoricalShaft([(0.05, 0.02), (0.02, 0.025), (0, 0.03)], (0, 1))
#roulement1 = shafts.TheoricalShaft([(0, 0.035), (0.03, 0.035), (0.03, 0.04)], (0, -1), [(0, 0.065), (0.03, 0.065)], (0, 1))
#pignon1 = shafts.TheoricalShaft([(0.03, 0.04), (0.08, 0.04),], (0, -1))
#circlip1 = shafts.TheoricalShaft([(0.08, 0.04), (0.08, 0.035), (0.085, 0.035), (0.085, 0.04)], (0, -1))
#synchro1 = shafts.TheoricalShaft([(0.085, 0.04), (0.115, 0.04)], (0, -1))
#pignon2 = shafts.TheoricalShaft([(0.115, 0.04), (0.165, 0.04)], (0, -1))
#circlip2 = shafts.TheoricalShaft([(0.165, 0.04), (0.165, 0.035), (0.17, 0.035), (0.17, 0.04)], (0, -1))
#pignon3 = shafts.TheoricalShaft([(0.17, 0.04), (0.22, 0.04)], (0, -1))
#circlip3 = shafts.TheoricalShaft([(0.22, 0.04), (0.22, 0.035), (0.225, 0.035), (0.225, 0.04)], (0, -1))
#synchro2 = shafts.TheoricalShaft([(0.225, 0.04), (0.25, 0.04)], (0, -1))
#pignon4 = shafts.TheoricalShaft([(0.25, 0.04), (0.30, 0.04), (0.30, 0.05)], (0, -1))
#
#pignon5 = shafts.TheoricalShaft([(0.30, 0.05) ,(0.30, 0.06), (0.33, 0.06), (0.33, 0.04)], (0, -1), montage='inbuilt')
#pignon6 = shafts.TheoricalShaft([(0.41, 0.04), (0.41, 0.05), (0.45, 0.05), (0.45, 0.045)], (0, -1), montage='inbuilt')
#
#roulement2 = shafts.TheoricalShaft([(0.45, 0.045), (0.45, 0.035), (0.48, 0.035)], (0, -1))
#joint1 = shafts.TheoricalShaft([(0.485, 0.03), (0.505, 0.03)], (0, -1))
#
#accessoires = [vis, roulement1, pignon1, circlip1, synchro1, pignon2,
#               circlip2, pignon3, circlip3, synchro2, pignon4, pignon5,
#               pignon6, roulement2, joint1]
#
##accessoires.reverse()
##accessoires.sort(key=lambda a: a.functional_points[0][1], reverse=True)
##accessoires = sorted(accessoires, key=lambda a: a.functional_points[0][0])
#
#functional_accessories = shafts.FunctionalAccessories(accessoires)
#shaft = shafts.Shaft(functional_accessories)
#
#print(shaft.ShaftCheck())
#shaft.Plot2()

###############################################################################

#plt.close('all')
#
#carter = ass.TheoricalShaft()
#shaft = ass.TheoricalShaft()
#vis = ass.TheoricalShaft()
#roulement1 = ass.TheoricalShaft()
#circlip1 = ass.TheoricalShaft()
#vis_test = ass.TheoricalShaft()
#roulement_test = ass.TheoricalShaft()
#vis_int = ass.TheoricalShaft()
#a = [carter, shaft, vis, roulement1, circlip1]
#
#s_circlip1 = ass.Surface.Instantiate([(0.08, 0.04), (0.08, 0.035), (0.085, 0.035), (0.085, 0.04)], circlip1, shaft)
#s_vis = ass.Surface.Instantiate([(0.05, 0.02), (0.02, 0.02), (0, 0.03)], vis, shaft)
#s_roulement1 = ass.Surface.Instantiate([(0, 0.035), (0.03, 0.035), (0.03, 0.04)], roulement1, shaft)
#s_roulement1_sup = ass.Surface.Instantiate([(0, 0.055), (0, 0.065), (0.03, 0.065)], carter, roulement1)
#s_vis_roulement = ass.Surface.Instantiate([(0, 0.035), (0, 0.045)], vis, roulement1)
#s_vis_test = ass.Surface.Instantiate([(0, 0.02), (0.02, 0.02), (0.05, 0.03)], shaft, vis_test)
#s_vis_test_2 = ass.Surface.Instantiate([(0.05, 0.035), (0.05, 0.045)], roulement_test, vis_test)

#s_vis_test = ass.Surface.Instantiate([(0, 0.02), (0.03, 0.02), (0.05, 0.03)], vis_test, shaft)
#s_vis_test_2 = ass.Surface.Instantiate([(0.05, 0.035), (0.05, 0.045)], vis_test, roulement_test) 
#
#s_vis_int = ass.Surface.Instantiate([(0.05, 0.02), (0.02, 0.02), (0, 0.03)], shaft, vis_int)
#s_vis_roulement_int = ass.Surface.Instantiate([(0, 0.035), (0, 0.045)], roulement1, vis_int)
#
#s_vis_int, s_vis_roulement_int, s_vis_test_2, s_vis_test, 

#s = [s_vis_roulement, s_vis, s_roulement1, s_roulement1_sup, s_circlip1]
#surfaces = ass.FunctionalSurfaces(s)
#surfaces.Plot()
#surfaces.PartPlot(roulement1)
#surfaces.PartPlot(vis)
#surfaces.SurfaceOrder(vis)
#surfaces.SurfaceOrder(roulement1)
#surfaces.SurfacePolygoneAngle(roulement1)
#print('\n')
#surfaces.PartPlot(vis_test)
#print('\n')
#surfaces.PartPlot(vis_int)

#surfaces.PartPlot(shaft)
#surfaces.SurfaceOrder(shaft)


#test = ass.TheoricalShaft()
#ext = ass.TheoricalShaft()
#s_test1 = ass.Surface.Instantiate([(-0.03, 0.01), (-0.02, 0.01)], ext, test)
#s_test2 = ass.Surface.Instantiate([(-0.005, 0.00), (0.005, 0.00)], ext, test)
#s_test3 = ass.Surface.Instantiate([(0.01, 0.01), (0.02, 0.01)], ext, test)
#s_test4 = ass.Surface.Instantiate([(-0.005, -0.02), (0.005, -0.02)], test, ext)
#
#surface = [s_test1, s_test2, s_test3, s_test4]
#fsurface = ass.FunctionalSurfaces(surface)
#fsurface.SurfaceOrder(test)



#test = ass.TheoricalShaft()
#ext = ass.TheoricalShaft()
#s_test1 = ass.Surface.Instantiate([(0,0.02), (0,0.03)], test, ext)
#s_test2 = ass.Surface.Instantiate([(0,0.01), (0.01,0)], ext, test)
#s_test3 = ass.Surface.Instantiate([(0.02,0), (0.025,0), (0.03,0)], ext, test)
#s_test4 = ass.Surface.Instantiate([(-0.005, -0.02), (0.005, -0.02)], test, ext)

#surface = [s_test1, s_test2, s_test3]
#fsurface = ass.FunctionalSurfaces(surface)
#lines = fsurface.SurfaceOrder(test)
#fsurface.PartPlot(test)
#print(lines)
#print()
#points = lines[0][3][1].points
#polygone = vm.Polygon2D(points)
#print(points)


###############################################################################


plt.close('all')

carter = ass.Part(name='carter')
shaft = ass.Part(name='shaft')
vis = ass.Part(name='vis')
roulement1 = ass.Part(name='roulement1')
circlip1 = ass.Part(montage='radial', name='cirlcip1')
vis_test = ass.Part(name='vis_test')
roulement_test = ass.Part(name='roulement_test')
vis_int = ass.Part(name='vis_int')
a = [carter, shaft, vis, roulement1, circlip1]

s_circlip1 = ass.Surface.Instantiate([(0.08, 0.04), (0.08, 0.035), (0.085, 0.035), (0.085, 0.04)], circlip1, shaft)
s_vis = ass.Surface.Instantiate([(0.05, 0.02), (0.02, 0.02), (0, 0.03)], vis, shaft)
s_roulement1 = ass.Surface.Instantiate([(0, 0.035), (0.03, 0.035), (0.03, 0.04)], roulement1, shaft)
s_roulement1_sup = ass.Surface.Instantiate([(0, 0.055), (0, 0.065), (0.03, 0.065)], carter, roulement1)
s_vis_roulement = ass.Surface.Instantiate([(0, 0.035), (0, 0.045)], vis, roulement1)


#polygon = vm.Polygon2D([(0,0), (1,0), (1,1), (0,1), (0,0)])
#print(polygon.SelfIntersect())



#shaft_carter = a[0].ViableShafts()
#shaft_shaft = a[1].ViableShafts()
#shaft_vis = a[2].ViableShafts()
#shaft_roulement1 = a[3].ViableShafts()
#shaft_circlip1 = a[4].ViableShafts()
#shafts = [shaft_shaft[0], shaft_vis[0], shaft_roulement1[0], shaft_circlip1[0], shaft_carter[0]]

shafts_product = ass.SurfaceRepertory([s_circlip1, s_roulement1, s_roulement1_sup, s_vis, s_vis_roulement]).shafts_product
#shafts = [shaft_shaft[0], shaft_vis[0], shaft_roulement1[0]]

assembly = ass.ShaftAssembly(shafts_product[0])

assemblies = assembly.viable_assembly_orders


#print(assemblies)
#for assemblyy in assemblies:
#    print('--')
#    assembly.GraphicMountability(assemblyy)
#    for shaft in assemblyy:
#        print(shaft.name)
#print(len(assemblies))
assembly.GraphicMountability(assemblies[0])































