#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:46:45 2019

@author: ringhausen
"""

#import MechanicalComponents.mechanical_components.shafts_assembly as ass
import mechanical_components.shafts_assembly as ass
import matplotlib.pyplot as plt
import math
import random
import volmdlr as vm

import time

###############################################################################

#
plt.close('all')
#
##carter = ass.Part(name='carter')
#shaft = ass.Part(name='shaft')
#vis = ass.Part(name='vis')
#roulement1 = ass.Part(name='roulement1')
#circlip1 = ass.Part(montage='radial', name='cirlcip1')
##circlip1 = ass.Part(montage='axial', name='cirlcip1')
#vis_test = ass.Part(name='vis_test')
#roulement_test = ass.Part(name='roulement_test')
#vis_int = ass.Part(name='vis_int')
##test = ass.Part(name='test')
##a = [carter, shaft, vis, roulement1, circlip1]
#
#s_circlip1 = ass.Surface.Instantiate([(0.08, 0.04), (0.08, 0.035), (0.085, 0.035), (0.085, 0.04)], circlip1, shaft)
#s_vis = ass.Surface.Instantiate([(0.05, 0.02), (0.02, 0.02), (0, 0.03)], vis, shaft)
#s_roulement1 = ass.Surface.Instantiate([(0, 0.035), (0.03, 0.035), (0.03, 0.04)], roulement1, shaft)
#s_roulement1_sup = ass.Surface.Instantiate([(0, 0.055), (0, 0.065), (0.03, 0.065)], None, roulement1)
#s_vis_roulement = ass.Surface.Instantiate([(0, 0.035), (0, 0.045)], vis, roulement1)
#
##contour_test = vm.primitives2D.RoundedLineSegments2D([(0, 0.055), (0, 0.065), (0.03, 0.065), (0.03, 0.08), (-0.01, 0.08), (-0.01, 0.055), (0, 0.055)], {})
#
##polygon = vm.Polygon2D([(0,0), (1,0), (1,1), (0,1), (0,0)])
##print(polygon.SelfIntersect())
#
#
##shaft_carter = a[0].ViableShafts()
##shaft_shaft = a[1].ViableShafts()
##shaft_vis = a[2].ViableShafts()
##shaft_roulement1 = a[3].ViableShafts()
##shaft_circlip1 = a[4].ViableShafts()
##shafts = [shaft_shaft[0], shaft_vis[0], shaft_roulement1[0], shaft_circlip1[0], shaft_carter[0]]
#
#
#shafts_product = ass.SurfaceRepertory([s_circlip1, s_roulement1, s_roulement1_sup, s_vis, s_vis_roulement]).shafts_product
##shafts_product = ass.SurfaceRepertory([s_roulement1, s_vis, s_vis_roulement]).shafts_product
##shafts_product = ass.SurfaceRepertory([s_circlip1]).shafts_product
#
#assembly = ass.ShaftAssembly(shafts_product[0])
#
#assemblies = assembly.viable_assembly_orders
#
#shafts_product[0][0].Wider(0.005)


# =============================================================================
# 
# =============================================================================

shaft = ass.Part(name='shaft')

#pts=[[-0.01,        0.23222834],
#[0.01 ,      0.23222834],
#[0.012  ,    0.23222834],
#[0.024   ,   0.23222834],
#[0.026  ,   0.2469398],
#[0.07067182 ,0.2469398 ],
#[0.19392024, 0.23222834],
#[0.20092024, 0.23222834],
#[0.20092024, 0.22722834],
#[0.19392024 ,0.22722834],
#[0.07067182, 0.2419398 ],
#[0.026 ,    0.2419398],
#[0.024 ,     0.22722834],
#[0.012  ,    0.22722834],
#[0.01   ,    0.22722834],
#[-0.01   ,     0.22722834]]

pts=[[0.014, 0.015],
[0.024, 0.015],
[0.08155265, 0.16405729],
[0.18008658, 0.16405729],
[0.18208658, 0.015     ],
[0.19408658, 0.015     ],
[0.19608658, 0.015     ],
[0.21608658, 0.015     ],
[0.21608658, 0.01      ],
[0.19608658, 0.01      ],
[0.19408658, 0.01      ],
[0.18208658, 0.01      ],
[0.18008658, 0.15905729],
[0.08155265, 0.15905729],
[0.024, 0.01 ],
[0.014, 0.01 ]]

points=[]
for pt in pts:
    points.append(vm.Point2D(tuple(pt)))
#points.reverse()
contour = vm.primitives2D.RoundedLineSegments2D(points, {}, closed=True)

s_shaft = ass.Surface.Instantiate(pts, None, shaft)
#s_shaft = ass.Surface.Instantiate(pts, shaft, None)
shaft_shaft = ass.Shaft(None, contour)
shaft_shaft.Plot()


debut = time.time()
shaft_shaft.Wider(0.01)
fin = time.time()
print('time = ', fin-debut)
print()
print(shaft_shaft.contour.points)



# =============================================================================
# 
# =============================================================================

#shaft = ass.Part(name='shaft')
#
#
##pts = [(0,0),
##       (0.05,0),
##       (0.05,0.05),
##       (0.07,0.05),
##       (0.07,0.06),
##       (0.03,0.06),
##       (0.03,0.03),
##       (0,0.03)]
#
#pts = [(0,0),
#       (0.05,0),
#       (0.05,0.05),
#       (0.07,0.05),
#       (0.07,0.06),
#       (0.04,0.06),
#       (0.03,0.03),
#       (0,0.03)]
#
#
#points=[]
#for pt in pts:
#    points.append(vm.Point2D(pt))
#contour = vm.primitives2D.RoundedLineSegments2D(points, {}, closed=True)
#s_shaft = ass.Surface.Instantiate(pts, shaft, None)
#shaft_shaft = ass.Shaft(None, contour)
#shaft_shaft.Plot()
#
#shaft_shaft.Wider(0.025)


















