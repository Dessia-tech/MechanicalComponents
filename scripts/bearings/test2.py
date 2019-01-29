#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.bearings as bearings
import numpy as npy
import volmdlr as vm
import copy
from mechanical_components.load import *

b0 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
b2 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
b3 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
b4 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = 1)
b5 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
list_bearing = [b1, b2, b3, b4]
BA = bearings.BearingCombination(list_bearing, radial_load_linkage = [True]*4, internal_pre_load = 0, 
                 connection_bi = ['n', 'p'], connection_be = ['p'], behavior_link = 'pn')

#k1 = 1e8
#k2 = 1e8
#matrice0 = npy.array([[k1, 0, -k1, 0], [0, k1+k2, -k2, -k1], [-k1, -k2, k1+k2, 0],
#                     [0, -k1, 0, k1]])
#mat = npy.zeros((6, 6))
#mat[0:4, 0:4] = mat[0:4, 0:4] + matrice0
#mat[2:6, 2:6] = mat[2:6, 2:6] + matrice0
#
#mat1 = copy.deepcopy(mat)
#mat1[0, 0:] = [1, 0, 0, 0, 0, 0]
#mat1[1, 0:] = [0, 1, 0, 0, 0, 0]
#b = npy.array([0, 0, 0, 0, -1, 0])
#x = npy.linalg.solve(mat1, b)
#load0 = npy.array([[k1, 0, -k1, 0], [0, k2, -k2, 0], [0, k1, 0, -k1]])
#ld = npy.zeros((6, 6))
#ld[0:3, 0:4] = ld[0:3, 0:4] + load0
#ld[3:6, 2:6] = ld[3:6, 2:6] + load0
#print(npy.dot(ld, x))
#print(npy.dot(mat, x))

k1 = 1e5
jm = 0.
k2 = 1e4
j2 = 0.
km = 1e6
f0 = 1e2
ka = 1e4

tab_pos = {0: 0,
           1: 0,
           2: 1,
           3: 1,
           4: 2,
           5: 2,
           6: 3, 7: 3,
           8: 4, 9: 4,
           10: 5, 11: 5,
           12: 6, 13: 6, 
           14: 7, 15: 7}

connection = {0: {2: [km, jm, 0.]},
              1: {3: [km, jm, 0.]},
              2: {0: [km, jm, 0.], 4: [k1, 0., 0.]},
              3: {1: [km, jm, 0.], 4: [k2, j2, 0.], 5: [k1, 0., 0.]},
              4: {3: [k2, j2, 0.], 2: [k1, 0., 0.], 6: [km, jm, 0.]},
              5: {3: [k1, 0., 0.], 7: [km, jm, 0.]},
              6: {4: [km, jm, 0.], 8: [k1, 0., 0.]},
              7: {5: [km, jm, 0.], 8: [k2, j2, 0.], 9: [k1, 0., 0.]},
              8: {7: [k2, j2, 0.], 6: [k1, 0., 0.], 10: [km, jm, 0.]},
              9: {7: [k1, 0., 0.]},
              10: {8: [km, jm, 0.], 12: [k1, 0., 0.], 13: [ka, 0., f0]},
              11: {13: [k1, 0., 0.]},
              12: {10: [k1, 0., 0.]},
              13: {11: [k1, 0., 0.], 10: [ka, 0., f0], 15: [km, jm, 0.]},
              15: {13: [km, jm, 0.]}
              }
connection_init = copy.deepcopy(connection)

#connection = {0: {2: [km, jm]},
#              2: {0: [km, jm]},
#              }

#vect = list(connection.keys())
#matrice = npy.zeros((len(vect), len(vect)))
#f = npy.zeros(len(vect))
#for node1, li_data in connection.items():
#    pos1 = vect.index(node1)
#    for node2, (k, j) in li_data.items():
#        pos2 = vect.index(node2)
#        matrice[pos1, pos1] = matrice[pos1, pos1] + k
#        matrice[pos1, pos2] = matrice[pos1, pos2] - k
#        if tab_pos[node1] < tab_pos[node2]:
#            f[pos1] = f[pos1] + k*j
#        else:
#            f[pos1] = f[pos1] - k*j
#        
#fixe = [0, 1]
#external_load = {4: -1e6}
#for i in fixe:
#    pos_i = vect.index(i)
#    matrice[pos_i, :] = [0]*len(vect)
#    matrice[pos_i, pos_i] = 1
#    f[pos_i] = 0
#    
#for node, l in external_load.items():
#    pos = vect.index(node)
#    f[pos] = f[pos] + l
#    
#x = npy.linalg.solve(matrice, f)

valid = False
fix = [0, 15]
load = {2: 110}

SM1 = StiffnessMatrix(connection, tab_pos, fix, load)

while valid:
    valid_int = True
    vect = list(connection.keys())
    matrice = npy.zeros((len(vect), len(vect)))
    f = npy.zeros(len(vect))
    for node1, li_data in connection.items():
        pos1 = vect.index(node1)
        for node2, (k, j, f0) in li_data.items():
            pos2 = vect.index(node2)
            matrice[pos1, pos1] = matrice[pos1, pos1] + k
            matrice[pos1, pos2] = matrice[pos1, pos2] - k
            if tab_pos[node1] < tab_pos[node2]:
                f[pos1] = f[pos1] + k*j
            else:
                f[pos1] = f[pos1] - k*j
            
    
    for i in fixe:
        if i in vect:
            pos_i = vect.index(i)
            matrice[pos_i, :] = [0]*len(vect)
            matrice[pos_i, pos_i] = 1
            f[pos_i] = 0
        
    for node, l in external_load.items():
        pos = vect.index(node)
        f[pos] = f[pos] + l
        
    x = npy.linalg.solve(matrice, f)
    list_del = []
    for node1, li_data in connection.items():
        pos1 = vect.index(node1)
        for node2, (k, j, f0) in li_data.items():
            pos2 = vect.index(node2)
            if tab_pos[node1] < tab_pos[node2]:
                load = k*(x[pos1] - x[pos2]) - k*j
                dist = x[pos1] - x[pos2]
            else:
                load = -k*(x[pos1] - x[pos2]) - k*j
                dist = -(x[pos1] - x[pos2])
#            if f0 == 0:
#                if load < 0:
#                    valid_int = False
#                    list_del.append((node1, node2))
            if f0 != 0:
                print(dist)
                external_m = copy.deepcopy(external_load)
                if tab_pos[node1] < tab_pos[node2]:
                    a = 1*(f0 - k*dist)
                    print(a)
                    external_load[node2] = round(a, 2)
                else:
                    a = 1*(-f0 + k*dist)
                    print(a)
                    external_load[node2] = round(a, 2)
                
                if external_m != external_load:
                    valid_int = False
            if tab_pos[node1] < tab_pos[node2]:
                print('internal load ({}, {}) {}'.format(node1, node2, load))
                
    list_add = []
    for node1, li_data in connection_init.items():
        for node2, (k, j, f0) in li_data.items():
            if f0 != 0:
                pos1 = vect.index(node1)
                pos2 = vect.index(node2)
                if tab_pos[node1] < tab_pos[node2]:
                    dist = x[pos1] - x[pos2]
                else:
                    dist = -(x[pos1] - x[pos2])
                if dist >= 0:
                    list_add.append((node1, node2))
            
    for nd1, nd2 in list_del:
        del connection[nd1][nd2]
        print('del ({}, {})'.format(nd1, nd2))
    list_del = []
    for node1, li_data in connection.items():
        if li_data == {}:
            list_del.append(node1)
    for nd in list_del:
        del connection[nd]
        print('del node {}'.format(nd))
        
    drap = True
    input_node = list(external_load.keys())
    for node1, li_data in connection.items():
        pos1 = vect.index(node1)
        for node2, (k, j, f0) in li_data.items():
            pos2 = vect.index(node2)
            if f0 != 0:
                input_node.extend([node1, node2])
    input_node = list(set(input_node))
    print(11, input_node)
    while drap:
        input_plus = copy.deepcopy(input_node)
        for nd in input_node:
            input_plus.extend(list(connection[nd].keys()))
        print(set(input_node), set(input_plus))
        if set(input_node) == set(input_plus):
            drap = False        
        input_node = list(set(input_plus))
    print(input_node)
    connection_plus = {}
    for nd in input_node:
        connection_plus[nd] = connection[nd]
    connection = copy.deepcopy(connection_plus)
    
    if valid_int == True:
        valid = False
    
        
print(npy.dot(matrice, x))
#for li_bg, nx_graph in zip(BA.graph, BA.nx_graph):
#    print(li_bg[0].load_bearing_results[0].sub_direction)
##    print(nx_graph[0].edges(data=True))
#    i1 = li_bg[0].load_bearing_results[0].list_node[3]
#    i2 = li_bg[-1].load_bearing_results[0].list_node[5]
##    print(i1, i2)
#    sol = nx.all_shortest_paths(nx_graph[0], source=i1,target=i2)
#    try:
#        print(11,len([p for p in sol]))
#    except:
#        pass

#BA.PlotGraph()
#fa = BA.SearchBestGraph()
#BA.BearingCombinationLoad(fr=1, fa=0)
#axial_load, axial_pre_load = BA.CheckViabilityAxialPath(list_bearing, 0)
#d = BA.Dict()
#obj = bearings.BearingCombination.DictToObject(d)
#d = obj.Dict()
#print(d['bearings'])
#obj.Plot(typ=None, box=False)
#print(obj.PlotData(typ='Load'))
#BA.Plot(box = False, typ = 'Load')

#BA.PlotGraph()

#d = BA.bearings[3].Dict()
#obj = bearings.RadialBearing.DictToObject(d)
#obj.Plot()
#import json
##print(json.dumps(d))
#
#sol = bearings.BearingCombination.DictToObject(d)
#sol.Plot(typ='Load', box=False)
#
#bg = BA.bearings_solution[0].Plot(typ=None)
#export = BA.bearings[2].PlotData()
#
#sol = BA.PlotData()
##print(BA.PlotD3())
##print(json.dumps(export))
#print(json.dumps(BA.PlotData()))
#
##ax = bg.MPLPlot(style='-ob')
##li = []
##for item in bg.basis_primitives:
##    if 'Arc2D' in str(item.__class__):
##        li.extend(item.Discret())
##c=vm.Contour2D(li)
##c.MPLPlot(style='ob')