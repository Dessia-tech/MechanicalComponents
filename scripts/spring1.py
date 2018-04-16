#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 18:01:54 2018

@author: jezequel
"""
import math
import mechanical_components.springs as springs
import numpy as npy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
from scipy.spatial import ConvexHull
from time import time
from bs4 import BeautifulSoup
import pandas as pd
from dessia_api_client import Client

F1 = 100
F2 = 500
stroke = 0.005
l1_max = 0.1
k_percent = 0.100
prod_volume = 50
r1 = 0.090
r2 = 0.120
n_springs = [i + 3 for i in range(8)]

# =============================================================================
# Catalog optimization
# =============================================================================
#co = springs.CatalogOptimizer(springs.ferroflex_catalog, F1, F2, stroke, k_percent, l1_max, r1, r2, n_springs, 'circular')
#dictionnary = {springs.ferroflex_catalog.name : co.opti_indices}
#cor = springs.CatalogOptimizationResults(dictionnary, springs.catalogs, 50)
#list_springs_assembly = cor.results

#products = [product.catalog.Price(product.product_index, 20)  for res in cor.results for product in res.products]
#products




# =============================================================================
# Assembly optimization
# =============================================================================
sao = springs.SpringAssemblyOptimizer(F1, F2, stroke, n_springs, r1, r2, l1_max, 'circular')

l_assemb = [assembly.n_springs for assembly in sao.assemblies]
k_assemb = [spring.Stiffness() for assembly in sao.assemblies for spring in assembly.springs]
D_assemb = [spring.D for assembly in sao.assemblies for spring in assembly.springs]
d_assemb = [spring.d for assembly in sao.assemblies for spring in assembly.springs]
m_assemb = [assembly.Mass() for assembly in sao.assemblies]
cost_assemb = [assembly.Cost() for assembly in sao.assemblies]
l0_assemb = [spring.Length(0) for assembly in sao.assemblies for i, spring in enumerate(assembly.springs) if i == 0]
lw_assemb = [spring.WireLength() for assembly in sao.assemblies for i, spring in enumerate(assembly.springs) if i == 0]
test_assemb = [spring.ResultingForce(spring.l_min)/spring.Stiffness() + spring.lc for assembly in sao.assemblies for i, spring in enumerate(assembly.springs) if i == 0]

m_assemb3 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 3]
m_assemb4 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 4]
m_assemb5 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 5]
m_assemb6 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 6]
m_assemb7 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 7]
m_assemb8 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 8]
m_assemb9 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 9]
m_assemb10 = [assembly.Mass() for assembly in sao.assemblies if assembly.n_springs == 10]
cost_assemb3 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 3]
cost_assemb4 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 4]
cost_assemb5 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 5]
cost_assemb6 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 6]
cost_assemb7 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 7]
cost_assemb8 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 8]
cost_assemb9 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 9]
cost_assemb10 = [assembly.Cost() for assembly in sao.assemblies if assembly.n_springs == 10]

materials = [springs.steel1,
             springs.steel2,
             springs.steel3,
             springs.steel4,
             springs.copper_tin_alloy,
             springs.copper_zinc_alloy]

n_spires = npy.linspace(2.5, 9.5, 8)
diametres_mm = [0.12, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
                0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.8, 3, 3.2, 3.4, 3.5, 3.6, 3.8,
                4, 4.2, 4.5, 4.70, 4.8, 5, 5.3, 5.5, 5.6,6, 6.3, 6.5, 6.70, 7, 7.5,
                8, 8.5, 9, 9.5, 10, 11, 12, 13, 14]

diametres_m = [d*10**-3 for d in diametres_mm]
security_factor = 0.80
springs_ok = []
all_springs = []
resultats = []
for d in diametres_m:
    for n in n_spires:
        
        for material in materials:
            sdo = springs.SpringDiscreteOptimizer(F1, F2, stroke, d, n, material)
            all_springs.append(springs.Spring(sdo.D, d, n, sdo.l0, material))
            
            if sdo.tau_k < material.tau_max\
            and d >= material.d_min and d <= material.d_max\
            and sdo.D/d > 5 and sdo.D/d < 13\
            and (sdo.D+d)/(sdo.D-d) > 1.4 and (sdo.D+d)/(sdo.D-d) < 2\
            and sdo.l1 < l1_max:
                spring = springs.Spring(sdo.D, d, n, sdo.l0, material)
                so = springs.SpringOptimizer(spring, {'l0' : (0.001, 5)}, sdo.F2)
                res = so.Optimize()
                
                if res['success']:
                    springs_ok.append(so.spring)
                    resultats.append(res)
                
D = [s.D for s in springs_ok]
d = [s.d for s in springs_ok]
#w = [s.D/s.d for s in springs_ok]
#a = [(s.D+s.d)/(s.D-s.d) for s in springs_ok]
#n = [s.n for s in springs_ok]
#k = [s.Stiffness() for s in springs_ok]
#i = [s.InclinationAngle()*180/math.pi for s in springs_ok]
#tan_i = [math.tan(s.InclinationAngle()) for s in springs_ok]
m = [s.Mass() for s in springs_ok]
#lw = [s.WireLength() for s in springs_ok]
#Sa = [s.Sa for s in springs_ok]
#l0 = [s.l0 for s in springs_ok]

Da = [s.D for s in all_springs]
da = [s.d for s in all_springs]

points = npy.array([[s.d, s.D] for s in springs_ok])
#
#fig = plt.figure()
#
## XC60
#Dc = [s.D for s in springs_ok if s.material.name == 'XC60']
#dc = [s.d for s in springs_ok if s.material.name == 'XC60']
#Dac = [s.D for s in all_springs if s.material.name == 'XC60']
#dac = [s.d for s in all_springs if s.material.name == 'XC60']
#pc = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == 'XC60'])
#
#ac = fig.add_subplot(321)
#ac.set_title('XC60')
#ac.set_xlabel('d (m)')
#ac.set_ylabel('D (m)')
#ac.plot(da, Da, '.', color = 'gray')
#ac.plot(dac, Dac, 'k.')
#ac.plot(dc, Dc, 'co')
#ac.plot(d_assemb, D_assemb, 'o', color = 'darkcyan')
#    
## XC95
#Db = [s.D for s in springs_ok if s.material.name == 'XC95']
#db = [s.d for s in springs_ok if s.material.name == 'XC95']
#Dab = [s.D for s in all_springs if s.material.name == 'XC95']
#dab = [s.d for s in all_springs if s.material.name == 'XC95']
#pb = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == 'XC95'])
#    
#ab = fig.add_subplot(322)
#ab.set_title('XC95')
#ab.set_xlabel('d (m)')
#ab.set_ylabel('D (m)')
#ab.plot(da, Da, '.', color = 'gray')
#ab.plot(dab, Dab, 'k.')
#ab.plot(db, Db, 'bo')
#ab.plot(d_assemb, D_assemb, 'o', color = 'navy')
#
## 50CV4
#Dg = [s.D for s in springs_ok if s.material.name == '50CV4']
#dg = [s.d for s in springs_ok if s.material.name == '50CV4']
#Dag = [s.D for s in all_springs if s.material.name == '50CV4']
#dag = [s.d for s in all_springs if s.material.name == '50CV4']
#pg = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == '50CV4'])
#    
#ag = fig.add_subplot(323)
#ag.set_title('50CV4')
#ag.set_xlabel('d (m)')
#ag.set_ylabel('D (m)')
#ag.plot(da, Da, '.', color = 'gray')
#ag.plot(dag, Dag, 'k.')
#ag.plot(dg, Dg, 'go')
#ag.plot(d_assemb, D_assemb, 'o', color = 'darkgreen')
#
## Z15CN17-03
#Dy = [s.D for s in springs_ok if s.material.name == 'Z15CN17-03']
#dy = [s.d for s in springs_ok if s.material.name == 'Z15CN17-03']
#Day = [s.D for s in all_springs if s.material.name == 'Z15CN17-03']
#day = [s.d for s in all_springs if s.material.name == 'Z15CN17-03']
#py = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == 'Z15CN17-03'])
#    
#ay = fig.add_subplot(324)
#ay.set_title('Z15CN17-03')
#ay.set_xlabel('d (m)')
#ay.set_ylabel('D (m)')
#ay.plot(da, Da, '.', color = 'gray')
#ay.plot(day, Day, 'k.')
#ay.plot(dy, Dy, 'o', color = 'yellow')
#ay.plot(d_assemb, D_assemb, 'o', color = 'y')
#
## CuSN6 R950
#Do = [s.D for s in springs_ok if s.material.name == 'CuSn6 R950']
#do = [s.d for s in springs_ok if s.material.name == 'CuSn6 R950']
#Dao = [s.D for s in all_springs if s.material.name == 'CuSn6 R950']
#dao = [s.d for s in all_springs if s.material.name == 'CuSn6 R950']
#po = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == 'CuSn6 R950'])
##hullo = ConvexHull(po)
#    
#ao = fig.add_subplot(325)
#ao.set_title('CuSn6 R950')
#ao.set_xlabel('d (m)')
#ao.set_ylabel('D (m)')
#ao.plot(da, Da, '.', color = 'gray')
#ao.plot(dao, Dao, 'k.')
#ao.plot(do, Do, 'o', color = 'orange')
#ao.plot(d_assemb, D_assemb, 'o', color = 'darkorange')
#
## CuZn36 R700
#Dr = [s.D for s in springs_ok if s.material.name == 'CuZn36 R700']
#dr = [s.d for s in springs_ok if s.material.name == 'CuZn36 R700']
#Dar = [s.D for s in all_springs if s.material.name == 'CuZn36 R700']
#dar = [s.d for s in all_springs if s.material.name == 'CuZn36 R700']
#pr = npy.array([[s.d, s.D] for s in springs_ok if s.material.name == 'CuZn36 R700'])
#    
#ar = fig.add_subplot(326)
#ar.set_title('CuSn36 R700')
#ar.set_xlabel('d (m)')
#ar.set_ylabel('D (m)')
#ar.plot(da, Da, '.', color = 'gray')
#ar.plot(dar, Dar, 'k.')
#ar.plot(dr, Dr, 'ro')
#ar.plot(d_assemb, D_assemb, 'o', color = 'darkred')

#c = springs.Catalog('ferroflex/catalog_SI', prod_volume)
#c.CorrectionDynParameters()

#catalogs = co.opti_catalogs

#fig = plt.figure()
#plt.plot(da, Da, '.', color = 'darkgrey')
#plt.plot(d, D, '.', color = 'grey')
#plt.plot(d_assemb, D_assemb, 'k.')
#[plt.plot(catalog.products['d'], catalog.products['D'], '.', color = colors.hsv_to_rgb((1/(i+1), 0.7, 0.7)), label = str(n_springs[i])) for i, catalog in enumerate(catalogs)]
##plt.plot(c.catalog_init['d'], c.catalog_init['R'], 'k.')
##plt.plot(co.opti_cat.catalog['d'], co.opti_cat.catalog['R'],'r.')
#plt.xlabel('d')
#plt.ylabel('D')
#plt.title('F1 = {0}, F2 = {1}, s = {2}'.format(F1, F2, stroke))
#plt.legend()
#fig.canvas.set_window_title('%k = ' + str(k_percent*100))


#fig2 = plt.figure()
#plt.plot(cost_assemb3, m_assemb3, 'c.', label = '3')
#plt.plot(cost_assemb4, m_assemb4, 'b.', label = '4')
#plt.plot(cost_assemb5, m_assemb5, 'g.', label = '5')
#plt.plot(cost_assemb6, m_assemb6, 'y.', label = '6')
#plt.plot(cost_assemb7, m_assemb7, '.', color = 'orange', label = '7')
#plt.plot(cost_assemb8, m_assemb8, 'r.', label = '8')
#plt.plot(cost_assemb9, m_assemb9, '.', color = 'gray', label = '9')
#plt.plot(cost_assemb10, m_assemb10, 'k.', label = '10')
#plt.xlabel('Cost')
#plt.ylabel('Mass')
#plt.legend()
#tf = time()
#print('t = ', tf - td, 's')

#def pareto_frontier_multi_max(myArray):
#    # Sort on first dimension (descending value) (weighted sum method using the sum of all elements)
#    myArray[:] = myArray[npy.array([sum(x) for x in myArray]).argsort()][::-1]
#    pareto_frontier = None
#    # Test next rows against the last row in pareto_frontier
#    while True:
#        # The top element is by definition Pareto (or it would have been removed by domination)
#        if(pareto_frontier == None):
#            pareto_frontier = myArray[0:1,:]
#        else:
#            pareto_frontier = npy.concatenate((pareto_frontier, myArray[0:1,:]))
# 
#        # remove the Pareto point that we've added to the Pareto frontier
#        myArray = npy.delete(myArray, 0, 0)
#        
#        rowNr = 0
#        while len(myArray) != 0 and rowNr < len(myArray):
#            row = myArray[rowNr:rowNr+1,:][0]
#            if sum([row[x] <= pareto_frontier[-1][x]
#                    for x in range(len(row))]) == len(row):
#                # If it is worse on all features remove the row from the array
#                myArray = npy.delete(myArray, rowNr, 0)
#            else:
#                rowNr += 1
#                
#        if len(myArray) == 0:
#            break
#    return pareto_frontier
#
#def pareto_frontier(Xs, Ys, maxX = True, maxY = True):
#    # Sort the list in either ascending or descending order of X
#    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
#
#    # Start the Pareto frontier with the first value in the sorted list
#    p_front = [myList[0]]    
#
#    # Loop through the sorted list
#    for pair in myList[1:]:
#        if maxY: 
#            if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
#                p_front.append(pair) # … and add them to the Pareto frontier
#        else:
#            if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
#                p_front.append(pair) # … and add them to the Pareto frontier
#
#    # Turn resulting pairs back into a list of Xs and Ys
#    p_frontX = [pair[0] for pair in p_front]
#    p_frontY = [pair[1] for pair in p_front]
#    index = [i for i, input_pair in enumerate(myList) if [Xs[i], Ys[i]] in p_front]
#    return p_frontX, p_frontY, index

#def yields(row, candidateRow):
#    return sum([row[x] <= candidateRow[x] for x in range(len(row))]) == len(row) 
#
#def simple_cull(inputPoints, yields):
#    paretoPoints = set()
#    candidateRowNr = 0
#    dominatedPoints = set()
#    while True:
#        candidateRow = inputPoints[candidateRowNr]
#        inputPoints.remove(candidateRow)
#        rowNr = 0
#        nonDominated = True
#        while len(inputPoints) != 0 and rowNr < len(inputPoints):
#            row = inputPoints[rowNr]
#            if dominates(candidateRow, row):
#                # If it is worse on all features remove the row from the array
#                inputPoints.remove(row)
#                dominatedPoints.add(tuple(row))
#            elif dominates(row, candidateRow):
#                nonDominated = False
#                dominatedPoints.add(tuple(candidateRow))
#                rowNr += 1
#            else:
#                rowNr += 1
#
#        if nonDominated:
#            # add the non-dominated point to the Pareto frontier
#            paretoPoints.add(tuple(candidateRow))
#
#        if len(inputPoints) == 0:
#            break
#    return paretoPoints, dominatedPoints

#array_test = [[cost_assemb[i], m_assemb[i]] for i in range(len(cost_assemb))]

#res_x, res_y, ind = pareto_frontier(cost_assemb, m_assemb, maxX = False, maxY = False)
#p_front, dominated = simple_cull(array_test, yields)
#res = list(p_front)
#res_x = [r[0] for r in res]
#res_y = [r[1] for r in res]
#plt.plot(cost_assemb, m_assemb, 'b.')
#plt.plot(res_x, res_y, 'r.')
#plt.plot(res_x, res_y, 'r')

#results = [sao.assemblies[i] for i in ind]


#plt.plot(cost_assemb, m_assemb, 'b.')
#plt.plot(res_x, res_y, 'r')