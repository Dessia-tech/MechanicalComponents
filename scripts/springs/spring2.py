#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 14:20:52 2018

@author: jezequel
"""
import mechanical_components.springs as springs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import json
import jsonpickle
#import pickle
from dessia_api_client import Client

input_data = [{'F1' : 100, 'F2' : 500, 'stroke' : 0.005,
               'l1_max' : 0.100, 'r1' : 0.090, 'r2' : 0.120,
               'n_springs1' : 3, 'n_springs2' : 10, 'pattern' : 'circular'},
              {'stiffness_precision' : 0.05, 'prod_volume' : 50}]

spring_spec = input_data[0]
catalog_spec = input_data[1]

n_springs = [i + spring_spec['n_springs1'] for i in range(spring_spec['n_springs2'] - spring_spec['n_springs1'] + 1)]

# =============================================================================
# Assembly optimization
# =============================================================================
sao = springs.SpringAssemblyOptimizer(spring_spec['F1'],
                                      spring_spec['F2'],
                                      spring_spec['stroke'],
                                      n_springs,
                                      spring_spec['r1'],
                                      spring_spec['r2'],
                                      spring_spec['l1_max'],
                                      spring_spec['pattern'].lower())

saor = springs.SpringAssemblyOptimizationResults(sao.assemblies, input_data)

#for r in saor.results:
#    mpas = r.matching_product_assemblies
#    print('--')
#    for m in mpas:
#        products = m.products
#        for p in products:
#            print(p.product_index)

#s2=pickle.dumps(saor)
    
saor_d=saor.Dict()
i = []
ps = []
for r in saor_d['results']:
    mpas = r['matching_product_assemblies']
    for m in mpas:
        products = m['products']
        if products[0] not in ps:
            ps.append(products[0])
        for p in products:
            i.append(p['product_index'])
                    
j=json.dumps(saor_d)
s=jsonpickle.dumps(saor)
#pickle.dumps(saor)

# =============================================================================
# Export FreeCAD
# =============================================================================
#sa = saor.results[-1]
#sa.CADExport('spring1',export_types=['fcstd'])

#c=Client()
#r=c.SubmitJob('mc_spring_assembly',input_data)
