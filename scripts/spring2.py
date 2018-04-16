#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 14:20:52 2018

@author: jezequel
"""
import mechanical_components.springs as springs
import matplotlib.pyplot as plt
import matplotlib.colors as colors

input_data = [{'F1' : 100, 'F2' : 500, 'stroke' : 0.005,
               'l1_max' : 0.100, 'r1' : 0.090, 'r2' : 0.120,
               'n_springs1' : 3, 'n_springs2' : 10, 'pattern' : 'circular'},
              {'k_precision' : 0.05, 'prod_volume' : 50}]

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

saor_d=saor.Dict()
print(saor_d)
import json
j=json.dumps(saor_d)

if catalog_spec:
    saor.CatalogStudy()
    
from dessia_api_client import Client

c=Client()
r=c.SubmitJob('mc_spring',input_data)