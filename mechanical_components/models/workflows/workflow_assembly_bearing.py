#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
from dessia_common import DessiaObject, dict_merge
from dessia_common import Evolution, CombinationEvolution
import dectree
from dessia_api_client import Client
from dessia_common import workflow as wf
import dessia_common as dc
from volmdlr import plot_data

import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
from mechanical_components.models.catalogs import schaeffler_catalog

from itertools import product
import networkx as nx
from random import random
from scipy.optimize import minimize, fsolve
from dataclasses import dataclass
from numpy import allclose
import matplotlib.pyplot as plt
from scipy import interpolate
import math
import json
import pkg_resources

blockA = wf.InstanciateModel(bearings_opt.BearingAssemblyOptimizer,  name='BearingAssemblyOptimizer')
optimizeA = wf.ModelMethod(bearings_opt.BearingAssemblyOptimizer, 'optimize', name='BearingAssemblyOptimizer-optimize')
attribute_selection1 = wf.ModelAttribute('bearing_assembly_simulations')

filters = [
          {'attribute' : 'bearing_assembly.overall_length', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly.mass', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly.cost', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly.number_bearing', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly.number_bearing_first_bc', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly.number_bearing_second_bc', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly_simulation_result.L10', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly_simulation_result.bearing_combination_first.max_axial_load', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly_simulation_result.bearing_combination_first.max_radial_load', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly_simulation_result.bearing_combination_second.max_axial_load', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_assembly_simulation_result.bearing_combination_second.max_radial_load', 'operator' : 'gt', 'bound' : -100},
           ]

filter_analyze= wf.Filter(filters)

input_values = {}
blocks = []

blocks.extend([blockA, optimizeA, attribute_selection1, 
                filter_analyze
                ])

pipes = [wf.Pipe(blockA.outputs[0], optimizeA.inputs[0]),
         wf.Pipe(optimizeA.outputs[1], attribute_selection1.inputs[0]),
          wf.Pipe(attribute_selection1.outputs[0], filter_analyze.inputs[0]),
         ]

workflow = wf.Workflow(blocks, pipes, filter_analyze.outputs[0])
    
input_values = {workflow.index(blockA.inputs[0]): [[[[0.1595, 0, 0], [0, -14000, 0], [0, 0, 0]]]],
                workflow.index(blockA.inputs[1]): [157.07],
                workflow.index(blockA.inputs[2]): [3600000],
                workflow.index(blockA.inputs[3]): [0.035, 0.035],
                workflow.index(blockA.inputs[4]): [0.072, 0.072],
                workflow.index(blockA.inputs[5]): [0, 0.3], 
                workflow.index(blockA.inputs[6]): [0.03, 0.03],
                workflow.index(blockA.inputs[7]): [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                                   bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                workflow.index(blockA.inputs[8]): [bearings.CombinationMounting([bearings.Mounting(), bearings.Mounting(left=True)])],
                workflow.index(blockA.inputs[9]): [[1], [1]],
                workflow.index(blockA.inputs[12]): schaeffler_catalog,
                workflow.index(optimizeA.inputs[1]): 10,
                }

a = workflow.to_dict()
obj = wf.Workflow.dict_to_object(a)
##
workflow_assembly_bearing = workflow.run(input_values)

    