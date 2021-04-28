#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 10:48:16 2021

@author: dasilva
"""
from dessia_api_client import Client
from dessia_common import workflow as wf

import plot_data
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
import mechanical_components
schaeffler_catalog = mechanical_components.models.schaeffler_catalog

block_optimizer = wf.InstanciateModel(bearings_opt.BearingAssemblyOptimizer,  name='BearingAssemblyOptimizer')
block_optimize = wf.ModelMethod(bearings_opt.BearingAssemblyOptimizer, 'optimize', name='BearingAssemblyOptimizer-optimize')
attribute_selection1 = wf.ModelAttribute('bearing_assembly_simulations') 

# filters = [
#           {'attribute' : 'bearing_assembly/overall_length', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly/mass', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly/cost', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly/number_bearing', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly/number_bearing_first_bc', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly/number_bearing_second_bc', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly_simulation_result/L10', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly_simulation_result/bearing_combination_first/max_axial_load', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly_simulation_result/bearing_combination_first/max_radial_load', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly_simulation_result/bearing_combination_second/max_axial_load', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_assembly_simulation_result/bearing_combination_second/max_radial_load', 'operator' : 'gt', 'bound' : -100},
#           ]

list_attribute = ['overall_length', 'mass', 'cost', 'number_bearing','number_bearing_first_bc','number_bearing_second_bc', 'first_max_axial_load','first_max_radial_load', 'second_max_axial_load','second_max_radial_load' ]


display = wf.MultiPlot(list_attribute,order = 1, name= 'Display')

# filter_analyze= wf.Filter(filters)
workflow_block = [block_optimizer, block_optimize, attribute_selection1, 
                  # filter_analyze,
                  display] 

workflow_pipe = [wf.Pipe(block_optimizer.outputs[0], block_optimize.inputs[0]),
                 wf.Pipe(block_optimize.outputs[1], attribute_selection1.inputs[0]),
                 # wf.Pipe(attribute_selection1.outputs[0], filter_analyze.inputs[0]), 
                 wf.Pipe(attribute_selection1.outputs[0], display.inputs[0])]

workflow = wf.Workflow(workflow_block, workflow_pipe, attribute_selection1.outputs[0])
workflow.plot_jointjs()

input_values = {workflow.index(block_optimizer.inputs[0]): [[[[0.1595, 0, 0], [0, -14000, 0], [0, 0, 0]]]],
                workflow.index(block_optimizer.inputs[1]): [157.07],
                workflow.index(block_optimizer.inputs[2]): [3600000],
                workflow.index(block_optimizer.inputs[3]): [0.035, 0.035],
                workflow.index(block_optimizer.inputs[4]): [0.072, 0.072],
                workflow.index(block_optimizer.inputs[5]): [0, 0.3],
                workflow.index(block_optimizer.inputs[6]): [0.1, 0.1],
                workflow.index(block_optimizer.inputs[7]): [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                                   bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                workflow.index(block_optimizer.inputs[8]): [bearings.CombinationMounting([bearings.Mounting(), bearings.Mounting(left=True)])],
                workflow.index(block_optimizer.inputs[9]): [[1, 2], [1, 2]],
                workflow.index(block_optimizer.inputs[12]): schaeffler_catalog,
                workflow.index(block_optimize.inputs[1]): 10,
                }

##
workflow_run = workflow.run(input_values)


