#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:27:41 2021
<<<<<<< HEAD

=======
>>>>>>> dev_meshes
@author: dasilva
"""

import mechanical_components.optimization.meshes as meshes_opt
import dessia_common.workflow as wf
import mechanical_components.meshes as me

from dessia_api_client import Client
import numpy as npy

block_optimizer = wf.InstanciateModel(meshes_opt.MeshAssemblyOptimizer, name = 'Mesh Assemby Optimizer')

block_optimize= wf.ModelMethod(meshes_opt.MeshAssemblyOptimizer, 'Optimize', name = 'Optimizer')

block_attributs = wf.ModelAttribute(attribute_name= 'solutions', name= 'Solutions Mesh Assembly')

# block_attributs_MesComb = wf.ModelAttribute(attribute_name= 'mesh_combinations', name= 'Mesh Combinations')


# attributes = []
# display = wf.MultiPlot(attributes=attributes, order = 1, name = 'Display')

# block_workflow = [block_optimizer, block_optimize, display]
# pipe_workflow = [wf.Pipe(block_optimizer.outputs[0], block_optimize.inputs[0]),
#                  wf.Pipe(block_optimize.outputs[0], display.inputs[0])]

block_workflow = [block_optimizer, block_optimize, 
                   block_attributs,
                   # block_attributs_MesComb
                  ]

pipe_workflow = [wf.Pipe(block_optimizer.outputs[0], block_optimize.inputs[0]),
                  wf.Pipe(block_optimize.outputs[1], block_attributs.inputs[0]), 
                  # wf.Pipe(block_attributs.outputs[0], block_attributs_MesComb.inputs[0])
                 ]


workflow = wf.Workflow(block_workflow, pipe_workflow, 
                        block_attributs.outputs[0]
                        
                       )
workflow.plot_jointjs()

connections = [(0, 1)]


gear_speeds = {0: (1878.1453579221634, 1974.460504482274),
               1: (449.8807309958231, 472.95153771355757)}
               
center_dists = [(0.11134984458664793, 0.1293457790652981)]

torques = {0: 'output', 1: 40}

rack = meshes_opt.RackOpti(transverse_pressure_angle_0=[20/180.*npy.pi,20/180.*npy.pi], module=[2*1e-3,2*1e-3],
             coeff_gear_addendum=[1,1],coeff_gear_dedendum=[1.25,1.25],coeff_root_radius=[0.38,0.38],
             coeff_circular_tooth_thickness=[0.5,0.5],helix_angle=[21,60])
        
meshoptis = []
for i, speed_input in enumerate(gear_speeds.values()):
    meshoptis.append(meshes_opt.MeshOpti(rack = rack, torque_input = torques[i], speed_input = speed_input))

center_distances = []
for i , center_distance in enumerate(center_dists):

    center_distances.append(meshes_opt.CenterDistanceOpti(center_distance = center_distance, meshes = [meshoptis[connections[i][0]], meshoptis[connections[i][1]]]))


cycles = {0: 1272321481513.054}

input_values = {workflow.index(block_optimizer.inputs[0]): center_distances,
                workflow.index(block_optimizer.inputs[1]):  {0: 1272321481513.054},
                workflow.index(block_optimize.inputs[1]): 5,
                workflow.index(block_optimize.inputs[3]): True
                }

workflow_run = workflow.run(input_values)

d1 = workflow_run.to_dict()
obj = wf.WorkflowRun.dict_to_object(d1)


import json

object1=json.dumps(d1)

object2=json.loads(object1)

# ta_class.dict_to_object(object2)

c = Client(api_url = 'https://api.demo.dessia.tech')
r = c.create_object_from_python_object(workflow_run)
