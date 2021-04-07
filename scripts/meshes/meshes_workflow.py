#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:27:41 2021

@author: dasilva
"""

import mechanical_components.optimization.meshes as meshes_opt
import dessia_common.workflow as wf

from dessia_api_client import Client
import numpy as np

block_optimizer = wf.InstanciateModel(meshes_opt.MeshAssemblyOptimizer, name = 'Mesh Assemby Optimizer')

block_optimize= wf.ModelMethod(meshes_opt.MeshAssemblyOptimizer, 'Optimize', name = 'Optimizer')

block_instantiate_meshes = wf.InstanciateModel(meshes_opt.Instanciate, name = ' Instanciate meshes')
block_instantiate = wf.ModelMethod(meshes_opt.Instanciate, 'instantiate', name = 'Instanciate' ) 


display = wf.MultiPlot(order=1, name = 'Display')

block_workflow = [block_optimizer, block_optimize, block_instantiate_meshes, block_instantiate, display]
pipe_workflow = [wf.Pipe(block_optimizer.outputs[0], block_optimize.inputs[0]),
                 wf.Pipe(block_instantiate.outputs[0], block_optimizer.inputs[0]),
                 wf.Pipe(block_instantiate_meshes.outputs[0], block_instantiate.inputs[0]),
                 wf.Pipe(block_optimize.outputs[0], display.inputs[1])]

workflow = wf.Workflow(block_workflow, pipe_workflow, block_optimize.outputs[0])

# block_rack = wf.InstantiateModel(meshes_opt.RackOpti, name = 'rack1')

# block_meshopti1 = wf.InstanciateModel(meshes_opt.MeshOpti, name = 'mesh opti1')
# block_meshopti2 = wf.InstanciateModel(meshes_opt.MeshOpti, name = 'mesh opti2')

# block_centerdistance = wf.InstanciateModel(meshes_opt.CenterDistanceOpti, name = 'Center Distance Optimization')

class Instanciate(DessiaObject):
    
    def __init__(self, gear_speeds: Dict, center_distances:List, torques: Dict):
        self.gear_speeds = gear_speeds
        self.center_distances = center_distances
        self.torques = torques 
        
        
    def instanciate(self):
        
        rack = RackOpti(transverse_pressure_angle_0=[20/180.*npy.pi,20/180.*npy.pi], module=[2*1e-3,2*1e-3],
             coeff_gear_addendum=[1,1],coeff_gear_dedendum=[1.25,1.25],coeff_root_radius=[0.38,0.38],
             coeff_circular_tooth_thickness=[0.5,0.5],helix_angle=[21,60])
        
        meshoptis = []
        for i, speed_input in enumerate(self.gear_speeds.values()):
            meshoptis.append(MeshOpti(rack = rack, torque_input = self.torques[i], speed_input = speed_input))
        
        center_distances = []
        for i, center_distance in enumerate(self.center_distances):
            center_distances.append([CenterDistanceOpti(center_distance = center_distance, meshes = [meshoptis[i], meshoptis[i+1]])])
            
        
        return center_distances
            
            
            
        

        