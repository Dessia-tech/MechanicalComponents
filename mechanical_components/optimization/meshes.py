#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""

"""
import numpy as npy
from dessia_common import DessiaObject
from typing import  List, Tuple, Dict
from mechanical_components.meshes import MeshAssembly, hardened_alloy_steel,\
        gear_graph_simple,Material

import math
from typing import  List, Tuple
import copy

try:
    _open_source = True
    import mechanical_components.optimization.meshes_protected as protected_module
except (ModuleNotFoundError, ImportError) as _:
    _open_source = False

#class ContinuousMeshesAssemblyOptimizer(protected_module.ProtectedContinuousMeshesAssemblyOptimizer if _open_source==True else object):


class RackOpti(DessiaObject):
     _standalone_in_db = True

     def __init__(self, transverse_pressure_angle_0: Tuple[float,float]=None, module: Tuple[float,float]=None,
                 coeff_gear_addendum : Tuple[float,float]=None, coeff_gear_dedendum: Tuple[float,float]=None,
                 coeff_root_radius: Tuple[float,float]=None, coeff_circular_tooth_thickness: Tuple[float,float]=None,
                 helix_angle: Tuple[float,float]=None ,
                 name : str=''):

         self.transverse_pressure_angle_0=transverse_pressure_angle_0
         self.module=module
         if helix_angle==None:
             helix_angle=[0,0]
         self.helix_angle=helix_angle
         
         self.coeff_gear_addendum=coeff_gear_addendum
         self.coeff_gear_dedendum=coeff_gear_dedendum
         self.coeff_root_radius=coeff_root_radius
         self.coeff_circular_tooth_thickness=coeff_circular_tooth_thickness
         self.list_gear=[]
        
         self.name=name
         
         DessiaObject.__init__(self, name=name)

class MeshOpti(DessiaObject):
    _standalone_in_db = True

    def __init__(self,torque_input: float,speed_input : Tuple[float,float],Z:int=0 ,rack: RackOpti=None,gearing_interior:str='False',
                 coefficient_profile_shift: Tuple[float,float]=None,material: Material=None,transverse_pressure_angle: Tuple[float,float]=None, name: str=''):
        self.rack=rack
        self.name=name
        self.torque_input=torque_input
        self.Z=Z
        self.material=material
        self.speed_input=speed_input
        self.gearing_interior=gearing_interior
        self.transverse_pressure_angle=transverse_pressure_angle
        
        if coefficient_profile_shift==None:
            coefficient_profile_shift=[0.2,0.8]
        self.coefficient_profile_shift=coefficient_profile_shift
        DessiaObject.__init__(self, name=name)

class CenterDistanceOpti(DessiaObject):
    _standalone_in_db = True

    def __init__(self,center_distance:Tuple[float,float],meshes:List[MeshOpti], name:str='', 
                 constraint_root_diameter:List[bool]=None, CA_min:List[float]=None, 
                 constraint_SAP_diameter:List[bool]=None, 
                 distance_SAP_root_diameter_active_min:List[float]=None,
                 total_contact_ratio_min: List[float]=None,axial_contact_ratio: List[float]=None,
                 transverse_contact_ratio_min: List[float]=None,):

        self.meshes = meshes
        self.center_distance = center_distance
        
        if constraint_root_diameter == None:
            
            self.constraint_root_diameter = [False]*(int(len(meshes)/2))
        else:
            self.constraint_root_diameter = constraint_root_diameter 
        
        if CA_min == None: 
            self.CA_min = [0]*(int(len(meshes)/2))
        else:
            self.CA_min = CA_min
            
        if constraint_SAP_diameter == None:
            self.constraint_SAP_diameter = [False]*(int(len(meshes)/2))
        else:
            self.constraint_SAP_diameter = constraint_SAP_diameter
        
        if distance_SAP_root_diameter_active_min == None:
            self.distance_SAP_root_diameter_active_min = [0]*(int(len(meshes)/2))   
        
        else:
            self.distance_SAP_root_diameter_active_min = distance_SAP_root_diameter_active_min
            
        
        if total_contact_ratio_min == None:
            self.total_contact_ratio_min = [0]*(int(len(meshes)/2))   
        
        else:
            self.total_contact_ratio_min = total_contact_ratio_min
            
        if axial_contact_ratio == None:
            self.axial_contact_ratio = [0]*(int(len(meshes)/2))   
        
        else:
            self.axial_contact_ratio = axial_contact_ratio
        
        if transverse_contact_ratio_min == None:
            self.transverse_contact_ratio_min = [0]*(int(len(meshes)/2))   
        
        else:
            self.transverse_contact_ratio_min = transverse_contact_ratio_min
            
            
            
        
            
        DessiaObject.__init__(self, name=name)


class MeshAssemblyOptimizer(protected_module.MeshAssemblyOptimizer if _open_source==True else DessiaObject):
    _standalone_in_db = True
    """
    Gear mesh assembly optimizer supervisor

    :param connections: List of tuples defining gear mesh connection [(node1,node2), (node2,node3)...]
    :param gear_speeds: Dictionary defining minimum and maximum speed for each gear mesh {node1: [speed_min,speed_max], node2: [speed_min,speed_max]...}
    :param center_distances: List of two elements define gear mesh connection [[mesh1_centerdistance_min, mesh1_centerdistance_max], [mesh2_centerdistance_min, mesh2_centerdistance_max] ...] with mesh1 the mesh between node1 and node2 ...
    :param transverse_pressure_angle: List of two elements define the transversal pressure angle interval for each mesh [[mesh1_transversepressure_min, mesh1_transversepressure_max], [mesh2_transversepressure_min, mesh2_transversepressure_max] ...]
    :param helix_angle: Dictionary to define for one mesh the minimum and maximum helix angle {node2: [mesh1_helixangle_min, mesh1_helixangle_max]}
    :param gear_width: Dictionary to define for each gear mesh the minimum and maximum gear width {node1: [node1_gearwidth_min,node1_gearwidth_max] ,node2: [node2_gearwidth_min,node2_gearwidth_max] ...}
    :param forbidden_frequencies: List of two elements defining unacceptable frequency interval [(freq1_min,freq1_max), (freq2_min,freq2_max)...]
    :param coefficient_profile_shift: Dictionary defining the minimum and maximum coefficient profile shift for each node {node1: [node1_coeffshift_min,node1_coeffshift_max],node2: [node2_coeffshift_min,node2_coeffshift_max] ...}
    :param rack_list: Dictionary define all admissible rack {rack1:mechanical_components.meshes.Rack, rack2:mechanical_components.meshes.Rack ...}
    :param rack_choice: Dictionary assign for each gear mesh a list of acceptable rack {node1:[rack1,rack2], node2:[rack1],node3:[rack4,rack5] ...}
    :param material: Dictionary defining material for each gear mesh {node1:mechanical_components.meshes.Material, node2:mechanical_components.meshes.Material ...}
    :param torques: Dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
    :param cycles: Dictionary defining the number of cycles for one node {node3: number_cycles3}
    :param safety_factor: Safety factor used for the ISO design

    >>> list_cd=[[0.117,0.117],[0.12,0.13]]
    >>> list_gear_set=[(1,0),(0,2)]
    >>> list_speed={1:[1000*npy.pi/30.,1500*npy.pi/30.],0:[4100*npy.pi/30.,
                   4300*npy.pi/30.],2:[200*npy.pi/30.,400*npy.pi/30.]}
    >>> GA = meshes_opt.MeshAssemblyOptimizer(connections = list_gear_set,
                                gear_speed = list_speed,
                                center_distance = list_cd)
    """


    def __init__(self,center_distances: List[CenterDistanceOpti],cycles : List[float],
                 rigid_links: List[Tuple[MeshOpti,MeshOpti]] =None,
                 safety_factor: int =1, verbose : int =False,name: str=''):

        list_gear=[]
        connections=[]
        cd = []
        DessiaObject.__init__(self, name=name)
        self.constraints_root_diameter={}
        self.list_CA_min={}
        self.constraints_SAP_diameter={}
        self.distances_SAP_root_diameter_active_min={}
        self.axial_contact_ratio={}
        self.total_contact_ratio_min={}
        self.transverse_contact_ratio_min={}
        for center_distance in center_distances:
            connections_plan=[]
            for gear in center_distance.meshes:
                if gear not in list_gear:
                    list_gear.append(gear)
            
            for i in range(int(len(center_distance.meshes)/2)):
                gear_1=center_distance.meshes[2*i]
                gear_2=center_distance.meshes[2*i+1]
                
                connections_plan.append((list_gear.index(gear_1),list_gear.index(gear_2)))
                self.constraints_root_diameter[connections_plan[-1]]=center_distance.constraint_root_diameter[i]
                self.list_CA_min[connections_plan[-1]]=center_distance.CA_min[i]
                self.distances_SAP_root_diameter_active_min[connections_plan[-1]]=center_distance.distance_SAP_root_diameter_active_min[i]
                self.constraints_SAP_diameter[connections_plan[-1]]=center_distance.constraint_SAP_diameter[i]
                
                self.axial_contact_ratio[connections_plan[-1]]=center_distance.axial_contact_ratio[i]
                self.total_contact_ratio_min[connections_plan[-1]]=center_distance.total_contact_ratio_min[i]
                self.transverse_contact_ratio_min[connections_plan[-1]]=center_distance.transverse_contact_ratio_min[i]
            cd.append(center_distance.center_distance)
            connections.append(connections_plan)

       
        rack_dict={}
        rack_list=[]
        if not rigid_links:
            rigid_links=[]
        gear_speeds={}
        external_torques={}
        Z={}
        rack_choice={}
        transverse_pressure_angle={}
        coefficient_profile_shift={}
        
        number_rack=0
        self.list_gearing_interior=[]
        material={}
        list_rack_gear={}
        rigid_links_init=[]
        for list_link in rigid_links:
           rigid_links_init.append((list_gear.index(list_link[0]),list_gear.index(list_link[1])))
        for i,gear in enumerate(list_gear):
            gear_speeds[i]=gear.speed_input
          
            external_torques[i]=gear.torque_input
            
            coefficient_profile_shift[i]=gear.coefficient_profile_shift
            if gear.Z:
                Z[i]=gear.Z
            if gear.gearing_interior=='True':
                self.list_gearing_interior.append(i)
            if gear.material:
                material[i]=gear.material
            
            if gear.transverse_pressure_angle:
                transverse_pressure_angle[i]=gear.transverse_pressure_angle
            if not gear.rack in rack_list:
                rack_dict[number_rack]=gear.rack
                rack_list.append(gear.rack)
                list_rack_gear[number_rack]=[i]
                
                number_rack+=1
            else:
                num_rack=rack_list.index(gear.rack)
                list_rack_gear[num_rack].append(i)
                
        
            
            rack_choice[i]=[rack_list.index(gear.rack)]
            
        for num_rack,rack in enumerate(rack_list):
            rack.list_gear=list_rack_gear[num_rack]
        a=0
        for num_gear in rack_dict:
            if rack_dict[num_gear]!=None:
                a=1
        if a==0:
            rack_dict=None
            rack_choice=None
        if not Z:
            Z=None

        if isinstance(cycles,list):
            cycles2={}
            for i,element in enumerate(cycles):
                cycles2[i]=element
            cycles=cycles2
            
       
        self.initialisation(connections=connections,gear_speeds=gear_speeds,center_distances=cd,
                            external_torques=external_torques,cycles=cycles, rigid_links=rigid_links_init,Z=Z,
                            rack_list=rack_dict,rack_choice=rack_choice,safety_factor=safety_factor,
                            material=material,
                            coefficient_profile_shift=coefficient_profile_shift,
                            transverse_pressure_angle=transverse_pressure_angle,verbose=verbose)




    def initialisation(self, connections, gear_speeds, center_distances, external_torques, cycles,
                       rigid_links=[], Z=None, transverse_pressure_angle={},
                       gear_width=None, forbidden_frequencies=[],
                       coefficient_profile_shift=None, rack_list=None,
                       rack_choice=None, material=None,
                       safety_factor=1, verbose=False):

        self._drap_gear_graph=False
        self._drap_list_gear=False
        self._drap_connections_kinematic_dfs=False
        self._drap_sub_graph_dfs=False

        list_gear=[] # list of all gears

        for gs in connections:
            for (eng1,eng2) in gs:
                if eng1 not in list_gear:
                    list_gear.append(eng1)
                if eng2 not in list_gear:
                    list_gear.append(eng2)
        number_mesh=0
        for gs in connections:
            for (eng1,eng2) in gs:
                number_mesh+=1

        # default parameters
       
        if len(transverse_pressure_angle.keys())<number_mesh:
            for num_mesh in list_gear:
                if num_mesh not in transverse_pressure_angle.keys():
                    transverse_pressure_angle[num_mesh]=[15/180.*math.pi,30/180.*math.pi]

        # if helix_angle==None:
        #     helix_angle={list_gear[0]:[15/180.*math.pi,25/180.*math.pi]}

        if gear_width==None:
            gear_width={list_gear[0]:[15*1e-3,25*1e-3]}
        gw_min = math.inf
        gw_max = -math.inf
        for ne in gear_width.keys():
            if gear_width[ne][0]<gw_min:
                gw_min=gear_width[ne][0]
            if gear_width[ne][1]>gw_max:
                gw_max=gear_width[ne][1]
        for ne in list_gear:
            if ne not in gear_width.keys():
                gear_width[ne]=[gw_min,gw_max]

        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne]=[-0.8,0.8]

        speed_min,speed_max=[math.inf,-math.inf]
        # definition min/max absolute speed
        for num_engr,(speed_interval_min,speed_interval_max) in gear_speeds.items():
            if speed_interval_min<speed_min:
                speed_min=speed_interval_min
            if speed_interval_max>speed_max:
                speed_max=speed_interval_max

        if rack_list==None:
#            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,2.5*1e-3],
#                          'transverse_pressure_angle_rack':[20*npy.pi/180.,20*npy.pi/180.],
#                          'coeff_gear_addendum':[1,1],
#                          'coeff_gear_dedendum':[1.25,1.25],
#                          'coeff_root_radius':[0.38,0.38],
#                          'coeff_circular_tooth_thickness':[0.5,0.5]}}
            rack_list={0:RackOpti()}
        for num_rack, rack in rack_list.items():
            if  not rack.module:
                rack_list[num_rack].module = [1.2*1e-3,2*1e-3]
            if  not rack.transverse_pressure_angle_0:
                rack_list[num_rack].transverse_pressure_angle_0 =[15/180.*math.pi,30/180.*math.pi]
            if  not  rack.coeff_gear_addendum:
                rack_list[num_rack].coeff_gear_addendum = [1.00,1.00]
            if  not  rack.coeff_gear_dedendum:
                rack_list[num_rack].coeff_gear_dedendum = [1.25,1.25]
            if  not  rack.coeff_root_radius:
                rack_list[num_rack].coeff_root_radius = [0.38,0.38]
            if not  rack.coeff_circular_tooth_thickness:
                rack_list[num_rack].coeff_circular_tooth_thickness = [0.5,0.5]
            
            
            if not  rack.helix_angle:
                rack_list[num_rack].helix_angle = [0,0]
            if  not  rack.name:
                rack_list[num_rack].name = 'Optim_Module'

       
        if rack_choice==None:
            rack_choice={list_gear[0]:[list(rack_list.keys())[0]]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne]=[list(rack_list.keys())[0]]

        if material==None:
            material={list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne]=hardened_alloy_steel

#        if torques==None:
#            torques=[{list_gear[0]:100,list_gear[1]:'output'}]
#
#        if cycle==None:
#            cycle={list_gear[0]:1e6}

        if Z==None:
            Z={}

        self.Z=Z
        self.connections = connections
        self.gear_speeds = gear_speeds
        self.forbidden_frequencies = forbidden_frequencies
        self.center_distances = center_distances
        self.transverse_pressure_angle = transverse_pressure_angle
        self.coefficient_profile_shift = coefficient_profile_shift
        self.rack_list = rack_list
        self.rack_choice = rack_choice
        self.material = material
        self.external_torques =external_torques
        self.cycles = cycles
        self.safety_factor = safety_factor
        self.rigid_links = rigid_links
        self.safety_factor = safety_factor
        
        self.nb_rack = len(self.rack_list.keys())
        self.check = True

        if self.Z == {}:
            var_Z = self.AnalyseZ()
            for num, li_z in var_Z.items():
                if li_z[0] > li_z[1]:
                    self.check = False
            self.Z=var_Z



        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]

# class Instanciate(DessiaObject):
    
#     def __init__(self, gear_speeds: Dict, center_distances:List, torques: Dict):
#         self.gear_speeds = gear_speeds
#         self.center_distances = center_distances
#         self.torques = torques 
        
        
#     def instanciate(self):
        
#         rack = RackOpti(transverse_pressure_angle_0=[20/180.*npy.pi,20/180.*npy.pi], module=[2*1e-3,2*1e-3],
#              coeff_gear_addendum=[1,1],coeff_gear_dedendum=[1.25,1.25],coeff_root_radius=[0.38,0.38],
#              coeff_circular_tooth_thickness=[0.5,0.5],helix_angle=[21,60])
        
#         meshoptis = []
#         for i, speed_input in enumerate(self.gear_speeds.values()):
#             meshoptis.append(MeshOpti(rack = rack, torque_input = self.torques[i], speed_input = speed_input))
        
#         center_distances = []
#         for i, center_distance in enumerate(self.center_distances):
#             if center_distance != self.center_distances[-1]: 
#                 center_distances.append([CenterDistanceOpti(center_distance = center_distance, meshes = [meshoptis[i], meshoptis[i+1]])])
#             else:
#                 center_distances.append([CenterDistanceOpti(center_distance = center_distance, meshes= [meshoptis[0], meshoptis[-1]])])
        
#         return center_distances
            
            
            
