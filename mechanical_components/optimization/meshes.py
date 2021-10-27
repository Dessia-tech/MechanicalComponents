#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""

"""
import numpy as npy
from dessia_common import DessiaObject
from typing import  List, Tuple, Dict
from mechanical_components.meshes import MeshAssembly, hardened_alloy_steel,\
        gear_graph_simple, Material

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
    
    
     """
    Rack Optimisation definition

    :param transverse_pressure_angle_0: Tuple of 2 float which define the limit of the transverse_pressure_angle of the rack (min,max)
    :type transverse_pressure_angle_0: Tuple[float, float]
    :param module: Tuple of 2 float which define the limit of the module of the rack (min,max)
    :type module: Tuple[float, float]
    :param coeff_gear_addendum: Tuple of 2 float which define the limit of the gear addendum coefficient (min,max) (gear_addendum = coeff_gear_addendum*module)
    :type coeff_gear_addendum: Tuple[float, float]
    :param coeff_gear_dedendum: Tuple of 2 float which define the limit of the gear dedendum coefficient (min,max) (gear_dedendum = coeff_gear_dedendum*module)
    :type coeff_gear_dedendum:  Tuple[float, float]
    :param coeff_root_radius:  Tuple of 2 float which define the limit of the root radius coefficient(min,max) (root_radius = coeff_root_radius*module)
    :type coeff_root_radius:  Tuple[float, float]
    :param coeff_circular_tooth_thickness: Tuple of 2 float which define the limit of  the circular tooth thickness coefficient (min,max) (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)
    :type coeff_circular_tooth_thickness: Tuple[float, float]
    :param helix_angle: Tuple of 2 float which define the limit of  the helix_angle of the rack
    :type helix_angle: Tuple[float, float]

     """
    
     _standalone_in_db = True

  

     def __init__(self, transverse_pressure_angle_0: Tuple[float, float] = (15/180.*math.pi, 30/180.*math.pi), module: Tuple[float, float] = (1.2*1e-3, 2*1e-3),
                  coeff_gear_addendum: Tuple[float, float] = (1.00, 1.00), coeff_gear_dedendum: Tuple[float, float] = (1.25, 1.25),
                  coeff_root_radius: Tuple[float, float] = (0.38, 0.38), coeff_circular_tooth_thickness: Tuple[float, float] = (0.5, 0.5),
                  helix_angle: Tuple[float, float] = (0, 0),
                  name: str = ''):


         self.transverse_pressure_angle_0 = transverse_pressure_angle_0
         self.module = module
         if helix_angle == None:
             helix_angle = [0, 0]
         self.helix_angle = helix_angle

         self.coeff_gear_addendum = coeff_gear_addendum
         self.coeff_gear_dedendum = coeff_gear_dedendum
         self.coeff_root_radius = coeff_root_radius
         self.coeff_circular_tooth_thickness = coeff_circular_tooth_thickness
         self.list_gear = []

         self.name = name

         DessiaObject.__init__(self, name=name)

class MeshOpti(DessiaObject):
    _standalone_in_db = True
    """
    Mesh Optimisation definition

    :param torque_input: float , represent the torque imposed on the gear
    :param speed_input: Tuple of 2 float which define the limit of the speed imposed on the gear (min,max)
    :param rack: RackOpti objet which define the specification for the rack of the mesh
    :param Z: int to impose the value of the number_tooth of the gear
    :param gearing_interior: string, equal 'True' if it's an interior gear, 'False' if not
    :param coefficient_profile_shift: Tuple of 2 float which define the limit of  the coefficient_profile_shift of the gear(min,max)
    :param material: Material Objet with define the material of the mesh
    :param transverse_pressure_angle: Tuple of 2 float which define the limit of  the coefficient_profile_shift of the gear(min,max)
    """

    def __init__(self, torque_input: float, speed_input: Tuple[float, float], Z: int = 0, rack: RackOpti = None, gearing_interior: str = 'False',
                 coefficient_profile_shift: Tuple[float, float] = (0.2, 0.8), material: Material = None, transverse_pressure_angle: Tuple[float, float] = None, name: str = ''):
        self.rack = rack
        self.name = name
        self.torque_input = torque_input
        self.Z = Z
        self.material = material
        self.speed_input = speed_input
        self.gearing_interior = gearing_interior
        self.transverse_pressure_angle = transverse_pressure_angle

        if coefficient_profile_shift == None:
            coefficient_profile_shift = (0.2, 0.8)
        self.coefficient_profile_shift = coefficient_profile_shift
        DessiaObject.__init__(self, name=name)

class CenterDistanceOpti(DessiaObject):
    """

    Objet for define the specification of the connexion between minimum 2 gears for the optimisation.

    :param center_distance: Tuple of 2 float which define the limit of the center_distance (min,max)
    :param meshes: List of MeshOpti which represent the meshes which have the same center_distances.
    [MeshOpti1,MeshOpti2,MeshOpti3,MeshOpti4] -> MeshOpti1 are connected to the MeshOpti2 and MeshOpti3 are connected to the MeshOpti4 but the two connections have the same center_distance

    :param constraint_root_diameter: List of bool which specify if the diameter take for the condition of No Collisions between gears is the Root Diameter (If it's False the diameter taken is the Root_diameter active)
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param CA_min: List of float which define the CA_min  for the connexions,  the CA is the space between the root_diameter of the first_gear and the tip_diameter of the other_gear
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param constraint_SAP_diameter: List of bool which activate a condition of distance betwen the SAP_diameter(Start of Active Profile) of the first gear and the root_diameter_active of the second.
     The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param distance_SAP_root_diameter_active_min: List of float which define the distance minimum between the SAP_diameter of the first gear and the root_diameter_active of the second
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param total_contact_ratio_min: List of float which define the condition total_contact_ratio_min (axial_contact_ratio+transverse_contact_ratio) between the gears
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param axial_contact_ratio: List of float which define the  axial_contact_ratio  imposed between the gears
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param transverse_contact_ratio_min: List of float which define the condition of transverse_contact_ratio_min between the gears
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param percentage_width_difference_pinion_gear: List of float which define the percentage of difference between the width of a pinion and the width of a gear (width_pinion = width_gears*(1+percentage))
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    :param max_width_difference_pinion_gear: List of float which define the max difference allowed between the width of a pinion and the width of a gear
    (width_pinion = width_gears*(1+percentage) but if width_gears*percentage>param max_width_difference_pinion_gear so width_pinion = width_gears + max_width_difference_pinion_gear)
    The first element of this list which specify the connexion between MeshOpti1 and MeshOpti2, the second the connexion between MeshOpti3 and MeshOpti4 and etc...

    """

    _standalone_in_db = True

    def __init__(self, center_distance: Tuple[float, float], meshes: List[MeshOpti], name: str = '',
                 constraint_root_diameter: List[bool] = None, CA_min: List[float] = None,
                 constraint_SAP_diameter: List[bool] = None,
                 distance_SAP_root_diameter_active_min: List[float] = None,
                 total_contact_ratio_min: List[float] = None, axial_contact_ratio: List[float] = None,
                 transverse_contact_ratio_min: List[float] = None,
                 percentage_width_difference_pinion_gear: List[float] = None,
                 max_width_difference_pinion_gear: List[float] = None):

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

        if percentage_width_difference_pinion_gear == None:
            self.percentage_width_difference_pinion_gear = [0]*(int(len(meshes)/2))

        else:
            self.percentage_width_difference_pinion_gear = percentage_width_difference_pinion_gear

        if max_width_difference_pinion_gear == None:
            self.max_width_difference_pinion_gear = [0]*(int(len(meshes)/2))

        else:
            self.max_width_difference_pinion_gear = max_width_difference_pinion_gear





        DessiaObject.__init__(self, name=name)


class MeshAssemblyOptimizer(protected_module.MeshAssemblyOptimizer if _open_source == True else DessiaObject):
    _standalone_in_db = True
    """
    Gear mesh assembly optimizer supervisor

    :param center_distances: List of CenterDistanceOpti which define all the connexion of the MeshAssembly
    :param cycles: List of float which define the cycle of hte gears
    :param rigid_links: List of Tuple of 2 MeshOpti which define if there have a rigid connexion between mesh


    """



    def __init__(self, center_distances: List[CenterDistanceOpti], cycles: List[float],
                 rigid_links: List[Tuple[MeshOpti, MeshOpti]] = None,
                 safety_factor: int = 1, verbose: int = False, name: str = ''):

        list_gear = []
        connections = []
        cd = []
        DessiaObject.__init__(self, name=name)
        self.constraints_root_diameter = {}
        self.list_CA_min = {}
        self.constraints_SAP_diameter = {}
        self.distances_SAP_root_diameter_active_min = {}
        self.axial_contact_ratio = {}
        self.total_contact_ratio_min = {}
        self.transverse_contact_ratio_min = {}
        self.percentage_width_difference_pinion_gear = {}
        self.max_width_difference_pinion_gear = {}
        for center_distance in center_distances:
            connections_plan = []
            for gear in center_distance.meshes:
                if gear not in list_gear:
                    list_gear.append(gear)

            for i in range(int(len(center_distance.meshes)/2)):
                gear_1 = center_distance.meshes[2*i]
                gear_2 = center_distance.meshes[2*i+1]

                connections_plan.append((list_gear.index(gear_1), list_gear.index(gear_2)))
                self.constraints_root_diameter[connections_plan[-1]] = center_distance.constraint_root_diameter[i]
                self.list_CA_min[connections_plan[-1]] = center_distance.CA_min[i]
                self.distances_SAP_root_diameter_active_min[connections_plan[-1]] = center_distance.distance_SAP_root_diameter_active_min[i]
                self.constraints_SAP_diameter[connections_plan[-1]] = center_distance.constraint_SAP_diameter[i]

                self.axial_contact_ratio[connections_plan[-1]] = center_distance.axial_contact_ratio[i]
                self.total_contact_ratio_min[connections_plan[-1]] = center_distance.total_contact_ratio_min[i]
                self.transverse_contact_ratio_min[connections_plan[-1]] = center_distance.transverse_contact_ratio_min[i]
                self.percentage_width_difference_pinion_gear[connections_plan[-1]] = center_distance.percentage_width_difference_pinion_gear[i]
                self.max_width_difference_pinion_gear[connections_plan[-1]] = center_distance.max_width_difference_pinion_gear[i]
            cd.append(center_distance.center_distance)
            connections.append(connections_plan)


        rack_dict = {}
        rack_list = []
        if not rigid_links:
            rigid_links = []
        gear_speeds = {}
        external_torques = {}
        Z = {}
        rack_choice = {}
        transverse_pressure_angle = {}
        coefficient_profile_shift = {}

        number_rack = 0
        self.list_gearing_interior = []
        
        self.list_mesh = list_gear
        material = {}
        list_rack_gear = {}
        rigid_links_init = []
        for list_link in rigid_links:
           rigid_links_init.append((list_gear.index(list_link[0]), list_gear.index(list_link[1])))
        for i, gear in enumerate(list_gear):
            gear_speeds[i] = gear.speed_input

            external_torques[i] = gear.torque_input

            coefficient_profile_shift[i] = gear.coefficient_profile_shift
            if gear.Z:
                Z[i] = gear.Z
            if gear.gearing_interior == 'True':
                self.list_gearing_interior.append(i)
            if gear.material:
                material[i] = gear.material

            if gear.transverse_pressure_angle:
                transverse_pressure_angle[i] = gear.transverse_pressure_angle
            if not gear.rack in rack_list:
                rack_dict[number_rack] = gear.rack
                rack_list.append(gear.rack)
                list_rack_gear[number_rack] = [i]

                number_rack += 1
            else:
                num_rack = rack_list.index(gear.rack)
                list_rack_gear[num_rack].append(i)



            rack_choice[i] = [rack_list.index(gear.rack)]

        for num_rack, rack in enumerate(rack_list):
            rack.list_gear = list_rack_gear[num_rack]
        a = 0
        for num_gear in rack_dict:
            if rack_dict[num_gear] != None:
                a = 1
        if a == 0:
            rack_dict = None
            rack_choice = None
        if not Z:
            Z = None

        if isinstance(cycles, list):
            cycles2 = {}
            for i, element in enumerate(cycles):
                cycles2[i] = element
            cycles = cycles2


        self.initialisation(connections=connections, gear_speeds=gear_speeds, center_distances=cd,
                            external_torques=external_torques, cycles=cycles, rigid_links=rigid_links_init, Z=Z,
                            rack_list=rack_dict, rack_choice=rack_choice, safety_factor=safety_factor,
                            material=material,
                            coefficient_profile_shift=coefficient_profile_shift,
                            transverse_pressure_angle=transverse_pressure_angle, verbose=verbose)




    def initialisation(self, connections, gear_speeds, center_distances, external_torques, cycles,
                       rigid_links=[], Z=None, transverse_pressure_angle={},
                       gear_width=None, forbidden_frequencies=[],
                       coefficient_profile_shift=None, rack_list=None,
                       rack_choice=None, material=None,
                       safety_factor=1, verbose=False):

        self._drap_gear_graph = False
        self._drap_list_gear = False
        self._drap_connections_kinematic_dfs = False
        self._drap_sub_graph_dfs = False

        list_gear = [] # list of all gears

        for gs in connections:
            for (eng1, eng2) in gs:
                if eng1 not in list_gear:
                    list_gear.append(eng1)
                if eng2 not in list_gear:
                    list_gear.append(eng2)
        number_mesh = 0
        for gs in connections:
            for (eng1, eng2) in gs:
                number_mesh += 1

        # default parameters

        if len(transverse_pressure_angle.keys()) < number_mesh:
            for num_mesh in list_gear:
                if num_mesh not in transverse_pressure_angle.keys():
                    transverse_pressure_angle[num_mesh] = [15/180.*math.pi, 30/180.*math.pi]

        # if helix_angle==None:
        #     helix_angle={list_gear[0]:[15/180.*math.pi,25/180.*math.pi]}

        if gear_width == None:
            gear_width = {list_gear[0]:[15*1e-3, 25*1e-3]}
        gw_min = math.inf
        gw_max = -math.inf
        for ne in gear_width.keys():
            if gear_width[ne][0] < gw_min:
                gw_min = gear_width[ne][0]
            if gear_width[ne][1] > gw_max:
                gw_max = gear_width[ne][1]
        for ne in list_gear:
            if ne not in gear_width.keys():
                gear_width[ne] = [gw_min, gw_max]

        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne] = [-0.8, 0.8]

        speed_min, speed_max = [math.inf, -math.inf]
        # definition min/max absolute speed
        for num_engr, (speed_interval_min, speed_interval_max) in gear_speeds.items():
            if speed_interval_min < speed_min:
                speed_min = speed_interval_min
            if speed_interval_max > speed_max:
                speed_max = speed_interval_max

        if rack_list == None:
#            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,2.5*1e-3],
#                          'transverse_pressure_angle_rack':[20*npy.pi/180.,20*npy.pi/180.],
#                          'coeff_gear_addendum':[1,1],
#                          'coeff_gear_dedendum':[1.25,1.25],
#                          'coeff_root_radius':[0.38,0.38],
#                          'coeff_circular_tooth_thickness':[0.5,0.5]}}
            rack_list = {0:RackOpti()}
        for num_rack, rack in rack_list.items():
            if  not rack.module:
                rack_list[num_rack].module = (1.2*1e-3, 2*1e-3)
            if  not rack.transverse_pressure_angle_0:
                rack_list[num_rack].transverse_pressure_angle_0 = (15/180.*math.pi, 30/180.*math.pi)
            if  not  rack.coeff_gear_addendum:
                rack_list[num_rack].coeff_gear_addendum = (1.00, 1.00)
            if  not  rack.coeff_gear_dedendum:
                rack_list[num_rack].coeff_gear_dedendum = (1.25, 1.25)
            if  not  rack.coeff_root_radius:
                rack_list[num_rack].coeff_root_radius = (0.38, 0.38)
            if not  rack.coeff_circular_tooth_thickness:
                rack_list[num_rack].coeff_circular_tooth_thickness = (0.5, 0.5)


            if not  rack.helix_angle:
                rack_list[num_rack].helix_angle = (0, 0)
            if  not  rack.name:
                rack_list[num_rack].name = 'Optim_Module'


        if rack_choice == None:
            rack_choice = {list_gear[0]:[list(rack_list.keys())[0]]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne] = [list(rack_list.keys())[0]]

        if material == None:
            material = {list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne] = hardened_alloy_steel

#        if torques==None:
#            torques=[{list_gear[0]:100,list_gear[1]:'output'}]
#
#        if cycle==None:
#            cycle={list_gear[0]:1e6}

        if Z == None:
            Z = {}

        self.Z = Z
        self.connections = connections
        self.gear_speeds = gear_speeds
        self.forbidden_frequencies = forbidden_frequencies
        self.center_distances = center_distances
        self.transverse_pressure_angle = transverse_pressure_angle
        self.coefficient_profile_shift = coefficient_profile_shift
        self.rack_list = rack_list
        self.rack_choice = rack_choice
        self.material = material
        self.external_torques = external_torques
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
            self.Z = var_Z



        self.solutions = []
        self.solutions_search = []
        self.analyse = []

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



