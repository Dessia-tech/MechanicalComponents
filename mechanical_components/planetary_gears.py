#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:56:02 2020

@author: launay
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as npy
import dectree as dt
from scipy.linalg import solve
import time
import math as m
import copy
import volmdlr as vm
import volmdlr.primitives3d as p3d
import volmdlr.primitives2d as p2d
import mechanical_components.meshes as meshes
from dessia_common import DessiaObject
from typing import Tuple, List, TypeVar
import numpy as np
import genmechanics as genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import plot_data as vmp
import scipy.optimize as op
from mechanical_components.meshes import hardened_alloy_steel
import mechanical_components.optimization.meshes_protected as mg
import mechanical_components.optimization.meshes as meshes_opt
import volmdlr.primitives3d as primitives3d
# from dessia_common.list_eq import list_eq
import genmechanics.geometry as gm_geo
import volmdlr.core_compiled as vm_compiled
import plot_data as pld
def list_eq(list_1, list_2):

    if len(list_1) == len(list_2):
        index_list = []
        for element_1 in list_1:
            for i, element_2 in enumerate(list_2):
                if element_1 == element_2:
                    index = i

                    if index not in index_list:
                        index_list.append(index)
                        break

        if len(index_list) == len(list_2):
            return True
    return False



class Gears(DessiaObject):


    def __init__(self, Z: int, name: str = ''):
        self.Z = Z



        DessiaObject.__init__(self, name=name)

    # def voldmr_volume_model(self):
    #     model = self.volume_model()
    #     return model

    # def volume_model(self, module, xy_position,z_position,lenght):
    #        self.module = module
    #        self.d = module*self.Z

    #        x = vm.Vector3D((1, 0, 0))
    #        y = vm.Vector3D((0, 1, 0))
    #        z = vm.Vector3D((0, 0, 1))
    #        rack = meshes.Rack(20*3.14/180, module)
    #        radius = m.cos(rack.transverse_pressure_angle)*module*self.Z
    #        meshes_1 = meshes.Mesh(self.Z, radius, 0.01, rack)
    #        Gears3D = {0:meshes_1.Contour(3)}

    #        export = []
    #        center_2 = (xy_position[0], xy_position[1])
    #        center = vm.Point2D(center_2)
    #        model_trans = Gears3D[0][0].Translation(center)

    #        model_trans_rot = model_trans.Rotation(center, 0.1)
    #        Gears3D_Rotate = [model_trans_rot]

    #        export = []

    #        for (i, center, k) in zip([Gears3D[0]], [center_2], [-1]):

    #                       model_export = []

    #                       for m2 in i:

    #                           center = vm.Point2D(center)
    #                           model_trans = m2.Translation(center)
    #                           model_trans_rot = model_trans.Rotation(center, k)
    #                           model_export.append(model_trans_rot)

    #                       export.append(model_export)

    #        Gears3D_Rotate = export
    #        vect_x = z_position*z
    #        extrusion_vector1 = lenght*z
    #        C1 = Gears3D_Rotate[0]
    #        L=[]
           # for gear in Gears3D_Rotate[0]:
           #     L.append(gear.plot_data('contour', stroke_width=8))

           # vmp.plot([L[0]])
           # t1 = p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1[0], [], vm.Vector3D(extrusion_vector1))

    #      if module==0:
    #          return None

    #      pos=vm.Point3D(self.position)
    #      axis=vm.vector3D((0,0,1))
    #      radius=(self.module*self.Z)/2
    #      cylinder=vm.cylinder(pos,axis,radius,self.length)
    #      return cylinder

class Planetary(Gears):
    _eq_is_data_eq = False

    '''

    Define a planetary

    :param Z: The number of tooth
    :type Z: int
    :param planetary_type: The type of the planetary:

        - ' Ring' for ring

        - 'Sun' for sun

    :type planetary_type: str
    :param name: Name
    :type name: str, optional


    '''

    def __init__(self, Z: int, planetary_type: str, name: str = '', speed_input: List[float] = None,
                 torque_input: List[float] = None, position: List[Tuple[float, float, float]] = None, module: float = 0):




        self.planetary_type = planetary_type
        self.p = 0
        self.speed = 0
        self.module = module
        self.d = 0
        self.speed_input = speed_input
        self.position = position
        Gears.__init__(self, Z=Z, name=name)
        self.length = 0.2
        self.Z = Z
        self.name = name
        self.torque_input = torque_input
        # self.torque_signe=1
        self.power = 0
        if planetary_type == 'Sun':

            self.p = 1

        else:

            self.p = -1

    def voldmr_volume_model(self):
        model = self.volume_model()
        return model

    def volume_model(self):
        position = self.position
        module = self.module

        if self.planetary_type == 'Sun':
            pos = vm.Point3D(position)
            axis = vm.Vector3D((1, 0, 0))
            radius = (module*self.Z)/2

            cylinder = p3d.Cylinder(pos, axis, radius, position[0]+self.length)
            return cylinder

        radius = (module*self.Z)/2
        p1 = vm.Point2D((radius, position[0]+self.length/2))
        p2 = vm.Point2D((radius+radius*0.1, position[0]+self.length/2))
        p3 = vm.Point2D((radius, position[0]-self.length/2))
        p4 = vm.Point2D((radius+radius*0.1, position[0]-self.length/2))
        points1 = [p1, p2, p3, p4]
        c1 = vm.Polygon2D(points1)
        vector_1 = vm.Point3D((0, 0, 0))

        profile1 = p3d.RevolvedProfile(vector_1, vm.Z3D, vm.X3D, c1, vm.O3D, vm.X3D)

        return profile1


class Planet(Gears):
    _eq_is_data_eq = False

    '''
    Define a planet

    :param Z: The number of tooth
    :type Z: int
    :param name: Name
    :type name: str, optional


    '''

    def __init__(self, Z: int, name: str = '', positions: List[Tuple[float, float, float]] = None, module: float = 0):



        self.length = 0.2
        self.speed = 0
        self.torque = 0
        self.module = module
        self.speed_input = [0, 0]
        self.Z = Z
        self.name = name
        self.positions = positions
        Gears.__init__(self, Z, name)

    def voldmr_volume_model(self):
        model = self.volume_model()
        return model

    def volume_model(self):

        positions = self.positions
        module = self.module
        model = []
        for position in positions:
            pos = vm.Point3D(position)

            axis = vm.Vector3D((1, 0, 0))
            radius = (module*self.Z)/2

            cylinder = p3d.Cylinder(pos, axis, radius, self.length)

            model.append(cylinder)

        return model




class PlanetCarrier(DessiaObject):
    _eq_is_data_eq = False

    '''
    Define a planet carrier

    :param name: Name
    :type name: str, optional



    '''

    def __init__(self, name: str = '', speed_input: List[float] = [0, 0], torque_input: List[float] = [0, 0]):


        self.speed = 0
        self.speed_input = speed_input

        self.torque_input = torque_input
        self.torque_signe = 1
        self.power = 0
        DessiaObject.__init__(self, name=name)



class Meshing(DessiaObject):
    _eq_is_data_eq = False



    def __init__(self, nodes: List[Gears], name: str = ''):

        self.nodes = nodes
        DessiaObject.__init__(self, name=name)


class MeshingPlanetary(Meshing):
    _eq_is_data_eq = False
    def __init__(self, nodes: List[Gears], name: str = ''):

        self.type = 'GI'
        Meshing.__init__(self, nodes, name)




    def speed_system_equations(self):

        for node in self.nodes:
            if isinstance(node, Planet):
                Z_planet = node.Z
            else:
                Z_planetary = node.p*node.Z

        matrix = npy.array([Z_planetary, Z_planet, -Z_planetary])
        rhs = npy.array([0])
        return matrix, rhs

    def torque_system_equations(self):
        for node in self.nodes:
            if isinstance(node, Planet):
                Z_planet = node.Z
            else:
                Z_planetary = node.p*node.Z
        matrix = npy.array([[-1/Z_planetary, 1/Z_planet],
                            [1/Z_planetary, 1/Z_planet]])
        rhs = npy.array([0, 0])
        return matrix, rhs


class MeshingPlanet(Meshing):
        _eq_is_data_eq = False
        def __init__(self, nodes: List[Gears], name: str = ''):

            self.type = 'GI'
            Meshing.__init__(self, nodes, name)
            self.Z_planets = []



        def speed_system_equations(self):

            matrix = npy.array([self.nodes[0].Z, self.nodes[1].Z])
            rhs = npy.array([0])
            return matrix, rhs

        def torque_system_equations(self):
            matrix = npy.array([-1/self.nodes[0].Z, 1/self.nodes[1].Z])
            rhs = npy.array([0])
            return matrix, rhs


class Pivot(DessiaObject):

    A = TypeVar('A', Planetary, PlanetCarrier)
    def __init__(self, nodes: A, name: str = ''):
        self.nodes = nodes

        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Fixed(DessiaObject):

    def __init__(self, nodes, name: str = ''):
        self.nodes = nodes

        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Double(DessiaObject):
    _eq_is_data_eq = False

    def __init__(self, nodes: List[Planet], name: str = ''):

        self.nodes = nodes

        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

    def voldmr_volume_model(self, axis=(1, 0, 0), center=(0, 0, 0)):
        model = self.volume_model(axis=axis, center=center)
        return model

    def volume_model(self, axis=(1, 0, 0), center=(0, 0, 0)):
         position_planet_1 = self.nodes[0].positions
         position_planet_2 = self.nodes[1].positions

         model = []
         axis = vm.Vector3D(axis[0], axis[1], axis[2])
         for i in range(len(position_planet_1)):

             if position_planet_2[i][0] > position_planet_1[i][0]:

                 if position_planet_2[i][0] > 0:
                     position = ((position_planet_2[i][0]-position_planet_1[i][0])/2, position_planet_1[i][1], position_planet_1[i][2])
                 else:
                     position = ((position_planet_1[i][0]-position_planet_2[i][0])/2, position_planet_1[i][1], position_planet_1[i][2])
             else:
                 if position_planet_1[i][0] > 0:
                     position = ((position_planet_1[i][0]-position_planet_2[i][0])/2, position_planet_1[i][1], position_planet_1[i][2])
                 else:
                     position = ((position_planet_2[i][0]-position_planet_1[i][0])/2, position_planet_1[i][1], position_planet_1[i][2])

             if not center == (0, 0, 0):

                 position = (position[0]+center[0], position[1]+center[1],
                             position[2]+center[2])
                 pos = vm.Point3D(position[0], position[1], position[2])
             if not axis == (1, 0, 0):
                 vector = npy.cross((axis[0], axis[1], axis[2]), (1, 0, 0))
                 axis_rotation = vm.Vector3D(vector[0], vector[1], vector[2])
                 axis_vector = vm.Vector3D(axis[0], axis[1], axis[2])
                 axis_origin = vm.Vector3D(1, 0, 0)

                 axis_vector_norme = copy.copy(axis_vector)
                 axis_vector_norme.normalize()

                 axis_origin_norme = copy.copy(axis_origin)
                 axis_origin_norme.normalize()

                 if axis_vector_norme.dot(axis_origin_norme) == 0:
                        angle = m.pi/2
                 else:
                        angle = m.acos(axis_vector.dot(axis_origin)/axis_vector_norme.dot(axis_origin_norme))
                 position2 = vm.Vector3D(position[0], position[1], position[2])

                 position2.rotation(center=vm.Vector3D(center[0], center[1], center[2]), axis=axis_rotation, angle=angle, copy=False)

                 pos = vm.Point3D(position2[0], position2[1], position2[2])


             else:
                 pos = vm.Point3D(position[0], position[1], position[2])


             cylinder = p3d.Cylinder(pos, axis, (self.nodes[0].Z*self.nodes[0].module)/10, abs(position_planet_1[i][0]-position_planet_2[i][0]))
             model.append(cylinder)
         return model

class ImposeSpeed(DessiaObject):
    A = TypeVar('A', Gears, PlanetCarrier)
    def __init__(self, node: A, input_speed: float, name: str = ''):
        self.input_speed = input_speed
        self.node = node

        DessiaObject.__init__(self, name=name)

    def speed_system_equations(self):
        matrix = npy.array([1])
        rhs = npy.array([self.input_speed])
        return matrix, rhs




class Connection(DessiaObject):
    _eq_is_data_eq = False



    '''
    Define a connection


    :param nodes: The 2 elements connected
    :type nodes: List[Planet,Planetary]
    :param connection_type: The type of the connection :

        -'D' is for Double

        -'GI' is when the first element of nodes meshing to the second inward of the planetary gear

        -'GE' is when the first element of nodes meshing to the second outward of the planetary gear


    :type connection_type: str

    :param name: Name
    :type name: str, optional



    '''

    def __init__(self, nodes: List[Gears], connection_type: str, name: str = ''):

        self.nodes = nodes
        self.connection_type = connection_type
        DessiaObject.__init__(self, name=name)

    # def __eq__(self,other):
    #     if isinstance(other,Connection):
    #         if self.connection_type==other.connection_type:
    #             nodeZ_1=[]
    #             nodeZ_2=[]
    #             node_name1=[]
    #             node_name2=[]
    #             for node in self.nodes:
    #                 nodeZ_1.append(node.Z)

    #             for node in other.nodes:
    #                 nodeZ_2.append(node.Z)


    #             if nodeZ_1==nodeZ_2:
    #                 return True
    #     return False

    # def __hash__(self):
    #     nodeZ_1=[]
    #     for node in self.nodes:
    #         nodeZ_1.append(node.Z)
    #     return node.Z





class PlanetaryGear(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['mech', 'mech_dict', 'max_length_meshing_chain', 'center', 'axis']
    _eq_is_data_eq = False



    '''
    Define a Planetary Gears

    :param planetaries: The planetaries of the planetary gear
    :type planetaries: List[Planetary]
    :param planets: The planets of the planetary gear
    :type planets: List[Planet]
    :param planet_carrier: The planet_carrer of the planetary gear
    :type planet_carrier: PlanetCarrier
    :param connections: List of the connection bettween element ( meshing and Double)
    :type connections: List[Connection]
    :param name: name
    :type name: str,optional
    '''
    def __init__(self, planetaries: List[Planetary], planets: List[Planet],
                 planet_carrier: PlanetCarrier, connections: List[Connection], number_branch_planet: int = 3, name: str = '',
                 number_group_solution_planet_structure: int = 0, number_group_solution_architecture: int = 0):

        self.number_group_solution_planet_structure = number_group_solution_planet_structure
        self.number_group_solution_architecture = number_group_solution_architecture
        self.length_double = 0.05
        self.length = 0
        self.number_branch_planet = number_branch_planet
        self.d_min = 0

        self.planetaries = planetaries
        self.planets = planets
        self.planet_carrier = planet_carrier
        self.elements = self.planetaries + self.planets + [self.planet_carrier]
        self.elements_name = []
        self.center = (0, 0, 0)
        self.axis = (1, 0, 0)

        for element in self.elements:
            self.elements_name.append(element.name)
        self.max_length_meshing_chain = []
        self.mech = 0
        self.mech_dict = 0
        self.sum_Z_planetary = 0

        self.max_Z_planetary = 0
        self.min_Z_planetary = 100000

        self.z_min = min([element.Z for element in self.planets+self.planetaries])
        for planetary in self.planetaries:
            self.sum_Z_planetary += planetary.Z

            if self.max_Z_planetary < planetary.Z:
                self.max_Z_planetary = planetary.Z

            if self.min_Z_planetary > planetary.Z:
                self.min_Z_planetary = planetary.Z







        self.connections = connections
        self.meshings = []
        self.doubles = []
        DessiaObject.__init__(self, name=name)
        self.position = False




        for i, connection in  enumerate(connections):


          ## Check to be sure that all the object in connection are in planetaries,
          ## planets, or planet_carrier ##

          if not connection.nodes[1] in self.elements:

                 if isinstance(connection.nodes[1], Planetary):

                     self.elements[self.elements_name.index(connection.nodes[1].name)].planetary_type = connection.nodes[1].planetary_type
                     self.elements[self.elements_name.index(connection.nodes[1].name)].p = connection.nodes[1].p

                 connection.nodes[1] = self.elements[self.elements_name.index(connection.nodes[1].name)]



          if not connection.nodes[0] in self.elements:
                 if isinstance(connection.nodes[0], Planetary):
                     self.elements[self.elements_name.index(connection.nodes[0].name)].planetary_type = connection.nodes[0].planetary_type
                     self.elements[self.elements_name.index(connection.nodes[0].name)].p = connection.nodes[0].p
                 connection.nodes[0] = self.elements[self.elements_name.index(connection.nodes[0].name)]




          if connection.connection_type != 'D':

              if isinstance(connection.nodes[0], Planet) and isinstance(connection.nodes[1], Planet):
                self.meshings.append(MeshingPlanet([connection.nodes[0], connection.nodes[1]], 'meshing'+str(i)))


              else:
                self.meshings.append(MeshingPlanetary([connection.nodes[0], connection.nodes[1]], 'meshing'+str(i),))


              self.meshings[-1].type = connection.connection_type


          else:
              self.doubles.append(Double([connection.nodes[0], connection.nodes[1]], 'Double'+str(i)))

        self.relations = self.meshings + self.doubles

        if self.planet_carrier.speed_input[0] == self.planet_carrier.speed_input[1] and self.planet_carrier.speed_input[0] == 0:
            self.speed_max_planet = 0
        else:
            self.speed_max_planet = self.speed_max_planets()

        self.number_ring = 0
        self.number_sun = 0

        for planetary in self.planetaries:
            if planetary.planetary_type == 'Ring':
                self.number_ring += 1
            else:
                self.number_sun += 1

    def __str__(self):

        Z_planets = {}

        for planet in self.planets:
            Z_planets[planet.name] = planet.Z

        Z_planetaries = {}
        number_ring = 0
        number_sun = 0

        for planetary in self.planetaries:
            Z_planetaries[planetary.name] = planetary.Z

            if planetary.planetary_type == 'Sun':
                number_sun += 1

            else:
                number_ring += 1
        connections_name = []
        for i in range(len(self.connections)):
            connections_name.append([self.connections[i].nodes[0].name, self.connections[i].nodes[1].name,
                                     self.connections[i].connection_type])

        return 'Name:' + self.name + '\n\n' + \
               'Planetary Number:' + str(len(self.planetaries)) + '\n' + \
               'Ring Number:'+ str(number_ring) + '\n' + \
               'Sun_Number:' + str(number_sun) + '\n' + \
               'Z_planetaries:' + str(Z_planetaries) + '\n\n' + \
               'Planets_Number:' + str(len(self.planets)) + '\n' + \
               'Planets_Double_Number:' + str(len(self.doubles)) + '\n' + \
               'Z_Planets:' + str(Z_planets) + '\n\n' + \
                str(connections_name) + '\n\n\n'

    def __eq__(self, other):

        if isinstance(other, PlanetaryGear):
            connection_Z = self.connections_Z()
            if list_eq(other.connections_Z(), connection_Z):
                return True
        return False

    def __hash__(self):
        Z = 0
        for planetary in self.planetaries:
            Z += planetary.Z
        for planet in self.planets:
            Z += planet.Z
        return Z


    def connections_Z(self):
        connections_Z = []
        for connection in self.connections:
            connections_Z.append([connection.nodes[0].Z, connection.nodes[1].Z, connection.connection_type])
        return connections_Z








    def matrix_position(self, element):
        '''Give the position of the element in the speed solve matrix and in the speed result liste

        :param element: the element whose position we want to know
        :type element: Planet,Planetary or PlanetCarrier

        :return: The position
        :rtype: int

        '''

        return self.elements.index(element)

    def graph(self):

        graph_planetary_gear = nx.Graph()

        for relation in self.relations:

            graph_planetary_gear.add_edge(str(relation.nodes[0]), str(relation))
            graph_planetary_gear.add_edge(str(relation), str(relation.nodes[1]))

            nx.set_node_attributes(graph_planetary_gear, relation.nodes[0], str(relation.nodes[0]))
            nx.set_node_attributes(graph_planetary_gear, relation.nodes[1], str(relation.nodes[1]))
            nx.set_node_attributes(graph_planetary_gear, relation, str(relation))



        for k, planet in enumerate(self.planets):

            graph_planetary_gear.add_edge(str(self.planet_carrier), 'Pv'+str(k))
            graph_planetary_gear.add_edge('Pv'+str(k), str(planet))
            nx.set_node_attributes(graph_planetary_gear, 'Pv'+str(k), 'Pv'+str(k))


        # plt.figure()
        # nx.draw_kamada_kawai(graph_planetary_gear, with_labels=True)
        return graph_planetary_gear

    # def plot_graph(self):



    #     graph_planetary_gears = self.graph()
    #     plt.figure()
    #     nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)



    def plot_kinematic_graph_gear(self, coordinate, lenght, diameter,
                                  diameter_pivot, lenght_pivot, color, plot_data):

        list_color = [pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN, pld.colors.RED, 'blue', 'green',
                      pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN]

        x = [coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
             coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2]

        y = [coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2,
             coordinate[1]+diameter/2, coordinate[1]+diameter/2, coordinate[1]+diameter/2]

        # plt.plot(x, y, list_color[color])

        for i in range(len(x)-1):

                point1 = vm.Point2D(x[i], y[i])
                point2 = vm.Point2D(x[i+1], y[i+1])
                line = vm.edges.LineSegment2D(point1, point2)
                edge_style = pld.EdgeStyle(line_width=2, color_stroke=list_color[color])
                plot_data.append(line.plot_data(edge_style=edge_style))

        x = [coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
             coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2]

        y = [coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2,
             coordinate[1]-diameter/2, coordinate[1]-diameter/2, coordinate[1]-diameter/2]

        # plt.plot(x, y, list_color[color])

        for i in range(len(x)-1):
                point1 = vm.Point2D(x[i], y[i])
                point2 = vm.Point2D(x[i+1], y[i+1])
                line = vm.edges.LineSegment2D(point1, point2)
                edge_style = pld.EdgeStyle(line_width=2, color_stroke=list_color[color])
                plot_data.append(line.plot_data(edge_style=edge_style))

    def plot_kinematic_graph_double(self, coordinate, diameter, lenght, color, plot_data):

        list_color = [pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN, pld.colors.RED, 'blue', 'green',
                      pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN]

        x = [coordinate[0], coordinate[0]+lenght]
        y = [coordinate[1]+diameter/2, coordinate[1]+diameter/2]

        # plt.plot(x, y, list_color[color])

        for i in range(len(x)-1):
                point1 = vm.Point2D(x[i], y[i])
                point2 = vm.Point2D(x[i+1], y[i+1])
                line = vm.edges.LineSegment2D(point1, point2)
                edge_style = pld.EdgeStyle(line_width=2, color_stroke=list_color[color])
                plot_data.append(line.plot_data(edge_style=edge_style))

        x = [coordinate[0], coordinate[0]+lenght]
        y = [coordinate[1]-diameter/2, coordinate[1]-diameter/2]

        for i in range(len(x)-1):
                point1 = vm.Point2D(x[i], y[i])
                point2 = vm.Point2D(x[i+1], y[i+1])
                line = vm.edges.LineSegment2D(point1, point2)
                edge_style = pld.EdgeStyle(line_width=2, color_stroke=list_color[color])
                plot_data.append(line.plot_data(edge_style=edge_style))

        # plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y, plot_data):

        coordinate_y_min = 0
        coordinate_y_max = 0
        coordinate_x_max = 0

        for coordinate in coordinates:

            if coordinate[0] > coordinate_x_max:
                coordinate_x_max = coordinate[0]

            if coordinate[1] < coordinate_y_min:
                coordinate_y_min = coordinate[1]

            if coordinate[1] > coordinate_y_max:
                coordinate_y_max = coordinate[1]

        coordinate_planet_carrier = [coordinate_x_max+planet_carrier_x, coordinate_y_min-planet_carrier_y]

        for coordinate in coordinates:
            x = [coordinate[0]-planet_carrier_x, coordinate_planet_carrier[0]]
            y = [coordinate[1], coordinate[1]]

            for i in range(len(x)-1):
                point1 = vm.Point2D(x[i], y[i])
                point2 = vm.Point2D(x[i+1], y[i+1])
                line = vm.edges.LineSegment2D(point1, point2)
                edge_style = pld.EdgeStyle(line_width=2, color_stroke=pld.colors.BLACK)
                plot_data.append(line.plot_data(edge_style=edge_style))

            # plt.plot(x, y, 'r')

        x = [coordinate_planet_carrier[0]+planet_carrier_x, coordinate_planet_carrier[0], coordinate_planet_carrier[0]]
        y = [coordinate_planet_carrier[1], coordinate_planet_carrier[1], coordinate_y_max]

        for i in range(len(x)-1):
            point1 = vm.Point2D(x[i], y[i])
            point2 = vm.Point2D(x[i+1], y[i+1])
            line = vm.edges.LineSegment2D(point1, point2)
            edge_style = pld.EdgeStyle(line_width=2, color_stroke=pld.colors.BLACK)
            plot_data.append(line.plot_data(edge_style=edge_style))
        # plt.plot(x, y, 'r')

        return coordinate_planet_carrier

    def plot_kinematic_graph_ring(self, coordinate, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color, plot_data):

        list_color = [pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN, pld.colors.RED, 'blue', 'green',
                      pld.colors.RED, pld.colors.BLUE, pld.colors.GREEN]

        x = [coordinate[0]-lenght_gear/2, coordinate[0]+lenght_gear/2, coordinate[0], coordinate[0],
             coordinate_planet_carrier[0]+lenght_ring, coordinate_planet_carrier[0]+lenght_ring]
        y = [coordinate[1], coordinate[1], coordinate[1], diameter_ring/2, diameter_ring/2, coordinate_planet_carrier[1]]

        for i in range(len(x)-1):
            point1 = vm.Point2D(x[i], y[i])
            point2 = vm.Point2D(x[i+1], y[i+1])
            line = vm.edges.LineSegment2D(point1, point2)
            edge_style = pld.EdgeStyle(line_width=2, color_stroke=list_color[color])
            plot_data.append(line.plot_data(edge_style=edge_style))

        # plt.plot(x, y, list_color[color])
        coordinate[1] -= (abs(coordinate[1]-coordinate_planet_carrier[1]))*2


    def plot_kinematic_graph(self, lenght_gear=0.1, diameter_gear=1, lenght_double=2, diameter_pivot=0.2, lenght_pivot=0.5,
                             planet_carrier_x=2, planet_carrier_y=2, diameter_ring_ini=10):
        '''
        Plot the kinematic graph of the planetary gear

        :param lenght_gear: The width of  the gears. The default is 0.1.
        :type lenght_gear: float, optional

        :param diameter_gear: The diameter of the gears. The default is 1
        :type diameter_gear: float, optional


        :param lenght_double: The lenght of the connections between 2 double planets. The default is 2
        :type lenght_double: float, optional

        :param diameter_pivot: The diameter of the representatives signs for pivot. The default is 0.2.
        :type diameter_pivot: float, optional

        :param lenght_pivot: The length of the representatives signs for pivot. The default is 0.5.
        :type lenght_pivot: float, optional

        :param planet_carrier_x: The parameter for the position of planet carrer in x. The default is 2.
        :type planet_carrier_x: float, optional

        :param planet_carrier_y: The parameter for the position of planet carrer in y. The default is 2.
        :type planet_carrier_y: float, optional

        :param diameter_ring_ini: The diameter of ring.  The default is 10.
        :type diameter_ring_ini: float, optional




        '''
        plot_data = []

        graph_path = self.path_planetary_to_planetary()

        plt.figure()

        previous_relation_double = []
        previous_relation_meshing = []

        previous_planet_meshing = []
        previous_planet_double = []
        inverse_relation_double = []

        coordinate_planet = [[0, 0]]
        coordinate = [0, 0]
        index_coordinate_planet = []
        flag_first_planet = 0
        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0, plot_data)
        for path in graph_path:



            previous_element = 0
            flag_way_inv_double = 0
            coordinate = [0, 0]

            color = 0

            for i, element in enumerate(path):
                if not flag_first_planet and isinstance(element, Planet):
                    if len(coordinate_planet) < 2:
                        index_coordinate_planet.append(element)
                        flag_first_planet = 1


                if isinstance(element, Double):



                    if element in  inverse_relation_double:

                        coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]






                    elif ((element.nodes[0] in previous_planet_double or element.nodes[1] in previous_planet_double) \
                    and not element  in previous_relation_double):

                        for double in previous_relation_double:

                            for node in double.nodes:

                                if element.nodes[0] == node or element.nodes[1] == node:

                                    if  not double == previous_element:

                                        if double in inverse_relation_double:
                                            flag_way_inv_double = 1
                                        else:
                                            flag_way_inv_double = 0

                                    else:
                                        if not double in inverse_relation_double:
                                            flag_way_inv_double = 1
                                        else:
                                            flag_way_inv_double = 0

                        if flag_way_inv_double:

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, lenght_double/(1+i*0.2), color, plot_data)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:
                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color, plot_data)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color, plot_data)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]

                            else:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color, plot_data)
                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]
                        else:

                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                    previous_relation_double.append(element)
                    previous_planet_double.extend([element.nodes[0], element.nodes[1]])

                elif isinstance(element, MeshingPlanet):

                    color += 1

                    if element.type == 'GI':

                        if previous_planet == element.nodes[0]:

                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]
                        else:
                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                    elif element.type == 'GE':

                        if previous_planet == element.nodes[0]:

                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]
                        else:
                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]





                    previous_relation_meshing.append(element)
                    previous_planet_meshing.extend([element.nodes[0], element.nodes[1]])

                if isinstance(element, Planet):
                    previous_planet = element

                if not isinstance(element, Planet) and not isinstance(element, Planetary) \
                and not isinstance(element, MeshingPlanetary):

                    if  not coordinate in coordinate_planet:
                        coordinate_planet.append(coordinate)

                        if not element.nodes[1] in index_coordinate_planet:
                            index_coordinate_planet.append(element.nodes[1])

                        else:

                            index_coordinate_planet.append(element.nodes[0])
                        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color, plot_data)


                    previous_element = element


        coordinate_planet_carrier = self.plot_kinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y, plot_data)
        lenght_ring_ini = 5
        for meshing in self.meshings:

            if isinstance(meshing, MeshingPlanetary):
                color += 1

                if (isinstance(meshing.nodes[0], Planetary) and meshing.nodes[0].planetary_type == 'Sun') \
                or (isinstance(meshing.nodes[1], Planetary) and meshing.nodes[1].planetary_type == 'Sun'):

                    if isinstance(meshing.nodes[0], Planetary):
                        index = index_coordinate_planet.index(meshing.nodes[1])
                    else:
                        index = index_coordinate_planet.index(meshing.nodes[0])

                    planetary_diameter = ((coordinate_planet[index][1]-diameter_gear/2)-coordinate_planet_carrier[1])*2

                    self.plot_kinematic_graph_gear([coordinate_planet[index][0], coordinate_planet_carrier[1]], lenght_gear,
                                                   planetary_diameter, diameter_pivot, lenght_pivot, color, plot_data)

                else:
                    if isinstance(meshing.nodes[0], Planetary):
                        index = index_coordinate_planet.index(meshing.nodes[1])
                    else:
                        index = index_coordinate_planet.index(meshing.nodes[0])

                    lenght_ring = lenght_ring_ini-(((coordinate_planet[index][0])*+100)/50)
                    diameter_ring = diameter_ring_ini-(((coordinate_planet[index][0])*10+100)/50)
                    coordinate_ring = [coordinate_planet[index][0], coordinate_planet[index][1]+diameter_gear/2]

                    self.plot_kinematic_graph_ring(coordinate_ring, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color, plot_data)


        return [pld.PrimitiveGroup(primitives=plot_data)]







    def volmdlr_primitives(self, frame=vm.OXYZ):
        
        axis = self.axis
        center = self.center
        components = self.doubles
        li_box = []
        li_box.extend(self.mesh_generation(axis=axis, center=center))

        for component in components:
            shell = component.volume_model(axis=axis, center=center)
            for shell_planet in shell:
                    li_box.append(shell_planet)

        return li_box

    def plot_data(self, positions_gearing=0):
        if self.planets[0].positions:
            plot_data = []
            primitive_2D = []
            meshing_chains = self.meshing_chain()
            list_color = [pld.colors.BLUE, pld.colors.RED, pld.colors.GREEN, pld.colors.BLACK]

            if self.d_min == 0:
                planetary_gear = PlanetaryGear(self.planetaries, self.planets, self.planet_carrier, self.connections)
                self.d_min = planetary_gear.d_min

            list_d = []
            list_element = self.planetaries+self.planets
            for element in list_element:
                list_d.append(element.module*element.Z)
            list_d_sorted = sorted(list_d)[::-1]

            for d in list_d_sorted:
                i = list_d.index(d)
                element = list_element[i]
                for i, meshing_chain in enumerate(meshing_chains):
                    if element in meshing_chain:
                        color = list_color[i]
                if element in self.planetaries:
                    position = vm.Point2D(0, 0)



                    circle = vm.wires.Circle2D(position, d/2)


                    edge_style = pld.EdgeStyle(line_width=2, color_stroke=color)
                    surface_style = pld.SurfaceStyle(color_fill=pld.colors.WHITE)
                    plot_data.append(circle.plot_data(edge_style=edge_style, surface_style=surface_style))
                else:
                    for position in element.positions:
                        position_2 = vm.Point2D(position[1], position[2])

                        contour = vm.wires.Circle2D(position_2, d/2)


                        edge_style = pld.EdgeStyle(line_width=2, color_stroke=color)
                        surface_style = pld.SurfaceStyle(color_fill=pld.colors.WHITE)
                        plot_data.append(contour.plot_data(edge_style=edge_style, surface_style=surface_style))


            if positions_gearing:

                for position in positions_gearing:
                    point = vm.Point2D(position[1], position[2])
                    plot_data.append(point.plot_data('.', size=5, color=color))

            return [pld.PrimitiveGroup(primitives=plot_data)]

        else:
            return self.plot_kinematic_graph()

    def mesh_generation(self, axis=(1, 0, 0), center=(0, 0, 0)):
        meshing_chains = self.meshing_chain()
        primitives = []
        cycles = {0: 1272321481513.054}
        material = {0:hardened_alloy_steel}
        max_length_meshing_chain = []
        for meshing_chain in meshing_chains:
            module = meshing_chain[0].module
            
           
            rack = meshes_opt.RackOpti(module=[module, module], transverse_pressure_angle_0=[20/180.*npy.pi, 20/180.*npy.pi],
                                       coeff_gear_addendum=[1, 1.25], coeff_gear_dedendum=[1, 1.25], coeff_root_radius=[0.35, 0.4],
                                       coeff_circular_tooth_thickness=[0.5, 0.5])
            list_rack = {0:rack}
            rack_choices = {}
            center_distances = []
            dbs = []
            coefficient_profile_shift = []
            transverse_pressure_angle = {}
            Z = {}
            torques = {0: 'output'}
            connections = []
            centers = {}
            number_primitive_planet = []
            planet = []
            rack_choice = {}

            for i, element in enumerate(meshing_chain):

                Z_2 = element.Z
                if isinstance(element, Planetary) and element.planetary_type == 'Ring':
                    Z_2 = -Z_2
                if i > 0:
                    if isinstance(element, Planetary):
                        torques[i] = element.torque_input[1]
                    else:
                        torques[i] = element.torque

                    connections.append([(i-1, i)])
                    center_distances.append([abs((Z_2+Z_1)*module/2)-0.001, abs((Z_2+Z_1)*module/2)+0.001])
                    transverse_pressure_angle[i-1] = [20/180.*npy.pi-0.00001, 20/180.*npy.pi]

                db = m.cos(20/180.*npy.pi)*module*abs(Z_2)
                dbs.append([db-0.00001, db])
                rack_choices[i] = 0
                Z[i] = Z_2
                
                coefficient_profile_shift.append([0.001, 0.6])

                if isinstance(element, Planetary):
                    centers[i] = element.position
                    if not center == (0, 0, 0):

                        centers[i] = (element.position[0]+center[0], element.position[1]+center[1],
                                      element.position[2]+center[2])
                    else:
                        centers[i] = element.position
                    if not axis == (1, 0, 0):

                        axis_rotation = vm.Vector3D(npy.cross(axis, (1, 0, 0)))

                        axis_vector = vm.Vector3D(axis)
                        axis_origin = vm.Vector3D((1, 0, 0))
                        axis_vector_norme = copy.copy(axis_vector)
                        axis_vector_norme.Normalize()
                        axis_origin_norme = copy.copy(axis_origin)
                        axis_origin_norme.Normalize()
                        if axis_vector_norme.Dot(axis_origin_norme) == 0:
                            angle = m.pi/2
                        else:
                            angle = m.acos(axis_vector.Dot(axis_origin)/axis_vector_norme.Dot(axis_origin_norme))

                        center2 = vm.Vector3D(centers[i])

                        center2.Rotation(center=vm.Vector3D(center), axis=axis_rotation, angle=angle, copy=False)
                        centers[i] = center2.vector



                else:
                    if not center == (0, 0, 0):

                        centers[i] = (element.positions[0][0]+center[0], element.positions[0][1]+center[1],
                                      element.positions[0][2]+center[2])
                    else:

                        centers[i] = element.positions[0]

                    if not axis == (1, 0, 0):
                        axis_rotation = vm.Vector3D(npy.cross(axis, (1, 0, 0)))
                        axis_vector = vm.Vector3D(axis)
                        axis_origin = vm.Vector3D((1, 0, 0))
                        axis_vector_norme = copy.copy(axis_vector)
                        axis_vector_norme.Normalize()
                        axis_origin_norme = copy.copy(axis_origin)
                        axis_origin_norme.Normalize()

                        if axis_vector_norme.Dot(axis_origin_norme) == 0:
                            angle = m.pi/2
                        else:
                            angle = m.acos(axis_vector.Dot(axis_origin)/axis_vector_norme.Dot(axis_origin_norme))

                        center2 = vm.Vector3D(centers[i])

                        center2.Rotation(center=vm.Vector3D(center), axis=axis_rotation, angle=angle, copy=False)
                        centers[i] = center2.vector

                    number_primitive_planet.append(i)
                    planet.append(element)
                previous_element = element
                Z_1 = Z_2

            mesh_optimizer = mg.ContinuousMeshesAssemblyOptimizer(Z=Z, center_distances=center_distances, connections=connections, rigid_links=[],
                                                                  transverse_pressure_angle=transverse_pressure_angle, rack_list=list_rack, rack_choice=rack_choices, material=material,
                                                                  external_torques=torques, cycles=cycles, safety_factor=1, db=dbs, coefficient_profile_shift=coefficient_profile_shift)

            mesh_optimizer.Optimize(verbose=True)
            if not mesh_optimizer.solutions:
                print('Convergence Problem')

            solution = mesh_optimizer.solutions[0].mesh_combinations[0]
            primitive = solution.volmdlr_primitives(axis=axis, centers=centers)
            gear_width_max = 0

            for gear_width in solution.gear_width:
                if gear_width > gear_width_max:
                    gear_width_max = gear_width

            max_length_meshing_chain.append(gear_width_max)
            for j, number in enumerate(number_primitive_planet):
                for n in range(self.number_branch_planet-1):
                    primitive_planet = copy.copy(primitive[number])

                    x = vm.Vector3D(axis[0], axis[1], axis[2])

                    # primitive_planet.Rotation(center=vm.Point3D((0,0,0)),axis=vm.Vector3D((1,0,0)),angle=(n+1)*2*m.pi/self.number_branch_planet)

                    vect_x = -0.5*solution.gear_width[number]*x +  x.dot(vm.Vector3D(centers[number][0], centers[number][1], centers[number][2]))*x

                    vect_center = vm.Vector3D(center[0], center[1], center[2])
                    C2 = primitive_planet.outer_contour2d.rotation(center=vm.Vector2D(vect_center.dot(primitive_planet.x), vect_center.dot(primitive_planet.y)),
                                                                   angle=(n+1)*2*m.pi/self.number_branch_planet)
                    primitive.append(primitives3d.ExtrudedProfile(vm.Vector3D(vect_x[0], vect_x[1], vect_x[2]), primitive_planet.x, primitive_planet.y,
                                                                  C2, [], primitive_planet.extrusion_vector))
            primitives.extend(primitive)

        self.max_length_meshing_chain = max_length_meshing_chain
        return primitives








    def path_planetary_to_planetary(self, planetaries=[]):
        '''
        A function which give all the path betwen the first planetary of the list planetaries (input) and the others

        The path includes the planets and the connections(meshing and Doubles)


        :param planetaries: The first planetary of the list is the beginning of all the path , the others planetaries of the list are the endings of the paths.
                            The default is the list of planetaries .

        :type planetaries: List[Planetary], optional

        :return: list_path
        :rtype: List[List[Planet,meshing,Double,Planetary]]

        '''
        if not planetaries:
            planetaries = self.planetaries

        graph_planetary_gears = self.graph()
        graph_planetary_gears.remove_node(str(self.planet_carrier))
        list_path = []

        for planetary in planetaries[1:]:


            list_path.append(nx.shortest_path(graph_planetary_gears,
                                              str(planetaries[0]), str(planetary)))

        for path in list_path:

            for i in range(len(path)):
                path[i] = nx.get_node_attributes(graph_planetary_gears, path[i])[path[i]]


        return list_path

    # def path_planetary_to_planetary_type(self):
    #     list_path = self.path_planetary_to_planetary()
    #     for i in range(len(list_path)):

    #         for j in range(len(list_path[i])):

    #             if isinstance(list_path[i][j], Planetary):

    #                 list_path[i][j] = list_path[i][j].planetary_type

    #             else:

    #                 list_path[i][j] = str(type(list_path[i][j]))
    #     return list_path



    def reason_abs(self, path):
        '''
        A function wich give the reason ( Willis relation) of a planetary gear

        :param path: The path betwen the two planetaries for which we want to calculate the reason (give by the method path_planetary_to_planetary)
        :type path: List[Planet, meshing, Double]

        :return: reason
        :rtype: float


        '''
        reason = 1
        for i, element in enumerate(path):

            if isinstance(element, Meshing):


                reason = reason*path[i-1].Z/path[i+1].Z


        return reason

    def reason(self, path):
        '''
        A function which give the reason ( Willis relation) of a planetary gear with the coefficient (-1)^n ( n = the number of meshing)

        :param path: The path betwen the two planetaries for which we want to calculate the reason
        :type path: List[Planet, meshing, Double]


        :return: reason
        :rtype: float


        '''
        reason = 1

        for i, element in enumerate(path):

            if isinstance(element, Meshing):

                reason = reason*-path[i-1].Z/path[i+1].Z

                if (isinstance(path[i-1], Planetary) and path[i-1].planetary_type == 'Ring')\
                or (isinstance(path[i+1], Planetary) and path[i+1].planetary_type == 'Ring'):

                    reason = -reason

        return reason

    def speed_range_simplex_intervalle_max(self, intervals_max, planetaries, range_planetary_max, range_planet_carrier,
                                           speed_min_max, reasons):

        c = [-3, -1, 0, 0]
        A = [[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, 1], [0, 1, 0, -1]]
        b = [range_planetary_max[1], -range_planetary_max[0],
             range_planet_carrier[1], -range_planet_carrier[0]]
        speed_diff_1_bound = (0, None)
        speed_diff_2_bound = (0, None)
        speed_1_bound = (None, None)
        speed_2_bound = (None, None)

        for i, planetary in enumerate(planetaries):
            speed_input_planetary = planetary.speed_input
            # reason=planetary_gear.reason(list_path[i][0])

            reason = reasons[i]

            if reason < 0:
                A.extend([[-reason, (1-reason), -reason, -(1-reason)], [-reason, 1-reason, reason, 1-reason]])
            else:

                if reason < 1:
                    A.extend([[reason, 1-reason, -reason, -(1-reason)], [reason, 1-reason, reason, 1-reason]])
                else:
                    A.extend([[reason, -(1-reason), -reason, -(1-reason)], [reason, -(1-reason), reason, (1-reason)]])

            b.extend([-speed_input_planetary[0], speed_input_planetary[1]])

        res = op.linprog(c, A_ub=A, b_ub=b, bounds=[speed_diff_1_bound, speed_diff_2_bound, speed_1_bound, speed_2_bound])

        if res.success:
            return [res.x[0]-res.x[2], res.x[0]+res.x[2]], [res.x[1]-res.x[3], res.x[1]+res.x[3]]
        else:

            return [], []













    def speed_range(self, element_1, element_2, list_planetary=[], generator=0, list_path=[]):
        '''


        A function which give the real speed_range of 2 planetaries ( or planet_carrier) which allow to fulfill
        the condition of input speed of all the other planetaries ( or planet_carrier)
        ( We need to have speed input into all planetaries and planet_carrer)

        :param input_1: The first input
        :type input_1: Planetary or PlanetCarrier

        :param input_2: The second input
        :type input_2: Planetary or PlanetCarrier

        :param list_planetary: The list of planetary that we want to check the input speed condition.
         The default is all the planetaries.
        :type list_planetary: List[Planetary], optional

        :return: A dictionary where all the planetary and planet_carrier are associated with their speed range
        :rtype: Dictionary

        '''
        if isinstance(element_1, PlanetCarrier):
            input_1 = element_2
            input_2 = element_1
        else:
            input_1 = element_1
            input_2 = element_2

        if list_planetary == []:
            list_planetary = copy.copy(self.planetaries)

        range_input_1 = [input_1.speed_input[0], input_1.speed_input[1]]
        range_input_2 = [input_2.speed_input[0], input_2.speed_input[1]]
        range_planet_carrier = copy.copy(self.planet_carrier.speed_input)
        range_planet_carrier = [range_planet_carrier[0], range_planet_carrier[1]]
        ranges_input_1 = []
        ranges_max_input_1 = []


        ranges_planet_carrier = []
        ranges_min_planet_carrier = []
        coeffs_input_1 = []
        coeffs_planet_carrier = []
        reasons = []
        speeds_min_max = []

        if not isinstance(input_1, PlanetCarrier) and not isinstance(input_2, PlanetCarrier):
            range_input_1_for = copy.copy(range_input_1)
            range_planet_carrier_for = copy.copy(range_planet_carrier)
            range_max_input_1 = copy.copy(range_input_1)


            path = self.path_planetary_to_planetary([input_1, input_2])
            reason = self.reason(path[0])
            index = self.matrix_position(self.planet_carrier)

            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))
            coeff_input_2 = 2*(range_input_2[1]-range_input_2[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))
            if reason < 0:
                reason_abs = abs(reason)
                speed_min = (reason*input_1.speed_input[0]-input_2.speed_input[0])/(reason-1)
                speed_max = (reason*input_1.speed_input[1]-input_2.speed_input[1])/(reason-1)


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1_for[0] += coeff_input_1*speed_diff

                    speed_diff_2_min = (self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)
                    speed_diff_1_max = 0
                    if speed_diff_2_min > input_2.speed_input[1]-input_2.speed_input[0]:
                        speed_diff_2_min = input_2.speed_input[1]-input_2.speed_input[0]
                        speed_diff_1_max = ((self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)-speed_diff_2_min)/reason_abs


                    range_max_input_1[0] += speed_diff_1_max


                    range_input_2[0] += coeff_input_2*speed_diff

                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier_for[0] = speed_min

                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1_for[1] -= coeff_input_1*speed_diff
                    speed_diff_2_min = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)
                    speed_diff_1_max = 0
                    if speed_diff_2_min > input_2.speed_input[1]-input_2.speed_input[0]:
                        speed_diff_2_min = input_2.speed_input[1]-input_2.speed_input[0]
                        speed_diff_1_max = ((speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)-speed_diff_2_min)/reason_abs

                    range_max_input_1[1] -= speed_diff_1_max

                    range_input_2[1] -= coeff_input_2*speed_diff

                elif speed_max > self.planet_carrier.speed_input[0]:
                    range_planet_carrier_for[1] = speed_max

                else:
                    return []


            else:


                if reason < 1:
                    speed_min = (reason*input_1.speed_input[1]-input_2.speed_input[0])/(reason-1)
                    speed_max = (reason*input_1.speed_input[0]-input_2.speed_input[1])/(reason-1)

                else:
                    speed_min = (reason*input_1.speed_input[0]-input_2.speed_input[1])/(reason-1)
                    speed_max = (reason*input_1.speed_input[1]-input_2.speed_input[0])/(reason-1)


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1-reason)/(coeff_input_2+coeff_input_1*reason)
                    speed_diff_2_min = (self.planet_carrier.speed_input[0]-speed_min)*(1-reason)
                    speed_diff_1_max = 0
                    if speed_diff_2_min > input_2.speed_input[1]-input_2.speed_input[0]:
                            speed_diff_2_min = input_2.speed_input[1]-input_2.speed_input[0]
                            speed_diff_1_max = ((self.planet_carrier.speed_input[0]-speed_min)*(1-reason)-speed_diff_2_min)/reason

                    if reason < 1:
                        range_input_1_for[1] -= coeff_input_1*speed_diff

                        range_max_input_1[1] -= speed_diff_1_max

                        range_input_2[0] += coeff_input_2*speed_diff

                    else:
                        range_input_1_for[0] -= coeff_input_1*speed_diff
                        range_max_input_1[0] -= speed_diff_1_max

                        range_input_2[1] += coeff_input_2*speed_diff

                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier_for[0] = speed_min

                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)/(coeff_input_2+coeff_input_1*reason)
                    speed_diff_2_min = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)
                    speed_diff_1_max = 0
                    if speed_diff_2_min > input_2.speed_input[1]-input_2.speed_input[0]:
                            speed_diff_2_min = input_2.speed_input[1]-input_2.speed_input[0]
                            speed_diff_1_max = ((speed_max-self.planet_carrier.speed_input[1])*(1-reason)-speed_diff_2_min)/reason

                    if reason < 1:
                        range_input_1_for[0] += coeff_input_1*speed_diff

                        range_max_input_1[0] += speed_diff_1_max

                        range_input_2[1] -= coeff_input_2*speed_diff

                    else:
                        range_input_1_for[1] += coeff_input_1*speed_diff

                        range_max_input_1[1] += speed_diff_1_max

                        range_input_2[0] -= coeff_input_2*speed_diff

                elif speed_max > self.planet_carrier.speed_input[0]:
                    range_planet_carrier_for[1] = speed_max

                else:
                    return []


            ranges_planet_carrier.append(range_planet_carrier_for)
            ranges_min_planet_carrier.append(range_planet_carrier_for)

            ranges_input_1.append(range_input_1_for)
            ranges_max_input_1.append(range_max_input_1)

            list_planetary.remove(input_1)
            list_planetary.remove(input_2)
            input_1_for = input_1

        else:
            if  isinstance(input_1, PlanetCarrier):
                list_planetary.remove(input_2)
                input_1_for = input_2

                range_input_1 = copy.copy(range_input_2)
            else:
                list_planetary.remove(input_1)
                input_1_for = input_1

        range_output = {}
        for i, planetary in enumerate(list_planetary):
            range_for_planetary_input_1 = copy.copy(range_input_1)
            range_for_planetary_planet_carrier = copy.copy(range_planet_carrier)
            range_for_planetary_planet_carrier = [range_for_planetary_planet_carrier[0], range_for_planetary_planet_carrier[1]]
            range_output[planetary] = [planetary.speed_input[0], planetary.speed_input[1]]

            range_max_input_1 = copy.copy(range_input_1)
            range_min_planet_carrier = copy.copy(range_planet_carrier)

            if not list_path:
                path = self.path_planetary_to_planetary([input_1_for, planetary])

                reason = self.reason(path[0])

            else:
                path = list_path[i]

                reason = self.reason(path[0])

            index = self.matrix_position(planetary)

            if range_input_1[1] == range_input_1[0] and range_planet_carrier[1] == range_planet_carrier[0]:
                return []

            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))
            coeff_input_planet_carrier = 2*(range_planet_carrier[1]-range_planet_carrier[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))

            if reason < 0:

                speed_min = reason*range_input_1[1]+(1-reason)*range_planet_carrier[0]
                speed_max = reason*range_input_1[0]+(1-reason)*range_planet_carrier[1]
                reason_abs = abs(reason)


                if speed_min < planetary.speed_input[0]:

                    speed_diff = (planetary.speed_input[0]-speed_min)/((1+reason_abs)*coeff_input_planet_carrier+reason_abs*coeff_input_1)

                    range_for_planetary_input_1[1] -= coeff_input_1*speed_diff
                    range_for_planetary_planet_carrier[0] += coeff_input_planet_carrier*speed_diff

                    speed_diff_2_min = (planetary.speed_input[0]-speed_min)/(1+reason_abs)

                    if speed_diff_2_min > self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]:
                        speed_diff_2_min = self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]
                        speed_diff_1_max = (planetary.speed_input[0]-speed_min-(speed_diff_2_min)*(1+reason_abs))/(reason_abs)
                        range_max_input_1[1] -= speed_diff_1_max

                    range_min_planet_carrier[0] += speed_diff_2_min




                elif speed_min < planetary.speed_input[1]:
                    range_output[planetary][0] = speed_min

                else:

                    return[]


                if speed_max > planetary.speed_input[1]:
                    speed_diff = (speed_max-planetary.speed_input[1])/((1+reason_abs)*coeff_input_planet_carrier+reason_abs*coeff_input_1)
                    range_for_planetary_input_1[0] += coeff_input_1*speed_diff
                    range_for_planetary_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff

                    speed_diff_2_min = (speed_max-planetary.speed_input[1])/(1+reason_abs)

                    if speed_diff_2_min > self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]:
                        speed_diff_2_min = self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]
                        speed_diff_1_max = (speed_max-planetary.speed_input[1]-(speed_diff_2_min)*(1+reason_abs))/(reason_abs)
                        range_max_input_1[0] += speed_diff_1_max

                    range_min_planet_carrier[1] -= speed_diff_2_min


                elif speed_max > planetary.speed_input[0]:
                    range_output[planetary][1] = speed_max
                else:

                    return []


            else:

                if reason < 1:

                    speed_min = reason*range_input_1[0]+(1-reason)*range_planet_carrier[0]
                    speed_max = reason*range_input_1[1]+(1-reason)*range_planet_carrier[1]

                else:
                    speed_min = reason*range_input_1[0]+(1-reason)*range_planet_carrier[1]
                    speed_max = reason*range_input_1[1]+(1-reason)*range_planet_carrier[0]

                if speed_min < planetary.speed_input[0]:

                    if reason < 1:
                        speed_diff = (planetary.speed_input[0]-speed_min)/(coeff_input_1*reason+coeff_input_planet_carrier*(1-reason))

                        range_for_planetary_input_1[0] += coeff_input_1*speed_diff
                        range_for_planetary_planet_carrier[0] += coeff_input_planet_carrier*speed_diff




                    else:
                        speed_diff = (planetary.speed_input[0]-speed_min)/(reason*coeff_input_1+coeff_input_planet_carrier*(reason-1))
                        range_for_planetary_input_1[0] += coeff_input_1*speed_diff
                        range_for_planetary_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff





                elif speed_min < planetary.speed_input[1]:
                    range_output[planetary][0] = speed_min

                else:

                    return []


                if speed_max > planetary.speed_input[1]:

                    if reason < 1:
                        speed_diff = (speed_max-planetary.speed_input[1])/(coeff_input_1*reason+coeff_input_planet_carrier*(1-reason))
                        range_for_planetary_input_1[1] -= coeff_input_1*speed_diff
                        range_for_planetary_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff







                    else:
                        speed_diff = (speed_max-planetary.speed_input[1])/(reason*coeff_input_1+coeff_input_planet_carrier*(reason-1))
                        range_for_planetary_input_1[1] -= coeff_input_1*speed_diff
                        range_for_planetary_planet_carrier[0] += coeff_input_planet_carrier*speed_diff








                elif speed_max > planetary.speed_input[0]:
                    range_output[planetary][1] = speed_max

                else:

                    return []

            ranges_input_1.append(range_for_planetary_input_1)
            ranges_planet_carrier.append(range_for_planetary_planet_carrier)

            ranges_min_planet_carrier.append(range_min_planet_carrier)

            ranges_max_input_1.append(range_max_input_1)
            speeds_min_max.append([speed_min, speed_max])
            reasons.append(reason)
            coeffs_input_1.append(coeff_input_1)
            coeffs_planet_carrier.append(coeff_input_planet_carrier)


        for i, range_input_1_for in enumerate(ranges_input_1):

            if range_input_1[0] > range_input_1_for[1] or range_input_1[1] < range_input_1_for[0]:
                if generator:
                    return 'simplex'

                range_input_1 = []
                break
            if range_input_1[1] > range_input_1_for[1]:
                range_input_1[1] = range_input_1_for[1]

            if range_input_1[0] < range_input_1_for[0]:
                range_input_1[0] = range_input_1_for[0]



        for i, range_planet_carrier_for in enumerate(ranges_planet_carrier):

            if range_planet_carrier[0] > range_planet_carrier_for[1] or range_planet_carrier[1] < range_planet_carrier_for[0]:
                range_input_1 = []
                if generator:
                    return 'simplex'
                break
            if range_planet_carrier[1] > range_planet_carrier_for[1]:
                range_planet_carrier[1] = range_planet_carrier_for[1]

            if range_planet_carrier[0] < range_planet_carrier_for[0]:
                range_planet_carrier[0] = range_planet_carrier_for[0]

        if not range_input_1:

            range_input_1, range_planet_carrier = self.speed_range_simplex_intervalle_max(ranges_max_input_1, list_planetary,
                                                                                          input_1.speed_input, self.planet_carrier.speed_input,
                                                                                          speeds_min_max, reasons)
            if not range_input_1:
                return []

        if not isinstance(input_1, PlanetCarrier) and not isinstance(input_2, PlanetCarrier):

            path = self.path_planetary_to_planetary([input_1, input_2])
            reason = self.reason(path[0])
            index = self.matrix_position(self.planet_carrier)


            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))
            coeff_input_2 = 2*(range_input_2[1]-range_input_2[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))
            if reason < 0:
                reason_abs = abs(reason)
                speed_min = (reason*input_1.speed_input[0]-input_2.speed_input[0])/(reason-1)
                speed_max = (reason*input_1.speed_input[1]-input_2.speed_input[1])/(reason-1)


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1[0] += coeff_input_1*speed_diff



                    range_input_2[0] += coeff_input_2*speed_diff


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1[1] -= coeff_input_1*speed_diff

                    range_input_2[1] -= coeff_input_2*speed_diff


            else:


                if reason < 1:
                    speed_min = (reason*input_1.speed_input[1]-input_2.speed_input[0])/(reason-1)
                    speed_max = (reason*input_1.speed_input[0]-input_2.speed_input[1])/(reason-1)

                else:
                    speed_min = (reason*input_1.speed_input[0]-input_2.speed_input[1])/(reason-1)
                    speed_max = (reason*input_1.speed_input[1]-input_2.speed_input[0])/(reason-1)


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1-reason)/(coeff_input_2+coeff_input_1*reason)

                    if reason < 1:
                        range_input_1[1] -= coeff_input_1*speed_diff
                        range_input_2[0] += coeff_input_2*speed_diff

                    else:
                        range_input_1[0] -= coeff_input_1*speed_diff
                        range_input_2[1] += coeff_input_2*speed_diff


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)/(coeff_input_2+coeff_input_1*reason)

                    if reason < 1:
                        range_input_1[0] += coeff_input_1*speed_diff

                        range_input_2[1] -= coeff_input_2*speed_diff

                    else:
                        range_input_1[1] += coeff_input_1*speed_diff
                        range_input_2[0] -= coeff_input_2*speed_diff





            range_output[input_1] = range_input_1
            range_output[input_2] = range_input_2
            range_output[self.planet_carrier] = range_planet_carrier

            return range_output


        range_output[input_1_for] = range_input_1
        range_output[self.planet_carrier] = range_planet_carrier
        return range_output


    def torque_min_max_signe_input(self, element1, element2):
        element1_signe = []
        element2_signe = []



        if isinstance(element2, PlanetCarrier):
            list_planetary = copy.copy(self.planetaries)
            list_planetary.remove(element1)
            list_path = self.path_planetary_to_planetary([element1]+list_planetary)
            for path in list_path:
                reason = self.reason(path)
                if reason < 0:
                    element1_signe.append(1)
                    element2_signe.append(-1)
                else:
                    element1_signe.append(-1)
                    if reason < 1:
                        element2_signe.append(-1)
                    else:
                        element2_signe.append(1)
        else:
            list_planetary = copy.copy(self.planetaries)
            list_planetary.remove(element1)
            list_planetary.remove(element2)
            list_path = self.path_planetary_to_planetary([element1, element2]+list_planetary)

            reason_1 = self.reason(list_path[0][0])
            list_path.remove(list_path[0])
            for path in list_path:
                reason = self.reason(path[0])

                if (reason+reason_1*(1-reason)/(reason_1-1)) < 0:
                    element1_signe.append(1)
                else:
                    element1_signe.append(-1)

                if  ((1-reason)/(reason_1-1)) < 0:
                    element2_signe.append(-1)
                else:
                    element2_signe.append(1)



            if reason_1/(reason_1-1) < 0:
                element1_signe.append(1)

            else:
                element1_signe.append(-1)

            if 1/(reason_1-1) < 0:
                element2_signe.append(-1)

            else:
                element2_signe.append(1)


        return element1_signe, element2_signe


    def torque_range(self, elements):

        element_list_2 = copy.copy(self.planetaries)
        for element in elements:
            if element in element_list_2:
                element_list_2.remove(element)

        element1 = element_list_2[0]

        if len(element_list_2) > 1:
            element2 = element_list_2[1]
        else:
            element2 = self.planet_carrier



        num_var = (len(self.planetaries)-1)*2
        num_eq = (len(self.planetaries)-1+2)*2
        A = np.zeros((num_eq, num_var))
        y = 0
        b = []
        c = [0]*num_var


        bounds = []



        if not isinstance(element2, PlanetCarrier):
            list_planetary = copy.copy(self.planetaries)
            list_planetary.remove(element1)
            list_planetary.remove(element2)
            list_path = self.path_planetary_to_planetary([element1, element2]+list_planetary)
            reason_second_planetary = self.reason(list_path[0][0])
            list_path.remove(list_path[0])

            for i in range(int(num_var/2)):
                A[y][2*i] = 1
                A[y][2*i+1] = 1
                y += 1
                A[y][2*i] = -1
                A[y][2*i+1] = 1
                y += 1
                c[2*i+1] = -1
                bounds.extend([(None, None), (0, None)])
                if i != num_var/2-1:
                    b.extend([list_planetary[i].torque_input[1], -list_planetary[i].torque_input[0]])
                else:
                    position_planet_carrier = 2*i
                    b.extend([self.planet_carrier.torque_input[1], -self.planet_carrier.torque_input[0]])


            b.extend([element1.torque_input[1], -element1.torque_input[0], element2.torque_input[1], -element2.torque_input[0]])

            for i, planetary in enumerate(list_planetary):
                reason_planetary = self.reason(list_path[i][0])

                coefficient_1 = -(reason_planetary+(1-reason_planetary)*reason_second_planetary/(reason_second_planetary-1))
                if coefficient_1 < 0:
                    A[-4][2*i] = coefficient_1
                    A[-4][2*i+1] = -coefficient_1
                    A[-3][2*i] = -coefficient_1
                    A[-3][2*i+1] = -coefficient_1
                else:
                    A[-4][2*i] = coefficient_1
                    A[-4][2*i+1] = coefficient_1
                    A[-3][2*i] = -coefficient_1
                    A[-3][2*i+1] = coefficient_1
                coefficient_2 = (1-reason_planetary)/(reason_second_planetary-1)

                if coefficient_2 < 0:
                    A[-2][2*i] = coefficient_2
                    A[-2][2*i+1] = -coefficient_2
                    A[-1][2*i] = -coefficient_2
                    A[-1][2*i+1] = -coefficient_2
                else:
                    A[-2][2*i] = coefficient_2
                    A[-2][2*i+1] = coefficient_2
                    A[-1][2*i] = -coefficient_2
                    A[-1][2*i+1] = coefficient_2


            coefficient_1 = -reason_second_planetary/(reason_second_planetary-1)

            if coefficient_1 < 0:
                A[-4][position_planet_carrier] = coefficient_1
                A[-4][position_planet_carrier+1] = -coefficient_1
                A[-3][position_planet_carrier] = -coefficient_1
                A[-3][position_planet_carrier+1] = -coefficient_1
            else:
                A[-4][position_planet_carrier] = coefficient_1
                A[-4][position_planet_carrier+1] = coefficient_1
                A[-3][position_planet_carrier] = -coefficient_1
                A[-3][position_planet_carrier+1] = coefficient_1

            coefficient_2 = 1/(reason_second_planetary-1)
            if coefficient_2 < 0:
                A[-2][position_planet_carrier] = coefficient_2
                A[-2][position_planet_carrier+1] = -coefficient_2
                A[-1][position_planet_carrier] = -coefficient_2
                A[-1][position_planet_carrier+1] = -coefficient_2
            else:
                A[-2][position_planet_carrier] = coefficient_2
                A[-2][position_planet_carrier+1] = coefficient_2
                A[-1][position_planet_carrier] = -coefficient_2
                A[-1][position_planet_carrier+1] = coefficient_2

        else:
            list_planetary = copy.copy(self.planetaries)
            list_planetary.remove(element1)
            list_path = self.path_planetary_to_planetary([element1]+list_planetary)

            for i in range(int(num_var/2)):
                A[y][2*i] = 1
                A[y][2*i+1] = 1
                y += 1
                A[y][2*i] = -1
                A[y][2*i+1] = 1
                y += 1
                c[2*i+1] = -1
                bounds.extend([(None, None), (0, None)])

                b.extend([list_planetary[i].torque_input[1], -list_planetary[i].torque_input[0]])





            b.extend([element1.torque_input[1], -element1.torque_input[0], element2.torque_input[1], -element2.torque_input[0]])

            for i, planetary in enumerate(list_planetary):

                reason_planetary = self.reason(list_path[i])

                coefficient_1 = -reason_planetary
                if coefficient_1 < 0:
                    A[-4][2*i] = coefficient_1
                    A[-4][2*i+1] = -coefficient_1
                    A[-3][2*i] = -coefficient_1
                    A[-3][2*i+1] = -coefficient_1
                else:
                    A[-4][2*i] = coefficient_1
                    A[-4][2*i+1] = coefficient_1
                    A[-3][2*i] = -coefficient_1
                    A[-3][2*i+1] = coefficient_1
                coefficient_2 = -(1-reason_planetary)

                if coefficient_2 < 0:
                    A[-2][2*i] = coefficient_2
                    A[-2][2*i+1] = -coefficient_2
                    A[-1][2*i] = -coefficient_2
                    A[-1][2*i+1] = -coefficient_2
                else:
                    A[-2][2*i] = coefficient_2
                    A[-2][2*i+1] = coefficient_2
                    A[-1][2*i] = -coefficient_2
                    A[-1][2*i+1] = coefficient_2




        res = op.linprog(c, A_ub=A, b_ub=b, bounds=bounds)

        if res.success:
            result = {}
            for i, planetary in enumerate(list_planetary):
                result[planetary] = [res.x[2*i]-res.x[2*i+1], res.x[2*i]+res.x[2*i+1]]


            return result
        else:
            return []





    def meshing_chain_recursive_function(self, number_meshing_chain, element, graph_planetary_gear, possibilities,
                                         list_possibilities, previous_relation, meshing_way, previous_meshing_chains):


        neighbors = nx.all_neighbors(graph_planetary_gear, str(element))
        meshing = 0
        neighbors_2 = []
        for neighbor in neighbors:
            neighbor = nx.get_node_attributes(graph_planetary_gear, neighbor)[neighbor]

            neighbors_2.append(neighbor)

        if not possibilities:
            possibilities.append(element)

        if previous_relation:

            neighbors_2.remove(previous_relation)




        for neighbor in neighbors_2:

            meshing = copy.copy(meshing)

            previous_relation = neighbor

            if isinstance(neighbor, Meshing):

                if meshing == 0 or meshing_way == 0:
                    if neighbor.nodes[0] == element:
                        possibilities.append(neighbor.nodes[1])
                        element_3 = neighbor.nodes[1]

                    else:
                        possibilities.append(neighbor.nodes[0])
                        element_3 = neighbor.nodes[0]

                    self.meshing_chain_recursive_function(number_meshing_chain, element_3, graph_planetary_gear, possibilities,
                                                          list_possibilities, previous_relation, 0, previous_meshing_chains)
                    meshing_way = -1

                if meshing == 1 or meshing_way == 1:

                    if neighbor.nodes[0] == element:
                        possibilities = [neighbor.nodes[1]] + possibilities
                        element_3 = neighbor.nodes[1]

                    else:
                        possibilities = [neighbor.nodes[0]] + possibilities

                        element_3 = neighbor.nodes[0]
                    self.meshing_chain_recursive_function(number_meshing_chain, element_3, graph_planetary_gear, possibilities,
                                                          list_possibilities, previous_relation, 1, previous_meshing_chains)
                meshing = 1


            elif isinstance(neighbor, Double):

                if neighbor.nodes[0] == element:

                    element_3 = neighbor.nodes[1]

                else:
                    element_3 = neighbor.nodes[0]

                self.meshing_chain_recursive_function(number_meshing_chain+1, element_3, graph_planetary_gear, [],
                                                      list_possibilities, previous_relation, 0, previous_meshing_chains)



        if not meshing:

                if not number_meshing_chain in previous_meshing_chains:

                    previous_meshing_chains.append(number_meshing_chain)
                    list_possibilities.append(possibilities)

                else:
                    index = previous_meshing_chains.index(number_meshing_chain)
                    list_possibilities[index] = possibilities




        return list_possibilities

    def meshing_chain(self):
        """
        A function wich return all the meshing chain in the planetary gear.
        A meshing chain is a list of planetaries and planets which meshed together

        :return: the list of meshing_chains
        :rtype: List[List[Planetary,Planet]]

        """
        graph_planetary_gear = self.graph()

        list_possibilities = self.meshing_chain_recursive_function(0, self.planetaries[0], graph_planetary_gear, [], [], 0, 0, [])

        return list_possibilities

    def meshing_chain_position_z(self, meshing_chains):
        z = [0]*len(meshing_chains)
        length = self.length_double
        z[0] = 0
        z_ini = 0

        doubles = copy.copy(self.doubles)
        orientation = [0]*len(meshing_chains)
        orientation[0] = 0

        for i, meshing_chain in enumerate(meshing_chains):
            number_double = 0

            for double in doubles:
                if double.nodes[0] in meshing_chain:
                     for j, meshing_chain in enumerate(meshing_chains):
                         if double.nodes[1] in meshing_chain:

                             if orientation[i] == 0:

                                 if number_double == 0:
                                     z[j] = z[i]+length
                                     number_double += 1
                                     orientation[j] = 1
                                 else:
                                     z[j] = z[i]-length
                                     orientation[j] = -1

                             else:
                                 z[j] = z[i]+orientation[i]*length
                                 orientation[j] = orientation[i]
                     doubles.remove(double)

                elif double.nodes[1] in meshing_chain:

                         for j, meshing_chain in enumerate(meshing_chains):
                                 if double.nodes[0] in meshing_chain:

                                     if orientation[i] == 0:

                                         if number_double == 0:
                                             z[j] = z[i]+length
                                             number_double += 1
                                             orientation[j] = 1
                                         else:
                                             z[j] = z[i]-length
                                             orientation[j] = -1

                                     else:
                                         z[j] = z[i]+orientation[i]*length
                                         orientation[j] = orientation[i]


                         doubles.remove(double)
        return z


    def speed_max_planets(self):
        speed_max = 0

        for planetary in self.planetaries:
            if planetary.planetary_type == 'Ring':

                speed_max_planetary = self.speed_solve({planetary:planetary.speed_input[1], self.planet_carrier:self.planet_carrier.speed_input[0]})
            else:
                speed_max_planetary = self.speed_solve({planetary:planetary.speed_input[0], self.planet_carrier:self.planet_carrier.speed_input[1]})

            for planet in self.planets:

                if abs(speed_max_planetary[planet]) > speed_max:

                    speed_max = abs(speed_max_planetary[planet])


        return speed_max

    def torque_max_planets(self):
        first_input = {}
        # for planetary in self.planetaries[1:]:
        #     first_input[planetary]=planetary.torque_input[1]

        # # first_input[self.planet_carrier]=self.planet_carrier.torque_input[1]
        # print(first_input)
        self.mech = 0
        if not self.mech:
            for i in range(len(self.planetaries)-1):
                first_input[self.planetaries[i]] = self.planetaries[i].torque_input[0]/self.number_branch_planet

            self.torque_resolution_PFS(first_input)

        torque_max_planet = 0
        speed_planetaries = [0]*len(self.planetaries)
        torque_planetaries = [0]*len(self.planetaries)
        for i, gearing_planetary in enumerate(self.mech_dict['gearings_planetary']):
            try:
                power = self.mech.TransmittedLinkagePower(gearing_planetary, 4)
            except:
                power = 1000000

            if gearing_planetary.part1 in self.mech_dict['part_planetaries']:
                index = self.mech_dict['part_planetaries'].index(gearing_planetary.part1)
            else:
                index = self.mech_dict['part_planetaries'].index(gearing_planetary.part2)
            speed = self.mech.kinematic_results[self.mech_dict['link_planetaries_ball'][index]][0]

            if speed == 0:
                torque_planetaries[index] = 0
            else:
                torque_planetaries[index] = (power/speed)
            speed_planetaries[index] = speed

        if self.mech_dict['flag_load_planet_carrier_unknow']:

            torque = self.mech.static_results[self.mech_dict['load_planet_carrier']][0]
        else:
            torque = self.mech_dict['input_torque_and_composant'][self.planet_carrier]


        # print(mech.kinematic_results)
        speed = self.mech.kinematic_results[self.mech_dict['planet_carrier_a']][0]
        planet_carrier_torque = torque
        planet_carrier_speed = speed

        powers = []
        power_element = []
        torques = []
        for i, planetary in enumerate(self.planetaries):
            if planetary.torque_input[1] < 0:
                powers.append(speed_planetaries[i]*planetary.torque_input[0])
            else:
                powers.append(speed_planetaries[i]*planetary.torque_input[1])
            power_element.append(planetary)

        if self.planet_carrier.torque_input[1] < 0:
                powers.append(planet_carrier_speed*self.planet_carrier.torque_input[0])
        else:
                powers.append(planet_carrier_speed*self.planet_carrier.torque_input[1])


        power_element.append(self.planet_carrier)
        speed_planetaries.append(planet_carrier_speed)

        planetary_input = []
        torques = []

        power_max = max(powers)
        power_tot = power_max
        planetary_input.append(power_element[powers.index(power_max)])
        if speed_planetaries[powers.index(power_max)] != 0:
            torques.append(power_max/speed_planetaries[powers.index(power_max)])
        else:
            torques.append(0)
        power_element.remove(power_element[powers.index(power_max)])
        speed_planetaries.remove(speed_planetaries[powers.index(power_max)])
        powers.remove(power_max)

        for i in range(len(self.planetaries)-2):
            power_max = 0
            for y, planetary in enumerate(power_element):
                if power_tot+powers[y] > power_max:
                    index_planetary = y
                    power_max = power_tot+powers[y]
            if power_max == 0:
                index_planetary = 0

            power_tot = power_max

            planetary_input.append(power_element[index_planetary])
            if speed_planetaries[index_planetary] != 0:
                torques.append(power_max/speed_planetaries[index_planetary])
            else:
                torques.append(0)

            power_element.remove(power_element[index_planetary])
            speed_planetaries.remove(speed_planetaries[index_planetary])
            powers.remove(powers[index_planetary])





        max_input = {}

        for i, planetary in enumerate(planetary_input):

            max_input[planetary] = planetary.torque_input[i]/self.number_branch_planet


        self.update_load_mech(max_input)

        for i, pivot in enumerate(self.mech_dict['pivot_planets']):
            try:
                power = self.mech.TransmittedLinkagePower(pivot, 1)
            except genmechanics.ModelError:
                power = 100000000000
            speed = self.mech.kinematic_results[pivot][0]
            if speed == 0:
                torque = 0
            else:
                torque = power/speed
            self.planets[i].torque = torque


            # print(torque)
            if abs(torque) > torque_max_planet:
                torque_max_planet = abs(torque)


        return  torque_max_planet


    def test_assembly_condition(self, number_planet, planetaries=[]):
        '''
        A function which test the assembly condition for the planetary gear

        :param number_planet: The number of planet which are arround the planetary gear ( exemple: 3,4 or 5)
        :type number_planet: Int

        :param planetaries: The list of the two planetary which we want to test the assembly condition. The default is all the planetary of the planetary gear.
        :type planetaries: List[Planetary], optional

        :return: The result of the test
        :rtype: Boolean


        '''
        if not planetaries:
            planetaries = self.planetaries
        valid = True
        list_path = self.path_planetary_to_planetary(planetaries)



        for path in list_path:
            if valid:

                basic_ratio_planet_1 = 1
                basic_ratio_planet_2 = 1
                planet_1 = 0
                planet_2 = 0

                for obj in path:

                    if isinstance(obj, Double):
                        planet_1 = obj.nodes[0]
                        planet_2 = obj.nodes[1]

                        break

                if not planet_1:

                    for obj in path:

                        if isinstance(obj, Planet):
                            planet_1 = obj

                            break

                list_nodes_2 = path

                position_planet = list_nodes_2.index(planet_1)
                inv_list_nodes = list_nodes_2[:position_planet+1]
                inv_list_nodes = inv_list_nodes[::-1]


                for j, node_2 in enumerate(inv_list_nodes):

                    if isinstance(node_2, (MeshingPlanetary, MeshingPlanet)):

                        basic_ratio_planet_1 = basic_ratio_planet_1*\
                        (-inv_list_nodes[j-1].Z/inv_list_nodes[j+1].Z)

                        if isinstance(inv_list_nodes[j+1], Planetary):

                            if inv_list_nodes[j+1].planetary_type == 'Ring':

                                basic_ratio_planet_1 = basic_ratio_planet_1 *-1



                for j, node_2 in enumerate(list_nodes_2[position_planet:]):

                    if isinstance(node_2, (MeshingPlanetary, MeshingPlanet)):

                        basic_ratio_planet_2 = basic_ratio_planet_2 * \
                        (-list_nodes_2[position_planet+j-1].Z / list_nodes_2[position_planet+j+1].Z)

                        if isinstance(list_nodes_2[position_planet+j+1], Planetary):

                            if list_nodes_2[position_planet+j+1].planetary_type == 'Ring':

                                basic_ratio_planet_2 = basic_ratio_planet_2 *-1


                if planet_2:
                    equation = (1/number_planet)*\
                        (1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z*planet_2.Z)

                else:
                    equation = (1/number_planet)*\
                        (1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z)

                valid = (int(equation) == equation)
        return valid



    def speed_system_equations(self):
        #initialize system matrix
        n_equations = len(self.relations)
        n_variables = len(self.elements)
        system_matrix = npy.zeros((n_equations, n_variables))
        rhs = npy.zeros(n_equations)

        for i, relation in enumerate(self.relations):
            matrix_relation, rhs_relation = relation.speed_system_equations()

            if isinstance(relation, MeshingPlanetary):

                if isinstance(relation.nodes[0], Planetary):
                    matrix_position_planetary = self.matrix_position(relation.nodes[0])
                    matrix_position_planet = self.matrix_position(relation.nodes[1])

                else:
                    matrix_position_planetary = self.matrix_position(relation.nodes[1])
                    matrix_position_planet = self.matrix_position(relation.nodes[0])

                system_matrix[i][matrix_position_planetary] = matrix_relation[0]
                system_matrix[i][matrix_position_planet] = matrix_relation[1]
                system_matrix[i][-1] = matrix_relation[2]

                rhs[i] = rhs_relation[0]

            else:

                matrix_position_planet_1 = self.matrix_position(relation.nodes[0])
                matrix_position_planet_2 = self.matrix_position(relation.nodes[1])
                system_matrix[i][matrix_position_planet_1] = matrix_relation[0]
                system_matrix[i][matrix_position_planet_2] = matrix_relation[1]

                rhs[i] = rhs_relation[0]


        return system_matrix, rhs

    def speed_solve(self, input_speeds_and_composants):
        '''
        A function which give the speed of all the elements(Planetary,Planet,Planet_carrier) of planetary gear
        whith 2 input elements and speeds

        :param input_speeds_and_composants: A dictionary where the element input are associated with their speed input
        :type input_speeds_and_composants: Dictionary{Planetary, Planet, PlanetCarrier : float}

        :return: A list where the first elements are the speeds of the planetaries, then, there are the speeds of the planets,
                and to finish the last element is the speed of the planet_carrier.
                We can know the position of an element in the list by using the function matrix position whith the element in input

        :rtype: List[float]


        '''

        system_matrix, vector_b = self.speed_system_equations()
        n_equations = len(self.relations)
        n_variables = len(self.elements)

        system_matrix_speed_solve_0 = npy.zeros((len(input_speeds_and_composants), n_variables))
        vector_b_speed_solve_0 = npy.zeros(len(input_speeds_and_composants))

        system_matrix = npy.concatenate((system_matrix, system_matrix_speed_solve_0), axis=0)
        vector_b = npy.concatenate((vector_b, vector_b_speed_solve_0), axis=0)
        impose_speeds = []

        for i, composant in enumerate(input_speeds_and_composants):
            if isinstance(input_speeds_and_composants[composant], PlanetCarrier) or isinstance(input_speeds_and_composants[composant], Planetary):
                position_element_1 = self.matrix_position(composant)
                position_element_2 = self.matrix_position(input_speeds_and_composants[composant])
                system_matrix[n_equations][position_element_1] = 1
                system_matrix[n_equations][position_element_2] = -1
                vector_b[n_equations] = 0
                n_equations += 1
            else:
                impose_speeds.append(ImposeSpeed(composant, input_speeds_and_composants[composant], 'ImposeSpeed'+str(i)))

        for impose_speed in impose_speeds:
            position_element = self.matrix_position(impose_speed.node)
            system_matrix[n_equations][position_element] = impose_speed.speed_system_equations()[0]
            vector_b[n_equations] = impose_speed.speed_system_equations()[1]
            n_equations += 1


        solution = solve(system_matrix, vector_b)


        element_association = {}
        for i in range(len(self.elements)):
            self.elements[i].speed = solution[i]
            element_association[self.elements[i]] = solution[i]

        return element_association


    # def path_planetary_to_double(self):
    #     graph_planetary_gears = self.graph()
    def torque_system_equation(self):
        n_meshing_planetary = 0
        for meshing in self.meshings:
            if isinstance(meshing, MeshingPlanetary):
                n_meshing_planetary += 1

        n_equations = (len(self.meshings)-n_meshing_planetary)+n_meshing_planetary*2
        n_variables = (len(self.meshings)-n_meshing_planetary)*2+n_meshing_planetary*3 + 1

        system_matrix = npy.zeros((n_equations, n_variables))
        rhs = npy.zeros(n_equations)
        num_element = 0
        num_equation = 0
        element_association = {}

        for element in self.elements:
           element_association[element] = []

        for meshing in self.meshings:
            matrix_relation, rhs_relation = meshing.torque_system_equations()

            if isinstance(meshing, MeshingPlanetary):

                system_matrix[num_equation][num_element] = matrix_relation[0][0]
                system_matrix[num_equation][num_element+1] = matrix_relation[0][1]
                system_matrix[num_equation+1][num_element+2] = matrix_relation[1][0]
                system_matrix[num_equation+1][num_element+1] = matrix_relation[1][1]

                rhs[num_equation] = rhs_relation[0]
                rhs[num_equation+1] = rhs_relation[1]

                if isinstance(meshing.nodes[0], Planetary):

                    element_association[meshing.nodes[0]].append(num_element)
                    element_association[meshing.nodes[1]].append(num_element+1)

                else:
                    element_association[meshing.nodes[1]].append(num_element)
                    element_association[meshing.nodes[0]].append(num_element+1)

                element_association[self.planet_carrier].append(num_element+2)

                num_element += 3
                num_equation += 2

            else:

                system_matrix[num_equation][num_element] = matrix_relation[0]
                system_matrix[num_equation][num_element+1] = matrix_relation[1]

                rhs[num_equation] = rhs_relation[0]

                element_association[meshing.nodes[0]].append(num_element)
                element_association[meshing.nodes[1]].append(num_element+1)


                num_element += 2
                num_equation += 1

        element_without_doubles = copy.copy(self.elements)

        for double in self.doubles:
            matrix_association_element = npy.zeros((1, system_matrix.shape[1]))

            for association in element_association[double.nodes[0]]:
                    matrix_association_element[0][association] = 1


            for association in element_association[double.nodes[1]]:
                    matrix_association_element[0][association] = 1


            system_matrix = npy.concatenate((system_matrix, matrix_association_element))
            rhs = npy.concatenate((rhs, [0]))

            element_without_doubles.remove(double.nodes[0])
            element_without_doubles.remove(double.nodes[1])



        for element in element_without_doubles:

            if len(element_association[element]) > 1:

                matrix_association_element = npy.zeros((1, system_matrix.shape[1]))

                for association in element_association[element]:
                    matrix_association_element[0][association] = -1

                if element == self.planet_carrier:
                    matrix_association_element[0][-1] = 1

                system_matrix = npy.concatenate((system_matrix, matrix_association_element))
                rhs = npy.concatenate((rhs, [0]))

        return system_matrix, rhs, element_association


    def torque_resolution_PFS(self, input_torque_and_composant, meshing_chains=[]):
        Ca = 0
        Cr = 0
        Cf = 0
        Cwb = 0# Speed coeff for bearings
        Cvgs = 0# Speed coeff for gear sets
        self.mech_dict = {}

        alpha_gs1 = 15/360*2*3.1415
        beta_gs1 = 0


        egs1 = gm_geo.Direction2Euler((0, 0, 1))

        part_planets = []
        part_planetaries = []
        part_planet_carrier = genmechanics.Part('planet_carrier')
        self.mech_dict['part_planet_carrier'] = part_planet_carrier
        ground = genmechanics.Part('ground')
        planet_carrier_a = linkages.FrictionlessBallLinkage(ground, part_planet_carrier, [0, 0, 0], [0, 0, 0], 'planet_carrier_a')
        planet_carrier_b = linkages.FrictionlessLinearAnnularLinkage(ground, part_planet_carrier, [0.1, 0, 0], [0, 0, 0], 'planet_carrier_b')
        link_planetaries_ball = []
        link_planetaries_linear_annular = []
        pivot_planets = []
        flag_double = 0
        for i, planet in enumerate(self.planets):

            for double in self.doubles:
                if planet in double.nodes:
                    pivot_planets.append(0)
                    part_planets.append(0)
                    flag_double = 1
                    break
            if not flag_double:
                part_planets.append(genmechanics.Part('planet'+str(i)))

                position = [planet.positions[0][0], planet.positions[0][1], planet.positions[0][2]]
                pivot_planets.append(linkages.FrictionlessRevoluteLinkage(part_planet_carrier, part_planets[-1],
                                                                          np.array(position), [0, 0, 0], 'pivot'+str(i)))

            flag_double = 0

        self.mech_dict['part_planets'] = part_planets
        previous_nodes = []

        for i, double in  enumerate(self.doubles):
            if double.nodes[0] in previous_nodes:
                planet_double = part_planets[self.planets.index(double.nodes[0])]
                link = pivot_planets[self.planets.index(double.nodes[0])]
                part_planets[self.planets.index(double.nodes[1])] = planet_double
                pivot_planets[self.planets.index(double.nodes[1])] = link


            elif double.nodes[1] in previous_nodes:
                planet_double = part_planets[self.planets.index(double.nodes[1])]
                link = pivot_planets[self.planets.index(double.nodes[1])]
                part_planets[self.planets.index(double.nodes[0])] = planet_double
                pivot_planets[self.planets.index(double.nodes[0])] = link



            else:
                 planet_double = genmechanics.Part('planet_double'+str(i))
                 position = [planet.positions[0][0], planet.positions[0][1], planet.positions[0][2]]
                 link = linkages.FrictionlessRevoluteLinkage(part_planet_carrier, planet_double,
                                                             np.array(position), [0, 0, 0], 'pivot_double'+str(i))

                 part_planets[self.planets.index(double.nodes[0])] = planet_double
                 part_planets[self.planets.index(double.nodes[1])] = planet_double
                 pivot_planets[self.planets.index(double.nodes[0])] = link
                 pivot_planets[self.planets.index(double.nodes[1])] = link

            previous_nodes.extend(double.nodes)




        for i, planetary in enumerate(self.planetaries):

            part_planetaries.append(genmechanics.Part('planetary'+str(i)))
            position = [planetary.position[0], planetary.position[1], planetary.position[2]]
            link_planetaries_ball.append(linkages.FrictionlessBallLinkage(ground, part_planetaries[-1], np.array(position),
                                                                          [0, 0, 0], 'planetary_ball'+str(i)))

            link_planetaries_linear_annular.append(linkages.FrictionlessLinearAnnularLinkage(ground, part_planetaries[-1],
                                                                                             np.array([position[0]+0.1, position[1], position[2]]),
                                                                                             [0, 0, 0], 'planetary_linear_angular'+str(i)))
        self.mech_dict['part_planetaries'] = part_planetaries
        self.mech_dict['pivot_planets'] = pivot_planets





        if not meshing_chains:
            meshing_chains = self.meshing_chain()

        gearings = []
        gearings_planetary = []
        position_gearings = []
        for j, meshing_chain in enumerate(meshing_chains):
            previous_element = meshing_chain[0]
            if isinstance(previous_element, Planet):
                new_position = [previous_element.positions[0][0], previous_element.positions[0][1], previous_element.positions[0][2]]
                previous_position = np.array(new_position)
                previous_part = part_planets[self.planets.index(previous_element)]
            else:
                new_position = [previous_element.position[0], previous_element.position[1], previous_element.position[2]]
                previous_position = np.array(new_position)
                previous_part = part_planetaries[self.planetaries.index(previous_element)]


            for i, element in enumerate(meshing_chain):
                if i > 0:
                    if isinstance(element, Planetary):
                        new_position2 = [element.position[0], element.position[1], element.position[2]]
                        position = np.array(new_position2)
                        part = part_planetaries[self.planetaries.index(element)]

                        if element.planetary_type == 'Ring':
                            orientation_gearing = previous_position-position
                            if orientation_gearing[1] == 0:
                                if orientation_gearing[2] == 0:
                                    angular = 0
                                else:
                                    angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                                signe = 1
                            else:
                                angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                                signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])



                        else:
                            orientation_gearing = (position-previous_position)
                            angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                            signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])






                    else:
                        new_position2 = [element.positions[0][0], element.positions[0][1], element.positions[0][2]]
                        position = np.array(new_position2)
                        part = part_planets[self.planets.index(element)]

                        if isinstance(previous_element, Planetary) and previous_element.planetary_type == 'Ring':
                            orientation_gearing = position-previous_position
                            angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                            signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([position[0],
                                                         position[1]+element.Z*element.module*m.cos(angular)*0.5*signe,
                                                         position[2]+element.Z*element.module*m.sin(angular)*0.5*signe])



                        else:

                            orientation_gearing = (position-previous_position)
                            if orientation_gearing[1] == 0:
                                if orientation_gearing[2] == 0:
                                    angular = 0
                                else:
                                    angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                                signe = 1
                            else:
                                angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                                signe = (orientation_gearing[1]/abs(orientation_gearing[1]))

                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])

                    position_gearings.append(position_gearing)

                    gearings.append(linkages.FrictionLessGearSetLinkage(previous_part, part, position_gearing, egs1, alpha_gs1, beta_gs1,
                                                                        'Gear set '+str(i) + str(j)))

                    if isinstance(previous_element, Planetary) or isinstance(element, Planetary):
                        gearings_planetary.append(gearings[-1])

                    previous_element = element

                    if isinstance(previous_element, Planet):
                        new_position = [previous_element.positions[0][0], previous_element.positions[0][1], previous_element.positions[0][2]]
                        previous_position = np.array(new_position)
                        previous_part = part_planets[self.planets.index(previous_element)]
                    else:
                        new_position = [previous_element.position[0], previous_element.position[1], previous_element.position[2]]
                        previous_position = np.array(new_position)
                        previous_part = part_planetaries[self.planetaries.index(previous_element)]



        # fig, ax = plt.subplots()
        # ax.set_xlim(-100, 100)
        # ax.set_ylim(-100, 100)
        # new_position_gearing=[]
        # for position in position_gearings:
        #     new_position_gearing.append([position[0],position[1],position[2]])
        # vmp.plot(self.plot_data(new_position_gearing),ax)
        list_all_input = self.planetaries+[self.planet_carrier]
        loads_known = []
        load_planet_carrier = 0
        for i, input_composant in enumerate(input_torque_and_composant):
            list_all_input.remove(input_composant)
            if input_composant == self.planet_carrier:
                loads_known.append(loads.KnownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [0, 0, 0],
                                                   [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
                load_planet_carrier = loads_known[-1]

            else:

                loads_known.append(loads.KnownLoad(part_planetaries[self.planetaries.index(input_composant)], [0, 0, 0], [0, 0, 0],
                                                   [0, 0, 0], [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
        loads_unknown = []
        flag_load_planet_carrier_unknow = 0
        for i, unknow_input in enumerate(list_all_input):
            if unknow_input == self.planet_carrier:
                loads_unknown.append(loads.SimpleUnknownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))
                load_planet_carrier = loads_unknown[-1]
                flag_load_planet_carrier_unknow = 1
            else:

                loads_unknown.append(loads.SimpleUnknownLoad(part_planetaries[self.planetaries.index(unknow_input)],
                                                             [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))

        pivot_planets_without_double = []
        for pivot in pivot_planets:
            if not pivot in pivot_planets_without_double:
                 pivot_planets_without_double.append(pivot)

        self.mech_dict['pivot_planets_without_double'] = pivot_planets_without_double
        self.mech_dict['gearings_planetary'] = gearings_planetary
        self.mech_dict['part_planetaries'] = part_planetaries
        self.mech_dict['load_planet_carrier'] = load_planet_carrier
        self.mech_dict['link_planetaries_ball'] = link_planetaries_ball
        self.mech_dict['planet_carrier_a'] = planet_carrier_a
        self.mech_dict['input_torque_and_composant'] = input_torque_and_composant
        self.mech_dict['flag_load_planet_carrier_unknow'] = flag_load_planet_carrier_unknow


        list_parts = gearings+pivot_planets_without_double+link_planetaries_ball+ link_planetaries_linear_annular+ [planet_carrier_a] + [planet_carrier_b]

        for i, planet in enumerate(self.planets):
            self.mech_dict[planet] = pivot_planets[i]
        self.mech_dict['gearing_end'] = list_parts.index(gearings[-1])





        imposed_speeds = [(link_planetaries_ball[0], 0, self.planetaries[0].speed_input[0]), (link_planetaries_ball[-1], 0, self.planetaries[-1].speed_input[0])]


        mech = genmechanics.Mechanism(list_parts, ground, imposed_speeds, loads_known, loads_unknown)

        self.mech = mech
        # for l, lv in mech.static_results.items():

        #     for d, v in lv.items():
        #         print(l.name, d, v)
        # torque_max_planet = 0

        # for i,gearing_planetary in enumerate(gearings_planetary):

        #     power = mech.TransmittedLinkagePower(gearing_planetary, 4)

        #     if gearing_planetary.part1 in part_planetaries:
        #         index=part_planetaries.index(gearing_planetary.part1)
        #     else:
        #         index=part_planetaries.index(gearing_planetary.part2)
        #     speed = mech.kinematic_results[link_planetaries_ball[index]][0]
        #     print(speed)
        #     print(power)
        #     self.planetaries[index].torque_signe = (power/speed)/abs(power/speed)
        #     self.planetaries[index].power = power

        # if flag_load_planet_carrier_unknow:

        #     torque = mech.static_results[load_planet_carrier][0]
        # else:
        #     torque = input_torque_and_composant[self.planet_carrier]
        # # print(mech.kinematic_results)
        # speed = mech.kinematic_results[planet_carrier_a][0]
        # self.planet_carrier.torque_signe=torque/abs(torque)
        # self.planet_carrier.power=torque*speed

        # for pivot in pivot_planets_without_double:
        #     power = mech.TransmittedLinkagePower(pivot, 1)
        #     speed = mech.kinematic_results[pivot][0]
        #     torque = power/speed

        #     # print(torque)
        #     if abs(torque) > torque_max_planet:
        #         torque_max_planet = abs(torque)


        # mech.GlobalSankey()

    def linkages(self, listshaftplanetary, shaft_planet_carrier, center):
        Ca = 0
        Cr = 0
        Cf = 0
        Cwb = 0# Speed coeff for bearings
        Cvgs = 0# Speed coeff for gear sets
        self.mech_dict = {}

        alpha_gs1 = 15/360*2*3.1415
        beta_gs1 = 0
        # egs1= gm_geo.Direction2Euler((0, 0, 1))


        part_planets = []
        part_planetaries = []
        part_planet_carrier = shaft_planet_carrier.part

        ground = genmechanics.Part('ground')
        # planet_carrier_a = linkages.FrictionlessBallLinkage(ground, part_planet_carrier, center, [0, 0, 0], 'planet_carrier_a')
        # planet_carrier_b = linkages.FrictionlessLinearAnnularLinkage(ground, part_planet_carrier, [center[0]+0.1, center[1], center[2]], [0, 0, 0], 'planet_carrier_b')
        link_planetaries_ball = []
        link_planetaries_linear_annular = []
        pivot_planets = []
        flag_double = 0
        for i, planet in enumerate(self.planets):

            for double in self.doubles:
                if planet in double.nodes:
                    pivot_planets.append(0)
                    part_planets.append(0)
                    flag_double = 1
                    break
            if not flag_double:
                part_planets.append(genmechanics.Part('planet'+str(i)))

                position = [planet.positions[0][0]+center[0], planet.positions[0][1]+center[1], planet.positions[0][2]+center[2]]
                pivot_planets.append(linkages.FrictionlessRevoluteLinkage(part_planet_carrier, part_planets[-1],
                                                                          np.array(position), [0, 0, 0], 'pivot'+str(i)))

            flag_double = 0


        previous_nodes = []

        for i, double in  enumerate(self.doubles):
            if double.nodes[0] in previous_nodes:
                planet_double = part_planets[self.planets.index(double.nodes[0])]
                link = pivot_planets[self.planets.index(double.nodes[0])]
                part_planets[self.planets.index(double.nodes[1])] = planet_double
                pivot_planets[self.planets.index(double.nodes[1])] = link


            elif double.nodes[1] in previous_nodes:
                planet_double = part_planets[self.planets.index(double.nodes[1])]
                link = pivot_planets[self.planets.index(double.nodes[1])]
                part_planets[self.planets.index(double.nodes[0])] = planet_double
                pivot_planets[self.planets.index(double.nodes[0])] = link



            else:
                 planet_double = genmechanics.Part('planet_double'+str(i))
                 position = [planet.positions[0][0]+center[0], planet.positions[0][1]+center[1], planet.positions[0][2]+center[2]]
                 link = linkages.FrictionlessRevoluteLinkage(part_planet_carrier, planet_double,
                                                             np.array(position), [0, 0, 0], 'pivot_double'+str(i))

                 part_planets[self.planets.index(double.nodes[0])] = planet_double
                 part_planets[self.planets.index(double.nodes[1])] = planet_double
                 pivot_planets[self.planets.index(double.nodes[0])] = link
                 pivot_planets[self.planets.index(double.nodes[1])] = link

            previous_nodes.extend(double.nodes)




        for i, planetary in enumerate(self.planetaries):

            part_planetaries.append(listshaftplanetary[i].part)
            # position=[planetary.position[0]+center[0],planetary.position[1]+center[1],planetary.position[2]+center[2]]
            # link_planetaries_ball.append(linkages.FrictionlessBallLinkage(ground, part_planetaries[-1], np.array(position),
            #                                                   [0, 0, 0], 'planetary_ball'+str(i)))

            # link_planetaries_linear_annular.append(linkages.FrictionlessLinearAnnularLinkage(ground, part_planetaries[-1],
            #                                                                      np.array([position[0]+0.1, position[1], position[2]]),
            #                                                                      [0, 0, 0], 'planetary_linear_angular'+str(i)))






        meshing_chains = self.meshing_chain()

        gearings = []
        gearings_planetary = []
        position_gearings = []
        for j, meshing_chain in enumerate(meshing_chains):
            previous_element = meshing_chain[0]
            if isinstance(previous_element, Planet):
                new_position = [previous_element.positions[0][0]+center[0], previous_element.positions[0][1]+center[1], previous_element.positions[0][2]+center[2]]
                previous_position = np.array(new_position)
                previous_part = part_planets[self.planets.index(previous_element)]
            else:
                new_position = [previous_element.position[0]+center[0], previous_element.position[1]+center[1], previous_element.position[2]+center[2]]
                previous_position = np.array(new_position)
                previous_part = part_planetaries[self.planetaries.index(previous_element)]


            for i, element in enumerate(meshing_chain):
                if i > 0:
                    if isinstance(element, Planetary):
                        new_position2 = [element.position[0]+center[0], element.position[1]+center[1], element.position[2]+center[2]]
                        position = np.array(new_position2)
                        part = part_planetaries[self.planetaries.index(element)]

                        if element.planetary_type == 'Ring':
                            orientation_gearing = previous_position-position
                            if orientation_gearing[1] == 0:
                                if orientation_gearing[2] == 0:
                                    angular = 0
                                else:
                                    angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                                signe = 1
                            else:
                                angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                                signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])



                        else:
                            orientation_gearing = (position-previous_position)
                            angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                            signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])






                    else:
                        new_position2 = [element.positions[0][0]+center[0], element.positions[0][1]+center[1], element.positions[0][2]+center[2]]
                        position = np.array(new_position2)
                        part = part_planets[self.planets.index(element)]

                        if isinstance(previous_element, Planetary) and previous_element.planetary_type == 'Ring':
                            orientation_gearing = position-previous_position
                            angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                            signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                            position_gearing = np.array([position[0],
                                                         position[1]+element.Z*element.module*m.cos(angular)*0.5*signe,
                                                         position[2]+element.Z*element.module*m.sin(angular)*0.5*signe])



                        else:

                            orientation_gearing = (position-previous_position)
                            if orientation_gearing[1] == 0:
                                if orientation_gearing[2] == 0:
                                    angular = 0
                                else:
                                    angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                                signe = 1
                            else:
                                angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                                signe = (orientation_gearing[1]/abs(orientation_gearing[1]))

                            position_gearing = np.array([previous_position[0],
                                                         previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                         previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])


                    position_gearings.append(position_gearing)

                    dgs = npy.cross(npy.array([0, position[1]-previous_position[1], position[2]-previous_position[2]]), (1, 0, 0))
                    egs1 = gm_geo.Direction2Euler(dgs, (1, 0, 0))

                    gearings.append(linkages.FrictionLessGearSetLinkage(previous_part, part, position_gearing, egs1, alpha_gs1, beta_gs1,
                                                                        'Gear set '+str(i) + str(j)))

                    if isinstance(previous_element, Planetary) or isinstance(element, Planetary):
                        gearings_planetary.append(gearings[-1])

                    previous_element = element

                    if isinstance(previous_element, Planet):
                        new_position = [previous_element.positions[0][0]+center[0], previous_element.positions[0][1]+center[1], previous_element.positions[0][2]+center[2]]
                        previous_position = np.array(new_position)
                        previous_part = part_planets[self.planets.index(previous_element)]
                    else:
                        new_position = [previous_element.position[0]+center[0], previous_element.position[1]+center[1], previous_element.position[2]+center[2]]
                        previous_position = np.array(new_position)
                        previous_part = part_planetaries[self.planetaries.index(previous_element)]



        # fig, ax = plt.subplots()
        # ax.set_xlim(-100, 100)
        # ax.set_ylim(-100, 100)
        # new_position_gearing=[]
        # for position in position_gearings:
        #     new_position_gearing.append([position[0],position[1],position[2]])
        # vmp.plot(self.plot_data(new_position_gearing),ax)
        # list_all_input = self.planetaries+[self.planet_carrier]
        # loads_known = []
        # load_planet_carrier = 0
        # for i, input_composant in enumerate(input_torque_and_composant):
        #     list_all_input.remove(input_composant)
        #     if input_composant == self.planet_carrier:
        #         loads_known.append(loads.KnownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [0, 0, 0],
        #                                            [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
        #         load_planet_carrier=loads_known[-1]

        #     else:

        #         loads_known.append(loads.KnownLoad(part_planetaries[self.planetaries.index(input_composant)], [0, 0, 0], [0, 0, 0],
        #                                            [0, 0, 0], [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
        # loads_unknown = []
        # flag_load_planet_carrier_unknow=0
        # for i, unknow_input in enumerate(list_all_input):
        #     if unknow_input == self.planet_carrier:
        #         loads_unknown.append(loads.SimpleUnknownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))
        #         load_planet_carrier=loads_unknown[-1]
        #         flag_load_planet_carrier_unknow=1
        #     else:

        #         loads_unknown.append(loads.SimpleUnknownLoad(part_planetaries[self.planetaries.index(unknow_input)],
        #                                                      [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))

        pivot_planets_without_double = []
        for pivot in pivot_planets:
            if not pivot in pivot_planets_without_double:
                 pivot_planets_without_double.append(pivot)




        list_parts = gearings+pivot_planets_without_double

        return list_parts


    def power_graph(self):
        labels = {}
        widths = []
        edges = []
        G = self.mech.DrawPowerGraph()
        plt.figure()
        for part in self.mech.parts:


            labels[part] = part.name

        for linkage in self.mech.linkages:

            labels[linkage] = linkage.name

            widths.append(abs(self.mech.TransmittedLinkagePower(linkage, 0)))

            edges.append((linkage, linkage.part1))

            widths.append(abs(self.mech.TransmittedLinkagePower(linkage, 1)))

            edges.append((linkage, linkage.part2))

        for load in self.mech.unknown_static_loads+self.mech.known_static_loads:

            widths.append(abs(self.mech.LoadPower(load)))
            labels[load] = load.name
            edges.append((load, load.part))

        max_widths = max(widths)
        widths = [6*w/max_widths for w in widths]

        # pos=nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos, nodelist=self.mech.linkages, node_color='grey')
        nx.draw_networkx_nodes(G, pos, nodelist=self.mech.parts)
        nx.draw_networkx_nodes(G, pos, nodelist=self.mech.unknown_static_loads, node_color='red')
        nx.draw_networkx_nodes(G, pos, nodelist=self.mech.known_static_loads, node_color='green')
        nx.draw_networkx_nodes(G, pos, nodelist=self.mech.parts, node_color='cyan')
        nx.draw_networkx_labels(G, pos, labels)
        nx.draw_networkx_edges(G, pos, edges, width=widths, edge_color='blue')
        nx.draw_networkx_edges(G, pos)
        G2 = G.copy()
        G2.remove_node(self.mech.ground)
        G3 = nx.find_cycle(G2, orientation='ignore')
        circles = [G3]



    def update_position_mech(self):
        if not self.mech:
            self.torque_max_planets()

        for planet in self.planets:


            self.mech.linkages[self.mech.linkages.index(self.mech_dict[planet])].position = npy.array(planet.positions[0])

        for i in range(self.mech_dict['gearing_end']):
            gearing = self.mech.linkages[i]
            part1 = gearing.part1
            part2 = gearing.part2

            if part1 in self.mech_dict['part_planetaries']:
                previous_element = self.planetaries[self.mech_dict['part_planetaries'].index(part1)]
                previous_position = previous_element.position
            else:
                previous_element = self.planets[self.mech_dict['part_planets'].index(part1)]
                previous_position = previous_element.positions[0]



            if part2 in self.mech_dict['part_planetaries']:


                element = self.planetaries[self.mech_dict['part_planetaries'].index(part2)]
                new_position2 = element.position
                position = np.array(new_position2)
                if element.planetary_type == 'Ring':
                    orientation_gearing = previous_position-position
                    if orientation_gearing[1] == 0:
                            if orientation_gearing[2] == 0:
                                angular = 0
                            else:
                                angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                            signe = 1
                    else:
                        angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                        signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                    position_gearing = np.array([previous_position[0],
                                                 previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                 previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])



                else:
                    orientation_gearing = (position-previous_position)
                    angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                    signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                    position_gearing = np.array([previous_position[0],
                                                 previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                 previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])






            else:


                element = self.planets[self.mech_dict['part_planets'].index(part2)]
                new_position2 = element.positions[0]
                position = np.array(new_position2)
                if isinstance(previous_element, Planetary) and previous_element.planetary_type == 'Ring':
                    orientation_gearing = position-previous_position
                    angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                    signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                    position_gearing = np.array([position[0],
                                                 position[1]+element.Z*element.module*m.cos(angular)*0.5*signe,
                                                 position[2]+element.Z*element.module*m.sin(angular)*0.5*signe])



                else:

                    orientation_gearing = (position-previous_position)
                    if orientation_gearing[1] == 0:
                            if orientation_gearing[2] == 0:
                                angular = 0
                            else:
                                angular = (orientation_gearing[2]/abs(orientation_gearing[2]))*m.pi/2
                            signe = 1
                    else:
                        angular = m.atan(orientation_gearing[2]/orientation_gearing[1])
                        signe = (orientation_gearing[1]/abs(orientation_gearing[1]))
                    position_gearing = np.array([previous_position[0],
                                                 previous_position[1]+previous_element.Z*previous_element.module*m.cos(angular)*0.5*signe,
                                                 previous_position[2]+previous_element.Z*previous_element.module*m.sin(angular)*0.5*signe])

            gearing.position = position_gearing



    def update_load_mech(self, input_torque_and_composant):
        list_all_input = self.planetaries+[self.planet_carrier]
        loads_known = []
        load_planet_carrier = 0
        part_planet_carrier = self.mech_dict['part_planet_carrier']
        part_planetaries = self.mech_dict['part_planetaries']
        for i, input_composant in enumerate(input_torque_and_composant):
            list_all_input.remove(input_composant)
            if input_composant == self.planet_carrier:
                loads_known.append(loads.KnownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [0, 0, 0],
                                                   [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
                load_planet_carrier = loads_known[-1]

            else:

                loads_known.append(loads.KnownLoad(part_planetaries[self.planetaries.index(input_composant)], [0, 0, 0], [0, 0, 0],
                                                   [0, 0, 0], [input_torque_and_composant[input_composant], 0, 0], 'input torque'+str(i)))
        loads_unknown = []
        flag_load_planet_carrier_unknow = 0
        for i, unknow_input in enumerate(list_all_input):
            if unknow_input == self.planet_carrier:
                loads_unknown.append(loads.SimpleUnknownLoad(part_planet_carrier, [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))
                load_planet_carrier = loads_unknown[-1]
                flag_load_planet_carrier_unknow = 1
            else:

                loads_unknown.append(loads.SimpleUnknownLoad(part_planetaries[self.planetaries.index(unknow_input)],
                                                             [0, 0, 0], [0, 0, 0], [], [0], 'output torque'+str(i)))

        self.mech_dict['load_planet_carrier'] = load_planet_carrier
        self.mech_dict['flag_load_planet_carrier_unknow'] = flag_load_planet_carrier_unknow
        self.mech_dict['input_torque_and_composant'] = input_torque_and_composant
        self.mech.ChangeLoads(loads_known, loads_unknown)
        self.mech.static_results
















    def torque_solve(self, input_torque_and_composant):
        '''
        A function which give the torque of all the elements(Planetary, Planet, PlanetCarrier) of planetary gear
        whith n-2 input elements and torques (whith n= number of planetary + planet carrier )

        :param input_torque_and_composant: A dictionary where the element input are associated with their torque input
        :type input_torque_and_composant: Dictionary{ Planetary, Planet, PlanetCarrier : float}

        :return:  A dictionary where all the element are associated with their torque calculated
        :rtype: Dictionary{ Planetary, Planet, PlanetCarrier : float}


        '''
        system_matrix, vector_b, element_association = self.torque_system_equation()

        for composant in input_torque_and_composant:
            matrix_input = npy.zeros((1, system_matrix.shape[1]))
            if isinstance(composant, PlanetCarrier):
                matrix_input[0][-1] = 1
            else:
                matrix_input[0][element_association[composant]] = 1
            if isinstance(input_torque_and_composant[composant], Planetary):
                 matrix_input[0][element_association[input_torque_and_composant[composant]]] = 1
                 vector_b = npy.concatenate((vector_b, [0]))
            elif isinstance(input_torque_and_composant[composant], PlanetCarrier):
                matrix_input[0][-1] = 1
                vector_b = npy.concatenate((vector_b, [0]))
            else:
                vector_b = npy.concatenate((vector_b, [input_torque_and_composant[composant]]))

            system_matrix = npy.concatenate((system_matrix, matrix_input))


        solution = solve(system_matrix, vector_b)

        torque_element_association = {}

        for element in element_association:
            max_solution = 0
            for position in element_association[element]:
                if abs(solution[position]) > abs(max_solution):

                    max_solution = solution[position]
            torque_element_association[element] = max_solution

        torque_element_association[self.planet_carrier] = solution[-1]
        return torque_element_association

    def recirculation_power_recursive_function(self, l, node, circle, circles):
        flag_node = 0
        flag_append_circle = 0
        for nodes in l:
            if nodes[0] == node:
                flag_node = 1
                circle2 = copy.copy(circle)

                if nodes[1] in circle2:
                    index = circle2.index(nodes[1])

                    if not circle2[index:] in circles:
                        flag_similiti = 0
                        for circle in circles:

                            if circle2[index:][0] in circle:
                                len_simili = 0

                                for node2 in circle2[index:]:
                                    if node2 in circle:
                                        len_simili += 1

                                if len_similiti == len(circle):
                                    flag_similiti = 1

                        if not flag_similiti:
                            flag_append_circle = 1
                            circles.append(circle2[index:])


                else:
                    circle2.append(nodes[1])

                    # l2=[]
                    # for li in circle2:

                    #     l2.append(li.name)

                    self.recirculation_power_recursive_function(l, nodes[1], circle2, circles)
        if not flag_node and not flag_append_circle:
            for node2 in circle:
                for nodes in l:
                    if node2 == nodes[1] and nodes[0] not in circle:
                        circle2 = copy.copy(circle)

                        index = circle2.index(nodes[1])

                        circle2 = circle2[:index+1]
                        circle2.append(nodes[0])
                        # print(node2.name)

                        self.recirculation_power_recursive_function_inverse(l, nodes[0], circle2, circles)
        return circles

    def recirculation_power_recursive_function_inverse(self, l, node, circle, circles):
        flag_node = 0
        flag_append_circle = 0
        for nodes in l:

            if nodes[1] == node:
                flag_node = 1
                circle2 = copy.copy(circle)

                if nodes[0] in circle2:
                    index = circle2.index(nodes[0])
                    if not circle2[index:] in circles:

                        flag_similiti = 0
                        for circle in circles:

                            if circle2[index:][0] in circle:
                                len_simili = 0
                                for node2 in circle2[index:]:
                                    if node2 in circle:
                                        len_simili += 1
                                if len_simili == len(circle):
                                    flag_similiti = 1

                        if not flag_similiti:
                            flag_append_circle = 1
                            circles.append(circle2[index:])

                else:
                    circle2.append(nodes[0])

                    # l2=[]
                    # for li in circle2:

                    #     l2.append(li.name)

                    self.recirculation_power_recursive_function_inverse(l, nodes[0], circle2, circles)
        # if not flag_node and not flag_append_circle:
        #     for node2 in circle:
        #         for nodes in l:
        #             if node2==nodes[0] and nodes[1] not in circle:
        #                 circle2=copy.copy(circle)

        #                 index=circle2.index(nodes[0])

        #                 circle2=circle2[:index+1]
        #                 circle2.append(nodes[1])
        #                 # print(node2.name)
        #                 if len(circle2)>3:
        #                     return circles
        #                 self.recirculation_power_recursive_function(l,nodes[1],circle2,circles)

        return circles




    def recirculation_power(self):
        first_input = {}
        # if not self.mech:
        #     # for i in range(len(self.planetaries)-1):
        #     #    first_input[self.planetaries[i]]=self.planetaries[i].torque_input[0]/self.number_branch_planet
        #     # self.torque_resolution_PFS(first_input)
        try:
            self.torque_max_planets()
            G = self.mech.DrawPowerGraph(return_graph=True)
        except genmechanics.ModelError:
            return [[0, 1000000]]

        G2 = G.copy()
        G2.remove_node(self.mech.ground)

        # l=list(nx.edge_bfs(G2,list(G2.nodes)[0]))
        # l_name=[]
        # for node in l:
        #     l_name.append([node[0].name,node[1].name])


        circles = []
        for node in list(G2.nodes):
            l = list(nx.edge_bfs(G2, node))
            self.recirculation_power_recursive_function(l, node, [node], circles)

        circle_name = []
        for circle in circles:
            circle_name = []
            for name in circle:
                circle_name.append(name.name)
            # print(circle_name)





        power_circles = []
        for circle in circles:
            power_list = []
            for node in circle:
                if isinstance(node, linkages.Linkage):

                    power_list.append(abs(self.mech.TransmittedLinkagePower(node, 0)))


            min_power = min(power_list)
            max_power = max(power_list)
            power_circles.append([min_power, max_power])


        power_input = []
        # for loads in self.mech.known_static_loads:

        #     torque = loads.torques[0]

        #     for linkage in self.mech.linkages:
        #         if (linkage.part1==loads.part and linkage.part2==self.mech.ground )or (linkage.part2==loads.part and linkage.part1==self.mech.ground):
        #             speed=self.mech.kinematic_results[linkage][0]
        #             print(speed)
        #             break



        #     if speed*torque>0:
        #         power_input+=speed*torque

        # for loads in self.mech.unknown_static_loads:

        #     torque = self.mech.static_results[loads][0]

        #     for linkage in self.mech.linkages:
        #         if (linkage.part1==loads.part and linkage.part2==self.mech.ground )or (linkage.part2==loads.part and linkage.part1==self.mech.ground):
        #             speed=self.mech.kinematic_results[linkage][0]

        #             break



        #     if speed*torque>0:
        #         power_input+=speed*torque

        for circle in circles:
            power_input.append(0)
            for node in circle:
                neighbors = list(G2.neighbors(node))
                for neighbor in neighbors:
                    if  not neighbor in circle:

                        if isinstance(neighbor, loads.KnownLoad)or isinstance(neighbor, loads.SimpleUnknownLoad):
                            power = self.mech.LoadPower(neighbor)

                        else:
                            if neighbor.part1 == node:
                                power = self.mech.TransmittedLinkagePower(neighbor, 0)
                            else:
                                power = self.mech.TransmittedLinkagePower(neighbor, 1)

                        if power > 0:
                            power_input[-1] += power





        power_difference = []
        for i, power in enumerate(power_circles):
            if power_input[i] != 0:
                power_difference.append([(power[0]/power_input[i])*100, (power[1]/power_input[i])*100])
            else:
                power_difference.append([(power[0]/0.00000001)*100, (power[1]/0.000000001)*100])


        return power_difference

    def update_length(self):
        self.mesh_generation()
        max_length_meshing_chain = self.max_length_meshing_chain
        max_length_meshing_chain2 = copy.copy(max_length_meshing_chain)

        length_max1 = max(max_length_meshing_chain)
        max_length_meshing_chain2.remove(length_max1)

        if max_length_meshing_chain:
            length_max2 = max(max_length_meshing_chain2)
            self.length_double = ((length_max1+length_max2)/2)*1.1
        else:
            length_max2 = 0
            self.length_double = 0

        meshing_chains = self.meshing_chain()

        z = self.meshing_chain_position_z(meshing_chains)
        z_max = 0
        index_min = 0
        z_min = np.inf
        index_max = 0
        for i, meshing_chain in enumerate(meshing_chains):

            if z[i] < z_min:
                z_min = z[i]
                index_min = i

            if z[i] > z_max:
                z_max = z[i]
                index_max = i

            for element in meshing_chain:
                if isinstance(element, Planetary):
                    element.position = (z[i], element.position[1], element.position[2])

                elif isinstance(element, Planet):
                    for j, position in enumerate(element.positions):
                        element.positions[j] = (z[i], position[1], position[2])


        self.length = z_max-z_min+(max_length_meshing_chain[index_min]+max_length_meshing_chain[index_max])/2


    def update_length_without_mesh_generation(self):
        meshing_chains = self.meshing_chain()

        z = self.meshing_chain_position_z(meshing_chains)
        z_max = 0

        z_min = np.inf

        for i, meshing_chain in enumerate(meshing_chains):

            if z[i] < z_min:
                z_min = z[i]
                i

            if z[i] > z_max:
                z_max = z[i]


            for element in meshing_chain:
                if isinstance(element, Planetary):
                    element.position = (z[i], element.position[1], element.position[2])

                elif isinstance(element, Planet):
                    for j, position in enumerate(element.positions):
                        element.positions[j] = (z[i], position[1], position[2])

    def ConvertCenterPlanetaryGears(self, z):
        meshing_chains = self.meshing_chain()
        z_meshing_chain = self.meshing_chain_position_z(meshing_chains)
        z_meshing_chain_2 = copy.copy(z_meshing_chain)
        z_meshing_chain_2.sort()
        z2 = z+(self.max_length_meshing_chain[z_meshing_chain.index(z_meshing_chain_2[0])]-self.length)/2
        for i, z in enumerate(z_meshing_chain_2):
            if z_meshing_chain.index(z) == 0:
                return z2
            else:
                z2 += self.length_double+(self.max_length_meshing_chain[z_meshing_chain.index(z)]+self.max_length_meshing_chain[z_meshing_chain.index(z_meshing_chain_2[i+1])])/2











class PositionMinMaxPlanetaryGear(DessiaObject):


     def __init__(self, planetary_gear: PlanetaryGear, name: str = '', positions_min_max: List[float] = '', modules_min_max: List[float] = ''):

         self.planetary_gear = planetary_gear
         self.positions_min_max = positions_min_max
         self.modules_min_max = modules_min_max

         if self.positions_min_max == '' and self.modules_min_max == '':
             self.positions_min_max = []
             self.modules_min_max = []
             element_list = planetary_gear.planets+ planetary_gear.planetaries
             for element in element_list:
                self.positions_min_max.append([0, 0])
                self.modules_min_max.append([0, 0])

         DessiaObject.__init__(self, name=planetary_gear.name+'PostionMinMax')



     def enter_position(self, position, planetary_gear, element, min_max):
        element_list = planetary_gear.planets+ planetary_gear.planetaries
        if min_max == 'Min':
            self.positions_min_max[element_list.index(element)][0] = position
        if min_max == 'Max':
            self.positions_min_max[element_list.index(element)][1] = position

     def enter_module(self, module, planetary_gear, element, min_max):
        element_list = planetary_gear.planets+ planetary_gear.planetaries
        if min_max == 'Min':
            self.modules_min_max[element_list.index(element)][0] = module
        if min_max == 'Max':
            self.modules_min_max[element_list.index(element)][1] = module

     def get_position(self, element, planetary_gear, min_max):
        element_list = planetary_gear.planets+ planetary_gear.planetaries

        if min_max == 'Min':
            return self.positions_min_max[element_list.index(element)][0]
        if min_max == 'Max':
            return self.positions_min_max[element_list.index(element)][1]

     def get_module(self, element, planetary_gear, min_max):
        element_list = planetary_gear.planets+ planetary_gear.planetaries
        if min_max == 'Min':
            return self.modules_min_max[element_list.index(element)][0]
        if min_max == 'Max':
            return self.modules_min_max[element_list.index(element)][1]





class PlanetaryGearResult(DessiaObject):
    _standalone_in_db = True

    _eq_is_data_eq = False
    _non_serializable_attributes = ['planetaries', 'planets', 'planet_carrier', 'connections', 'doubles']
    def __init__(self, planetary_gear: PlanetaryGear, position_min_max: PositionMinMaxPlanetaryGear, geometry_min_max: str = 'Min', recycle_power: int = 0):
        self.planetary_gear = planetary_gear
        self.geometry_min_max = geometry_min_max
        self.position_min_max = position_min_max
        self.planetaries = planetary_gear.planetaries
        self.planets = planetary_gear.planets
        self.planet_carrier = planetary_gear.planet_carrier
        self.connections = planetary_gear.connections
        self.doubles = planetary_gear.doubles
        self.recycle_power = recycle_power
        # self.update_geometry()

        # if not self.recycle_power:
        #     planetary_gear_recirculation_power = self.planetary_gear.recirculation_power()
        #     max_recirculation_branch = []
        #     for recirculation_branch in planetary_gear_recirculation_power:
        #         max_recirculation_branch.append(recirculation_branch[1])
        #     if max_recirculation_branch:
        #         self.recycle_power = max(max_recirculation_branch)

        if not self.planetary_gear.speed_max_planet:
            self.speed_max_planet = self.planetary_gear.speed_max_planets()
            self.planetary_gear.speed_max_planet = self.speed_max_planet

        else:
            self.speed_max_planet = self.planetary_gear.speed_max_planet

        # self.torque_max_planet=0

        self.D_train = self.planetary_gear.d_min
        # self.D_mini=
        self.sum_Z_planetary = self.planetary_gear.sum_Z_planetary
        # self.sum_speed_planetary = self.planetary_gear.sum_speed_planetary

        self.max_Z_planetary = self.planetary_gear.max_Z_planetary
        self.min_Z_planetary = self.planetary_gear.min_Z_planetary
        # self.speed_planet_carrer = self.planetary_gear.speed_planet_carrer

        DessiaObject.__init__(self, name=self.planetary_gear.name+'Result')
        self.update_geometry()
        for planetary in self.planetaries:
            d = planetary.module*planetary.Z
            if planetary.planetary_type == 'Ring':
                d = d*1.3
            if d > self.D_train:
                 self.D_train = d

        for planet in self.planets:
            if planet.positions:
                d = planet.module*planet.Z+ 2*((planet.positions[0][2])**2+(planet.positions[0][1])**2)**0.5

            if d > self.D_train:
                self.D_train = d

    @classmethod
    def dict_to_object(cls, dict_):
        planetary_gear = PlanetaryGear.dict_to_object(dict_['planetary_gear'])
        position_min_max = PositionMinMaxPlanetaryGear.dict_to_object(dict_['position_min_max'])
        recycle_power = dict_['recycle_power']
        obj = cls(planetary_gear=planetary_gear, position_min_max=position_min_max,
                  recycle_power=recycle_power)
        obj.planetary_gear.length = dict_['planetary_gear']['length']
        obj.planetary_gear.length_double = dict_['planetary_gear']['length_double']
        obj.planetary_gear.update_length_without_mesh_generation()
        return obj


    def __str__(self):

        Z_planets = {}

        for planet in self.planets:
            Z_planets[planet.name] = planet.Z

        Z_planetaries = {}
        number_ring = 0
        number_sun = 0

        for planetary in self.planetaries:
            Z_planetaries[planetary.name] = planetary.Z

            if planetary.planetary_type == 'Sun':
                number_sun += 1

            else:
                number_ring += 1
        connections_name = []
        for i in range(len(self.connections)):
            connections_name.append([self.connections[i].nodes[0].name, self.connections[i].nodes[1].name,
                                     self.connections[i].connection_type])

        return 'Name:' + self.name + '\n\n' + \
               'Planetary Number:' + str(len(self.planetaries)) + '\n' + \
               'Ring Number:'+ str(number_ring) + '\n' + \
               'Sun_Number:' + str(number_sun) + '\n' + \
               'Z_planetaries:' + str(Z_planetaries) + '\n\n' + \
               'Planets_Number:' + str(len(self.planets)) + '\n' + \
               'Planets_Double_Number:' + str(len(self.doubles)) + '\n' + \
               'Z_Planets:' + str(Z_planets) + '\n\n' + \
                str(connections_name) + '\n\n\n'



    def update_geometry(self):
        if self.geometry_min_max:
            for planet in self.planetary_gear.planets:
                planet.positions = self.position_min_max.get_position(planet, self.planetary_gear, self.geometry_min_max)
                planet.module = self.position_min_max.get_module(planet, self.planetary_gear, self.geometry_min_max)

            for planetary in self.planetary_gear.planetaries:
                planetary.position = self.position_min_max.get_position(planetary, self.planetary_gear, self.geometry_min_max)
                planetary.module = self.position_min_max.get_module(planetary, self.planetary_gear, self.geometry_min_max)
        self.planetary_gear.update_length_without_mesh_generation()



    def update_speed(self, input_speeds_and_composants):
        list_composant = []
        for composant in input_speeds_and_composants:
            composant.speed_input = input_speeds_and_composants[composant]
            list_composant.append(composant)

        composant_1 = list_composant[0]
        composant_2 = list_composant[1]
        speed_1 = self.planetary_gear.speed_solve({composant_1:composant_1.speed_input[0], composant_2:composant_2.speed_input[0]})
        speed_2 = self.planetary_gear.speed_solve({composant_1:composant_1.speed_input[0], composant_2:composant_2.speed_input[1]})
        speed_3 = self.planetary_gear.speed_solve({composant_1:composant_1.speed_input[1], composant_2:composant_2.speed_input[0]})
        speed_4 = self.planetary_gear.speed_solve({composant_1:composant_1.speed_input[1], composant_2:composant_2.speed_input[1]})

        for planetary in self.planetary_gear.planetaries:
            speed_max = 0
            speed_min = np.inf
            for speed in [speed_1, speed_2, speed_3, speed_4]:
                if speed[planetary] > speed_max:
                    speed_max = speed[planetary]
                if speed[planetary] < speed_min:
                    speed_min = speed[planetary]

            planetary.speed_input = [speed_min, speed_max]
        speed_max = 0

        for planet in self.planetary_gear.planets:
            for speed in [speed_1, speed_2, speed_3, speed_4]:
                if speed[planet] > speed_max:
                    speed_max = speed[planetary]


        self.speed_max_planet = speed_max


    def update_torque(self, input_torques_and_composants):
        element_list_2 = copy.copy(self.planetaries)
        element_list_3 = copy.copy(self.planetaries)

        for element in input_torques_and_composants:

            element.torque_input = input_torques_and_composants[element]
            if element in element_list_2:
                element_list_2.remove(element)


        element_1 = element_list_2[0]
        element_list_3.remove(element_1)

        if len(element_list_2) > 1:
            element_2 = element_list_2[1]
            element_list_3.remove(element_2)
        else:
            element_2 = self.planet_carrier

        element1_signe, element2_signe = self.planetary_gear.torque_min_max_signe_input(element_1, element_2)
        max_input_1 = {}
        min_input_1 = {}
        max_input_2 = {}
        min_input_2 = {}
        for i, element in enumerate(element_list_3):
            if element1_signe[i] == 1:
               max_input_1[element] = element.torque_input[1]
               min_input_1[element] = element.torque_input[0]

            else:

                max_input_1[element] = element.torque_input[0]
                min_input_1[element] = element.torque_input[1]

            if element2_signe[i] == 1:
               max_input_2[element] = element.torque_input[1]
               min_input_2[element] = element.torque_input[0]

            else:
                max_input_2[element] = element.torque_input[0]
                min_input_2[element] = element.torque_input[1]

        if len(element_list_3) == len(self.planetaries)-2:
            if element1_signe[-1] == 1:
               max_input_1[self.planet_carrier] = self.planet_carrier.torque_input[1]
               min_input_1[self.planet_carrier] = self.planet_carrier.torque_input[0]

            else:
                max_input_1[self.planet_carrier] = self.planet_carrier.torque_input[0]
                min_input_1[self.planet_carrier] = self.planet_carrier.torque_input[1]

            if element2_signe[i] == 1:
               max_input_2[self.planet_carrier] = self.planet_carrier.torque_input[1]
               min_input_2[self.planet_carrier] = self.planet_carrier.torque_input[0]

            else:
                max_input_2[self.planet_carrier] = self.planet_carrier.torque_input[0]
                min_input_2[self.planet_carrier] = self.planet_carrier.torque_input[1]



        res_min_1 = self.planetary_gear.torque_solve(min_input_1)

        res_max_1 = self.planetary_gear.torque_solve(max_input_1)

        res_min_2 = self.planetary_gear.torque_solve(min_input_2)

        res_max_2 = self.planetary_gear.torque_solve(max_input_2)

        element_1.torque_input = [res_min_1[element_1], res_max_1[element_1]]


        element_2.torque_input = [res_min_2[element_2], res_max_2[element_2]]

        self.update_torque_max()



    def speed_range(self, element_1, element_2):

        return self.planetary_gear.speed_range(element_1, element_2)

    def torque_range(self, elements):

        return self.planetary_gear.torque_range(elements)

    def update_d_train(self):
        self.D_train = 0
        for planetary in self.planetaries:
            d = planetary.module*planetary.Z

            if planetary.planetary_type == 'Ring':
                d = d*1.3
            if d > self.D_train:
                 self.D_train = d

        for planet in self.planets:
            if planet.positions:
                d = planet.module*planet.Z+ 2*((planet.positions[0][1])**2+(planet.positions[0][2])**2)**0.5

            if d > self.D_train:
                self.D_train = d




    def update_torque_max(self):
          self.update_geometry()

          self.torque_max_planet = self.planetary_gear.torque_max_planets()



    def volmdlr_primitives(self, frame=vm.OXYZ):

        self.update_geometry()
        self.planetary_gear.update_length()
        li_box = self.planetary_gear.volmdlr_primitives()
        return li_box

    def plot_data(self):
        self.update_geometry()

        plot_data = self.planetary_gear.plot_data()
        return plot_data


    def Volume(self, primitive_volmdlr):
        z = primitive_volmdlr.x.cross(primitive_volmdlr.y)
        z.normalize()
        coeff = primitive_volmdlr.extrusion_vector.dot(z)

        area = primitive_volmdlr.outer_contour2d.area()
        if primitive_volmdlr.inner_contours2d:
            for inner_contour in primitive_volmdlr.inner_contours2d:
                area -= inner_contour.area()
        return area*coeff

    def mass(self):
        volumes = self.volmdlr_primitives()
        mass = 0

        for volume in volumes:
            if volume.__class__.__name__ == 'ExtrudedProfile':

                mass += self.Volume(volume) *hardened_alloy_steel.volumic_mass
            else:
                mass += volume.volume() *hardened_alloy_steel.volumic_mass

        return mass








class PlanetsStructure(DessiaObject):
    '''
    Define a PlanetsStructure (A planetary gears without planetaries)

    :param planets: The list of all the planets of the PlanetStructure
    :type planets: List[Planet]

    :param connections : List of the connection bettween Planet( meshing and Double)
    :type connections:List[Connection]

    :param name : Name
    :type name: str, optional


    '''
    _standalone_in_db = True

    _eq_is_data_eq = False
    def __init__(self, planets: List[Planet], connections: List[Connection], name: str = '', number_group_solution_planet_structure: int = 0):
        self.number_group_solution_planet_structure = number_group_solution_planet_structure
        self.planets = planets
        self.connections = connections

        self.meshings = []
        self.doubles = []
        DessiaObject.__init__(self, name=name)

        for i, connection in  enumerate(self.connections):

          if connection.connection_type != 'D':

             self.meshings.append(MeshingPlanet([connection.nodes[0], connection.nodes[1]], 'meshing'+str(i)))

          else:
             self.doubles.append(Double([connection.nodes[0], connection.nodes[1]], 'Double'+str(i)))

        self.relations = self.meshings + self.doubles

        self.number_double = len(self.doubles)
        self.number_meshing = len(self.meshings)

    def graph(self):

        graph_planetary_gear = nx.Graph()

        for relation in self.relations:

            graph_planetary_gear.add_edge(str(relation.nodes[0]), str(relation))
            graph_planetary_gear.add_edge(str(relation), str(relation.nodes[1]))

            nx.set_node_attributes(graph_planetary_gear, relation.nodes[0], str(relation.nodes[0]))
            nx.set_node_attributes(graph_planetary_gear, relation.nodes[1], str(relation.nodes[1]))
            nx.set_node_attributes(graph_planetary_gear, relation, str(relation))


        return graph_planetary_gear

    # def plot(self):

    #     graph_planetary_gears = self.graph()
    #     plt.figure()
    #     nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)

    def path_planet_to_planet(self):
        '''
        A function which give all the path betwen the first planet of the list planets(input of PlanetStructure) and the other.
        The path includes the planets and the connections(meshing and Doubles)

        :return: list_path
        :rtype: List[List[Planet,meshing,Double]]


        '''
        graph_planetary_gears = self.graph()
        list_path = []

        for planet in self.planets[1:]:
            list_path.append(nx.shortest_path(graph_planetary_gears,
                                              str(self.planets[0]), str(planet)))

        for path in list_path:

            for i in range(len(path)):
                path[i] = nx.get_node_attributes(graph_planetary_gears, path[i])[path[i]]
        return list_path

    def meshing_chain_recursive_function(self, n, planet, graph_planetary_gear, possibilities, list_possibilities, previous_relation):

        planet_2 = copy.copy(planet)
        neighbors = nx.all_neighbors(graph_planetary_gear, str(planet))
        meshing = 0
        neighbors_2 = []
        for neighbor in neighbors:
            neighbor = nx.get_node_attributes(graph_planetary_gear, neighbor)[neighbor]
            neighbors_2.append(neighbor)
        neighbors = neighbors_2
        if not possibilities:
            possibilities.append(planet)
        n += 1
        if previous_relation:

            neighbors.remove(previous_relation)




        for neighbor in neighbors:


            possibilities_2 = copy.copy(possibilities)
            previous_relation = neighbor

            if isinstance(neighbor, MeshingPlanet):
                meshing = 1

                if neighbor.nodes[0] == planet:

                    possibilities_2.append(neighbor.nodes[1])
                    planet_3 = neighbor.nodes[1]

                else:
                    possibilities_2.append(neighbor.nodes[0])
                    planet_3 = neighbor.nodes[0]
                self.meshing_chain_recursive_function(n, planet_3, graph_planetary_gear, possibilities_2,
                                                      list_possibilities, previous_relation)

            elif isinstance(neighbor, Double):

                if neighbor.nodes[0] == planet:

                    planet_3 = neighbor.nodes[1]

                else:
                    planet_3 = neighbor.nodes[0]

                self.meshing_chain_recursive_function(n, planet_3, graph_planetary_gear, [],
                                                      list_possibilities, previous_relation)



        if not meshing:
            list_possibilities.append(possibilities)
            return list_possibilities

        return list_possibilities

    def meshing_chain(self):
        '''
        A function wich return all the meshing chain in the planetary gear.
        A meshing chain is a list of planets which meshed together

        :return: List of meshing chains
        :rtype: List[Planet]


        '''
        graph_planetary_gear = self.graph()
        list_possibilities = self.meshing_chain_recursive_function(0, self.planets[0], graph_planetary_gear, [], [], 0)
        return list_possibilities



    def plot_kinematic_graph_gear(self, coordinate, lenght, diameter, diameter_pivot, lenght_pivot, color, plot_data):
        list_color = ['red', 'blue', 'green', 'red', 'blue', 'green',
                      'red', 'blue', 'green']

        x = [coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
             coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2]

        y = [coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2,
             coordinate[1]+diameter/2, coordinate[1]+diameter/2, coordinate[1]+diameter/2]

        for i in range(len(x)-1):
            point1 = vm.Point2D((x[i], y[i]))
            point2 = vm.Point2D((x[i+1], y[i+1]))
            line = vm.LineSegment2D(point1, point2)
            plot_data.append(line.plot_data('line', color=list_color[color]))

        # plt.plot(x, y, list_color[color])

        x = [coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0], coordinate[0],
             coordinate[0]+lenght/2, coordinate[0]-lenght/2]

        y = [coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2,
             coordinate[1]-diameter/2, coordinate[1]-diameter/2, coordinate[1]-diameter/2]

        for i in range(len(x)-1):
            point1 = vm.Point2D((x[i], y[i]))
            point2 = vm.Point2D((x[i+1], y[i+1]))
            line = vm.LineSegment2D(pont1, point2)
            plot_data.append(line.plot_data('line', color=list_color[color]))

        # plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_double(self, coordinate, diameter, lenght, color, plot_data):
        list_color = ['red', 'blue', 'green', 'red', 'blue', 'green',
                      'red', 'blue', 'green']
        line = []
        x = [coordinate[0], coordinate[0]+lenght]
        y = [coordinate[1]+diameter/2, coordinate[1]+diameter/2]
        for i in range(len(x)-1):
            point1 = vm.Point2D((x[i], y[i]))
            point2 = vm.Point2D((x[i+1], y[i+1]))
            line = vm.LineSegment2D(point1, point2)
            plot_data.append(line.plot_data('line', color=list_color[color]))


        # plt.plot(x, y, list_color[color])

        x = [coordinate[0], coordinate[0]+lenght]
        y = [coordinate[1]-diameter/2, coordinate[1]-diameter/2]
        for i in range(len(x)-1):
            point1 = vm.Point2D((x[i], y[i]))
            point2 = vm.Point2D((x[i+1], y[i+1]))
            line = vm.LineSegment2D(point1, point2)
            plot_data.append(line.plot_data('line', color=list_color[color]))

        # plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y, plot_data):
        coordinate_y_min = 0
        coordinate_y_max = 0
        coordinate_x_max = 0
        for coordinate in coordinates:
            if coordinate[0] > coordinate_x_max:
                coordinate_x_max = coordinate[0]
            if coordinate[1] < coordinate_y_min:
                coordinate_y_min = coordinate[1]
            if coordinate[1] > coordinate_y_max:
                coordinate_y_max = coordinate[1]

        coordinate_planet_carrier = [coordinate_x_max+planet_carrier_x, coordinate_y_min-planet_carrier_y]

        for coordinate in coordinates:
            x = [coordinate[0]-planet_carrier_x, coordinate_planet_carrier[0]]
            y = [coordinate[1], coordinate[1]]
            for i in range(len(x)-1):
                point1 = vm.Point2D((x[i], y[i]))
                point2 = vm.Point2D((x[i+1], y[i+1]))
                line = vm.LineSegment2D(point1, point2)
                plot_data.append(line.plot_data('line', color=list_color[color]))
            # plt.plot(x, y, 'r')

        x = [coordinate_planet_carrier[0]+planet_carrier_x, coordinate_planet_carrier[0], coordinate_planet_carrier[0]]
        y = [coordinate_planet_carrier[1], coordinate_planet_carrier[1], coordinate_y_max]

        for i in range(len(x)-1):
            point1 = vm.Point2D((x[i], y[i]))
            point2 = vm.Point2D((x[i+1], y[i+1]))
            line = vm.LineSegment2D(point1, point2)
            plot_data.append(line.plot_data('line', color=list_color[color]))
        # plt.plot(x, y, 'r')





    def plot_data(self, lenght_gear=0.1, diameter_gear=1, lenght_double=2, diameter_pivot=0.2,
                  lenght_pivot=0.5, planet_carrier_x=2, planet_carrier_y=2):
        '''


        Plot the kinematic graph of the planetary gear

        :param lenght_gear: The width of  the gears. The default is 0.1.
        :type length_gear: float, optional

        :param diameter_gear: The diameter of the gears. The default is 1.
        :type diameter_gear: float, optional

        :param lenght_double: The lenght of the connections betwen 2 double planets. The default is 2.
        :type length_double: float, optional

        :param diameter_pivot: The diameter of the representatives signs of pivot. The default is 0.2.
        :type diameter_pivot: float, optional

        :param lenght_pivot: The length of the representatives signs of pivot. The default is 0.5.
        :type lenght_pivot: float, optional

        :param planet_carrier_x: float, optional The parameter for the position of planet carrer in x. The default is 2.
        :type planet_carrer_x: float, optional

        :param planet_carrier_y: The parameter for the position of planet carrer in y. The default is 2.
        :type planet_carrier_y: float, optional

        '''

        graph_path = self.path_planet_to_planet()

        plt.figure()

        previous_relation_double = []
        previous_relation_meshing = []

        previous_planet_meshing = []
        previous_planet_double = []
        inverse_relation_double = []
        inverse_relation_meshing = []
        coordinate_planet = [[0, 0]]
        coordinate = [0, 0]

        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0, plot_data)
        for path in graph_path:

            plot_data = []
            flag_way_inv_meshing = 0
            flag_way_inv_double = 0
            coordinate = [0, 0]

            color = 0

            for i, element in enumerate(path):

                if isinstance(element, Double):



                    if element in  inverse_relation_double:

                        coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]

                    elif ((element.nodes[0] in previous_planet_double or element.nodes[1] in previous_planet_double) \
                    and not element  in previous_relation_double):

                        for double in previous_relation_double:

                            for node in double.nodes:

                                if element.nodes[0] == node or element.nodes[1] == node:

                                    if  not double == previous_element:

                                        if double in inverse_relation_double:
                                            flag_way_inv_double = 1
                                        else:
                                            flag_way_inv_double = 0

                                    else:

                                        if not double in inverse_relation_double:
                                            flag_way_inv_double = 1
                                        else:
                                            flag_way_inv_double = 0

                        if flag_way_inv_double:

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color, plot_data)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color, plot_data)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:

                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color, plot_data)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]


                            else:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color, plot_data)
                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]
                        else:

                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                    previous_relation_double.append(element)
                    previous_planet_double.extend([element.nodes[0], element.nodes[1]])

                elif isinstance(element, MeshingPlanet):
                    color += 1
                    if element in  inverse_relation_meshing:

                        coordinate = [coordinate[0], coordinate[1]-diameter_gear]


                    elif ((element.nodes[0] in previous_planet_meshing or element.nodes[1] in previous_planet_meshing) \
                    and not element  in previous_relation_meshing):

                        for meshing in previous_relation_meshing:

                            for node in meshing.nodes:

                                if element.nodes[0] == node or element.nodes[1] == node:

                                    if  not meshing == previous_element:

                                        if meshing in inverse_relation_meshing:

                                            flag_way_inv_meshing = 1
                                        else:
                                            flag_way_inv_meshing = 0
                                    else:

                                        if not meshing in inverse_relation_meshing:
                                            flag_way_inv_meshing = 1
                                        else:
                                            flag_way_inv_meshing = 0

                        if flag_way_inv_meshing:
                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                        else:

                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]
                            inverse_relation_meshing.append(element)




                    else:


                        coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                    previous_relation_meshing.append(element)
                    previous_planet_meshing.extend([element.nodes[0], element.nodes[1]])

                if not isinstance(element, Planet):

                    if  not coordinate in coordinate_planet:
                        coordinate_planet.append(coordinate)
                        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color, plot_data)

                    previous_element = element


        self.plot_kinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y, plot_data)

        return plot_data





