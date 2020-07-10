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
import volmdlr.primitives3D as p3d
import volmdlr.primitives2D as p2d
import mechanical_components.meshes as meshes
from dessia_common import DessiaObject
from typing import Tuple, List, TypeVar
import numpy as np
import mechanical_components.genmechanics2 as genmechanics
import mechanical_components.genmechanics2.linkages as linkages
import mechanical_components.genmechanics2.loads as loads


class Gears(DessiaObject):


    def __init__(self, Z: int, name: str = ''):
        self.Z = Z
        
       
        
        DessiaObject.__init__(self, name=name)
        
    # def voldmr_volume_model(self):
    #     model = self.volume_model()
    #     return model

    # def volume_model(self):
    #      # self.module = module
    #      # self.d = module*self.Z
    #      # radius = self.Z*module
    #      # x = vm.Vector3D((1, 0, 0))
    #      # y = vm.Vector3D((0, 1, 0))
    #      # z = vm.Vector3D((0, 0, 1))
    #      # rack = meshes.Rack(0.34, module)
    #      # meshes_1 = meshes.Mesh(self.Z, radius, 0.01, rack)
    #      # Gears3D = {0:meshes_1.Contour(3)}
    #      # export = []
    #      # center_2 = (xy_position[0], xy_position[1])
    #      # center = vm.Point2D(center_2)
    #      # model_trans = Gears3D[0][0].Translation(center)
    #      # model_trans_rot = model_trans.Rotation(center, 0.1)
    #      # Gears3D_Rotate = [model_trans_rot]

    #      # export = []

    #      # for (i, center, k) in zip([Gears3D[0]], [center_2], [-1]):

    #      #                model_export = []

    #      #                for m in i:

    #      #                    center = vm.Point2D(center)
    #      #                    model_trans = m.Translation(center)
    #      #                    model_trans_rot = model_trans.Rotation(center, k)
    #      #                    model_export.append(model_trans_rot)

    #      #                export.append(model_export)

    #      # Gears3D_Rotate = export
    #      # vect_x = z_position*z
    #      # extrusion_vector1 = lenght*z
    #      # C1 = vm.Contour2D(Gears3D_Rotate[0])
    #      # t1 = p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1, [], vm.Vector3D(extrusion_vector1))
    #      if module==0:
    #          return None
         
    #      pos=vm.Point3D(self.position)
    #      axis=vm.vector3D((0,0,1))
    #      radius=(self.module*self.Z)/2
    #      cylinder=vm.cylinder(pos,axis,radius,self.length)
    #      return cylinder

class Planetary(Gears):
    # _standalone_in_db = True

    # _generic_eq = True
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

    def __init__(self, Z: int, planetary_type: str, name: str = '', speed_input : List[float] ='',position : List[float] = '', module:float=0):
        



        self.planetary_type = planetary_type
        self.p = 0
        self.speed = 0
        self.module = module
        self.d = 0
        self.speed_input = speed_input
        self.position=position
        Gears.__init__(self,Z=Z, name=name)
        self.length=20
        self.Z=Z
        self.name=name

        if planetary_type == 'Sun':
            
            self.p = 1

        else:
            
            self.p = -1
    
    def voldmr_volume_model(self):
        model = self.volume_model()
        return model

    def volume_model(self):
        position=self.position
        module=self.module
            
        if self.planetary_type=='Sun':
            pos=vm.Point3D(position)
            axis=vm.Vector3D((0,0,1))
            radius=(module*self.Z)/2
            
            cylinder=p3d.Cylinder(pos,axis,radius,self.length)
            return cylinder
        
        radius=(module*self.Z )/2
        p1 = vm.Point2D((radius,position[2]+self.length/2))
        p2 = vm.Point2D((radius+radius*0.1, position[2]+self.length/2))
        p3 = vm.Point2D((radius,position[2]-self.length/2))
        p4 = vm.Point2D((radius+radius*0.1, position[2]-self.length/2))
        points1 = [p1, p2, p3, p4 ]
        c1 = vm.Polygon2D(points1)
        vector_1=vm.Point3D((0,0,0))
        
        profile1 = p3d.RevolvedProfile(vector_1, vm.Y3D, vm.Z3D, c1, vm.O3D, vm.Z3D)
        
        return profile1
         

class Planet(Gears):
    # _standalone_in_db = True

    # _generic_eq = True
    '''
    Define a planet

    :param Z: The number of tooth 
    :type Z: int
    :param name: Name 
    :type name: str, optional


    '''

    def __init__(self, Z: int,name: str = '',positions: List[float] ='',module:float =0):
        

       
        self.length=20
        self.speed = 0
        self.module = module
        self.speed_input = [0, 0]
        self.Z=Z
        self.name=name
        self.positions=positions
        Gears.__init__(self, Z, name)

    def voldmr_volume_model(self):
        model = self.volume_model()
        return model

    def volume_model(self):

        positions=self.positions
        module=self.module
        model=[]
        for position in positions:
            pos=vm.Point3D(position)
            
            axis=vm.Vector3D((0,0,1))
            radius=(module*self.Z)/2
            
            cylinder=p3d.Cylinder(pos,axis,radius,self.length)
            
            model.append(cylinder)
        
        return model
    
         
        
        
class PlanetCarrier(DessiaObject):
    '''
    Define a planet carrier

    :param name: Name
    :type name: str, optional

 

    '''

    def __init__(self, name: str = '', speed_input : List[float]=[0,0]):
        

        self.speed = 0
        self.speed_input = speed_input
        DessiaObject.__init__(self, name=name)



class Meshing(DessiaObject):


    def __init__(self, nodes: List[Gears], name: str = ''):

        self.nodes = nodes
        DessiaObject.__init__(self, name=name)


class MeshingPlanetary(Meshing):

    def __init__(self, nodes: List[Gears], name: str = ''):

        self.type = 'GI'
        Meshing.__init__(self, nodes, name)

        for node in self.nodes:
            if isinstance(node, Planet):
                self.Z_planet = node.Z
            else:
                self.Z_planetary = node.p*node.Z


    def speed_system_equations(self):

        matrix = npy.array([self.Z_planetary, self.Z_planet, -self.Z_planetary])
        rhs = npy.array([0])
        return matrix, rhs

    def torque_system_equations(self):
        # for node in self.nodes:
        #     if isinstance(node, Planet):
        #         Z_planet = node.Z
        #     else:
        #         Z_planetary = node.Z
        matrix = npy.array([[-1/self.Z_planetary, 1/self.Z_planet],
                            [1/self.Z_planetary, 1/self.Z_planet]])
        rhs = npy.array([0, 0])
        return matrix, rhs


class MeshingPlanet(Meshing):

        def __init__(self, nodes: List[Gears], name: str = ''):

            self.type = 'GI'
            Meshing.__init__(self, nodes, name)
            self.Z_planets = []

            for node in self.nodes:
                     self.Z_planets.append(node.Z)

        def speed_system_equations(self):

            matrix = npy.array([self.Z_planets[0], self.Z_planets[1]])
            rhs = npy.array([0])
            return matrix, rhs

        def torque_system_equations(self):
            matrix = npy.array([-1/self.Z_planets[0], 1/self.Z_planets[1]])
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
    # _standalone_in_db = True

    # _generic_eq = True
    def __init__(self, nodes: List[Planet], name: str = ''):

        self.nodes = nodes
        
        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

    def voldmr_volume_model(self):
        model = self.volume_model()
        return model

    def volume_model(self):
         position_planet_1=self.nodes[0].positions
         position_planet_2=self.nodes[1].positions
         model=[]
         axis=vm.Vector3D((0,0,1))
         for i in range(len(position_planet_1)):
             if position_planet_2[i][2]>position_planet_1[i][2]:
                 if position_planet_2[i][2]>0:
                     position=(position_planet_1[i][0],position_planet_1[i][1],(position_planet_2[i][2]-position_planet_1[i][2])/2)
                 else:
                     position=(position_planet_1[i][0],position_planet_1[i][1],(position_planet_1[i][2]-position_planet_2[i][2])/2)
             else:
                 if position_planet_1[i][2]>0:
                     position=(position_planet_1[i][0],position_planet_1[i][1],(position_planet_1[i][2]-position_planet_2[i][2])/2)
                 else:
                     position=(position_planet_1[i][0],position_planet_1[i][1],(position_planet_2[i][2]-position_planet_1[i][2])/2)
             pos=vm.Point3D(position)
             cylinder = p3d.Cylinder(pos, axis, (self.nodes[0].Z*self.nodes[0].module)/10, abs(position_planet_1[i][2]-position_planet_2[i][2]))
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
    # _standalone_in_db = True

    # _generic_eq = True
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

class PlanetaryGear(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True
    # _non_serializable_attributes =
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
                 planet_carrier: PlanetCarrier, connections: List[Connection], name: str = ''):

 
        self.d_min=0
        self.planetaries = planetaries
        self.planets = planets
        self.planet_carrier = planet_carrier
        self.elements = self.planetaries + self.planets + [self.planet_carrier]
        self.elements_name = []
        for element in self.elements:
            self.elements_name.append(element.name)
            
       
        self.sum_Z_planetary=0
        self.sum_speed_planetary=0
        self.max_Z_planetary=0
        self.min_Z_planetary=100000
        for planetary in self.planetaries:
            self.sum_Z_planetary+=planetary.Z
            d=planetary.module*planetary.Z
            if d>self.d_min:
                 self.d_min=d
            if self.max_Z_planetary<planetary.Z:
                self.max_Z_planetary=planetary.Z
                
            if self.min_Z_planetary>planetary.Z:
                self.min_Z_planetary=planetary.Z
            if planetary.speed_input:
                self.sum_speed_planetary+=planetary.speed_input[0]
        if planet_carrier.speed_input:
            self.speed_planet_carrer= planet_carrier.speed_input[0]
            
        for planet in self.planets:
            if planet.positions:
                d=planet.module*planet.Z+ 2*((planet.positions[0][0])**2+(planet.positions[0][1])**2)**0.5
                       
            if d>self.d_min:
                self.d_min=d
        self.connections = connections
        self.meshings = []
        self.doubles = []
        DessiaObject.__init__(self, name=name)
        self.position=False


            
        
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


        return graph_planetary_gear

    # def plot_graph(self):
 

       
    #     graph_planetary_gears = self.graph()
    #     plt.figure()
    #     nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)



    def plot_kinematic_graph_gear(self, coordinate, lenght, diameter,
                                  diameter_pivot, lenght_pivot, color):

        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
                       coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2])

        y = npy.array([coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2,
                       coordinate[1]+diameter/2, coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
                       coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2])

        y = npy.array([coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2,
                       coordinate[1]-diameter/2, coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_double(self, coordinate, diameter, lenght, color):

        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y):

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
            plt.plot(x, y, 'r')

        x = [coordinate_planet_carrier[0]+planet_carrier_x, coordinate_planet_carrier[0], coordinate_planet_carrier[0]]
        y = [coordinate_planet_carrier[1], coordinate_planet_carrier[1], coordinate_y_max]
        plt.plot(x, y, 'r')

        return coordinate_planet_carrier

    def plot_kinematic_graph_ring(self, coordinate, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color):

        list_color = ['steelblue', 'orchid', 'darkorange', 'palegreen', 'steelblue', 'orchid', 'darkorange', 'palegreen',
                      'steelblue', 'orchid', 'darkorange', 'palegreen', 'steelblue', 'orchid', 'darkorange', 'palegreen']

        x = [coordinate[0]-lenght_gear/2, coordinate[0]+lenght_gear/2, coordinate[0], coordinate[0],
             coordinate_planet_carrier[0]+lenght_ring, coordinate_planet_carrier[0]+lenght_ring]
        y = [coordinate[1], coordinate[1], coordinate[1], diameter_ring/2, diameter_ring/2, coordinate_planet_carrier[1]]

        plt.plot(x, y, list_color[color])
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
        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0)
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

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:
                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]

                            else:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
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
                        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color)


                    previous_element = element


        coordinate_planet_carrier = self.plot_kinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y)
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
                                                   planetary_diameter, diameter_pivot, lenght_pivot, color)

                else:
                    if isinstance(meshing.nodes[0], Planetary):
                        index = index_coordinate_planet.index(meshing.nodes[1])
                    else:
                        index = index_coordinate_planet.index(meshing.nodes[0])

                    lenght_ring = lenght_ring_ini-(((coordinate_planet[index][0])*+100)/50)
                    diameter_ring = diameter_ring_ini-(((coordinate_planet[index][0])*10+100)/50)
                    coordinate_ring = [coordinate_planet[index][0], coordinate_planet[index][1]+diameter_gear/2]

                    self.plot_kinematic_graph_ring(coordinate_ring, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color)






    
    def volmdlr_primitives(self,frame=vm.OXYZ):
        print('A')
        components=self.planetaries+self.planets+self.doubles
        li_box = []
        for component in components:
            shell=component.volume_model()
            if isinstance(component,Planet) or isinstance(component,Double):
                for shell_planet in shell:
                    li_box.append(shell_planet)
            else:
                li_box.append(shell)
               
        return li_box
    
    def plot_data(self):
        plot_data=[]
        primitive_2D=[]
        meshing_chains=self.meshing_chain()
        list_color = ['blue', 'red', 'green', 'black']
        if self.d_min==0:
            planetary_gear=PlanetaryGear(self.planetaries,self.planets,self.planet_carrier,self.connections)
            self.d_min=planetary_gear.d_min
        for planetary in self.planetaries:
            for i,meshing_chain in enumerate(meshing_chains):
                if planetary in meshing_chain:
                    color=list_color[i]
            position=vm.Point2D((0, 0))
            
            d=planetary.module*planetary.Z
            
            circle=vm.Circle2D(position,d/2)
            
            contour=vm.Contour2D([circle],True)        
            plot_data.append(contour.plot_data('contour',stroke_width=self.d_min,color=color)) 
        for planet in self.planets:
            d=planet.module*planet.Z
            # for i,meshing_chain in enumerate(meshing_chains):
            #     if planet in meshing_chain:
            #         color=list_color[i]
            for position in planet.positions:
                position_2=vm.Point2D((position[0], position[1]))
                
                circle=vm.Circle2D(position_2,d/2)

                contour=vm.Contour2D([circle],True)        
                plot_data.append(contour.plot_data('contour',stroke_width=self.d_min,color=color))        
        print(plot_data) 
        return plot_data
                
    


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
    
    def speed_range_test_intervalle_max(self,intervals_max,planetaries,range_planetary_max,range_planet_carrier,speed_min_max,reasons,coeffs_input_1,coeffs_planet_carrier,precision,):
        
        interval_max=intervals_max[0]
        for interval_max_2 in intervals_max:
            if interval_max[0]>interval_max_2 [1] or interval_max[1]<interval_max_2 [0]:
                return [],[]
               
            
            if interval_max[0]<interval_max_2 [0]:
                interval_max[0]=interval_max_2 [0]
                
            if interval_max[1]>interval_max_2 [1]:
                interval_max[1]=interval_max_2 [1]
        

        
        ranges_planet_carrier_speed=[]
        speeds_planetary=[]

        for speed in npy.arange(interval_max[0],interval_max[1]-precision,precision):
            ranges_min_planet_carrier=[]
            for i,planetary in enumerate(planetaries):
                range_min_planet_carrier=copy.copy(range_planet_carrier)
                speed_min=speed_min_max[i][0]
                speed_max=speed_min_max[i][1]
                reason=reasons[i]
                coeff_input_1=coeffs_input_1[i]
                coeff_input_planet_carrier=coeffs_planet_carrier[i]
                
                if reason < 0:
                    reason_abs = abs(reason)
                    
                    if speed_min < planetary.speed_input[0]:
    
                        speed_diff_1_max=range_planetary_max[1]-speed
                        speed_diff_2_min=(planetary.speed_input[0]-speed_min-(speed_diff_1_max)*reason_abs)/(1+reason_abs)
                        range_min_planet_carrier[0] += speed_diff_2_min
                        
                    if speed_max > planetary.speed_input[1]:
                        speed_diff_1_max=-range_planetary_max[0]+speed
                        speed_diff_2_min=(speed_max-planetary.speed_input[1]-(speed_diff_1_max)*reason_abs)/(1+reason_abs)
                        range_min_planet_carrier[1] -= speed_diff_2_min
                else:
                    
                    if reason<1:
                        
                        if speed_min < planetary.speed_input[0]:
                            speed_diff_1_max=-range_planetary_max[0]+speed
                            speed_diff_2_min=(planetary.speed_input[0]-speed_min-(speed_diff_1_max)*reason)/(1-reason)
                            range_min_planet_carrier[0] += speed_diff_2_min

                            
                        if speed_max > planetary.speed_input[1]:
                            speed_diff_1_max=range_planetary_max[1]-speed
                            speed_diff_2_min=(speed_max-planetary.speed_input[1]-speed_min-(speed_diff_1_max)*reason)/(1-reason)
                            range_min_planet_carrier[1] -= speed_diff_2_min
                            
                    else:
                        
                        if speed_min < planetary.speed_input[0]:
                            speed_diff_1_max=-range_planetary_max[0]+speed
                            speed_diff_2_min=(planetary.speed_input[0]-speed_min-(speed_diff_1_max)*reason)/(reason-1)

                            range_min_planet_carrier[1] -= speed_diff_2_min 
                        
                            
                        if speed_max > planetary.speed_input[1]: 
                            speed_diff_1_max=range_planetary_max[1]-speed
                            speed_diff_2_min=(speed_max-planetary.speed_input[1]-(speed_diff_1_max)*reason)/(reason-1)
                          
                            range_min_planet_carrier[0] += speed_diff_2_min 
                            
  
                if range_min_planet_carrier[1]-range_min_planet_carrier[0]<precision:
                    
                    break
                else:
                    ranges_min_planet_carrier.append(range_min_planet_carrier)
            
            if ranges_min_planet_carrier:
                
                range_planet_carrier_speed=ranges_min_planet_carrier[0]
                
                for range_min_planet_carrier in ranges_min_planet_carrier:
                    if range_planet_carrier_speed[0]>range_min_planet_carrier[1] or range_planet_carrier_speed[1]<range_min_planet_carrier[0]:
                        range_planet_carrier_speed=[]
                        break 
                       
                    
                    if range_planet_carrier_speed[0]<range_min_planet_carrier[0]:
                        range_planet_carrier_speed[0]=range_min_planet_carrier[0]
                        
                    if range_planet_carrier_speed[1]>range_min_planet_carrier[1]:
                        range_planet_carrier_speed[1]=range_min_planet_carrier[1]
                        
                if range_planet_carrier_speed:
                    
                    if (range_planet_carrier_speed[1]-range_planet_carrier_speed[0])>precision:
                        ranges_planet_carrier_speed.append(range_planet_carrier_speed)
                        speeds_planetary.append(speed)
        
        
                     
        if not ranges_planet_carrier_speed:
            return [],[]
        speed_diff_2=ranges_planet_carrier_speed[0][1]-ranges_planet_carrier_speed[0][0]
        number_range=0
        for i,range_planet_carrier_speed_2 in  enumerate(ranges_planet_carrier_speed):
            if (range_planet_carrier_speed_2[1]-range_planet_carrier_speed_2[0])>speed_diff_2:
                number_range=i
                speed_diff_2=range_planet_carrier_speed_2[1]-range_planet_carrier_speed_2[0]
                
            
        range_planet_carrier_speed=ranges_planet_carrier_speed[number_range] 
        a=0      
       
        for i,range_speed in enumerate(ranges_planet_carrier_speed):
            flag=0

            if range_planet_carrier_speed[0]>range_speed[1] or range_planet_carrier_speed[1]<range_speed[0]:
                    speeds_planetary.remove(speeds_planetary[i+a])
                    a-=1
                
            if range_planet_carrier_speed[0]>range_speed[0]:
                if range_planet_carrier_speed[1]-range_speed[0]<precision:
                   speeds_planetary.remove(speeds_planetary[i+a])
                   a-=1
                   flag=1
                else:
                    range_planet_carrier_speed[0]=range_speed[0]
                
            if range_planet_carrier_speed[1]>range_speed[1]:
                if -range_planet_carrier_speed[0]+range_speed[1]<precision:
                    if not flag:
                        speeds_planetary.remove(speeds_planetary[i+a])
                        a-=1
                else:        
                    range_planet_carrier_speed[1]=range_speed[1]
                
        
        speed_planetary_max=max(speeds_planetary)
        speed_planetary_min=min(speeds_planetary)

        
        if speed_planetary_min==speed_planetary_max or (range_planet_carrier_speed[1]-range_planet_carrier_speed[0])<precision:
            return [],[]
        
        else:
            return [speed_planetary_min,speed_planetary_max], range_planet_carrier
        
                
            
            
                
                
                        
                        

                

    def speed_range(self, element_1, element_2,precision,list_planetary=[],generator=0,list_path=[]):
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
        if isinstance(element_1,PlanetCarrier):
            input_1=element_2
            input_2=element_1
        else:
            input_1=element_1
            input_2=element_2
            
        if list_planetary == []:
            list_planetary = copy.copy(self.planetaries)

        range_input_1 = [input_1.speed_input[0], input_1.speed_input[1]]
        range_input_2 = [input_2.speed_input[0], input_2.speed_input[1]]
        range_planet_carrier = copy.copy(self.planet_carrier.speed_input)
        ranges_input_1=[]
        ranges_max_input_1=[]
        
        
        ranges_planet_carrier=[]
        ranges_min_planet_carrier=[]
        coeffs_input_1=[]
        coeffs_planet_carrier=[]
        reasons=[]
        speeds_min_max=[]
        
        if not isinstance(input_1, PlanetCarrier) and not isinstance(input_2, PlanetCarrier):
            range_input_1_for=copy.copy(range_input_1)
            range_planet_carrier_for=copy.copy(range_planet_carrier)
            range_max_input_1=copy.copy(range_input_1)


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
                    speed_diff_1_max=0
                    if speed_diff_2_min>input_2.speed_input[1]-input_2.speed_input[0]-precision:
                        speed_diff_2_min=input_2.speed_input[1]-input_2.speed_input[0]-precision
                        speed_diff_1_max=((self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)-speed_diff_2_min)/reason_abs
    
                    
                    range_max_input_1[0]+=speed_diff_1_max
                    
                    
                    range_input_2[0] += coeff_input_2*speed_diff
                 
                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier_for[0] = speed_min
                    
                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1_for[1] -= coeff_input_1*speed_diff
                    speed_diff_2_min = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)
                    speed_diff_1_max=0
                    if speed_diff_2_min>input_2.speed_input[1]-input_2.speed_input[0]-precision:
                        speed_diff_2_min=input_2.speed_input[1]-input_2.speed_input[0]-precision
                        speed_diff_1_max=((speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)-speed_diff_2_min)/reason_abs
                    
                    range_max_input_1[1]-=speed_diff_1_max
                    
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
                    speed_diff_1_max=0
                    if speed_diff_2_min>input_2.speed_input[1]-input_2.speed_input[0]-precision:
                            speed_diff_2_min=input_2.speed_input[1]-input_2.speed_input[0]-precision
                            speed_diff_1_max= ((self.planet_carrier.speed_input[0]-speed_min)*(1-reason)-speed_diff_2_min)/reason
                            
                    if reason < 1:
                        range_input_1_for[1] -= coeff_input_1*speed_diff  
                  
                        range_max_input_1[1]-=speed_diff_1_max
                        
                        range_input_2[0] += coeff_input_2*speed_diff

                    else:
                        range_input_1_for[0] -= coeff_input_1*speed_diff
                        range_max_input_1[0]-=speed_diff_1_max
                        
                        range_input_2[1] += coeff_input_2*speed_diff

                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier_for[0] = speed_min

                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)/(coeff_input_2+coeff_input_1*reason)
                    speed_diff_2_min = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)
                    speed_diff_1_max=0
                    if speed_diff_2_min>input_2.speed_input[1]-input_2.speed_input[0]-precision:
                            speed_diff_2_min=input_2.speed_input[1]-input_2.speed_input[0]-precision
                            speed_diff_1_max= ((speed_max-self.planet_carrier.speed_input[1])*(1-reason)-speed_diff_2_min)/reason
                            
                    if reason < 1:
                        range_input_1_for[0] += coeff_input_1*speed_diff
                        
                        range_max_input_1[0]+=speed_diff_1_max
                        
                        range_input_2[1] -= coeff_input_2*speed_diff

                    else:
                        range_input_1_for[1] += coeff_input_1*speed_diff
                        
                        range_max_input_1[1]+=speed_diff_1_max
                        
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
            input_1_for=input_1
           
        else:
            if  isinstance(input_1, PlanetCarrier):
                list_planetary.remove(input_2)
                input_1_for = input_2
                
                range_input_1= copy.copy(range_input_2)
            else:
                list_planetary.remove(input_1)
                input_1_for = input_1

        range_output = {}
        for i,planetary in enumerate(list_planetary):
            range_for_planetary_input_1=copy.copy(range_input_1)
            range_for_planetary_planet_carrier=copy.copy(range_planet_carrier)
            range_output[planetary] = [planetary.speed_input[0], planetary.speed_input[1]]
            
            range_max_input_1=copy.copy(range_input_1)
            range_min_planet_carrier=copy.copy(range_planet_carrier)
            
            if not list_path:
                path = self.path_planetary_to_planetary([input_1_for, planetary])
        
                reason = self.reason(path[0])
                
            else:
                path =list_path[i]
        
                reason = self.reason(path[0])

            index = self.matrix_position(planetary)
            
            if range_input_1[1] == range_input_1[0] and range_planet_carrier[1] == range_planet_carrier[0]:
                return []

            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))
            coeff_input_planet_carrier = 2*(range_planet_carrier[1]-range_planet_carrier[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))
            
            if reason < 0:

                speed_min = reason*range_input_1[1]+(1-reason)*range_planet_carrier[0] 
                speed_max =  reason*range_input_1[0]+(1-reason)*range_planet_carrier[1] 
                reason_abs = abs(reason)


                if speed_min < planetary.speed_input[0]:

                    speed_diff = (planetary.speed_input[0]-speed_min)/((1+reason_abs)*coeff_input_planet_carrier+reason_abs*coeff_input_1)

                    range_for_planetary_input_1[1] -= coeff_input_1*speed_diff
                    range_for_planetary_planet_carrier[0] += coeff_input_planet_carrier*speed_diff

                    speed_diff_2_min = (planetary.speed_input[0]-speed_min)/(1+reason_abs)
                    if speed_diff_2_min>self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]-precision:
                        speed_diff_2_min=self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]-precision
                        speed_diff_1_max=(planetary.speed_input[0]-speed_min-(speed_diff_2_min)*(1+reason_abs))/(reason_abs)
                        range_max_input_1[1]-=speed_diff_1_max
                        
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
                    
                    if speed_diff_2_min>self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]-precision:
                        speed_diff_2_min=self.planet_carrier.speed_input[1]-self.planet_carrier.speed_input[0]-precision
                        speed_diff_1_max=(speed_max-planetary.speed_input[1]-(speed_diff_2_min)*(1+reason_abs))/(reason_abs)
                        range_max_input_1[0]+=speed_diff_1_max
                    
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
            speeds_min_max.append([speed_min,speed_max])
            reasons.append(reason)
            coeffs_input_1.append(coeff_input_1)
            coeffs_planet_carrier.append(coeff_input_planet_carrier)
            
            
        for i,range_input_1_for in enumerate(ranges_input_1):
            
            if range_input_1[0]>range_input_1_for[1] or range_input_1[1]<range_input_1_for[0]:
                if generator:
                    return 'simplex'
                
                range_input_1=[]
                break
            if range_input_1[1]>range_input_1_for[1]:
                range_input_1[1]=range_input_1_for[1]
                
            if range_input_1[0]<range_input_1_for[0]:
                range_input_1[0]=range_input_1_for[0]
                

            
        for i,range_planet_carrier_for in enumerate(ranges_planet_carrier):
            
            if range_planet_carrier[0]>range_planet_carrier_for[1] or range_planet_carrier[1]<range_planet_carrier_for[0]:
                range_input_1=[]
                if generator:
                    return 'simplex'
                break        
            if range_planet_carrier[1]>range_planet_carrier_for[1]:
                range_planet_carrier[1]=range_planet_carrier_for[1]
                
            if range_planet_carrier[0]<range_planet_carrier_for[0]:
                range_planet_carrier[0]=range_planet_carrier_for[0]
                
        if not range_input_1:
         
            range_input_1,range_planet_carrier=self.speed_range_test_intervalle_max(ranges_max_input_1,list_planetary,
                                                                                    input_1.speed_input,self.planet_carrier.speed_input,
                                                                                    speeds_min_max,reasons, coeffs_input_1,coeffs_planet_carrier,precision)
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
               
                
                
             
            if range_input_1[1]-range_input_1[0]<precision or range_input_2[1]-range_input_2[0]<precision:
                return []
            else:
                range_output[input_1] = range_input_1
                range_output[input_2] = range_input_2
                range_output[self.planet_carrier] = range_planet_carrier

            return range_output
        if range_input_1[1]-range_input_1[0]<precision or range_planet_carrier[1]-range_planet_carrier[0]<precision:
            return[]
        else:
            range_output[input_1_for] = range_input_1
            range_output[self.planet_carrier] = range_planet_carrier
            return range_output



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
    
    def meshing_chain_position_z(self,meshing_chains):
        z=[0]*len(meshing_chains)
        length=30
        z[0]=0
        z_ini=0
        doubles=copy.copy(self.doubles)
        orientation=[0]*len(meshing_chains)
        orientation[0]=0
        for i,meshing_chain in enumerate(meshing_chains):
            number_double=0
            
            for double in doubles:
                if double.nodes[0] in meshing_chain :
                     for j,meshing_chain in enumerate(meshing_chains):
                         if double.nodes[1] in meshing_chain:
                             
                             if orientation[i]==0:
                                 
                                 if number_double==0:
                                     z[j]=z[i]+length
                                     number_double+=1
                                     orientation[j]=1
                                 else:
                                     z[j]=z[i]-length
                                     orientation[j]=-1
                                     
                             else:
                                 z[j]=z[i]+orientation[i]*length
                                 orientation[j]=orientation[i]
                     doubles.remove(double)
                                 
                elif double.nodes[1] in meshing_chain :
                            
                         for j,meshing_chain in enumerate(meshing_chains):
                                 if double.nodes[0] in meshing_chain:
                                     
                                     if orientation[i]==0:
                                         
                                         if number_double==0:
                                             z[j]=z[i]+length
                                             number_double+=1
                                             orientation[j]=1
                                         else:
                                             z[j]=z[i]-length
                                             orientation[j]=-1
                                             
                                     else:
                                         z[j]=z[i]+orientation[i]*length
                                         orientation[j]=orientation[i]
                         
                         doubles.remove(double)    
        return z    
        
    
    def speed_max_planets(self):
        speed_max=0
        for planetary in self.planetaries:
            if planetary.planetary_type=='Ring':
                speed_max_planetary=self.speed_solve({planetary:planetary.speed_input[1],self.planet_carrier:self.planet_carrier.speed_input[0]})
            else:
                speed_max_planetary=self.speed_solve({planetary:planetary.speed_input[0],self.planet_carrier:self.planet_carrier.speed_input[1]})
            for planet in self.planets:

                if abs(speed_max_planetary[planet])>speed_max:
                    speed_max=abs(speed_max_planetary[planet])
        

        return speed_max



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
            if isinstance(input_speeds_and_composants[composant],PlanetCarrier) or isinstance(input_speeds_and_composants[composant],Planetary):
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
       
        element_association={}
        for i in range(len(self.elements)):
            self.elements[i].speed = solution[i]
            element_association[self.elements[i]]=solution[i]
            
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
        print(system_matrix)
        return system_matrix, rhs, element_association
        

    def torque_resolution_PFS(self,input_torque_and_composant,meshing_chains=[]):
        Ca=0
        Cr=0
        Cf=0
        Cwb=0# Speed coeff for bearings
        Cvgs=0# Speed coeff for gear sets
        
        alpha_gs1=20/360*2*3.1415
        beta_gs1=25/360*2*3.1415
        
        egs1=[ 0,1.57079633 ,-1.57 ]
        
        
        
        part_planets=[]
        part_planetaries=[]
        part_planet_carrier=genmechanics.Part('planet_carrier')
        
        ground=genmechanics.Part('ground')
        planet_carrier_a=linkages.BallLinkage(ground,part_planet_carrier,[0,0,0],[0,0,0],Ca,Cr,Cwb,'planet_carrier_a')
        planet_carrier_b=linkages.LinearAnnularLinkage(ground,part_planet_carrier,[0,0,0.1],[0,0,0],Cr,Cwb,'planet_carrier_b')
        link_planetaries_ball=[]
        link_planetaries_linear_annular=[]
        pivot_planets=[]
        flag_double=0
        for i,planet in enumerate(self.planets):
            for double in self.doubles:
                if planet in double.nodes:
                    pivot_planets.append(0)
                    part_planets.append(0)
                    flag_double=1
                    break
            if not flag_double:
                part_planets.append(genmechanics.Part('planet'+str(i)))
                pivot_planets.append(linkages.FrictionlessRevoluteLinkage(part_planet_carrier,part_planets[-1],np.array(planet.positions[0]),[0,0,0],'pivot'+str(i)))
            flag_double=0
                
            
        previous_nodes=[]
     
        for i,double in  enumerate(self.doubles):
            if double.nodes[0] in previous_nodes:
                planet_double=part_planets[self.planets.index(double.nodes[0])]
                link=pivot_planets[self.planets.index(double.nodes[0])]
                part_planets[self.planets.index(double.nodes[1])]=planet_double
                pivot_planets[self.planets.index(double.nodes[1])]=link
                
            
            elif double.nodes[1] in previous_nodes:
                planet_double=part_planets[self.planets.index(double.nodes[1])]
                link=pivot_planets[self.planets.index(double.nodes[1])]
                part_planets[self.planets.index(double.nodes[0])]=planet_double
                pivot_planets[self.planets.index(double.nodes[0])]=link
                
                
            
            else:
                 planet_double=  genmechanics.Part('planet_double'+str(i))
                 link= linkages.FrictionlessRevoluteLinkage(part_planet_carrier,planet_double,np.array(planet.positions[0]),[0,0,0],'pivot'+str(i))
                 part_planets[self.planets.index(double.nodes[0])]=planet_double
                 part_planets[self.planets.index(double.nodes[1])]=planet_double
                 pivot_planets[self.planets.index(double.nodes[0])]=link
                 pivot_planets[self.planets.index(double.nodes[1])]=link
               
            previous_nodes.extend(double.nodes)
     
        
        
            
        for i,planetary in enumerate(self.planetaries):
            part_planetaries.append(genmechanics.Part('planetary'+str(i)))
            link_planetaries_ball.append(linkages.BallLinkage(ground,part_planetaries[-1],np.array(planetary.position),[0,0,0],Ca,Cr,Cwb,'planetary_ball'+str(i)))
            link_planetaries_linear_annular.append(linkages.LinearAnnularLinkage(ground,part_planetaries[-1],
                                                                                 np.array([planetary.position[0],planetary.position[1],planetary.position[2]+0.1]),
                                                                                 [0,0,0],Cr,Cwb,'planetary_linear_angular'+str(i)))
            
        
        
        if not meshing_chains:
            meshing_chains=self.meshing_chain()
        gearings=[]
        for j,meshing_chain in enumerate(meshing_chains):
            previous_element=meshing_chain[0]
            if isinstance(previous_element,Planet):
                previous_position=np.array(previous_element.positions[0])
                previous_part=part_planets[self.planets.index(previous_element)]
            else:
                previous_position=np.array(previous_element.position[0])
                previous_part=part_planetaries[self.planetaries.index(previous_element)]
            
            for i,element in enumerate(meshing_chain):
                if i>0:
                    if isinstance(element,Planetary) :
                        position=np.array(element.position)
                        part=part_planetaries[self.planetaries.index(element)]
                        
                        if element.planetary_type=='Ring':
                            orientation_gearing=previous_position-position
                            angular=m.atan(orientation_gearing[1]/orientation_gearing[0])
                            position_gearing=np.array([previous_position[0]+previous_element.Z*previous_element.module*m.cos(angular),
                                                       previous_position[1]+previous_element.Z*previous_element.module*m.sin(angular),previous_position[2]])
                            
                            
                            
                        else:
                            
                            position_gearing=(position-previous_position)/2 +previous_position
                            
                        
                            
                            
                        
                        
                    else:
                        position=np.array(element.positions[0])
                        part=part_planets[self.planets.index(element)]
                        
                        if isinstance(previous_element,Planetary) and previous_element.planetary_type=='Ring':
                            orientation_gearing=position-previous_position
                            angular=m.atan(orientation_gearing[1]/orientation_gearing[0])
                            position_gearing=np.array([position[0]+element.Z*element.module*m.cos(angular),
                                                       position[1]+element.Z*element.module*m.sin(angular),position[2]])
                            
                            
                        else:
                           
                            position_gearing=(position-previous_position)/2 +previous_position
                           
                        
                        
                        
                    gearings.append(linkages.GearSetLinkage(previous_part,part,position_gearing,egs1,alpha_gs1,beta_gs1,Cf,Cvgs,'Gear set '+str(i) + str(j)))
                    previous_element=element
                    if isinstance(previous_element,Planet):
                        previous_position=np.array(previous_element.positions[0])
                        previous_part=part_planets[self.planets.index(previous_element)]
                    else:
                        previous_position=np.array(previous_element.position[0])
                        previous_part=part_planetaries[self.planetaries.index(previous_element)]
                        
                        
        list_all_input=self.planetaries+[self.planet_carrier]
        loads_known=[]
        for i,input_composant in enumerate(input_torque_and_composant):
            list_all_input.remove(input_composant)
            if input_composant==self.planet_carrier:
                loads_known.append(loads.KnownLoad(part_planet_carrier,[0,0,0],[0,0,0],[0,0,0],[0,0,input_torque_and_composant[input_composant]],'input torque'+str(i)))
            
            else:
                
                loads_known.append(loads.KnownLoad(part_planetaries[self.planetaries.index(input_composant)],[0,0,0],[0,0,0],[0,0,0],[0,0,input_torque_and_composant[input_composant]],'input torque'+str(i)))
        loads_unknown=[]
        for i, unknow_input in enumerate(list_all_input):
            if unknow_input==self.planet_carrier:
                loads_unknown.append(loads.SimpleUnknownLoad(part_planet_carrier,[0,0,0],[0,0,0],[],[0],'output torque'+str(i)))
            else:
                loads_unknown.append(loads.SimpleUnknownLoad(part_planetaries[self.planetaries.index(unknow_input)],[0,0,0],[0,0,0],[],[0],'output torque'+str(i)))
        
        pivot_planets_without_double=[]
        for pivot in pivot_planets:
            if not pivot in pivot_planets_without_double:
                 pivot_planets_without_double.append(pivot)
        
        list_parts=gearings+pivot_planets_without_double+link_planetaries_ball+ link_planetaries_linear_annular+ [planet_carrier_a] + [planet_carrier_b]
        print(list_parts)
        imposed_speeds=[(link_planetaries_ball[0],0,200),(planet_carrier_a,0,200)]
        print(loads_known)
        print(loads_unknown)
        mech=genmechanics.Mechanism(list_parts,ground,imposed_speeds,loads_known,loads_unknown)
        
        for l,lv in mech.static_results.items():
    
            for d,v in lv.items():
                print(l.name,d,v)

         
        
        
                
                        
                    
                    
                    
                    
                    
        
    
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
            if isinstance(composant,PlanetCarrier):
                matrix_input[0][-1] = 1
            else:
                matrix_input[0][element_association[composant]] = 1
            if isinstance(input_torque_and_composant[composant],Planetary):
                 matrix_input[0][element_association[input_torque_and_composant[composant]]] = 1
                 vector_b = npy.concatenate((vector_b, [0]))
            elif isinstance(input_torque_and_composant[composant],PlanetCarrier):
                matrix_input[0][-1] = 1
                vector_b = npy.concatenate((vector_b, [0]))
            else:
                vector_b = npy.concatenate((vector_b, [input_torque_and_composant[composant]]))
                
            system_matrix = npy.concatenate((system_matrix, matrix_input))

        
        solution = solve(system_matrix, vector_b)
        print(solution)
        print(element_association)
        torque_element_association = {}

        for element in element_association:
            max_solution = 0
            for position in element_association[element]:
                if abs(solution[position]) > abs(max_solution):

                    max_solution = solution[position]
            torque_element_association[element] = max_solution

        torque_element_association[self.planet_carrier] = solution[-1]
        return torque_element_association


    def torque_solve_2(self, input_torque_and_composant):
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
            matrix_input[0][element_association[composant]] = 1
            if isinstance(input_torque_and_composant[composant],Planetary) or isinstance(input_torque_and_composant[composant],PlanetCarrier):
                 matrix_input[0][element_association[input_torque_and_composant[composant]]] = 1
                 vector_b = npy.concatenate((vector_b, [0]))
            else:
                vector_b = npy.concatenate((vector_b, [input_torque_and_composant[composant]]))
                
            system_matrix = npy.concatenate((system_matrix, matrix_input))
            
        

          
        solution = solve(system_matrix, vector_b)

        torque_element_association = {}
        print(solution)
        for element in element_association:
            torque_element_association[element] = solution[element_association[element]]

        torque_element_association[self.planet_carrier] = solution[-1]
        return torque_element_association



class PositionMinMaxPlanetaryGear(DessiaObject):
     def __init__(self,planetary_gear : PlanetaryGear,name : str='',positions_min_max:List[float]= '',modules_min_max:List[float]=''):
        
         self.planetary_gear=planetary_gear
         self.positions_min_max=positions_min_max
         self.modules_min_max=modules_min_max
         if self.positions_min_max=='' and self.modules_min_max=='':
             self.positions_min_max=[]
             self.modules_min_max=[]
         DessiaObject.__init__(self, name=self.planetary_gear.name+'PostionMinMax')
         element_list=planetary_gear.planets+ planetary_gear.planetaries
         for element in element_list:
             self.positions_min_max.append([0,0])
             self.modules_min_max.append([0,0])
         
     def enter_position(self,position,element,min_max):
        element_list=self.planetary_gear.planets+ self.planetary_gear.planetaries
        if min_max=='Min':
            self.positions_min_max[element_list.index(element)][0]=position
        if min_max=='Max':
            self.positions_min_max[element_list.index(element)][1]=position
            
     def enter_module(self,module,element,min_max):
        element_list=self.planetary_gear.planets+ self.planetary_gear.planetaries
        if min_max=='Min':
            self.modules_min_max[element_list.index(element)][0]=module
        if min_max=='Max':
            self.modules_min_max[element_list.index(element)][1]=module
            
     def get_position(self,element,planetary_gear,min_max):
        element_list=planetary_gear.planets+ planetary_gear.planetaries
        if min_max=='Min':
            return self.positions_min_max[element_list.index(element)][0]
        if min_max=='Max':
            return self.positions_min_max[element_list.index(element)][1]
        
     def get_module(self,element,planetary_gear,min_max):
        element_list=planetary_gear.planets+ planetary_gear.planetaries
        if min_max=='Min':
            return self.modules_min_max[element_list.index(element)][0]
        if min_max=='Max':
            return self.modules_min_max[element_list.index(element)][1]
        
        
        


class PlanetaryGearResult(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True
    
    def __init__(self,planetary_gear: PlanetaryGear,position_min_max : PositionMinMaxPlanetaryGear, geometry_min_max : str = 'Max'):
        self.planetary_gear=planetary_gear
        self.geometry_min_max=geometry_min_max
        self.position_min_max=position_min_max
        
        
        self.speed_max_planet=self.planetary_gear.speed_max_planets()
        
        self.d_min=self.planetary_gear.d_min
        self.sum_Z_planetary=self.planetary_gear.sum_Z_planetary
        self.sum_speed_planetary=self.planetary_gear.sum_speed_planetary
        
        self.max_Z_planetary=self.planetary_gear.max_Z_planetary
        self.min_Z_planetary=self.planetary_gear.min_Z_planetary
        self.speed_planet_carrer=self.planetary_gear.speed_planet_carrer
        DessiaObject.__init__(self, name=self.planetary_gear.name+'Result')
    
    def volmdlr_primitives(self,frame=vm.OXYZ):
        print(self.position_min_max.positions_min_max)
        if self.geometry_min_max:
            for planet in self.planetary_gear.planets:
                planet.positions=self.position_min_max.get_position(planet,self.planetary_gear,self.geometry_min_max)
                planet.module=self.position_min_max.get_module(planet,self.planetary_gear,self.geometry_min_max)
            for planetary in self.planetary_gear.planetaries:
                planetary.position=self.position_min_max.get_position(planetary,self.planetary_gear,self.geometry_min_max)
                planetary.module=self.position_min_max.get_module(planetary,self.planetary_gear,self.geometry_min_max)
              
        li_box=self.planetary_gear.volmdlr_primitives()
        return li_box
    
    def plot_data(self):
        # if self.geometry_min_max:
        #     for planet in self.planetary_gear.planets:
        #         planet.positions=self.position_min_max.get_position(planet,self.geometry_min_max)
        #         planet.module=self.position_min_max.get_module(planet,self.geometry_min_max)
        #     for planetary in self.planetary_gear.planetaries:
        #         planetary.position=self.position_min_max.get_position(planetary,self.geometry_min_max)
        #         planetary.module=self.position_min_max.module(planetary,self.geometry_min_max)
        
        plot_data=self.planetary_gear.plot_data()
        return plot_data




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

    def __init__(self, planets: List[Planet], connections: List[Connection], name: str = ''):
        
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



    def plot_kinematic_graph_gear(self, coordinate, lenght, diameter, diameter_pivot, lenght_pivot, color):
        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0],
                       coordinate[0], coordinate[0]+lenght/2, coordinate[0]-lenght/2])

        y = npy.array([coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2, coordinate[1]+diameter_pivot/2,
                       coordinate[1]+diameter/2, coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0]+lenght_pivot/2, coordinate[0]-lenght_pivot/2, coordinate[0], coordinate[0],
                       coordinate[0]+lenght/2, coordinate[0]-lenght/2])

        y = npy.array([coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2, coordinate[1]-diameter_pivot/2,
                       coordinate[1]-diameter/2, coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_double(self, coordinate, diameter, lenght, color):
        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_kinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y):
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
            plt.plot(x, y, 'r')

        x = [coordinate_planet_carrier[0]+planet_carrier_x, coordinate_planet_carrier[0], coordinate_planet_carrier[0]]
        y = [coordinate_planet_carrier[1], coordinate_planet_carrier[1], coordinate_y_max]
        plt.plot(x, y, 'r')





    def plot_kinematic_graph(self, lenght_gear=0.1, diameter_gear=1, lenght_double=2, diameter_pivot=0.2,
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

        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0)
        for path in graph_path:


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

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:

                            self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:

                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]


                            else:
                                self.plot_kinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
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
                        self.plot_kinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color)

                    previous_element = element


        self.plot_kinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y)



