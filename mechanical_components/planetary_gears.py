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
from typing import TypeVar, List

class Gears(DessiaObject):


    def __init__(self, name: str, Z: int):
        self.Z = Z

        DessiaObject.__init__(self, name=name)

    def volume_plot(self, xy_position, z_position, module, lenght):
         self.module = module
         self.d = module*self.Z
         radius = self.Z*module
         x = vm.Vector3D((1, 0, 0))
         y = vm.Vector3D((0, 1, 0))
         z = vm.Vector3D((0, 0, 1))
         rack = meshes.Rack(0.34, module)
         meshes_1 = meshes.Mesh(self.Z, radius, 0.01, rack)
         Gears3D = {0:meshes_1.Contour(3)}
         export = []
         center_2 = (xy_position[0], xy_position[1])
         center = vm.Point2D(center_2)
         model_trans = Gears3D[0][0].Translation(center)
         model_trans_rot = model_trans.Rotation(center, 0.1)
         Gears3D_Rotate = [model_trans_rot]

         export = []

         for (i, center, k) in zip([Gears3D[0]], [center_2], [-1]):

                        model_export = []

                        for m in i:

                            center = vm.Point2D(center)
                            model_trans = m.Translation(center)
                            model_trans_rot = model_trans.Rotation(center, k)
                            model_export.append(model_trans_rot)

                        export.append(model_export)

         Gears3D_Rotate = export
         vect_x = z_position*z
         extrusion_vector1 = lenght*z
         C1 = vm.Contour2D(Gears3D_Rotate[0])
         t1 = p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1, [], vm.Vector3D(extrusion_vector1))

         return t1

class Planetary(Gears):

    def __init__(self, name: str, Z: int, planetary_type: str):



        self.planetary_type = planetary_type
        self.p = 0
        self.speed = 0
        self.module = 0
        self.d = 0
        self.speed_input = [0, 0]
        Gears.__init__(self, name, Z)

        if planetary_type == 'Sun':
            self.p = 1

        else:
            self.p = -1

class Planet(Gears):

    def __init__(self, name: str, planet_type: str, Z: int):

        self.planet_type = planet_type

        self.speed = 0
        self.module = 0
        self.speed_input = [0, 0]

        Gears.__init__(self, name, Z)

class PlanetCarrier():


    def __init__(self, name):
        self.name = name
        self.speed = 0
        self.speed_input = [0, 0]
        # DessiaObject.__init__(self, name=name)



class Gearing(DessiaObject):


    def __init__(self, name, nodes):

        self.nodes = nodes
        DessiaObject.__init__(self, name=name)


class GearingPlanetary(Gearing):

    def __init__(self, name: str, nodes: List[Gears]):

        self.type = 'GI'
        Gearing.__init__(self, name, nodes)

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

        matrix = npy.array([[-1/self.Z_planetary, 1/self.Z_planet],
                            [1/self.Z_planetary, 1/self.Z_planet]])
        rhs = npy.array([0, 0])
        return matrix, rhs


class GearingPlanet(Gearing):

        def __init__(self, name: str, nodes: List[Gears]):

            self.type = 'GI'
            Gearing.__init__(self, name, nodes)
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
    def __init__(self, name: str, nodes: A):
        self.nodes = nodes

        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Fixed(DessiaObject):

    def __init__(self, name, nodes):
        self.nodes = nodes

        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Double(DessiaObject):

    def __init__(self, name: str, nodes: List[Planet]):

        self.nodes = nodes
        DessiaObject.__init__(self, name=name)
    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

    def volume_plot(self, xy_position, z_position, radius, lenght):

         pos = vm.Point3D((xy_position[0], xy_position[1], z_position))
         axis = vm.Vector3D((0, 0, 1))
         cylinder = p3d.Cylinder(pos, axis, radius, lenght)
         return cylinder

class ImposeSpeed(DessiaObject):
    A = TypeVar('A', Gears, PlanetCarrier)
    def __init__(self, name: str, node: A, input_speed: float):
        self.input_speed = input_speed
        self.node = node

        DessiaObject.__init__(self, name=name)

    def speed_system_equations(self):
        matrix = npy.array([1])
        rhs = npy.array([self.input_speed])
        return matrix, rhs


class PlanetaryGear(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True
    A = TypeVar('A', Gears, str)
    def __init__(self, name: str, planetaries: List[Planetary], planets: List[Planet],
                 planet_carrier: PlanetCarrier, connexions: List[List[A]]):
        '''
        Define a PlanetaryGears

        Parameters
        ----------
        name : String
            DESCRIPTION.
        planetaries : TYPE
            DESCRIPTION.
        planets : TYPE
            DESCRIPTION.
        planet_carrier : TYPE
            DESCRIPTION.
        connexions : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        self.planetaries = planetaries
        self.planets = planets
        self.planet_carrier = planet_carrier
        self.elements = self.planetaries + self.planets + [self.planet_carrier]
        self.elements_name = []
        for element in self.elements:
            self.elements_name.append(element.name)

        self.connexions = connexions
        self.gearings = []
        self.doubles = []
        DessiaObject.__init__(self, name=name)
        for i, connexion in  enumerate(connexions):

          ## Check to be sure that all the object in connexion are in planetaries,
          ## planets, or planet_carrier ##
          if not connexion[1] in self.elements:

                if isinstance(connexion[1], Planetary):

                    self.elements[self.elements_name.index(connexion[1].name)].planetary_type = connexion[1].planetary_type
                    self.elements[self.elements_name.index(connexion[1].name)].p = connexion[1].p

                connexion[1] = self.elements[self.elements_name.index(connexion[1].name)]

          if not connexion[0] in self.elements:
                if isinstance(connexion[0], Planetary):
                    self.elements[self.elements_name.index(connexion[0].name)].planetary_type = connexion[0].planetary_type
                    self.elements[self.elements_name.index(connexion[0].name)].p = connexion[0].p
                connexion[0] = self.elements[self.elements_name.index(connexion[0].name)]




          if connexion[2] != 'D':

              if isinstance(connexion[0], Planet) and isinstance(connexion[1], Planet):
                self.gearings.append(GearingPlanet('Gearing'+str(i), [connexion[0], connexion[1]]))


              else:
                self.gearings.append(GearingPlanetary('Gearing'+str(i),
                                                      [connexion[0], connexion[1]]))


              self.gearings[-1].type = connexion[2]


          else:
             self.doubles.append(Double('Double'+str(i), [connexion[0], connexion[1]]))

        self.relations = self.gearings + self.doubles

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
        connexions_name = copy.deepcopy(self.connexions)
        for i in range(len(connexions_name)):

            connexions_name[i][0] = connexions_name[i][0].name
            connexions_name[i][1] = connexions_name[i][1].name

        return 'Name:' + self.name + '\n\n' + \
               'Planetary Number:' + str(len(self.planetaries)) + '\n' + \
               'Ring Number:'+ str(number_ring) + '\n' + \
               'Sun_Number:' + str(number_sun) + '\n' + \
               'Z_planetaries:' + str(Z_planetaries) + '\n\n' + \
               'Planets_Number:' + str(len(self.planets)) + '\n' + \
               'Planets_Double_Number:' + str(len(self.doubles)) + '\n' + \
               'Z_Planets:' + str(Z_planets) + '\n\n' + \
                str(connexions_name) + '\n\n\n'


    def matrix_position(self, element):

        return self.elements_name.index(element.name)

    def graph(self):

        graph_planetary_gear = nx.Graph()

        for relation in self.relations:

            graph_planetary_gear.add_edge(relation.nodes[0].name, relation.name)
            graph_planetary_gear.add_edge(relation.name, relation.nodes[1].name)

            nx.set_node_attributes(graph_planetary_gear, relation.nodes[0], relation.nodes[0].name)
            nx.set_node_attributes(graph_planetary_gear, relation.nodes[1], relation.nodes[1].name)
            nx.set_node_attributes(graph_planetary_gear, relation, relation.name)



        for k, planet in enumerate(self.planets):

            graph_planetary_gear.add_edge(self.planet_carrier.name, 'Pv'+str(k))
            graph_planetary_gear.add_edge('Pv'+str(k), planet.name)
            nx.set_node_attributes(graph_planetary_gear, 'Pv'+str(k), 'Pv'+str(k))


        return graph_planetary_gear

    def plot(self):

        graph_planetary_gears = self.graph()
        plt.figure()
        nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)



    def plot_cinematic_graph_gear(self, coordinate, lenght, diameter,
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

    def plot_cinematic_graph_double(self, coordinate, diameter, lenght, color):

        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_cinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y):

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

    def plot_cinematic_graph_ring(self, coordinate, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color):

        list_color = ['steelblue', 'orchid', 'darkorange', 'palegreen', 'steelblue', 'orchid', 'darkorange', 'palegreen',
                      'steelblue', 'orchid', 'darkorange', 'palegreen', 'steelblue', 'orchid', 'darkorange', 'palegreen']

        x = [coordinate[0]-lenght_gear/2, coordinate[0]+lenght_gear/2, coordinate[0], coordinate[0],
             coordinate_planet_carrier[0]+lenght_ring, coordinate_planet_carrier[0]+lenght_ring]
        y = [coordinate[1], coordinate[1], coordinate[1], diameter_ring/2, diameter_ring/2, coordinate_planet_carrier[1]]

        plt.plot(x, y, list_color[color])
        coordinate[1] -= (abs(coordinate[1]-coordinate_planet_carrier[1]))*2


    def plot_cinematic_graph(self, lenght_gear=0.1, diameter_gear=1, lenght_double=2, diameter_pivot=0.2, lenght_pivot=0.5,
                             planet_carrier_x=2, planet_carrier_y=2, diameter_ring_ini=10):

        graph_path = self.path_planetary_to_planetary()

        plt.figure()

        previous_relation_double = []
        previous_relation_gearing = []

        previous_planet_gearing = []
        previous_planet_double = []
        inverse_relation_double = []

        coordinate_planet = [[0, 0]]
        coordinate = [0, 0]
        index_coordinate_planet = []
        flag_first_planet = 0
        self.plot_cinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0)
        for path in graph_path:



            previous_element = 0
            flag_way_inv_double = 0
            coordinate = [0, 0]

            color = 0

            for i, element in enumerate(path):
                if not flag_first_planet and isinstance(element, Planet):
                    if len(coordinate_planet) < 2:
                        index_coordinate_planet.append(element.name)
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

                            self.plot_cinematic_graph_double(coordinate, diameter_pivot, lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:
                            self.plot_cinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:
                                self.plot_cinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]

                            else:
                                self.plot_cinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]
                        else:

                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                    previous_relation_double.append(element)
                    previous_planet_double.extend([element.nodes[0], element.nodes[1]])

                elif isinstance(element, GearingPlanet):

                    color += 1

                    if element.type == 'GI':

                        if previous_planet.name == element.nodes[0].name:

                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]
                        else:
                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                    elif element.type == 'GE':

                        if previous_planet.name == element.nodes[0].name:

                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]
                        else:
                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]





                    previous_relation_gearing.append(element)
                    previous_planet_gearing.extend([element.nodes[0], element.nodes[1]])

                if isinstance(element, Planet):
                    previous_planet = element

                if not isinstance(element, Planet) and not isinstance(element, Planetary) \
                and not isinstance(element, GearingPlanetary):

                    if  not coordinate in coordinate_planet:
                        coordinate_planet.append(coordinate)

                        if not element.nodes[1].name in index_coordinate_planet:
                            index_coordinate_planet.append(element.nodes[1].name)

                        else:

                            index_coordinate_planet.append(element.nodes[0].name)
                        self.plot_cinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color)


                    previous_element = element


        coordinate_planet_carrier = self.plot_cinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y)
        lenght_ring_ini = 5
        for gearing in self.gearings:

            if isinstance(gearing, GearingPlanetary):
                color += 1

                if (isinstance(gearing.nodes[0], Planetary) and gearing.nodes[0].planetary_type == 'Sun') \
                or (isinstance(gearing.nodes[1], Planetary) and gearing.nodes[1].planetary_type == 'Sun'):

                    if isinstance(gearing.nodes[0], Planetary):
                        index = index_coordinate_planet.index(gearing.nodes[1].name)
                    else:
                        index = index_coordinate_planet.index(gearing.nodes[0].name)

                    planetary_diameter = ((coordinate_planet[index][1]-diameter_gear/2)-coordinate_planet_carrier[1])*2

                    self.plot_cinematic_graph_gear([coordinate_planet[index][0], coordinate_planet_carrier[1]], lenght_gear,
                                                   planetary_diameter, diameter_pivot, lenght_pivot, color)

                else:
                    if isinstance(gearing.nodes[0], Planetary):
                        index = index_coordinate_planet.index(gearing.nodes[1].name)
                    else:
                        index = index_coordinate_planet.index(gearing.nodes[0].name)

                    lenght_ring = lenght_ring_ini-(((coordinate_planet[index][0])*+100)/50)
                    diameter_ring = diameter_ring_ini-(((coordinate_planet[index][0])*10+100)/50)
                    coordinate_ring = [coordinate_planet[index][0], coordinate_planet[index][1]+diameter_gear/2]

                    self.plot_cinematic_graph_ring(coordinate_ring, lenght_gear, coordinate_planet_carrier, diameter_ring, lenght_ring, color)








    def path_planetary_to_planetary(self, planetaries=[]):
        if not planetaries:
            planetaries = self.planetaries
        graph_planetary_gears = self.graph()
        graph_planetary_gears.remove_node(self.planet_carrier.name)
        list_path = []

        for planetary in planetaries[1:]:
            list_path.append(nx.shortest_path(graph_planetary_gears,
                                              planetaries[0].name, planetary.name))

        for path in list_path:

            for i in range(len(path)):
                path[i] = nx.get_node_attributes(graph_planetary_gears, path[i])[path[i]]


        return list_path

    def path_planetary_to_planetary_type(self):
        list_path = self.path_planetary_to_planetary()
        for i in range(len(list_path)):

            for j in range(len(list_path[i])):

                if isinstance(list_path[i][j], Planetary):

                    list_path[i][j] = list_path[i][j].planetary_type

                else:

                    list_path[i][j] = str(type(list_path[i][j]))
        return list_path



    def reason_abs(self, path):
        reason = 1
        for i, element in enumerate(path):

            if isinstance(element, Gearing):


                reason = reason*path[i-1].Z/path[i+1].Z


        return reason

    def reason(self, path):
        reason = 1

        for i, element in enumerate(path):

            if isinstance(element, Gearing):

                reason = reason*-path[i-1].Z/path[i+1].Z

                if (isinstance(path[i-1], Planetary) and path[i-1].planetary_type == 'Ring')\
                or (isinstance(path[i+1], Planetary) and path[i+1].planetary_type == 'Ring'):

                    reason = -reason

        return reason

    def speed_range(self, input_1, input_2, list_planetary=[]):

        if list_planetary == []:
            list_planetary = copy.copy(self.planetaries)

        range_input_1 = [input_1.speed_input[0], input_1.speed_input[1]]
        range_input_2 = [input_2.speed_input[0], input_2.speed_input[1]]
        range_planet_carrier = copy.copy(self.planet_carrier.speed_input)

        if not isinstance(input_1, PlanetCarrier) and not isinstance(input_2, PlanetCarrier):

            path = self.path_planetary_to_planetary([input_1, input_2])
            reason = self.reason(path[0])
            index = self.matrix_position(self.planet_carrier)

            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))
            coeff_input_2 = 2*(range_input_2[1]-range_input_2[0])/((range_input_1[1]-range_input_1[0])+(range_input_2[1]-range_input_2[0]))

            if reason < 0:
                reason_abs = abs(reason)
                speed_min = self.speed_solve({input_1 : input_1.speed_input[0], input_2 : input_2.speed_input[0]})[index]
                speed_max = self.speed_solve({input_1 : input_1.speed_input[1], input_2 : input_2.speed_input[1]})[index]


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1[0] += coeff_input_1*speed_diff
                    range_input_2[0] += coeff_input_2*speed_diff

                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier[0] = speed_min

                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1+reason_abs)/(coeff_input_2+coeff_input_1*reason_abs)

                    range_input_1[1] -= coeff_input_1*speed_diff
                    range_input_2[1] -= coeff_input_2*speed_diff

                elif speed_max > self.planet_carrier.speed_input[0]:
                    range_planet_carrier[1] = speed_max

                else:
                    return []


            else:


                if reason < 1:
                    speed_min = self.speed_solve({input_1 : input_1.speed_input[1], input_2 : input_2.speed_input[0]})[index]
                    speed_max = self.speed_solve({input_1 : input_1.speed_input[0], input_2 : input_2.speed_input[1]})[index]

                else:
                    speed_min = self.speed_solve({input_1 : input_1.speed_input[0], input_2 : input_2.speed_input[1]})[index]
                    speed_max = self.speed_solve({input_1 : input_1.speed_input[1], input_2 : input_2.speed_input[0]})[index]


                if speed_min < self.planet_carrier.speed_input[0]:

                    speed_diff = (self.planet_carrier.speed_input[0]-speed_min)*(1-reason)/(coeff_input_2+coeff_input_1*reason)

                    if reason < 1:
                        range_input_1[1] -= coeff_input_1*speed_diff
                        range_input_2[0] += coeff_input_2*speed_diff

                    else:
                        range_input_1[0] -= coeff_input_1*speed_diff
                        range_input_2[1] += coeff_input_2*speed_diff

                elif speed_min < self.planet_carrier.speed_input[1]:
                    range_planet_carrier[0] = speed_min

                else:
                    return []


                if speed_max > self.planet_carrier.speed_input[1]:

                    speed_diff = (speed_max-self.planet_carrier.speed_input[1])*(1-reason)/(coeff_input_2+coeff_input_1*reason)

                    if reason < 1:
                        range_input_1[0] += coeff_input_1*speed_diff
                        range_input_2[1] -= coeff_input_2*speed_diff

                    else:
                        range_input_1[1] += coeff_input_1*speed_diff
                        range_input_2[0] -= coeff_input_2*speed_diff

                elif speed_max > self.planet_carrier.speed_input[0]:
                    range_planet_carrier[1] = speed_max

                else:
                    return []

            # print(range_input_1)
            # print(range_input_2)
            # print(range_planet_carrier)
            list_planetary.remove(input_1)
            list_planetary.remove(input_2)

        else:
            if  isinstance(input_1, PlanetCarrier):
                list_planetary.remove(input_2)
                input_1 = input_2
            else:
                list_planetary.remove(input_1)


        range_output = {}

        for planetary in list_planetary:

            range_output[planetary.name] = [planetary.speed_input[0], planetary.speed_input[1]]
            path = self.path_planetary_to_planetary([input_1, planetary])

            reason = self.reason(path[0])

            index = self.matrix_position(planetary)
            c = ((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))
            if range_input_1[1] == range_input_1[0] and range_planet_carrier[1] == range_planet_carrier[0]:
                return []

            coeff_input_1 = 2*(range_input_1[1]-range_input_1[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))
            coeff_input_planet_carrier = 2*(range_planet_carrier[1]-range_planet_carrier[0])/((range_input_1[1]-range_input_1[0])+(range_planet_carrier[1]-range_planet_carrier[0]))

            if reason < 0:

                speed_min = self.speed_solve({input_1 : range_input_1[1], self.planet_carrier : range_planet_carrier[0]})[index]
                speed_max = self.speed_solve({input_1 : range_input_1[0], self.planet_carrier : range_planet_carrier[1]})[index]
                reason_abs = abs(reason)


                if speed_min < planetary.speed_input[0]:

                    speed_diff = (planetary.speed_input[0]-speed_min)/((1+reason_abs)*coeff_input_planet_carrier+reason_abs*coeff_input_1)

                    range_input_1[1] -= coeff_input_1*speed_diff
                    range_planet_carrier[0] += coeff_input_planet_carrier*speed_diff

                elif speed_min < planetary.speed_input[1]:
                    range_output[planetary.name][0] = speed_min

                else:
                    return[]


                if speed_max > planetary.speed_input[1]:
                    speed_diff = (speed_max-planetary.speed_input[1])/((1+reason_abs)*coeff_input_planet_carrier+reason_abs*coeff_input_1)
                    range_input_1[0] += coeff_input_1*speed_diff
                    range_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff

                elif speed_max > planetary.speed_input[0]:
                    range_output[planetary.name][1] = speed_max
                else:
                    return []


            else:

                if reason < 1:

                    speed_min = self.speed_solve({input_1 : range_input_1[0], self.planet_carrier : range_planet_carrier[0]})[index]
                    speed_max = self.speed_solve({input_1 : range_input_1[1], self.planet_carrier : range_planet_carrier[1]})[index]

                else:
                    speed_min = self.speed_solve({input_1 : range_input_1[0], self.planet_carrier : range_planet_carrier[1]})[index]
                    speed_max = self.speed_solve({input_1 : range_input_1[1], self.planet_carrier : range_planet_carrier[0]})[index]

                if speed_min < planetary.speed_input[0]:

                    if reason < 1:
                        speed_diff = (planetary.speed_input[0]-speed_min)/(coeff_input_1*reason+coeff_input_planet_carrier*(1-reason))
                        range_input_1[0] += coeff_input_1*speed_diff
                        range_planet_carrier[0] += coeff_input_planet_carrier*speed_diff


                    else:
                        speed_diff = (planetary.speed_input[0]-speed_min)/(reason*coeff_input_1-coeff_input_planet_carrier*(reason-1))
                        range_input_1[0] += coeff_input_1*speed_diff
                        range_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff

                elif speed_min < planetary.speed_input[1]:
                    range_output[planetary.name][0] = speed_min

                else:
                    return []


                if speed_max > planetary.speed_input[1]:

                    if reason < 1:
                        speed_diff = (speed_max-planetary.speed_input[1])/(coeff_input_1*reason+coeff_input_planet_carrier*(1-reason))
                        range_input_1[1] -= coeff_input_1*speed_diff
                        range_planet_carrier[1] -= coeff_input_planet_carrier*speed_diff

                    else:
                        speed_diff = (speed_max-planetary.speed_input[1])/(reason*coeff_input_1-coeff_input_planet_carrier*(reason-1))
                        range_input_1[1] -= coeff_input_1*speed_diff
                        range_planet_carrier[0] += coeff_input_planet_carrier*speed_diff


                elif speed_max > planetary.speed_input[0]:
                    range_output[planetary.name][1] = speed_max

                else:
                    return []


        if not isinstance(input_1, PlanetCarrier) and not isinstance(input_2, PlanetCarrier):

            path = self.path_planetary_to_planetary([input_1, input_2])
            reason = self.reason(path[0])
            index = self.matrix_position(self.planet_carrier)

            if reason < 0:

                speed_min = self.speed_solve({input_1:range_input_1[0], input_2:range_input_2[0]})[index]
                speed_max = self.speed_solve({input_1:range_input_1[1], input_2:range_input_2[1]})[index]
                reason_abs = abs(reason)



                speed_diff_min = (range_planet_carrier[0]-speed_min)*(1+reason_abs)
                range_input_2[0] += speed_diff_min

                speed_diff_max = (speed_max-range_planet_carrier[1])*(1+reason_abs)
                range_input_2[1] -= speed_diff_max


            else:
                if reason < 1:

                    speed_min = self.speed_solve({input_1 : range_input_1[1], input_2 : range_input_2[0]})[index]
                    speed_max = self.speed_solve({input_1 : range_input_1[0], input_2 : range_input_2[1]})[index]

                else:
                    speed_min = self.speed_solve({input_1 : range_input_1[0], input_2 : range_input_2[1]})[index]
                    speed_max = self.speed_solve({input_1 : range_input_1[1], input_2 : range_input_2[0]})[index]



                speed_diff_min = (range_planet_carrier[0]-speed_min)*(1-reason)

                if reason < 1:
                    range_input_2[0] += speed_diff_min

                else:
                    range_input_2[1] += speed_diff_min


                speed_diff_max = (speed_max-range_planet_carrier[1])*(1-reason)

                if reason < 1:
                    range_input_2[1] -= speed_diff_max

                else:
                    range_input_2[0] -= speed_diff_max




            range_output[input_1.name] = range_input_1
            range_output[input_2.name] = range_input_2
            range_output[self.planet_carrier.name] = range_planet_carrier
            return range_output

        range_output[input_1.name] = range_input_1
        range_output[self.planet_carrier.name] = range_planet_carrier
        return range_output



    def gearing_chain_recursive_function(self, number_gearing_chain, element, graph_planetary_gear, possibilities,
                                         list_possibilities, previous_relation, gearing_way, previous_gearing_chains):


        neighbors = nx.all_neighbors(graph_planetary_gear, element.name)
        gearing = 0
        neighbors_2 = []
        for neighbor in neighbors:
            neighbor = nx.get_node_attributes(graph_planetary_gear, neighbor)[neighbor]

            neighbors_2.append(neighbor)

        if not possibilities:
            possibilities.append(element)

        if previous_relation:

            neighbors_2.remove(previous_relation)




        for neighbor in neighbors_2:

            gearing = copy.copy(gearing)

            previous_relation = neighbor

            if isinstance(neighbor, Gearing):

                if gearing == 0 or gearing_way == 0:
                    if neighbor.nodes[0].name == element.name:
                        possibilities.append(neighbor.nodes[1])
                        element_3 = neighbor.nodes[1]

                    else:
                        possibilities.append(neighbor.nodes[0])
                        element_3 = neighbor.nodes[0]

                    self.gearing_chain_recursive_function(number_gearing_chain, element_3, graph_planetary_gear, possibilities,
                                                          list_possibilities, previous_relation, 0, previous_gearing_chains)
                    gearing_way = -1

                if gearing == 1 or gearing_way == 1:

                    if neighbor.nodes[0].name == element.name:
                        possibilities = [neighbor.nodes[1]] + possibilities
                        element_3 = neighbor.nodes[1]

                    else:
                        possibilities = [neighbor.nodes[0]] + possibilities

                        element_3 = neighbor.nodes[0]
                    self.gearing_chain_recursive_function(number_gearing_chain, element_3, graph_planetary_gear, possibilities,
                                                          list_possibilities, previous_relation, 1, previous_gearing_chains)
                gearing = 1


            elif isinstance(neighbor, Double):

                if neighbor.nodes[0].name == element.name:

                    element_3 = neighbor.nodes[1]

                else:
                    element_3 = neighbor.nodes[0]

                self.gearing_chain_recursive_function(number_gearing_chain+1, element_3, graph_planetary_gear, [],
                                                      list_possibilities, previous_relation, 0, previous_gearing_chains)



        if not gearing:

                if not number_gearing_chain in previous_gearing_chains:

                    previous_gearing_chains.append(number_gearing_chain)
                    list_possibilities.append(possibilities)

                else:
                    index = previous_gearing_chains.index(number_gearing_chain)
                    list_possibilities[index] = possibilities




        return list_possibilities

    def gearing_chain(self):
        graph_planetary_gear = self.graph()
        list_possibilities = self.gearing_chain_recursive_function(0, self.planetaries[0], graph_planetary_gear, [], [], 0, 0, [])

        return list_possibilities



    def test_assembly_condition(self, number_planet, planetaries=[]):
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

                    if isinstance(node_2, (GearingPlanetary, GearingPlanet)):

                        basic_ratio_planet_1 = basic_ratio_planet_1*\
                        (-inv_list_nodes[j-1].Z/inv_list_nodes[j+1].Z)

                        if isinstance(inv_list_nodes[j+1], Planetary):

                            if inv_list_nodes[j+1].planetary_type == 'Ring':

                                basic_ratio_planet_1 = basic_ratio_planet_1 *-1



                for j, node_2 in enumerate(list_nodes_2[position_planet:]):

                    if isinstance(node_2, (GearingPlanetary, GearingPlanet)):

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

            if isinstance(relation, GearingPlanetary):

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

        system_matrix, vector_b = self.speed_system_equations()
        n_equations = len(self.relations)
        n_variables = len(self.elements)

        system_matrix_speed_solve_0 = npy.zeros((2, n_variables))
        vector_b_speed_solve_0 = npy.zeros(2)

        system_matrix = npy.concatenate((system_matrix, system_matrix_speed_solve_0), axis=0)
        vector_b = npy.concatenate((vector_b, vector_b_speed_solve_0), axis=0)
        impose_speeds = []

        for i, composant in enumerate(input_speeds_and_composants):
            impose_speeds.append(ImposeSpeed('ImposeSpeed'+str(i), composant,
                                             input_speeds_and_composants[composant]))

        for impose_speed in impose_speeds:
            position_element = self.matrix_position(impose_speed.node)
            system_matrix[n_equations][position_element] = impose_speed.speed_system_equations()[0]
            vector_b[n_equations] = impose_speed.speed_system_equations()[1]
            n_equations += 1


        solution = solve(system_matrix, vector_b)

        for i in range(len(self.elements)):
            self.elements[i].speed = solution[i]

        return solution

    def torque_system_equation(self):
        n_gearing_planetary = 0
        for gearing in self.gearings:
            if isinstance(gearing, GearingPlanetary):
                n_gearing_planetary += 1

        n_equations = (len(self.gearings)-n_gearing_planetary)+n_gearing_planetary*2
        n_variables = (len(self.gearings)-n_gearing_planetary)*2+n_gearing_planetary*3 + 1

        system_matrix = npy.zeros((n_equations, n_variables))
        rhs = npy.zeros(n_equations)
        num_element = 0
        num_equation = 0
        element_association = {}

        for element in self.elements:
           element_association[element] = []

        for gearing in self.gearings:
            matrix_relation, rhs_relation = gearing.torque_system_equations()

            if isinstance(gearing, GearingPlanetary):

                system_matrix[num_equation][num_element] = matrix_relation[0][0]
                system_matrix[num_equation][num_element+1] = matrix_relation[0][1]
                system_matrix[num_equation+1][num_element+2] = matrix_relation[1][0]
                system_matrix[num_equation+1][num_element+1] = matrix_relation[1][1]

                rhs[num_equation] = rhs_relation[0]
                rhs[num_equation+1] = rhs_relation[1]

                if isinstance(gearing.nodes[0], Planetary):

                    element_association[gearing.nodes[0]].append(num_element)
                    element_association[gearing.nodes[1]].append(num_element+1)

                else:
                    element_association[gearing.nodes[1]].append(num_element)
                    element_association[gearing.nodes[0]].append(num_element+1)

                element_association[self.planet_carrier].append(num_element+2)

                num_element += 3
                num_equation += 2

            else:

                system_matrix[num_equation][num_element] = matrix_relation[0]
                system_matrix[num_equation][num_element+1] = matrix_relation[1]

                rhs[num_equation] = rhs_relation[0]

                element_association[gearing.nodes[0]].append(num_element)
                element_association[gearing.nodes[1]].append(num_element+1)


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

    def torque_solve(self, input_torque_and_composant):

        system_matrix, vector_b, element_association = self.torque_system_equation()

        for composant in input_torque_and_composant:
            matrix_input = npy.zeros((1, system_matrix.shape[1]))
            matrix_input[0][element_association[composant]] = 1
            system_matrix = npy.concatenate((system_matrix, matrix_input))
            vector_b = npy.concatenate((vector_b, [input_torque_and_composant[composant]]))


        solution = solve(system_matrix, vector_b)

        torque_element_association = {}

        for element in element_association:
            max_solution = 0
            for position in element_association[element]:
                if abs(solution[position]) > abs(max_solution):

                    max_solution = solution[position]
            torque_element_association[element.name] = max_solution

        torque_element_association[self.planet_carrier.name] = solution[-1]
        return torque_element_association



class PlanetsStructure(DessiaObject):
    A = TypeVar('A', int, str)
    def __init__(self, name: str, architecture: List[List[str]], branch_connexions: List[List[A]]):
        self.architecture = architecture
        self.branch_connexions = branch_connexions
        self.planets = []
        self.connexions = []
        number_planet = 0
        self.gearings = []
        self.doubles = []
        DessiaObject.__init__(self, name=name)
        first_last_composant_branch = []
        for branch in self.architecture:

            flag_first_planet_branch = 0

            for planet in branch:
                number_planet += 1
                self.planets.append(Planet('Pl'+str(number_planet), planet, 7))

                if flag_first_planet_branch:

                    if planet == 'Double':
                        self.connexions.append([self.planets[-2], self.planets[-1], 'D'])
                    else:
                        self.connexions.append([self.planets[-2], self.planets[-1], 'GI'])

                else:
                    first_composant_branch = self.planets[-1]
                    flag_first_planet_branch = 1

            last_composant_branch = self.planets[-1]
            first_last_composant_branch.append([first_composant_branch, last_composant_branch])

        for branch_connexion in self.branch_connexions:

            if branch_connexion[2] == 'Simple':
                self.connexions.append([first_last_composant_branch[branch_connexion[0]-1][1],
                                        first_last_composant_branch[branch_connexion[1]-1][0], 'GI'])

            else:
                self.connexions.append([first_last_composant_branch[branch_connexion[0]-1][1],
                                        first_last_composant_branch[branch_connexion[1]-1][0], 'D'])


        for i, connexion in  enumerate(self.connexions):

          if connexion[2] != 'D':

             self.gearings.append(GearingPlanet('Gearing'+str(i), [connexion[0], connexion[1]]))

          else:
             self.doubles.append(Double('Double'+str(i), [connexion[0], connexion[1]]))

        self.relations = self.gearings + self.doubles

    def graph(self):

        graph_planetary_gear = nx.Graph()

        for relation in self.relations:

            graph_planetary_gear.add_edge(relation.nodes[0].name, relation.name)
            graph_planetary_gear.add_edge(relation.name, relation.nodes[1].name)

            nx.set_node_attributes(graph_planetary_gear, relation.nodes[0], relation.nodes[0].name)
            nx.set_node_attributes(graph_planetary_gear, relation.nodes[1], relation.nodes[1].name)
            nx.set_node_attributes(graph_planetary_gear, relation, relation.name)


        return graph_planetary_gear

    def plot(self):

        graph_planetary_gears = self.graph()
        plt.figure()
        nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)

    def path_planet_to_planet(self):
        graph_planetary_gears = self.graph()
        list_path = []

        for planet in self.planets[1:]:
            list_path.append(nx.shortest_path(graph_planetary_gears,
                                              self.planets[0].name, planet.name))

        for path in list_path:

            for i in range(len(path)):
                path[i] = nx.get_node_attributes(graph_planetary_gears, path[i])[path[i]]
        return list_path

    def gearing_chain_recursive_function(self, n, planet, graph_planetary_gear, possibilities, list_possibilities, previous_relation):

        planet_2 = copy.copy(planet)
        neighbors = nx.all_neighbors(graph_planetary_gear, planet.name)
        gearing = 0
        neighbors_2 = []
        for neighbor in neighbors:
            neighbor = nx.get_node_attributes(graph_planetary_gear, neighbor)[neighbor]
            neighbors_2.append(neighbor)
        neighbors = neighbors_2
        if not possibilities:
            possibilities.append(planet_2)
        n += 1
        if previous_relation:

            neighbors.remove(previous_relation)




        for neighbor in neighbors:


            possibilities_2 = copy.copy(possibilities)
            previous_relation = neighbor

            if isinstance(neighbor, GearingPlanet):
                gearing = 1

                if neighbor.nodes[0].name == planet_2.name:

                    possibilities_2.append(neighbor.nodes[1])
                    planet_3 = neighbor.nodes[1]

                else:
                    possibilities_2.append(neighbor.nodes[0])
                    planet_3 = neighbor.nodes[0]
                self.gearing_chain_recursive_function(n, planet_3, graph_planetary_gear, possibilities_2,
                                                      list_possibilities, previous_relation)

            elif isinstance(neighbor, Double):

                if neighbor.nodes[0].name == planet_2.name:

                    planet_3 = neighbor.nodes[1]

                else:
                    planet_3 = neighbor.nodes[0]

                self.gearing_chain_recursive_function(n, planet_3, graph_planetary_gear, [],
                                                      list_possibilities, previous_relation)



        if not gearing:
            list_possibilities.append(possibilities)
            return list_possibilities

        return list_possibilities

    def gearing_chain(self):
        graph_planetary_gear = self.graph()
        list_possibilities = self.gearing_chain_recursive_function(0, self.planets[0], graph_planetary_gear, [], [], 0)
        return list_possibilities



    def plot_cinematic_graph_gear(self, coordinate, lenght, diameter, diameter_pivot, lenght_pivot, color):
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

    def plot_cinematic_graph_double(self, coordinate, diameter, lenght, color):
        list_color = ['mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k',
                      'mediumblue', 'purple', 'green', 'k', 'mediumblue', 'purple', 'green', 'k']

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]+diameter/2, coordinate[1]+diameter/2])

        plt.plot(x, y, list_color[color])

        x = npy.array([coordinate[0], coordinate[0]+lenght])
        y = npy.array([coordinate[1]-diameter/2, coordinate[1]-diameter/2])

        plt.plot(x, y, list_color[color])

    def plot_cinematic_graph_planet_carrier(self, coordinates, planet_carrier_x, planet_carrier_y):
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





    def plot_cinematic_graph(self, lenght_gear=0.1, diameter_gear=1, lenght_double=2, diameter_pivot=0.2,
                             lenght_pivot=0.5, planet_carrier_x=2, planet_carrier_y=2):

        graph_path = self.path_planet_to_planet()

        plt.figure()

        previous_relation_double = []
        previous_relation_gearing = []

        previous_planet_gearing = []
        previous_planet_double = []
        inverse_relation_double = []
        inverse_relation_gearing = []
        coordinate_planet = [[0, 0]]
        coordinate = [0, 0]

        self.plot_cinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, 0)
        for path in graph_path:


            flag_way_inv_gearing = 0
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

                            self.plot_cinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                        else:

                            self.plot_cinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                            coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]
                            inverse_relation_double.append(element)




                    else:

                        if not element in previous_relation_double:

                            if previous_relation_double and previous_relation_double[-1] in inverse_relation_double:

                                self.plot_cinematic_graph_double(coordinate, diameter_pivot, -lenght_double/(1+i*0.2), color)
                                inverse_relation_double.append(element)
                                coordinate = [coordinate[0]-lenght_double/(1+i*0.2), coordinate[1]]


                            else:
                                self.plot_cinematic_graph_double(coordinate, diameter_pivot, +lenght_double/(1+i*0.2), color)
                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]
                        else:

                                coordinate = [coordinate[0]+lenght_double/(1+i*0.2), coordinate[1]]

                    previous_relation_double.append(element)
                    previous_planet_double.extend([element.nodes[0], element.nodes[1]])

                elif isinstance(element, GearingPlanet):
                    color += 1
                    if element in  inverse_relation_gearing:

                        coordinate = [coordinate[0], coordinate[1]-diameter_gear]


                    elif ((element.nodes[0] in previous_planet_gearing or element.nodes[1] in previous_planet_gearing) \
                    and not element  in previous_relation_gearing):

                        for gearing in previous_relation_gearing:

                            for node in gearing.nodes:

                                if element.nodes[0] == node or element.nodes[1] == node:

                                    if  not gearing == previous_element:

                                        if gearing in inverse_relation_gearing:

                                            flag_way_inv_gearing = 1
                                        else:
                                            flag_way_inv_gearing = 0
                                    else:

                                        if not gearing in inverse_relation_gearing:
                                            flag_way_inv_gearing = 1
                                        else:
                                            flag_way_inv_gearing = 0

                        if flag_way_inv_gearing:
                            coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                        else:

                            coordinate = [coordinate[0], coordinate[1]-diameter_gear]
                            inverse_relation_gearing.append(element)




                    else:


                        coordinate = [coordinate[0], coordinate[1]+diameter_gear]

                    previous_relation_gearing.append(element)
                    previous_planet_gearing.extend([element.nodes[0], element.nodes[1]])

                if not isinstance(element, Planet):

                    if  not coordinate in coordinate_planet:
                        coordinate_planet.append(coordinate)
                        self.plot_cinematic_graph_gear(coordinate, lenght_gear, diameter_gear, diameter_pivot, lenght_pivot, color)

                    previous_element = element


        self.plot_cinematic_graph_planet_carrier(coordinate_planet, planet_carrier_x, planet_carrier_y)



class OptimizerPlanetStructure(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True

    def __init__(self, number_max_planet: int, number_junction: int, number_max_junction_by_planet: int, min_planet_branch: int):
        self.number_max_planet = number_max_planet
        self.number_junction = number_junction
        self.number_max_junction_by_planet = number_max_junction_by_planet
        self.min_planet_branch = min_planet_branch
        DessiaObject.__init__(self, name='Optimizer')

    ## Recursive Function which give all the possibilities  of planet_type in a branch for a Planet number fixed ##
    ## Exemple Input:(0,[],[],4) -> Output:[['Simple', 'Double', 'Double', 'Simple'], ['Simple', 'Double', 'Double', 'Double'], ['Simple', 'Simple', 'Double', 'Double'],
    ##['Simple', 'Simple', 'Simple', 'Simple'], ['Double', 'Double', 'Simple', 'Simple'], ['Double', 'Double', 'Double', 'Simple'], ['Double', 'Double', 'Double', 'Double']] ##
    def list_possibilities_planets_type(self, n, list_planet_type, planet_type, number_max_planet):

        if n == number_max_planet:
            planet_type_2 = copy.copy(planet_type)
            list_planet_type.append(planet_type_2)

            return list_planet_type

        if not planet_type:
            planet_type_1 = copy.copy(planet_type)
            planet_type_1.append('Simple')

            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_1, number_max_planet)
        else:

            planet_type_1 = copy.copy(planet_type)

            planet_type_2 = copy.copy(planet_type)
            planet_type_1.append('Simple')

            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_1, number_max_planet)

            planet_type_2.append('Double')
            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_2, number_max_planet)

        return list_planet_type


    ## Recursive Function which give all the possibilities for a junction number fixed ##
    ## Exemple Input:([],[0,0,0,0],2,2,[2,2,2,2]) -> Output:[[2, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 2, 0, 0],
    ## [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 2, 0], [0, 0, 1, 1], [0, 0, 0, 2]] ##
    def list_possibilities_junction(self, list_possibilities, possibilitie, number_junction_recursive_fonction,
                                    sum_number_junction_max_by_planet):

        number_junction_recursive_fonction_2 = copy.copy(number_junction_recursive_fonction)
        possibilitie_2 = copy.copy(possibilitie)
        sum_number_junction_max_by_planet_2 = copy.copy(sum_number_junction_max_by_planet)
        number_junction_max = self.number_junction
        if number_junction_recursive_fonction == number_junction_max:

            if not possibilitie_2  in list_possibilities:
                list_possibilities.append(possibilitie_2)

            return list_possibilities


        for i in range(1, len(possibilitie)):

            if i > 1 and flag:
                possibilitie_2[i-1] -= 1
                sum_number_junction_max_by_planet_2[i] = sum_number_junction_max_by_planet[i]

            flag = 1

            if possibilitie_2[i] < sum_number_junction_max_by_planet_2[i]:

                if i+1 < len(possibilitie)-1:
                    sum_number_junction_max_by_planet_2[i+1] += self.number_junction_max_by_planet

                possibilitie_2[i] += 1

                self.list_possibilities_junction(list_possibilities, possibilitie_2, number_junction_recursive_fonction_2+1,
                                                 sum_number_junction_max_by_planet_2)

            else:
                flag = 0

        return list_possibilities


    ## Recursive Function which give all the possibilities for a limited number of connexion( number_max_connextion) in a list_connexion ##
    ## Exemple Input:(0 ,[], [[0, 0], [0, 0]] ,[[2, 4], [2, 5], [3, 6]],2,2) -> Output:[[[2, 4], [2, 5]], [[2, 4], [3, 6]], [[2, 5], [3, 6]]] ##
    def list_possibilities_connexion_branch_step_1(self, n, list_possibilities, possibilitie,
                                                   list_connexion_branch, number_max_connexion):

        possibilitie_2 = copy.copy(possibilitie)
        if  n == number_max_connexion:

            if not possibilitie in list_possibilities:
                flag_similaritie = 0

                for possibilitie_3 in list_possibilities:
                    similaritie = 0

                    for connexion in possibilitie_2:

                        if connexion in possibilitie_3:
                            similaritie += 1

                    if similaritie == number_max_connexion:
                        flag_similaritie = 1
                        break

                if not flag_similaritie:
                    list_possibilities.append(possibilitie_2)


        else:
            n += 1
            for element_1 in list_connexion_branch:

                number_connexion_by_branch = 0
                flag_branch = 1

                for element_2 in possibilitie_2:

                    if element_2[0] == element_1[0]:
                        number_connexion_by_branch += 1

                    if element_2[1] == element_1[1]:
                        flag_branch = 0

                if number_connexion_by_branch <= self.number_max_junction_by_planet and flag_branch:
                    possibilitie_2[n-1] = element_1

                    self.list_possibilities_connexion_branch_step_1(n, list_possibilities, possibilitie_2,
                                                                    list_connexion_branch, number_max_connexion)
                    possibilitie_2[n-1] = possibilitie[n-1]
        return list_possibilities









    ## Recursive Function which give all the possibilities  of connexion  for a number_junction fixed and a number branch fixed ##
    ## Exemple Input:(0 ,[0,0],[],[],[2,3,4,5,6],[1],[2,2],2) -> Output:[[[[1, 2], [1, 3]], [[2, 4], [2, 5]]], [[[1, 2], [1, 3]], [[2, 4], [3, 5]]], [[[1, 2], [1, 3]], [[2, 5], [3, 4]]], [[[1, 2], [1, 3]], [[3, 4], [3, 5]]]] ##
    def list_possibilities_connexion_branch_step_2(self, n, possibilitie, list_possibilities,
                                                   list_previous_branch, remaning_branch, previous_branch,
                                                   number_junction_branch):

        possibilitie_2 = copy.copy(possibilitie)
        remaning_branch_2 = copy.copy(remaning_branch)
        previous_branch_2 = copy.copy(previous_branch)

        if n == len(number_junction_branch):
            list_possibilities.append(possibilitie_2)
            list_previous_branch.append(previous_branch_2)
            return None

        number_connexion = number_junction_branch[n]
        n += 1

        list_possibilities_connexion = []
        for element in previous_branch_2:
            for j in range(number_connexion):
                list_possibilities_connexion.append([element, remaning_branch_2[j]])

        list_possibilities_connexion_2 = []
        possibilitie_3 = []

        for i in range(number_connexion):

            possibilitie_3.append([0, 0])

        self.list_possibilities_connexion_branch_step_1(0, list_possibilities_connexion_2, possibilitie_3,
                                                        list_possibilities_connexion, number_connexion)



        for element_1 in list_possibilities_connexion_2:
            remaning_branch_3 = copy.copy(remaning_branch_2)
            previous_branch_3 = copy.copy(previous_branch_2)

            for i in range(number_connexion):
                previous_branch_3.append(remaning_branch_2[i])
                remaning_branch_3.remove(remaning_branch_2[i])

            possibilitie_2[n-1] = element_1

            for element_2 in element_1:

                if element_2[0] in previous_branch_3:
                    previous_branch_3.remove(element_2[0])

            self.list_possibilities_connexion_branch_step_2(n, possibilitie_2, list_possibilities, list_previous_branch,
                                                            remaning_branch_3, previous_branch_3, number_junction_branch)

        return list_possibilities


    ## Recursive Function which give all the possibilities to distribute a number planet fixed in a number_branch fixed  ##
    ## Exemple Input:([0,0,0],[],6,0,3,0,1) -> Output:[[1, 1, 4], [1, 2, 3], [1, 3, 2], [1, 4, 1], [2, 1, 3], [2, 2, 2], [2, 3, 1], [3, 1, 2], [3, 2, 1], [4, 1, 1]] ##
    def list_possibilities_planets_by_branch_step_1(self, possibilities, list_possibilities,
                                                    number_planet, number_branch, number_branch_max,
                                                    number_other_planet):


        possibilities_2 = copy.copy(possibilities)
        number_branch += 1

        if number_branch == number_branch_max:

            if number_planet-number_other_planet > 0:
                possibilities_2[number_branch-1] = number_planet-number_other_planet
                list_possibilities.append(possibilities_2)

            return list_possibilities


        for i in range(self.min_planet_branch, number_planet-number_other_planet-(number_branch_max-number_branch-1)):

            possibilities_2[number_branch-1] = i
            self.list_possibilities_planets_by_branch_step_1(possibilities_2, list_possibilities,
                                                             number_planet, number_branch, number_branch_max, number_other_planet+i)



    ##  Function which give all the possibilities of branch for a possibility junction  ##
    ## Exemple Input:([0,1,1],7,3,1) -> Output:[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]],
    ##[[2, 1, 1, 2, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]], [[2, 1, 2, 1, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]]etc... ##
    def list_possibilities_planets_by_branch_step_2(self, junction):

        list_number_planet_branch = []
        number_planet_branch = 0
        list_number_junctions = []
        number_branch = 1

        for i in range(len(junction)):

            number_planet_branch += 1

            if junction[i]:

                number_branch += junction[i]+1
                list_number_planet_branch.append(number_planet_branch)
                number_planet_branch = 0
                list_number_junctions.append(junction[i]+1)

        list_connexion = []
        connexion = []
        remaning_branch = []
        previous_branch = [1]
        list_number_planet = [0]

        for i in range(number_branch-1):
            remaning_branch.append(i+2)
            list_number_planet.append(0)

        for i in range(len(list_number_junctions)):
            connexion.append(0)

        list_previous_branch = []

        self.list_possibilities_connexion_branch_step_2(0, connexion, list_connexion, list_previous_branch, remaning_branch,
                                                        previous_branch, list_number_junctions)

        list_architecture = []
        for i, connexion_2 in enumerate(list_connexion):

            number_planet_know = 0
            previous_branch_final = list_previous_branch[i]
            previous_connexion_by_branch = []

            for i, connexion_by_branch in enumerate(connexion_2):

                for element in connexion_by_branch:

                    if not element[0] in previous_connexion_by_branch:

                        list_number_planet[element[0]-1] = list_number_planet_branch[i]
                        number_planet_know += list_number_planet_branch[i]
                        previous_connexion_by_branch.append(element[0])



            number_planet_unknow = self.number_max_planet-number_planet_know
            possibilitie = []

            for i in range(len(previous_branch_final)):
                possibilitie.append(0)

            list_planet_by_branch = []

            self.list_possibilities_planets_by_branch_step_1(possibilitie, list_planet_by_branch, number_planet_unknow, 0,
                                                             len(previous_branch_final), 0)

            for planet_by_branch in list_planet_by_branch:

                for i, branch in enumerate(previous_branch_final):

                    list_number_planet[branch-1] = planet_by_branch[i]

                list_architecture.append([copy.copy(list_number_planet), connexion_2])

        return list_architecture, number_branch



    ##  Recursive Function which give all the possibilities of architecure (branch and planet_type) for a possibility junction  ##
    ## Exemple Input:(1,5,[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]],[0,0,0,0,0],[],[],[]]) -> Output: list of PlanetsStructure objet ##
    def list_possibilities_architecture_planet(self, branch, branch_maximum, architecture, possibilities,
                                               list_possibilities, possibilities_connexion, list_possibilities_connexion):



        possibilities_connexion_2 = copy.copy(possibilities_connexion)
        number_planet = architecture[0][branch-1]

        list_planet_possibilities = []

        if number_planet == 1:

            for branch_connexion in possibilities_connexion_2:
                if branch_connexion[1] == branch:
                    list_planet_possibilities = [[branch_connexion[2]]]
                    break

        else:
            self.list_possibilities_planets_type(0, list_planet_possibilities, [], number_planet)


        number_connexion = 0
        for connexion in architecture[1]:
            for branch_2 in connexion:
                if branch_2[0] == branch:
                    number_connexion += 1

                    possibilities_connexion_2.append([branch, branch_2[1], 'Simple'])



        for planet_possibilities in list_planet_possibilities:

            possibilities_2 = copy.copy(possibilities)
            possibilities_2[branch-1] = planet_possibilities

            if branch == branch_maximum:

                possibilities_connexion_3 = copy.deepcopy(possibilities_connexion_2)
                list_possibilities.append(PlanetsStructure('PlanetsStructure', possibilities_2, possibilities_connexion_3))

            else:

                if number_connexion == 0:
                    self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                list_possibilities, possibilities_connexion_2,
                                                                list_possibilities_connexion)



                elif number_connexion == 2:

                    possibilities_connexion_2[-1][2] = 'Double'
                    possibilities_connexion_2[-2][2] = 'Simple'

                    self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                list_possibilities, possibilities_connexion_2,
                                                                list_possibilities_connexion)

                    possibilities_connexion_2[-1][2] = 'Simple'
                    possibilities_connexion_2[-2][2] = 'Double'
                    self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                list_possibilities, possibilities_connexion_2,
                                                                list_possibilities_connexion)


                    if planet_possibilities[-1] == 'Simple':

                        possibilities_connexion_2[-1][2] = 'Double'
                        possibilities_connexion_2[-2][2] = 'Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                    elif planet_possibilities[-1] == 'Double':

                        possibilities_connexion_2[-1][2] = 'Simple'
                        possibilities_connexion_2[-2][2] = 'Simple'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)



                elif number_connexion == 3:

                    if planet_possibilities[-1] == 'Simple':

                        possibilities_connexion_2[-1][2] = 'Double'
                        possibilities_connexion_2[-2][2] = 'Double'
                        possibilities_connexion_2[-3][2] = 'Simple'

                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                        possibilities_connexion_2[-1][2] = 'Simple'
                        possibilities_connexion_2[-2][2] = 'Double'
                        possibilities_connexion_2[-3][2] = 'Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                        possibilities_connexion_2[-1][2] = 'Double'
                        possibilities_connexion_2[-2][2] = 'Simple'
                        possibilities_connexion_2[-3][2] = 'Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                    elif planet_possibilities[-1] == 'Double':

                        possibilities_connexion_2[-1][2] = 'Simple'
                        possibilities_connexion_2[-2][2] = 'Simple'
                        possibilities_connexion_2[-3][2] = 'Double'

                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                        possibilities_connexion_2[-1][2] = 'Simple'
                        possibilities_connexion_2[-2][2] = 'Double'
                        possibilities_connexion_2[-3][2] = 'Simple'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)

                        possibilities_connexion_2[-1][2] = 'Double'
                        possibilities_connexion_2[-2][2] = 'Simple'
                        possibilities_connexion_2[-3][2] = 'Simple'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2,
                                                                    list_possibilities_connexion)


    def solution_sort_recursive_function(self,new_first_node,first_node_check,new_graph,graph_check,previous_node,previous_node_check):
        list_number_false=[]
        valid=False

        for new_neighbor in nx.neighbors(new_graph,new_first_node):
            if new_neighbor != previous_node:
                number_false=0
                if len(list( nx.neighbors(graph_check,first_node_check))) != len(list( nx.neighbors(new_graph,new_first_node))):
                     
                    valid=True
                    return valid
                
                for neighbor_check in nx.neighbors(graph_check,first_node_check):
                    if len(previous_node_check)<2 or neighbor_check != previous_node_check[-2]:
                        if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor]) == type(nx.get_node_attributes(graph_check,neighbor_check)[neighbor_check]):
                            valid=self.solution_sort_recursive_function(new_neighbor,neighbor_check,new_graph,graph_check,new_first_node,previous_node_check + [neighbor_check])
                            
                            
                        else:
                            number_false+=1
                        
                        
                                
                if number_false == len(list(nx.neighbors(graph_check,first_node_check)))-(len(previous_node_check)>1):
                    
                    valid=True
                    return valid
                
                else:
                    list_number_false.append([number_false,type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor])])

        sum_number_false=0
        list_previous_type_false=[]
        for number_false in list_number_false:
            if not number_false[1] in list_previous_type_false:
                sum_number_false += number_false[0]
                list_previous_type_false.append(number_false[1])

        if sum_number_false!=0 and sum_number_false!=len(list_number_false):
            valid=True
            
            return valid
        
        return valid
                


    def solution_sort(self, new_planet_structure, planet_structures_check):
        
        new_graph = new_planet_structure.graph()
        
        for node in nx.nodes(new_graph):
            
            if len(list(nx.neighbors(new_graph,node))) == 1:
                first_node=node
                break
       
        for planet_structure in planet_structures_check:
            
            
            graph_check= planet_structure.graph()
            possible_first_node_check=[]
            for node in nx.nodes(graph_check):
                
                if len(list(nx.neighbors(graph_check,node))) == 1 and type(nx.get_node_attributes(graph_check,node)[node]) == type(nx.get_node_attributes(new_graph,first_node)[first_node]):
                    possible_first_node_check.append(node)
                    
                    
            for node in possible_first_node_check:
                valid=self.solution_sort_recursive_function(first_node,node,new_graph,graph_check,first_node,[node])
                
                if  not valid:
                    return False
                
        return True           
            
        
    
    def decision_tree(self):
        tree = dt.DecisionTree()


        list_possibilities_junction = []
        list_planet = []
        sum_number_max_junction_by_planet = []
        list_solution = []
        for i in range(self.number_max_planet-2):
            list_planet.append(0)
            sum_number_max_junction_by_planet.append(self.number_max_junction_by_planet)

        self.list_possibilities_junction(list_possibilities_junction, list_planet, 0, sum_number_max_junction_by_planet)

        tree.SetCurrentNodeNumberPossibilities(len(list_possibilities_junction))
        node = tree.NextNode(True)
        

        while not tree.finished:


            if len(node) == 1:

                list_junction = list_possibilities_junction[node[0]]

                list_global_architecture, number_branch = self.list_possibilities_planets_by_branch_step_2(list_junction)

                tree.SetCurrentNodeNumberPossibilities(len(list_global_architecture))



            if len(node) == 2:

                global_architecture = list_global_architecture[node[1]]
                list_branch = []

                for i in range(number_branch):
                    list_branch.append(0)


                list_connexion = []
                list_possibilities = []
                self.list_possibilities_architecture_planet(1, number_branch, global_architecture, list_branch,
                                                            list_possibilities, [], list_connexion)

               
                tree.SetCurrentNodeNumberPossibilities(len(list_possibilities))
            if len(node)==3:
                planet_structure=list_possibilities[node[2]]
                # planet_structure.plot_cinematic_graph()
                if self.solution_sort(planet_structure,list_solution):
                    list_solution.append(planet_structure)
                tree.SetCurrentNodeNumberPossibilities(0)
            node = tree.NextNode(True)

        return list_solution

class OptimizerPlanetaryGearsArchitecture(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True
    def __init__(self, planet_structures: PlanetsStructure, input_speeds: List[List[float]]):


        self.planet_structures = planet_structures
        self.number_input = len(input_speeds)
        DessiaObject.__init__(self, name='Optimizer')

    def list_possibilities_planetaries_recursive_function(self, n, planetaries, possibilitie, list_possibilities, gearing_chains,
                                                          gearing_chains_planetary_type, gearing_chains_occupation,
                                                          gearing_chains_planet_index):


        possibilitie_2 = copy.copy(possibilitie)
        list_planetary_type = [['Ring', 'GI', -1], ['Sun', 'GE', 1]]
        gearing_chains_2 = copy.copy(gearing_chains)
        number = 0



        if n == len(planetaries):

            if not 0 in gearing_chains_occupation:
                list_possibilities.append(possibilitie_2)

            return list_possibilities

        planetary = planetaries[n]

        for i, gearing_chain in enumerate(gearing_chains_2):
            number = copy.copy(i)

            if gearing_chains_occupation[i] < 2:

                gearing_chains_occupation_2 = copy.copy(gearing_chains_occupation)
                gearing_chains_occupation_2[i] += 1



                if gearing_chains_planet_index[i] == 0:

                    gearing_chains_planet_index_2 = copy.copy(gearing_chains_planet_index)
                    gearing_chains_planet_index_2[i] = 1

                    for planetary_type in list_planetary_type:
                        if not gearing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            gearing_chain_2 = copy.copy(gearing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            gearing_chains_planetary_type_2 = copy.copy(gearing_chains_planetary_type)
                            gearing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, gearing_chain_2[0], planetary_type[1], number]
                            self.list_possibilities_planetaries_recursive_function(n+1, planetaries, possibilitie_2,
                                                                                   list_possibilities, gearing_chains,
                                                                                   gearing_chains_planetary_type_2,
                                                                                   gearing_chains_occupation_2,
                                                                                   gearing_chains_planet_index_2)


                    gearing_chains_planet_index_2[i] = -1


                    for planetary_type in list_planetary_type:

                        if not gearing_chains_planetary_type[i] == planetary_type[0]:
                            planetary_2 = copy.copy(planetary)
                            gearing_chain_2 = copy.copy(gearing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            gearing_chains_planetary_type_2 = copy.copy(gearing_chains_planetary_type)
                            gearing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, gearing_chain_2[-1], planetary_type[1], number]
                            self.list_possibilities_planetaries_recursive_function(n+1, planetaries, possibilitie_2,
                                                                                   list_possibilities, gearing_chains,
                                                                                   gearing_chains_planetary_type_2,
                                                                                   gearing_chains_occupation_2,
                                                                                   gearing_chains_planet_index_2)



                elif gearing_chains_planet_index[i] == 1:

                    gearing_chains_planet_index_2 = copy.copy(gearing_chains_planet_index)
                    gearing_chains_planet_index_2[i] = -1

                    for planetary_type in list_planetary_type:

                        if not gearing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            gearing_chain_2 = copy.copy(gearing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            gearing_chains_planetary_type_2 = copy.copy(gearing_chains_planetary_type)
                            gearing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, gearing_chain_2[-1], planetary_type[1], number]
                            self.list_possibilities_planetaries_recursive_function(n+1, planetaries, possibilitie_2,
                                                                                   list_possibilities, gearing_chains,
                                                                                   gearing_chains_planetary_type_2,
                                                                                   gearing_chains_occupation_2,
                                                                                   gearing_chains_planet_index_2)


                else:
                    gearing_chains_planet_index_2 = copy.copy(gearing_chains_planet_index)
                    gearing_chains_planet_index_2[i] = -1

                    for planetary_type in list_planetary_type:
                        if not gearing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            gearing_chain_2 = copy.copy(gearing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            gearing_chains_planetary_type_2 = copy.copy(gearing_chains_planetary_type)
                            gearing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, gearing_chain_2[0], planetary_type[1], number]
                            self.list_possibilities_planetaries_recursive_function(n+1, planetaries, possibilitie_2,
                                                                                   list_possibilities, gearing_chains,
                                                                                   gearing_chains_planetary_type_2,
                                                                                   gearing_chains_occupation_2,
                                                                                   gearing_chains_planet_index_2)

        return list_possibilities



    def list_possibilities_planetaries(self, planetaries, planets_structure, planet_carrier):
        gearing_chains = planets_structure.gearing_chain()



        if len(planetaries) < len(gearing_chains):
            return []

        possibilitie = []
        connexion = []
        gearing_chains_planetary_type = []
        gearing_chains_planet_index = []
        gearing_chains_occupation = []

        for i in enumerate(planetaries):
            possibilitie.append(0)

        for i in enumerate(gearing_chains):
            gearing_chains_planetary_type.append(0)
            gearing_chains_occupation.append(0)
            gearing_chains_planet_index.append(0)

        list_possibilities = self.list_possibilities_planetaries_recursive_function(0, planetaries, possibilitie, [],
                                                                                    gearing_chains, gearing_chains_planetary_type,
                                                                                    gearing_chains_occupation,
                                                                                    gearing_chains_planet_index)


        for double in planets_structure.doubles:

            connexion.append([double.nodes[0], double.nodes[1], 'D'])

        list_solution = []

        for i, possibilitie_2 in enumerate(list_possibilities):

            connexion_2 = copy.copy(connexion)
            previous_planetary = []

            for planetary_connexion in possibilitie_2:
                connexion_2.append(planetary_connexion[:-1])

                if not planetary_connexion[3]  in previous_planetary:

                    gearing_chain = gearing_chains[planetary_connexion[3]]
                    previous_planetary.append(planetary_connexion[3])

                    if len(gearing_chain) > 1:

                        if gearing_chain[0] == planetary_connexion[1]:
                            for planet_1, planet_2 in zip(gearing_chain[:-1], gearing_chain[1:]):

                                connexion_2.append([copy.deepcopy(planet_1), copy.deepcopy(planet_2), planetary_connexion[2]])


                        elif gearing_chain[-1] == planetary_connexion[1]:

                           gearing_chain_2 = gearing_chain[::-1]

                           for planet_1, planet_2 in zip(gearing_chain_2[:-1], gearing_chain_2[1:]):

                                connexion_2.append([copy.deepcopy(planet_1), copy.deepcopy(planet_2), planetary_connexion[2]])


            planetary_gear = PlanetaryGear('PlanetaryGear'+str(i), planetaries, planets_structure.planets, planet_carrier, connexion_2)
            list_path = planetary_gear.path_planetary_to_planetary()
            list_planets = []

            for planet in planetary_gear.planets:
                list_planets.append(planet.name)

            for path in list_path:
                for element in path:

                    if element.name in list_planets:
                        list_planets.remove(element.name)


            if not list_planets:

                list_solution.append(PlanetaryGear('PlanetaryGear'+str(i), copy.deepcopy(planetaries),
                                                   copy.deepcopy(planets_structure.planets), copy.copy(planet_carrier), connexion_2))



        return list_solution


    def solution_sort_recursive_function(self,new_first_node,first_node_check,new_graph,graph_check,previous_node,previous_node_check):
        list_number_false=[]
        valid=False

        for new_neighbor in nx.neighbors(new_graph,new_first_node):
            if new_neighbor != previous_node:
                number_false=0
                if len(list( nx.neighbors(graph_check,first_node_check))) != len(list( nx.neighbors(new_graph,new_first_node))):
                     
                    valid=True
                    return valid
                
                for neighbor_check in nx.neighbors(graph_check,first_node_check):
                    if len(previous_node_check)<2 or neighbor_check != previous_node_check[-2]:
                        if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor]) == type(nx.get_node_attributes(graph_check,neighbor_check)[neighbor_check]):
                            
                            if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor]) == Planetary:
                                
                                if nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor].planetary_type == nx.get_node_attributes(graph_check,neighbor_check)[neighbor_check].planetary_type:
                                    
                                    valid=self.solution_sort_recursive_function(new_neighbor,neighbor_check,new_graph,graph_check,new_first_node,previous_node_check + [neighbor_check])
                                
                                else:
                                    
                                    number_false+=1
                            
                            else:
                                valid=self.solution_sort_recursive_function(new_neighbor,neighbor_check,new_graph,graph_check,new_first_node,previous_node_check + [neighbor_check])
                            
                            
                        else:
                            number_false+=1
                        
                        
                        if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor]) == Planetary:
                            print(number_false)
                if number_false == len(list(nx.neighbors(graph_check,first_node_check)))-(len(previous_node_check)>1):
                    if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor]) == Planetary:
                            print('a')
                    valid=True
                    return valid
                
                else:
                    if type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor])== Planetary:
                            list_number_false.append([number_false,nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor].planetary_type])
                    else:
                        list_number_false.append([number_false,type(nx.get_node_attributes(new_graph,new_neighbor)[new_neighbor])])
        # print(list( nx.neighbors(graph_check,first_node_check)))
        # print(list( nx.neighbors(new_graph,new_first_node)))
        # print(list_number_false)
         
        sum_number_false=0
        list_previous_type_false=[]
        for number_false in list_number_false:
            if number_false[1]=='Sun':
                print(number_false[0])
            if not number_false[1] in list_previous_type_false:
                sum_number_false += number_false[0]
                list_previous_type_false.append(number_false[1])

        if sum_number_false!=0 and sum_number_false!=len(list_number_false):
            valid=True
            print(valid)
            return valid
        
        return valid
                


    def solution_sort(self, new_planetary_gear, planetary_gears_check):
        
        new_graph = new_planetary_gear.graph()
        new_graph.remove_node(new_planetary_gear.planet_carrier.name)
        for i,planet in enumerate(new_planetary_gear.planets):
            new_graph.remove_node('Pv'+str(i))
        for node in nx.nodes(new_graph):
            
            if len(list(nx.neighbors(new_graph,node))) == 1:
                first_node=node
                break
       
        for planetary_gear in planetary_gears_check:
            
            
            graph_check= planetary_gear.graph()
            graph_check.remove_node(planetary_gear.planet_carrier.name)
            for i,planet in enumerate(planetary_gear.planets):
                graph_check.remove_node('Pv'+str(i))
            possible_first_node_check=[]
            for node in nx.nodes(graph_check):
                
                if len(list(nx.neighbors(graph_check,node))) == 1 and type(nx.get_node_attributes(graph_check,node)[node]) == type(nx.get_node_attributes(new_graph,first_node)[first_node]):
                    possible_first_node_check.append(node)
                    
                    
            for node in possible_first_node_check:
                valid=self.solution_sort_recursive_function(first_node,node,new_graph,graph_check,first_node,[node])
                
                if  not valid:
                    return False
                
        return True           

    def decision_tree(self):
        tree = dt.DecisionTree()


        planet_carrier = PlanetCarrier('PlanetCarrier')
        planetaries = []
        n = 0
        list_solution = []
        list_paths_type = []

        for i in range(self.number_input-1):
            planetaries.append(Planetary('Planetary_'+str(i), 7, 'Sun'))

        tree.SetCurrentNodeNumberPossibilities(len(self.planet_structures))
        node = tree.NextNode(True)


        while not tree.finished:
            valid = True

            if len(node) == 1:
             planet_architecture = self.planet_structures[node[0]]
             list_planetary_gears = self.list_possibilities_planetaries(planetaries, planet_architecture, planet_carrier)
             tree.SetCurrentNodeNumberPossibilities(len(list_planetary_gears))

            if len(node) == 2:

                planetary_gear = list_planetary_gears[node[1]]


                if self.solution_sort(planetary_gear,list_solution):
                    list_solution.append(planetary_gear)

                tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(valid)

        return list_solution

class OptimizerPlanetaryGearsZNumber(DessiaObject):
    _standalone_in_db = True

    _generic_eq = True

    def __init__(self, planetary_gear: PlanetaryGear, input_speeds: List[List[float]], diameter_cylinder: float, lenght_cylinder: float, Z_range_sun: List[int], Z_range_ring: List[int], number_planet: int):
        self.planetary_gear = planetary_gear
        self.input_speeds = input_speeds
        self.diameter = diameter_cylinder
        self.lenght = lenght_cylinder
        self.number_input = len(input_speeds)
        self.Z_range_sun = Z_range_sun
        self.Z_range_ring = Z_range_ring
        self.number_planet = number_planet
        DessiaObject.__init__(self, name='Optimizer')

    def multiplication_possibility_speed(self, list_1, n, element_multiplication, list_multiplication):


        for i in range(len(list_1)):

            element_multiplication_2 = copy.copy(element_multiplication)

            if  not list_1[i] in element_multiplication_2:
                element_multiplication_2.append(list_1[i])

                if n != len(list_1)-1:

                    self.multiplication_possibility_speed(list_1, n+1, element_multiplication_2, list_multiplication)

                else:
                    if not element_multiplication_2 in list_multiplication:

                        list_multiplication.append(element_multiplication_2)

        return list_multiplication


    def decision_tree_architecture(self, number_max_planet, number_junction, number_max_junction_by_planet, min_planet_branch):
        tree = dt.DecisionTree()


        planet_carrier = PlanetCarrier('PlanetCarrier')
        planetaries = []
        n = 0
        list_solution = []
        list_paths_type = []

        for i in range(self.number_input-1):
            planetaries.append(Planetary('Planetary_'+str(i), 7, 'Sun'))

        list_possibilities_junction = []
        list_planet = []
        sum_number_max_junction_by_planet = []

        for i in range(number_max_planet-2):
            list_planet.append(0)
            sum_number_max_junction_by_planet.append(number_max_junction_by_planet)

        self.list_possibilities_junction(list_possibilities_junction, list_planet, 0, number_junction,
                                         number_max_junction_by_planet, sum_number_max_junction_by_planet)

        tree.SetCurrentNodeNumberPossibilities(len(list_possibilities_junction))
        node = tree.NextNode(True)
        valid = True

        while not tree.finished:

            valid = True
            if len(node) == 1:

                list_junction = list_possibilities_junction[node[0]]

                list_global_architecture, number_branch = self.list_possibilities_planets_by_branch_step_2(list_junction,
                                                                                                           number_max_planet,
                                                                                                           number_max_junction_by_planet,
                                                                                                           min_planet_branch)

                tree.SetCurrentNodeNumberPossibilities(len(list_global_architecture))



            if len(node) == 2:

                global_architecture = list_global_architecture[node[1]]
                list_branch = []

                for i in range(number_branch):
                    list_branch.append(0)

                list_possibilitie_planet_architecture = []
                list_connexion = []
                self.list_possibilities_architecture_planet(1, number_branch, global_architecture, list_branch,
                                                            list_possibilitie_planet_architecture, [], list_connexion)

                tree.SetCurrentNodeNumberPossibilities(len(list_possibilitie_planet_architecture))



            if len(node) == 3:

                planet_architecture = list_possibilitie_planet_architecture[node[2]]

                list_possibilities = self.list_possibilities_planetaries(planetaries, planet_architecture, planet_carrier)
                tree.SetCurrentNodeNumberPossibilities(len(list_possibilities))


            if len(node) == 4:
                valid = True
                planetary_gear = list_possibilities[node[3]]

                paths_solution = planetary_gear.path_planetary_to_planetary_type()

                for paths in list_paths_type:

                    number_similaritie = 0

                    for path in paths:

                        for path_solution in paths_solution:
                            if path_solution == path:
                             number_similaritie += 1

                    if number_similaritie == len(paths):
                        valid = False
                        break


                if valid:
                    list_paths_type.append(paths_solution)

                tree.SetCurrentNodeNumberPossibilities(self.number_input)


            if len(node) == 5:

                planetary_gear.planet_carrier.speed_input = self.input_speeds[node[4]]

                input_speeds_2 = copy.copy(self.input_speeds)
                input_speeds_2.remove(self.input_speeds[node[4]])

                list_possibility_speed = []
                self.multiplication_possibility_speed(input_speeds_2, 0, [], list_possibility_speed)

                tree.SetCurrentNodeNumberPossibilities(len(list_possibility_speed))


            if len(node) == 6:

                possibility_speed = list_possibility_speed[node[5]]

                for i, planetarie in enumerate(planetary_gear.planetaries):
                    planetarie.speed_input = possibility_speed[i]

                list_solution.append(copy.deepcopy(planetary_gear))
                tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(valid)
        return list_solution


    def test_GCD(self, Z_1, Z_2):

        if m.gcd(Z_1, Z_2) != 1:

                     return False

        return True


    def test_vitesse_and_assembly_condition(self, planetary_gear, begin_gearing_chain, end_gearing_chain,
                                            list_previous_planetary):

        range_speed = planetary_gear.speed_range(begin_gearing_chain, planetary_gear.planet_carrier, list_previous_planetary)

        if not range_speed:
            return False


        if range_speed[begin_gearing_chain.name][0] > range_speed[begin_gearing_chain.name][1]:
            return False

        if range_speed[planetary_gear.planet_carrier.name][0] > range_speed[planetary_gear.planet_carrier.name][1]:
            return False

        if not planetary_gear.test_assembly_condition(self.number_planet, [begin_gearing_chain, end_gearing_chain]):
            return False

        return True

    def z_range_mini_max(self, planetary_gear, element, begin_gearing_chain, end_gearing_chain, path):

        if not element in path[0]:

            return []

        reason = planetary_gear.reason_abs(path[0])

        reason_1 = abs((begin_gearing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1])/
                       (end_gearing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1]))

        reason_2 = abs((begin_gearing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0])/
                       (end_gearing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0]))
        
        reason_3 = abs((begin_gearing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1])/
                       (end_gearing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1]))
        
        reason_4 = abs((begin_gearing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0])/
                       (end_gearing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0]))
        
        reason_min=min(reason_1,reason_2,reason_3,reason_4)
        reason_max=max(reason_1,reason_2,reason_3,reason_4)

        # if reason_1 < reason_2:
        #     reason_min = reason_1
        #     reason_max = reason_2

        # else:
        #     reason_min = reason_2
        #     reason_max = reason_1

        if reason_min and reason_max:
            Z_min = int(reason*element.Z/reason_max)-1
            Z_max = int(reason*element.Z/reason_min)+1

        elif not reason_min:
             Z_min = int(reason*element.Z/reason_max)-1
             Z_max = self.Z_range_sun[1]
        else:
            Z_min = 0
            Z_max = int(reason*element.Z/reason_min)+1
        Z_range_mini_maxi = [Z_min, Z_max]
       
        return Z_range_mini_maxi

    def decision_tree_speed_possibilities(self):
        tree = dt.DecisionTree()
        tree.SetCurrentNodeNumberPossibilities(self.number_input)
        node = tree.NextNode(True)
        list_solution = []
        while not tree.finished:


            if len(node) == 1:

                    self.planetary_gear.planet_carrier.speed_input = self.input_speeds[node[0]]

                    input_speeds_2 = copy.copy(self.input_speeds)
                    input_speeds_2.remove(self.input_speeds[node[0]])

                    list_possibility_speed = []
                    self.multiplication_possibility_speed(input_speeds_2, 0, [], list_possibility_speed)

                    tree.SetCurrentNodeNumberPossibilities(len(list_possibility_speed))


            if len(node) == 2:

                    possibility_speed = list_possibility_speed[node[1]]
                    for i, planetarie in enumerate(self.planetary_gear.planetaries):
                        planetarie.speed_input = possibility_speed[i]

                    list_solution.append(copy.deepcopy(self.planetary_gear))
                    tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(True)

        return list_solution
    def decision_tree(self):
        list_planetary_gears_speed = self.decision_tree_speed_possibilities()


        for i, planetary_gear in enumerate(list_planetary_gears_speed):
            print(i)
            planet_double = []


            for double in planetary_gear.doubles:

                if not double.nodes[0] in planet_double:
                    planet_double.append(double.nodes[0])

                if not double.nodes[1] in planet_double:
                    planet_double.append(double.nodes[1])

            list_planet_remove = []

            for planet in planetary_gear.planets:
                if not planet in planet_double:
                    list_planet_remove.append(planet)


            list_tree = []
            debut = time.time()
            list_node_range_data = []
            gearing_chains_modif = planetary_gear.gearing_chain()
            gearing_chains = copy.copy(gearing_chains_modif)
            number_element_gearing_chain = []
            numbers_gearing_chain = []
            number_gearing_chain = 0
            totals_element_previous_gearing_chain = []
            total_element_previous_gearing_chain = 0
            flags_gearing_change = []
            list_solution = []
            flag_gcd = []

            for i, gearing_chain in enumerate(gearing_chains_modif):
                if isinstance(gearing_chain[-1], Planetary):

                    if gearing_chain[-1].planetary_type == 'Ring':

                        gearing_chains_modif[i] = gearing_chain[::-1]
                        gearing_chains[i] = gearing_chain[::-1]

                    if  not isinstance(gearing_chain[0], Planetary) and gearing_chain[-1].planetary_type == 'Sun':

                        gearing_chains_modif[i] = gearing_chain[::-1]
                        gearing_chains[i] = gearing_chain[::-1]
                gearing_chain_2 = copy.copy(gearing_chain)

                for element in gearing_chain_2:

                    if element in list_planet_remove:

                        gearing_chains_modif[i].remove(element)


            for gearing_chain in gearing_chains_modif:

                number_element_gearing_chain.append(len(gearing_chain))
                flags_gearing_change.append(1)

                for i, element in enumerate(gearing_chain):
                    flag_gcd.append(2)

                    if i != 0:
                        flags_gearing_change.append(0)

                    totals_element_previous_gearing_chain.append(total_element_previous_gearing_chain)
                    numbers_gearing_chain.append(number_gearing_chain)

                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':

                        list_tree.append(self.Z_range_ring[1]-self.Z_range_ring[0])
                        list_node_range_data.append(self.Z_range_ring[0])

                    else:
                        list_tree.append(self.Z_range_sun[1]-self.Z_range_sun[0])
                        list_node_range_data.append(self.Z_range_sun[0])

                number_gearing_chain += 1
                total_element_previous_gearing_chain += len(gearing_chain)

            list_planet_remove_neighbour = []


            for i, planet in enumerate(list_planet_remove):
                planet.Z = 1
                list_planet_remove_neighbour.append([planet])

                for gearing in planetary_gear.gearings:

                    if gearing.nodes[0] == planet:
                        list_planet_remove_neighbour[i].append(gearing.nodes[1])

                    if gearing.nodes[1] == planet:
                        list_planet_remove_neighbour[i].append(gearing.nodes[0])


            tree = dt.RegularDecisionTree(list_tree)

            Z_range_mini_maxi = []
            Z_range_mini_maxi_2 = []
            flag_gearing_change = 0
            flag_Z_range_mini_maxi = 0
            flag_Z_range_mini_maxi_2 = 0
            number_max_z_planet = self.Z_range_sun[1]
            list_planetaries_Z_range_mini_maxi = []
            list_path = []
            while not tree.finished:

                valid = True
                node = tree.current_node

                number_gearing_chain = numbers_gearing_chain[len(node)-1]
                flag_gearing_change = flags_gearing_change[len(node)-1]
                total_element_previous_gearing_chain = totals_element_previous_gearing_chain[len(node)-1]

                element = gearing_chains_modif[number_gearing_chain][len(node)-total_element_previous_gearing_chain-1]
                element.Z = list_node_range_data[len(node)-1]+ node[len(node)-1]


                if len(node) == 1:


                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':
                        number_max_z_planet = element.Z



                elif not flag_gearing_change:

                    previous_element = gearing_chains_modif[number_gearing_chain][len(node)-total_element_previous_gearing_chain-2]

                    if len(gearing_chains_modif[number_gearing_chain]) > 2 \
                    and gearing_chains_modif[number_gearing_chain][0].planetary_type == 'Ring':

                        number_max_z_previous_planets = 0
                        for i in range(len(node)-total_element_previous_gearing_chain-2):
                            previous_planet = gearing_chains_modif[number_gearing_chain][i+1]


                            number_max_z_previous_planets += previous_planet.Z


                        if isinstance(element, Planetary):
                            flag_impose_z_number = True

                            for planet in list_planet_remove:
                              if planet in gearing_chains[number_gearing_chain]:
                                  flag_impose_z_number = False
                                  break

                            if flag_impose_z_number:

                                if element.Z != number_max_z_planet-number_max_z_previous_planets:

                                        valid = False


                            else:
                                if element.Z > number_max_z_planet-number_max_z_previous_planets:
                                        valid = False
                        else:
                            if element.Z > number_max_z_planet-number_max_z_previous_planets:

                                valid = False
                    else:
                        if element.Z > number_max_z_planet:

                                valid = False

                    if flag_gcd[len(node)-1] == 2 and valid:

                        for relation in planetary_gear.relations:

                            if relation.nodes[0] == previous_element and relation.nodes[1] == element:
                                flag_gcd[len(node)-1] = 1
                                break

                            if relation.nodes[1] == previous_element and relation.nodes[0] == element:
                                flag_gcd[len(node)-1] = 1
                                break


                        if flag_gcd[len(node)-1] == 2:
                           flag_gcd[len(node)-1] = 0


                    if flag_gcd[len(node)-1]:
                        if not self.test_GCD(previous_element.Z, element.Z):
                            valid = False




                else:

                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':
                        number_max_z_planet = element.Z


                    else:
                        number_max_z_planet = self.Z_range_sun[1]


                if len(node) == number_element_gearing_chain[number_gearing_chain]+total_element_previous_gearing_chain and valid:

                    begin_gearing_chain = gearing_chains_modif[number_gearing_chain][0]
                    end_gearing_chain = gearing_chains_modif[number_gearing_chain][-1]

                    planetary_gear = PlanetaryGear(planetary_gear.name, planetary_gear.planetaries, planetary_gear.planets,
                                                   planetary_gear.planet_carrier, planetary_gear.connexions)

                    if Z_range_mini_maxi:

                        if element.Z < Z_range_mini_maxi[0] or element.Z > Z_range_mini_maxi[1]:
                            valid = False

                    if Z_range_mini_maxi_2:

                        if element.Z < Z_range_mini_maxi_2[0]  or element.Z > Z_range_mini_maxi_2[1]:
                            valid = False

                    if valid:

                        if numbers_gearing_chain[len(node)-1] == 0:

                                first_planetary = begin_gearing_chain

                                if not first_planetary in list_planetaries_Z_range_mini_maxi:
                                    list_planetaries_Z_range_mini_maxi.append(first_planetary)


                                if isinstance(begin_gearing_chain, Planetary) and isinstance(end_gearing_chain, Planetary):

                                    if not end_gearing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_gearing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([begin_gearing_chain, end_gearing_chain]))

                                    if not flag_Z_range_mini_maxi:


                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element, begin_gearing_chain,
                                                                                  end_gearing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)-1])

                                        flag_Z_range_mini_maxi = 1



                                        if element.Z < Z_range_mini_maxi[0] or element.Z < Z_range_mini_maxi[1]:

                                            valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_gearing_chain,
                                                                                             end_gearing_chain,
                                                                                             list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)+1])

                                        else:
                                            valid = False
                                    else:

                                        valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_gearing_chain,
                                                                                         end_gearing_chain,
                                                                                         list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)+1])



                        else:

                            if not begin_gearing_chain in list_planetaries_Z_range_mini_maxi:


                                        list_planetaries_Z_range_mini_maxi.append(begin_gearing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, begin_gearing_chain]))

                            if not flag_Z_range_mini_maxi:

                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element,
                                                                                  first_planetary, begin_gearing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(begin_gearing_chain)-1])

                                        flag_Z_range_mini_maxi = 1

                                        if Z_range_mini_maxi:
                                            if element.Z < Z_range_mini_maxi[0] or element.Z < Z_range_mini_maxi[1]:
                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                                 begin_gearing_chain,
                                                                                                 list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(begin_gearing_chain)+1])

                                            else:
                                                valid = False

                            else:

                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                 begin_gearing_chain,
                                                                                 list_planetaries_Z_range_mini_maxi[:list_planetaries_Z_range_mini_maxi.index(begin_gearing_chain)+1])




                            if isinstance(end_gearing_chain, Planetary) and valid:


                                if not end_gearing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_gearing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, end_gearing_chain]))


                                if not flag_Z_range_mini_maxi_2:

                                        Z_range_mini_maxi_2 = self.z_range_mini_max(planetary_gear, element, first_planetary, end_gearing_chain,
                                                                                    list_path[list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)-1])
                                        flag_Z_range_mini_maxi_2 = 1
                                        if Z_range_mini_maxi_2:
                                            if element.Z < Z_range_mini_maxi_2[0] or element.Z < Z_range_mini_maxi_2[1]:

                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_gearing_chain,
                                                                                                 end_gearing_chain, number_planet,
                                                                                                 list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)+1])

                                            else:
                                                valid = False
                                else:

                                    valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                     end_gearing_chain, number_planet,
                                                                                     list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_gearing_chain)+1])






                        if number_gearing_chain == len(gearing_chains_modif)-1 and valid:

                            list_tree_planetary = []





                            if list_planet_remove:

                                for i in range(len(list_planet_remove)):
                                    list_planet_remove_2  = copy.copy(list_planet_remove)
                                    list_planet_remove_neighbour_2 = copy.copy(list_planet_remove_neighbour)
                                    list_planet_remove[i].Z = 1
                                    valid_planet = True

                                    for gearing_chain in gearing_chains:
                                        if list_planet_remove[i] in gearing_chain and valid_planet:

                                            if gearing_chain[0].planetary_type == 'Ring' :

                                                number_z_max = 0
                                                for element_2 in gearing_chain:
                                                    if  element_2 != list_planet_remove[i] and element_2 != gearing_chain[0]:
                                                        number_z_max += element_2.Z

                                                if  isinstance(gearing_chain[-1], Planetary):

                                                    list_planet_remove[i].Z = gearing_chain[0].Z-number_z_max
                                                    if list_planet_remove[i].Z < self.Z_range_sun[0]:
                                                        list_planet_remove[i].Z = self.Z_range_sun[0]
                                                        valid_planet = False
                                                    neighbour = list_planet_remove_neighbour[i]
                                                    if not self.test_GCD(list_planet_remove[i].Z, neighbour[1].Z):
                                                        valid_planet = False

                                                        break
                                                    if not self.test_GCD(list_planet_remove[i].Z, neighbour[2].Z):
                                                        valid_planet = False
                                                        break

                                                    list_planet_remove_2.remove(list_planet_remove[i])
                                                    list_planet_remove_neighbour_2.remove(neighbour)

                                                else:
                                                    if gearing_chain[0].Z-number_z_max - self.Z_range_sun[0] > 0:
                                                        list_tree_planetary.append(gearing_chain[0].Z-number_z_max - self.Z_range_sun[0])

                                                    else:
                                                        valid_planet = False
                                                        break

                                            else:
                                                list_tree_planetary.append(self.Z_range_sun[1]-self.Z_range_sun[0])
                                            break
                                if list_tree_planetary and valid_planet:
                                    tree_planet = dt.RegularDecisionTree(list_tree_planetary)

                                    while not tree_planet.finished:

                                        valid_planet = True
                                        node_planet = tree_planet.current_node
                                        element_planet = list_planet_remove_2[len(node_planet)-1]
                                        neighbour = list_planet_remove_neighbour_2[len(node_planet)-1]
                                        element_planet.Z = node_planet[len(node_planet)-1]+self.Z_range_sun[0]

                                        if not self.test_GCD(element_planet.Z, neighbour[1].Z):
                                            valid_planet = False

                                        if not self.test_GCD(element_planet.Z, neighbour[2].Z):
                                            valid_planet = False

                                        if valid_planet and len(node_planet) == len(list_tree_planetary):
                                            list_solution.append(planetary_gear)

                                            print(planetary_gear)
                                            # if list_solution:
                                            #     return list_solution

                                        tree_planet.NextNode(valid_planet)

                            else:
                                list_solution.append(planetary_gear)
                                print(planetary_gear)

                        # if len(list_solution) > 30:
                        #     return list_solution

                else:

                    Z_range_mini_maxi = []
                    Z_range_mini_maxi_2 = []
                    flag_Z_range_mini_maxi = 0
                    flag_Z_range_mini_maxi_2 = 0

                tree.NextNode(valid)

            fin = time.time()

            print(fin-debut)
        return list_solution




















