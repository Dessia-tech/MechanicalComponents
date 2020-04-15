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

class Gears():

    def __init__(self, name, Z):
        self.name = name
        self.Z = Z

    def volume_plot(self, xy_position, z_position, module, length):
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
         extrusion_vector1 = length*z
         C1 = vm.Contour2D(Gears3D_Rotate[0])
         t1 = p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1, [], vm.Vector3D(extrusion_vector1))

         return t1

class Planetary(Gears):

    def __init__(self, name, Z, planetary_type):

        self.name = name
        self.Z = Z
        self.planetary_type = planetary_type
        self.p = 0
        self.speed = 0
        self.module = 0
        self.d = 0
        Gears.__init__(self, self.name, self.Z)

        if planetary_type == 'Sun':
            self.p = 1

        else:
            self.p = -1

class Planet(Gears):

    def __init__(self, name, planet_type, Z):
        self.name = name
        self.planet_type = planet_type
        self.Z = Z
        self.speed = 0
        self.module = 0

        Gears.__init__(self, self.name, self.Z)

class PlanetCarrier():

    def __init__(self, name):
        self.name = name
        self.speed = 0


class Gearing():

    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes


class GearingPlanetary(Gearing):

    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes
        Gearing.__init__(self, self.name, self.nodes)

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

        def __init__(self, name, nodes):
            self.name = name
            self.nodes = nodes
            Gearing.__init__(self, self.name, self.nodes)
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


class Pivot():

    def __init__(self, name, nodes):
        self.nodes = nodes
        self.name = name

    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Fixed():

    def __init__(self, name, nodes):
        self.nodes = nodes
        self.name = name

    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

class Double():

    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes

    def speed_system_equations(self):
        matrix = npy.array([1, -1])
        rhs = npy.array([0])
        return matrix, rhs

    def volume_plot(self, xy_position, z_position, radius, length):

         pos = vm.Point3D((xy_position[0], xy_position[1], z_position))
         axis = vm.Vector3D((0, 0, 1))
         cylinder = p3d.Cylinder(pos, axis, radius, length)
         return cylinder

class ImposeSpeed():

    def __init__(self, name, node, input_speed):
        self.input_speed = input_speed
        self.node = node
        self.name = name

    def speed_system_equations(self):
        matrix = npy.array([1])
        rhs = npy.array([self.input_speed])
        return matrix, rhs


class PlanetaryGears():

    def __init__(self, name, planetaries, planets, planet_carrier, connexions):
        self.name = name
        self.planetaries = planetaries
        self.planets = planets
        self.planet_carrier = planet_carrier
        self.elements = self.planetaries + self.planets + [self.planet_carrier]
        self.connexions = connexions
        self.gearings = []
        self.doubles = []

        for i, connexion in  enumerate(connexions):

          if connexion[2] != 'D':

              if isinstance(connexion[0], Planet) and isinstance(connexion[1], Planet):
                self.gearings.append(GearingPlanet('Gearing'+str(i), [connexion[0], connexion[1]]))

              else:
                self.gearings.append(GearingPlanetary('Gearing'+str(i),
                                                      [connexion[0], connexion[1]]))

          else:
             self.doubles.append(Double('Double'+str(i), [connexion[0], connexion[1]]))

        self.relations = self.gearings + self.doubles

    def __str__(self):

        Z_planets = []

        for planet in self.planets:
            Z_planets.append(planet.Z)

        Z_planetaries = []
        number_ring = 0
        number_sun = 0

        for planetary in self.planetaries:
            Z_planetaries.append(planetary.Z)

            if planetary.planetary_type == 'Sun':
                number_sun += 1

            else:
                number_ring += 1

        return 'Name:' + self.name + '\n\n' + \
               'Planetary Number:' + str(len(self.planetaries)) + '\n' + \
               'Ring Number:'+ str(number_ring) + '\n' + \
               'Sun_Number:' + str(number_sun) + '\n' + \
               'Z_planetaries:' + str(Z_planetaries) + '\n\n' + \
               'Planets_Number:' + str(len(self.planets)) + '\n' + \
               'Planets_Double_Number:' + str(len(self.doubles)) + '\n' + \
               'Z_Planets:' + str(Z_planets) + '\n\n\n'


    def matrix_position(self, element):
        return self.elements.index(element)

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


        return graph_planetary_gear

    def plot(self):
        graph_planetary_gears = self.graph()

        nx.draw_kamada_kawai(graph_planetary_gears, with_labels=True)

    def path_planetary_to_planetary(self):
        graph_planetary_gears = self.graph()
        graph_planetary_gears.remove_node(self.planet_carrier.name)
        list_path = []

        for planetary in self.planetaries[1:]:
            list_path.append(nx.shortest_path(graph_planetary_gears,
                                              self.planetaries[0].name, planetary.name))

        for path in list_path:

            for i in range(len(path)):
                path[i] = nx.get_node_attributes(graph_planetary_gears, path[i])[path[i]]


        return list_path

    def test_assembly_condition(self, number_planet):
        valid = True
        list_path = self.path_planetary_to_planetary()



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

        element_without_doubles = self.elements

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
            
        torque_element_association[self.planet_carrier.name] =solution[-1]    
        return torque_element_association

