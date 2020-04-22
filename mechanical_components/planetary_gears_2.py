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


class PlanetaryGear():

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
    
 
    
class PlanetsStructure():
    def __init__(self,architecture,branch_connexions):
        self.architecture= architecture
        self.branch_connexions=branch_connexions
        self.planets=[]
        self.connexions=[]
        number_planet=0
        self.gearings = []
        self.doubles = []
        
        first_last_composant_branch=[]
        for branch in self.architecture:
            
            flag_first_planet_branch=0
            for planet in branch:
                number_planet+=1
                self.planets.append(Planet('Pl'+str(number_planet),planet,7))
                
                if flag_first_planet_branch :
                    if planet=='Double' and self.planets[-2].planet_type=='Double':
                        self.connexions.append([self.planets[-2], self.planets[-1],'D'])
                    else:
                        self.connexions.append([self.planets[-2], self.planets[-1],'GI'])
                else:
                    first_composant_branch=self.planets[-1]
                    flag_first_planet_branch=1
            last_composant_branch=self.planets[-1]
            first_last_composant_branch.append([first_composant_branch,last_composant_branch])
            
        for branch_connexion in self.branch_connexions:
            
            if branch_connexion[2]=='Simple':
                self.connexions.append([first_last_composant_branch[branch_connexion[0]-1][1], first_last_composant_branch[branch_connexion[1]-1][0],'GI'])
            else:
                self.connexions.append([first_last_composant_branch[branch_connexion[0]-1][1], first_last_composant_branch[branch_connexion[1]-1][0],'D'])     
               
        for i, connexion in  enumerate(self.connexions):

          if connexion[2] != 'D':

             self.gearings.append(GearingPlanet('Gearing'+str(i),
                                                      [connexion[0], connexion[1]]))

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
    

    def plot_cinematic_graph_gear(self, coordinate,length,diameter):
        x=npy.array([coordinate[0],coordinate[0]+length/2,coordinate[0]-length/2,coordinate[0],coordinate[0],coordinate[0]+length/2,coordinate[0]-length/2])
        y=npy.array([coordinate[1]+diameter/2,coordinate[1]+diameter/2,coordinate[1]+diameter/2,coordinate[1]+diameter/2,coordinate[1]-diameter/2,coordinate[1]-diameter/2,coordinate[1]-diameter/2])
        plt.plot(x,y)
        
    def plot_cinematic_graph(self,length_gear,diameter_gear,length_double):
        graph_path=self.path_planet_to_planet()
        plt.figure()
        
        previous_relation_double=[]
        previous_relation_gearing=[]
        
        previous_planet_gearing=[]
        previous_planet_double=[]
        inverse_relation_double=[]
        inverse_relation_gearing=[]
        
        for path in graph_path:
            
            flag_inv_gearing=0
            flag_way_inv_gearing=0
            flag_inv_double=0
            flag_way_inv_double=0
            coordinate=[0,0]
            self.plot_cinematic_graph_gear(coordinate,length_gear,diameter_gear)
            for i,element in enumerate(path):
                if isinstance(element,Double):
                            
                    
                    
                    if element in  inverse_relation_double:
                        x=[coordinate[0],coordinate[0]-length_double/(1+i*0.2)]
                        y=[coordinate[1],coordinate[1]]
                        coordinate=[coordinate[0]-length_double/(1+i*0.2),coordinate[1]]
                        


                        
                    
                    
                    elif ((element.nodes[0] in previous_planet_double or element.nodes[1] in previous_planet_double) and not element  in previous_relation_double ):
                        
                        for double in previous_relation_double:
                            for node in double.nodes:
                                if element.nodes[0]==node or element.nodes[1]==node:
                                    if  not double==previous_element:
                                        if double in inverse_relation_double:
                                            flag_way_inv_double=1
                                        else:
                                            flag_way_inv_double=0
                                    else:
                                        if not double in inverse_relation_double:
                                            flag_way_inv_double=1
                                        else:
                                            flag_way_inv_double=0
                                
                        if flag_way_inv_double:
                            x=[coordinate[0],coordinate[0]+length_double/(1+i*0.2)]
                            y=[coordinate[1],coordinate[1]]
                            coordinate=[coordinate[0]+length_double/(1+i*0.2),coordinate[1]]
                        else:
                            x=[coordinate[0],coordinate[0]-length_double/(1+i*0.2)]
                            y=[coordinate[1],coordinate[1]]
                            coordinate=[coordinate[0]-length_double/(1+i*0.2),coordinate[1]]
                            inverse_relation_double.append(element)
                            
                        

                            
                    else:    
                        x=[coordinate[0],coordinate[0]+length_double/(1+i*0.2)]
                        y=[coordinate[1],coordinate[1]]
                        coordinate=[coordinate[0]+length_double/(1+i*0.2),coordinate[1]]
                        
                    plt.plot(x,y)    
                    previous_relation_double.append(element)
                    previous_planet_double.extend([element.nodes[0],element.nodes[1]])
                
                elif isinstance(element,GearingPlanet)  :
                    if element in  inverse_relation_gearing:
                        coordinate=[coordinate[0],coordinate[1]-diameter_gear]

                    
                    if ((element.nodes[0] in previous_planet_gearing or element.nodes[1] in previous_planet_gearing) and not element  in previous_relation_gearing ):
                        
                        for gearing in previous_relation_gearing:
                            for node in gearing.nodes:
                                if element.nodes[0]==node or element.nodes[1]==node:

                                    if  not gearing==previous_element:
                                        if gearing in inverse_relation_gearing:
                                            
                                            flag_way_inv_gearing=1
                                        else:
                                            flag_way_inv_gearing=0
                                    else:
                                        if not gearing in inverse_relation_gearing:
                                            flag_way_inv_gearing=1
                                        else:
                                            flag_way_inv_gearing=0
                                            
                        if flag_way_inv_gearing:
                            coordinate=[coordinate[0],coordinate[1]+diameter_gear]
                            inverse_relation_gearing.append(element)
                        else:
                            coordinate=[coordinate[0],coordinate[1]-diameter_gear]
                      

                            
                    else:    

                        
                        coordinate=[coordinate[0],coordinate[1]+diameter_gear]
 
                    previous_relation_gearing.append(element)
                    previous_planet_gearing.extend([element.nodes[0],element.nodes[1]])
                
                if not isinstance(element,Planet) :   
                    self.plot_cinematic_graph_gear(coordinate,length_gear,diameter_gear)
                    previous_element=element
                    
        
        
          
                
            
    
    
class OptimizerPlanetaryGears():
    def __init__(self,input_speeds,diameter_cylinder,length_cylinder):
        self.input_speeds=input_speeds
        self.diameter=diameter_cylinder
        self.length=length_cylinder
        self.number_input=len(input_speeds)
        
    ## Recursive Function which give all the possibilities  of planet_type in a branch for a Planet number fixed ##
    ## Exemple Input:(0,[],[],4) -> Output:[['Simple', 'Double', 'Double', 'Simple'], ['Simple', 'Double', 'Double', 'Double'], ['Simple', 'Simple', 'Double', 'Double'], 
    ##['Simple', 'Simple', 'Simple', 'Simple'], ['Double', 'Double', 'Simple', 'Simple'], ['Double', 'Double', 'Double', 'Simple'], ['Double', 'Double', 'Double', 'Double']] ##    
    def list_possibilities_planets_type(self, n, list_planet_type, planet_type, number_max_planet):

        if n==number_max_planet:
            planet_type_2=copy.copy(planet_type)
            list_planet_type.append(planet_type_2)
        
            return list_planet_type
        
        if not planet_type:
            planet_type_1=copy.copy(planet_type)
            planet_type_2=copy.copy(planet_type)
            planet_type_1.append('Simple')
              
            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_1, number_max_planet)
            if number_max_planet>1:  
                planet_type_2.extend(['Double','Double'])
                self.list_possibilities_planets_type(n+2, list_planet_type, planet_type_2, number_max_planet)

            
        elif planet_type[-1]=='Simple':
            
            if n<number_max_planet-1:
                planet_type_2=copy.copy(planet_type)
                planet_type_2.extend(['Double','Double'])
                self.list_possibilities_planets_type(n+2, list_planet_type, planet_type_2, number_max_planet)
            planet_type_1=copy.copy(planet_type)
            planet_type_1.append('Simple')
            
            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_1, number_max_planet)
        elif planet_type[-1]=='Double':
            
            planet_type_1=copy.copy(planet_type)
            
            planet_type_2=copy.copy(planet_type)
            planet_type_1.append('Simple')
            
            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_1, number_max_planet)
            
            planet_type_2.append('Double')
            self.list_possibilities_planets_type(n+1, list_planet_type, planet_type_2, number_max_planet)
            
    
    
    ## Recursive Function which give all the possibilities for a junction number fixed ##
    ## Exemple Input:([],[0,0,0,0],2,2,[2,2,2,2]) -> Output:[[2, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 2, 0, 0], 
    ## [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 2, 0], [0, 0, 1, 1], [0, 0, 0, 2]] ##
    def list_possibilities_junction(self, list_possibilities, possibilitie, number_junction, 
                                    number_junction_max, number_junction_max_by_planet, sum_number_junction_max_by_planet):
        
        number_junction_2=copy.copy(number_junction)
        possibilitie_2=copy.copy(possibilitie)
        sum_number_junction_max_by_planet_2=copy.copy(sum_number_junction_max_by_planet)
        
        if number_junction==number_junction_max:
            
            if not possibilitie_2  in list_possibilities:
                list_possibilities.append(possibilitie_2)
                
            return list_possibilities


        for i in range(1,len(possibilitie)):

            if i>1 and flag:
                possibilitie_2[i-1]-=1
                sum_number_junction_max_by_planet_2[i]=sum_number_junction_max_by_planet[i]
                
            flag=1
            
            if possibilitie_2[i]<sum_number_junction_max_by_planet_2[i]:
                
                if i+1<len(possibilitie)-1:
                    sum_number_junction_max_by_planet_2[i+1]+=number_junction_max_by_planet
                    
                possibilitie_2[i]+=1
                
                self.list_possibilities_junction(list_possibilities, possibilitie_2, number_junction_2+1, 
                                                 number_junction_max, number_junction_max_by_planet, sum_number_junction_max_by_planet_2)
            
            else:
                flag=0
    
    
    ## Recursive Function which give all the possibilities for a limited number of connexion( number_max_connextion) in a list_connexion ##
    ## Exemple Input:(0 ,[], [[0, 0], [0, 0]] ,[[2, 4], [2, 5], [3, 6]],2,2) -> Output:[[[2, 4], [2, 5]], [[2, 4], [3, 6]], [[2, 5], [3, 6]]] ##
    def list_possibilities_connexion_branch_step_1(self, n, list_possibilities, possibilitie, 
                                                   list_connexion_branch, number_max_connexion, number_max_junction_by_branch):
        
        possibilitie_2=copy.copy(possibilitie)
        if  n==number_max_connexion:
            
            if not possibilitie in list_possibilities:
                flag_similaritie=0
                
                for possibilitie_3 in list_possibilities: 
                    similaritie=0
                    
                    for connexion in possibilitie_2:
                        
                        if connexion in possibilitie_3:
                            similaritie+=1
         
                    if similaritie==number_max_connexion:
                        flag_similaritie=1
                        break
                    
                if not flag_similaritie:
                    list_possibilities.append(possibilitie_2)
                    
                    
        else:
            n+=1
            for element_1 in list_connexion_branch:
                
                number_connexion_by_branch=0
                flag_branch=1
                
                for element_2 in possibilitie_2:
                    
                    if element_2[0]==element_1[0]:
                        number_connexion_by_branch+=1
                        
                    if element_2[1]==element_1[1]:
                        flag_branch=0
                        
                if number_connexion_by_branch<=number_max_junction_by_branch and flag_branch:          
                    possibilitie_2[n-1]=element_1
                    
                    self.list_possibilities_connexion_branch_step_1(n, list_possibilities, possibilitie_2, 
                                                                    list_connexion_branch,number_max_connexion,number_max_junction_by_branch)
                    possibilitie_2[n-1]=possibilitie[n-1]
            
            
            
            
            
                
                
    
    
    
    ## Recursive Function which give all the possibilities  of connexion  for a number_junction fixed and a number branch fixed ##
    ## Exemple Input:(0 ,[0,0],[],[],[2,3,4,5,6],[1],[2,2],2) -> Output:[[[[1, 2], [1, 3]], [[2, 4], [2, 5]]], [[[1, 2], [1, 3]], [[2, 4], [3, 5]]], [[[1, 2], [1, 3]], [[2, 5], [3, 4]]], [[[1, 2], [1, 3]], [[3, 4], [3, 5]]]] ##
    def list_possibilities_connexion_branch_step_2(self, n, possibilitie, list_possibilities, 
                                                   list_previous_branch, remaning_branch, previous_branch, 
                                                   number_junction_branch, number_max_junction_by_branch):
        
        possibilitie_2=copy.copy(possibilitie)
        remaning_branch_2=copy.copy(remaning_branch)
        previous_branch_2=copy.copy(previous_branch)
        
        if n==len(number_junction_branch):
            list_possibilities.append(possibilitie_2)
            list_previous_branch.append(previous_branch_2)
            return None

        number_connexion=number_junction_branch[n]
        n+=1

        list_possibilities_connexion=[]
        for element in previous_branch_2:
            for j in range(number_connexion):
                list_possibilities_connexion.append([element,remaning_branch_2[j]])

        list_possibilities_connexion_2=[]
        possibilitie_3=[]
        
        for i in range(number_connexion):

            possibilitie_3.append([0,0])
        
        self.list_possibilities_connexion_branch_step_1(0, list_possibilities_connexion_2, possibilitie_3, 
                                                        list_possibilities_connexion, number_connexion, number_max_junction_by_branch ) 
        
        
            
        for element_1 in list_possibilities_connexion_2:
            remaning_branch_3=copy.copy(remaning_branch_2)
            previous_branch_3=copy.copy(previous_branch_2)

            for i in range(number_connexion):
                previous_branch_3.append(remaning_branch_2[i])
                remaning_branch_3.remove(remaning_branch_2[i])
                
            possibilitie_2[n-1]=element_1
            
            for element_2 in element_1:
                
                if element_2[0] in previous_branch_3:
                    previous_branch_3.remove(element_2[0])
                    
            self.list_possibilities_connexion_branch_step_2(n, possibilitie_2, list_possibilities, list_previous_branch, 
                                                            remaning_branch_3, previous_branch_3, number_junction_branch, number_max_junction_by_branch)
            
        

            
    ## Recursive Function which give all the possibilities to distribute a number planet fixed in a number_branch fixed  ##
    ## Exemple Input:([0,0,0],[],6,0,3,0,1) -> Output:[[1, 1, 4], [1, 2, 3], [1, 3, 2], [1, 4, 1], [2, 1, 3], [2, 2, 2], [2, 3, 1], [3, 1, 2], [3, 2, 1], [4, 1, 1]] ##    
    def list_possibilities_planets_by_branch_step_1(self, possibilities, list_possibilities, 
                                                    number_planet, number_branch, number_branch_max, 
                                                    number_other_planet, min_planet_branch):
        
        possibilities_2=copy.copy(possibilities)
        number_branch+=1
        
        if number_branch==number_branch_max:
            
            if number_planet-number_other_planet>0:
                possibilities_2[number_branch-1]=number_planet-number_other_planet
                list_possibilities.append(possibilities_2)
                
            return list_possibilities
        

        for i in range (min_planet_branch,number_planet-number_other_planet-(number_branch_max-number_branch-1)):
            
            possibilities_2[number_branch-1]=i
            self.list_possibilities_planets_by_branch_step_1(possibilities_2, list_possibilities, 
                                                             number_planet, number_branch, number_branch_max, number_other_planet+i, 1)
    

    
    ##  Function which give all the possibilities of branch for a possibility junction  ##
    ## Exemple Input:([0,1,1],7,3,1) -> Output:[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]], 
    ##[[2, 1, 1, 2, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]], [[2, 1, 2, 1, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]]etc... ##               
    def list_possibilities_planets_by_branch_step_2(self, junction, number_planet, number_max_junction_by_branch, min_planet_branch):
        list_number_planet_branch=[]
        number_planet_branch=0
        list_number_junctions=[]
        number_branch=1
        
        for i in range(len(junction)):
            number_planet_branch+=1
            if junction[i]:
                number_branch+=junction[i]+1
                list_number_planet_branch.append(number_planet_branch)
                number_planet_branch=0
                list_number_junctions.append(junction[i]+1)
                
        list_connexion=[]
        connexion=[]
        remaning_branch=[]
        previous_branch=[1]
        list_number_planet=[0]
        
        for i in range(number_branch-1):
            remaning_branch.append(i+2)
            list_number_planet.append(0)
            
        for i in range(len(list_number_junctions)):
            connexion.append(0)
            
        list_previous_branch=[]
        
        self.list_possibilities_connexion_branch_step_2(0, connexion, list_connexion, list_previous_branch, 
                                                        remaning_branch, previous_branch, list_number_junctions, number_max_junction_by_branch)
        
        list_architecture=[]
        for i,connexion_2 in enumerate(list_connexion):
            
            number_planet_know=0
            previous_branch_final=list_previous_branch[i]
            previous_connexion_by_branch=[]
            
            for i,connexion_by_branch in enumerate(connexion_2):
                
                for element in connexion_by_branch:
                    
                    if not element[0] in previous_connexion_by_branch:
                        
                        list_number_planet[element[0]-1]=list_number_planet_branch[i]
                        number_planet_know+=list_number_planet_branch[i]
                        previous_connexion_by_branch.append(element[0])
                    
                    
       
            number_planet_unknow=number_planet-number_planet_know
            possibilitie=[]
            
            for i in range(len(previous_branch_final)):
                possibilitie.append(0)
                
            list_planet_by_branch=[]
            
            self.list_possibilities_planets_by_branch_step_1(possibilitie, list_planet_by_branch, number_planet_unknow, 0, 
                                                             len(previous_branch_final), 0, min_planet_branch)
            
            for planet_by_branch in list_planet_by_branch:
                
                for i,branch in enumerate(previous_branch_final):
                    
                    list_number_planet[branch-1]=planet_by_branch[i]
                    
                list_architecture.append([copy.copy(list_number_planet),connexion_2])
        
        return list_architecture,number_branch
    
    
    
    ##  Recursive Function which give all the possibilities of architecure (branch and planet_type) for a possibility junction  ##
    ## Exemple Input:(1,5,[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]],[0,0,0,0,0],[],[],[]]) -> Output: list of PlanetsStructure objet ##  
    def list_possibilities_architecture_planet(self, branch, branch_maximum,architecture, possibilities, 
                                               list_possibilities, possibilities_connexion, list_possibilities_connexion):
        

        
        possibilities_connexion_2=copy.copy(possibilities_connexion)
        number_planet=architecture[0][branch-1]
        
        list_planet_possibilities=[]                
        if number_planet==1:
            for branch_connexion in possibilities_connexion_2:
                if branch_connexion[1]==branch:
                    list_planet_possibilities= [[branch_connexion[2]]]
                    break
            
        else:   
            self.list_possibilities_planets_type(0, list_planet_possibilities, [], number_planet)
            
            
        number_connexion=0
        for connexion in architecture[1]:
            for branch_2 in connexion:
                if branch_2[0]==branch:
                    number_connexion+=1

                    possibilities_connexion_2.append([branch,branch_2[1],'Simple'])
                    

        
        for planet_possibilities in list_planet_possibilities:
            
            possibilities_2=copy.copy(possibilities)
            possibilities_2[branch-1]=planet_possibilities

            if branch==branch_maximum:
                possibilities_connexion_3=copy.deepcopy(possibilities_connexion_2) 
                list_possibilities.append(PlanetsStructure(possibilities_2,possibilities_connexion_3))
                
            else:
            
                if number_connexion==0:
                    self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                
                
                
                elif number_connexion==2:

                    possibilities_connexion_2[-1][2]='Double'
                    possibilities_connexion_2[-2][2]='Simple'
                    
                    self.list_possibilities_architecture_planet(branch+1,branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                    
                    possibilities_connexion_2[-1][2]='Simple'
                    possibilities_connexion_2[-2][2]='Double'
                    self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                    
                    
                    if planet_possibilities[-1]=='Simple':
                        
                        possibilities_connexion_2[-1][2]='Double'
                        possibilities_connexion_2[-2][2]='Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                    elif planet_possibilities[-1]=='Double':
                        
                        possibilities_connexion_2[-1][2]='Simple'
                        possibilities_connexion_2[-2][2]='Simple'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                
                
                
                elif number_connexion==3:
                    
                    if planet_possibilities[-1]=='Simple':
                        
                        possibilities_connexion_2[-1][2]='Double'
                        possibilities_connexion_2[-2][2]='Double'
                        possibilities_connexion_2[-3][2]='Simple'
                        
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                        possibilities_connexion_2[-1][2]='Simple'
                        possibilities_connexion_2[-2][2]='Double'
                        possibilities_connexion_2[-3][2]='Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                    list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                        possibilities_connexion_2[-1][2]='Double'
                        possibilities_connexion_2[-2][2]='Simple'
                        possibilities_connexion_2[-3][2]='Double'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2,
                                                                    list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                    elif planet_possibilities[-1]=='Double':
                        
                        possibilities_connexion_2[-1][2]='Simple'
                        possibilities_connexion_2[-2][2]='Simple'
                        possibilities_connexion_2[-3][2]='Double'
                       
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                    list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                        possibilities_connexion_2[-1][2]='Simple'
                        possibilities_connexion_2[-2][2]='Double'
                        possibilities_connexion_2[-3][2]='Simple'
                        self.list_possibilities_architecture_planet(branch+1, branch_maximum, architecture, possibilities_2, 
                                                                    list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                        
                        possibilities_connexion_2[-1][2]='Double'
                        possibilities_connexion_2[-2][2]='Simple'
                        possibilities_connexion_2[-3][2]='Simple'
                        self.list_possibilities_architecture_planet(branch+1,branch_maximum, architecture, possibilities_2, 
                                                                    list_possibilities, possibilities_connexion_2, list_possibilities_connexion)
                            
        
            
        
    def decission_tree(self, number_max_planet, number_junction, number_max_junction_by_planet, min_planet_branch):
        tree=dt.DecisionTree()
        tree.SetCurrentNodeNumberPossibilities(self.number_input)
        node=tree.NextNode(True)
        planet_carrier=PlanetCarrier('PlanetCarrier')
        planetaries=[]
        n=0
        list_solution=[]
        for i in range (self.number_input-1):
            planetaries.append(Planetary('Planetary_'+str(i),7,'Sun'))
            
        while not tree.finished:
            
            if len(node)==1:
                planet_carrier.speed=self.input_speeds[node[0]]
                num_planetary=0
                
                for i in range(self.number_input):
                    if i != node[0] :
                        planetaries[num_planetary] = self.input_speeds[i]
                        num_planetary+=1
                        
                list_possibilities_junction=[]
                list_planet=[]
                sum_number_max_junction_by_planet=[]
                
                for i in range(number_max_planet-2):
                    list_planet.append(0)
                    sum_number_max_junction_by_planet.append(number_max_junction_by_planet)
                
                self.list_possibilities_junction(list_possibilities_junction, list_planet, 0, number_junction, 
                                                 number_max_junction_by_planet, sum_number_max_junction_by_planet)
                
                tree.SetCurrentNodeNumberPossibilities(len(list_possibilities_junction))
                print(list_possibilities_junction)
                
                
                
                
                
                
                        
            if len(node)==2:
                
                list_junction=list_possibilities_junction[node[1]]
                
                list_global_architecture,number_branch=self.list_possibilities_planets_by_branch_step_2(list_junction, number_max_planet, 
                                                                                                        number_max_junction_by_planet, min_planet_branch)
                tree.SetCurrentNodeNumberPossibilities(len(list_global_architecture))
                print(list_global_architecture)
                
                
            if len(node)==3:
                global_architecture=list_global_architecture[node[2]]
                list_branch=[]
                for i in range(number_branch):
                    list_branch.append(0)
                    
                list_possibilitie_planet_architecture=[]
                list_connexion=[]
                self.list_possibilities_architecture_planet(1, number_branch, global_architecture, list_branch, 
                                                            list_possibilitie_planet_architecture, [], list_connexion)
                
                tree.SetCurrentNodeNumberPossibilities(len(list_possibilitie_planet_architecture))

                
            
            if len(node)==4:
                planet_architecture=list_possibilitie_planet_architecture[node[3]]
                #planet_architecture.plot_cinematic_graph(0.1,1,4)
                n+=1
                list_solution.append(planet_architecture)
                #planet_architecture.plot()
                tree.SetCurrentNodeNumberPossibilities(0)
                
            node=tree.NextNode(True)
        
        return list_solution
                
            
                
                                                       
                
                
                
                
                
                
                
                
        
        

