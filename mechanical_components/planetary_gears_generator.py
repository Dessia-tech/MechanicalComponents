#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 17:34:03 2020

@author: launay
"""

from dessia_common import DessiaObject
from typing import  List, Tuple
import networkx as nx

from mechanical_components.planetary_gears import Planetary, Planet, PlanetCarrier,  \
Connection, PlanetaryGear, PlanetsStructure, PositionMinMaxPlanetaryGear, PlanetaryGearResult

import dectree as dt
import time
import math as m
import copy
import scipy.optimize as op
import numpy as np
import matplotlib.pyplot as plt
import random as random
import cma
import pyDOE
from dessia_common.vectored_objects import Catalog, Objective, ParetoSettings, ObjectiveSettings, from_csv
import os
import plot_data
class GeneratorPlanetsStructure(DessiaObject):
    '''
    A geanerator of planet_structure

    :param number_max_planet: The number of planet in the planet structure
    :type number_max_planet: int

    :param number_junction: The number of junction in the planet structure
                             (a junction is for example when a planet is connected to 3 element (Planetary or Planets))
                             (when a planet is connected to 4 elements we consider that is equal to 2 junctions)
    :type number_junction: int

    :param number_max_junction_by_planet: The maximum of junction that we can have on 1 planet.
    :type number_max_junction_by_planet: int

    :param min_planet_branch: The minimum of planet that we want in a branch.
                               (a branch begining with of planetary or a junction and ending with a planetary or a junction)
                               (when there are a junction, 1 branch ending and 2 begining)

    :type min_planet_branch: int

    :param name: Name
    :type name: str, optional


    '''
    _standalone_in_db = True

    _eq_is_data_eq = False

    def __init__(self, number_max_planet: int, number_junction: int, number_max_junction_by_planet: int, min_planet_branch: int,
                 number_max_meshing_plan: int, name: str = ''):

        self.number_max_meshing_plan = number_max_meshing_plan
        self.number_max_planet = number_max_planet
        self.number_junction = number_junction
        self.number_max_junction_by_planet = number_max_junction_by_planet
        self.min_planet_branch = min_planet_branch
        DessiaObject.__init__(self, name=name)

    ## Recursive Function which give all the possibilities  of planet_type in a branch for a Planet number fixed ##
    ## Exemple Input:(0,[],[],4) -> Output:[['Simple', 'Double', 'Double', 'Simple'], ['Simple', 'Double', 'Double', 'Double'], ['Simple', 'Simple', 'Double', 'Double'],
    ##['Simple', 'Simple', 'Simple', 'Simple'], ['Double', 'Double', 'Simple', 'Simple'], ['Double', 'Double', 'Double', 'Simple'], ['Double', 'Double', 'Double', 'Double']] ##
    def planets_type_possibilities(self, n, list_planet_type, planet_type, number_max_planet):

        if n == number_max_planet:
            planet_type_2 = copy.copy(planet_type)
            list_planet_type.append(planet_type_2)

            return list_planet_type

        if not planet_type:
            planet_type_1 = copy.copy(planet_type)
            planet_type_1.append('Simple')

            self.planets_type_possibilities(n+1, list_planet_type, planet_type_1, number_max_planet)
        else:

            planet_type_1 = copy.copy(planet_type)

            planet_type_2 = copy.copy(planet_type)
            planet_type_1.append('Simple')

            self.planets_type_possibilities(n+1, list_planet_type, planet_type_1, number_max_planet)

            planet_type_2.append('Double')
            self.planets_type_possibilities(n+1, list_planet_type, planet_type_2, number_max_planet)

        return list_planet_type


    ## Recursive Function which give all the possibilities for a junction number fixed ##
    ## Exemple Input:([],[0,0,0,0],2,2,[2,2,2,2]) -> Output:[[2, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 2, 0, 0],
    ## [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 2, 0], [0, 0, 1, 1], [0, 0, 0, 2]] ##
    def junction_possibilities(self, list_possibilities, possibilitie, number_junction_recursive_fonction,
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
                    sum_number_junction_max_by_planet_2[i+1] += self.number_max_junction_by_planet

                possibilitie_2[i] += 1

                self.junction_possibilities(list_possibilities, possibilitie_2, number_junction_recursive_fonction_2+1,
                                            sum_number_junction_max_by_planet_2)

            else:
                flag = 0

        return list_possibilities


    ## Recursive Function which give all the possibilities for a limited number of connection( number_max_connextion) in a list_connection ##
    ## Exemple Input:(0 ,[], [[0, 0], [0, 0]] ,[[2, 4], [2, 5], [3, 6]],2,2) -> Output:[[[2, 4], [2, 5]], [[2, 4], [3, 6]], [[2, 5], [3, 6]]] ##
    def connection_branch_possibilities_step_1(self, n, list_possibilities, possibilitie,
                                               list_connection_branch, number_max_connection):

        possibilitie_2 = copy.copy(possibilitie)
        if  n == number_max_connection:

            if not possibilitie in list_possibilities:
                flag_similaritie = 0

                for possibilitie_3 in list_possibilities:
                    similaritie = 0

                    for connection in possibilitie_2:

                        if connection in possibilitie_3:
                            similaritie += 1

                    if similaritie == number_max_connection:
                        flag_similaritie = 1
                        break

                if not flag_similaritie:
                    list_possibilities.append(possibilitie_2)


        else:
            n += 1
            for element_1 in list_connection_branch:

                number_connection_by_branch = 0
                flag_branch = 1

                for element_2 in possibilitie_2:

                    if element_2[0] == element_1[0]:
                        number_connection_by_branch += 1

                    if element_2[1] == element_1[1]:
                        flag_branch = 0

                if number_connection_by_branch <= self.number_max_junction_by_planet and flag_branch:
                    possibilitie_2[n-1] = element_1

                    self.connection_branch_possibilities_step_1(n, list_possibilities, possibilitie_2,
                                                                list_connection_branch, number_max_connection)
                    possibilitie_2[n-1] = possibilitie[n-1]
        return list_possibilities









    ## Recursive Function which give all the possibilities  of connection  for a number_junction fixed and a number branch fixed ##
    ## Exemple Input:(0 ,[0,0],[],[],[2,3,4,5,6],[1],[2,2],2) -> Output:[[[[1, 2], [1, 3]], [[2, 4], [2, 5]]], [[[1, 2], [1, 3]], [[2, 4], [3, 5]]], [[[1, 2], [1, 3]], [[2, 5], [3, 4]]], [[[1, 2], [1, 3]], [[3, 4], [3, 5]]]] ##
    def connection_branch_possibilities_step_2(self, n, possibilitie, list_possibilities,
                                               list_previous_branch, remaning_branch, previous_branch,
                                               number_junction_branch):

        possibilitie_2 = copy.copy(possibilitie)
        remaning_branch_2 = copy.copy(remaning_branch)
        previous_branch_2 = copy.copy(previous_branch)

        if n == len(number_junction_branch):
            list_possibilities.append(possibilitie_2)
            list_previous_branch.append(previous_branch_2)
            return None

        number_connection = number_junction_branch[n]
        n += 1

        list_possibilities_connection = []
        for element in previous_branch_2:
            for j in range(number_connection):
                list_possibilities_connection.append([element, remaning_branch_2[j]])

        list_possibilities_connection_2 = []
        possibilitie_3 = []

        for i in range(number_connection):

            possibilitie_3.append([0, 0])

        self.connection_branch_possibilities_step_1(0, list_possibilities_connection_2, possibilitie_3,
                                                    list_possibilities_connection, number_connection)



        for element_1 in list_possibilities_connection_2:
            remaning_branch_3 = copy.copy(remaning_branch_2)
            previous_branch_3 = copy.copy(previous_branch_2)

            for i in range(number_connection):
                previous_branch_3.append(remaning_branch_2[i])
                remaning_branch_3.remove(remaning_branch_2[i])

            possibilitie_2[n-1] = element_1

            for element_2 in element_1:

                if element_2[0] in previous_branch_3:
                    previous_branch_3.remove(element_2[0])

            self.connection_branch_possibilities_step_2(n, possibilitie_2, list_possibilities, list_previous_branch,
                                                        remaning_branch_3, previous_branch_3, number_junction_branch)

        return list_possibilities


    ## Recursive Function which give all the possibilities to distribute a number planet fixed in a number_branch fixed  ##
    ## Exemple Input:([0,0,0],[],6,0,3,0,1) -> Output:[[1, 1, 4], [1, 2, 3], [1, 3, 2], [1, 4, 1], [2, 1, 3], [2, 2, 2], [2, 3, 1], [3, 1, 2], [3, 2, 1], [4, 1, 1]] ##
    def planets_by_branch_possibilities_step_1(self, possibilities, list_possibilities,
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
            self.planets_by_branch_possibilities_step_1(possibilities_2, list_possibilities,
                                                        number_planet, number_branch, number_branch_max, number_other_planet+i)



    ##  Function which give all the possibilities of branch for a possibility junction  ##
    ## Exemple Input:([0,1,1],7,3,1) -> Output:[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]],
    ##[[2, 1, 1, 2, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]], [[2, 1, 2, 1, 1], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]]etc... ##
    def planets_by_branch_possibilities_step_2(self, junction):

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

        list_connection = []
        connection = []
        remaning_branch = []
        previous_branch = [1]
        list_number_planet = [0]

        for i in range(number_branch-1):
            remaning_branch.append(i+2)
            list_number_planet.append(0)

        for i in range(len(list_number_junctions)):
            connection.append(0)

        list_previous_branch = []

        self.connection_branch_possibilities_step_2(0, connection, list_connection, list_previous_branch, remaning_branch,
                                                    previous_branch, list_number_junctions)

        list_architecture = []
        for i, connection_2 in enumerate(list_connection):

            number_planet_know = 0
            previous_branch_final = list_previous_branch[i]
            previous_connection_by_branch = []

            for i, connection_by_branch in enumerate(connection_2):

                for element in connection_by_branch:

                    if not element[0] in previous_connection_by_branch:

                        list_number_planet[element[0]-1] = list_number_planet_branch[i]
                        number_planet_know += list_number_planet_branch[i]
                        previous_connection_by_branch.append(element[0])



            number_planet_unknow = self.number_max_planet-number_planet_know
            possibilitie = []

            for i in range(len(previous_branch_final)):
                possibilitie.append(0)

            list_planet_by_branch = []

            self.planets_by_branch_possibilities_step_1(possibilitie, list_planet_by_branch, number_planet_unknow, 0,
                                                        len(previous_branch_final), 0)

            for planet_by_branch in list_planet_by_branch:

                for i, branch in enumerate(previous_branch_final):

                    list_number_planet[branch-1] = planet_by_branch[i]

                list_architecture.append([copy.copy(list_number_planet), connection_2])

        return list_architecture, number_branch



    ##  Recursive Function which give all the possibilities of architecure (branch and planet_type) for a possibility junction  ##
    ## Exemple Input:(1,5,[[[2, 1, 1, 1, 2], [[[1, 2], [1, 3]], [[2, 4], [2, 5]]]],[0,0,0,0,0],[],[],[]]) -> Output: list of PlanetsStructure objet ##
    def architecture_planet_possibilities(self, branch, branch_maximum, architecture, possibilities,
                                          list_possibilities, possibilities_connection, list_possibilities_connection):



        possibilities_connection_2 = copy.copy(possibilities_connection)
        number_planet = architecture[0][branch-1]

        list_planet_possibilities = []

        if number_planet == 1:

            for branch_connection in possibilities_connection_2:
                if branch_connection[1] == branch:
                    list_planet_possibilities = [[branch_connection[2]]]
                    break

        else:
            self.planets_type_possibilities(0, list_planet_possibilities, [], number_planet)


        number_connection = 0
        for connection in architecture[1]:
            for branch_2 in connection:
                if branch_2[0] == branch:
                    number_connection += 1

                    possibilities_connection_2.append([branch, branch_2[1], 'Simple'])



        for planet_possibilities in list_planet_possibilities:

            possibilities_2 = copy.copy(possibilities)
            possibilities_2[branch-1] = planet_possibilities

            if branch == branch_maximum:

                possibilities_connection_3 = copy.deepcopy(possibilities_connection_2)

                planets = []
                connections = []
                first_last_composant_branch = []
                Z = 7

                for branch_2 in possibilities_2:

                    flag_first_planet_branch = 0

                    for planet in branch_2:
                        number_planet += 1
                        planets.append(Planet(Z, 'Pl'+str(number_planet)))
                        Z += 1

                        if flag_first_planet_branch:

                            if planet == 'Double':
                                connections.append(Connection([planets[-2], planets[-1]], 'D'))
                            else:
                                connections.append(Connection([planets[-2], planets[-1]], 'GI'))

                        else:
                            first_composant_branch = planets[-1]
                            flag_first_planet_branch = 1

                    last_composant_branch = planets[-1]
                    first_last_composant_branch.append([first_composant_branch, last_composant_branch])

                for branch_connection in possibilities_connection_3:

                    if branch_connection[2] == 'Simple':
                        connections.append(Connection([first_last_composant_branch[branch_connection[0]-1][1],
                                                       first_last_composant_branch[branch_connection[1]-1][0]], 'GI'))

                    else:
                        connections.append(Connection([first_last_composant_branch[branch_connection[0]-1][1],
                                                       first_last_composant_branch[branch_connection[1]-1][0]], 'D'))

                list_possibilities.append(PlanetsStructure(planets, connections, 'PlanetsStructure'))

            else:

                if number_connection == 0:
                    self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                           list_possibilities, possibilities_connection_2,
                                                           list_possibilities_connection)



                elif number_connection == 2:

                    possibilities_connection_2[-1][2] = 'Double'
                    possibilities_connection_2[-2][2] = 'Simple'

                    self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                           list_possibilities, possibilities_connection_2,
                                                           list_possibilities_connection)

                    possibilities_connection_2[-1][2] = 'Simple'
                    possibilities_connection_2[-2][2] = 'Double'
                    self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                           list_possibilities, possibilities_connection_2,
                                                           list_possibilities_connection)


                    if planet_possibilities[-1] == 'Simple':

                        possibilities_connection_2[-1][2] = 'Double'
                        possibilities_connection_2[-2][2] = 'Double'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                    elif planet_possibilities[-1] == 'Double':

                        possibilities_connection_2[-1][2] = 'Simple'
                        possibilities_connection_2[-2][2] = 'Simple'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)



                elif number_connection == 3:

                    if planet_possibilities[-1] == 'Simple':

                        possibilities_connection_2[-1][2] = 'Double'
                        possibilities_connection_2[-2][2] = 'Double'
                        possibilities_connection_2[-3][2] = 'Simple'

                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                        possibilities_connection_2[-1][2] = 'Simple'
                        possibilities_connection_2[-2][2] = 'Double'
                        possibilities_connection_2[-3][2] = 'Double'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                        possibilities_connection_2[-1][2] = 'Double'
                        possibilities_connection_2[-2][2] = 'Simple'
                        possibilities_connection_2[-3][2] = 'Double'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                    elif planet_possibilities[-1] == 'Double':

                        possibilities_connection_2[-1][2] = 'Simple'
                        possibilities_connection_2[-2][2] = 'Simple'
                        possibilities_connection_2[-3][2] = 'Double'

                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                        possibilities_connection_2[-1][2] = 'Simple'
                        possibilities_connection_2[-2][2] = 'Double'
                        possibilities_connection_2[-3][2] = 'Simple'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)

                        possibilities_connection_2[-1][2] = 'Double'
                        possibilities_connection_2[-2][2] = 'Simple'
                        possibilities_connection_2[-3][2] = 'Simple'
                        self.architecture_planet_possibilities(branch+1, branch_maximum, architecture, possibilities_2,
                                                               list_possibilities, possibilities_connection_2,
                                                               list_possibilities_connection)


    def solution_sort_recursive_function(self, new_first_node, first_node_check, new_graph, graph_check, previous_node,
                                         previous_node_check):
        list_number_false = []
        valid = False

        for new_neighbor in nx.neighbors(new_graph, new_first_node):

            if new_neighbor != previous_node:

                number_false = 0
                if len(list(nx.neighbors(graph_check, first_node_check))) != len(list(nx.neighbors(new_graph, new_first_node))):

                    valid = True
                    return valid

                for neighbor_check in nx.neighbors(graph_check, first_node_check):

                    if len(previous_node_check) < 2 or neighbor_check != previous_node_check[-2]:

                        if type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor]) == \
                           type(nx.get_node_attributes(graph_check, neighbor_check)[neighbor_check]):

                            valid = self.solution_sort_recursive_function(new_neighbor, neighbor_check, new_graph,
                                                                          graph_check, new_first_node,
                                                                          previous_node_check + [neighbor_check])


                        else:
                            number_false += 1



                if number_false == len(list(nx.neighbors(graph_check, first_node_check)))-(len(previous_node_check) > 1):

                    valid = True
                    return valid

                else:
                    list_number_false.append([number_false, type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor])])

        sum_number_false = 0
        list_previous_type_false = []
        for number_false in list_number_false:

            if not number_false[1] in list_previous_type_false:
                sum_number_false += number_false[0]
                list_previous_type_false.append(number_false[1])

        if sum_number_false != 0 and sum_number_false != len(list_number_false):
            valid = True

            return valid

        return valid



    def solution_sort(self, new_planet_structure, planet_structures_check):

        new_graph = new_planet_structure.graph()
        if len(new_planet_structure.doubles)+1 != self.number_max_meshing_plan:
            return False
        for node in nx.nodes(new_graph):

            if len(list(nx.neighbors(new_graph, node))) == 1:
                first_node = node
                break

        for planet_structure in planet_structures_check:


            graph_check = planet_structure.graph()
            possible_first_node_check = []
            for node in nx.nodes(graph_check):

                if len(list(nx.neighbors(graph_check, node))) == 1 and \
                   type(nx.get_node_attributes(graph_check, node)[node]) == \
                   type(nx.get_node_attributes(new_graph, first_node)[first_node]):

                    possible_first_node_check.append(node)


            for node in possible_first_node_check:
                valid = self.solution_sort_recursive_function(first_node, node, new_graph, graph_check, first_node, [node])

                if  not valid:
                    return False

        return True



    def decision_tree(self) -> List[PlanetsStructure]:
        tree = dt.DecisionTree()


        list_possibilities_junction = []
        list_planet = []
        sum_number_max_junction_by_planet = []
        list_solution = []
        for i in range(self.number_max_planet-2):
            list_planet.append(0)
            sum_number_max_junction_by_planet.append(self.number_max_junction_by_planet)

        self.junction_possibilities(list_possibilities_junction, list_planet, 0, sum_number_max_junction_by_planet)

        tree.SetCurrentNodeNumberPossibilities(len(list_possibilities_junction))
        node = tree.NextNode(True)


        while not tree.finished:


            if len(node) == 1:

                list_junction = list_possibilities_junction[node[0]]

                list_global_architecture, number_branch = self.planets_by_branch_possibilities_step_2(list_junction)

                tree.SetCurrentNodeNumberPossibilities(len(list_global_architecture))



            if len(node) == 2:

                global_architecture = list_global_architecture[node[1]]
                list_branch = []

                for i in range(number_branch):
                    list_branch.append(0)


                list_connection = []
                list_possibilities = []
                self.architecture_planet_possibilities(1, number_branch, global_architecture, list_branch,
                                                       list_possibilities, [], list_connection)


                tree.SetCurrentNodeNumberPossibilities(len(list_possibilities))

            if len(node) == 3:

                planet_structure = list_possibilities[node[2]]
                # planet_structure.plot_kinematic_graph()
                if self.solution_sort(planet_structure, list_solution):
                    planet_structure.number_group_solution_planet_structure = len(list_solution)
                    list_solution.append(planet_structure)


                tree.SetCurrentNodeNumberPossibilities(0)


            node = tree.NextNode(True)

        return list_solution

class GeneratorPlanetaryGearsArchitecture(DessiaObject):
    '''
    A generator of architectures of planetary gears

    :param planet_structures: The list of Planets structure with which we want to generate planetary gears architrectures
    :type planet_structures: List[PlanetsStructure]

    :param input_speeds: The list of speed range input
    :type input_speeds: List[List[float]]

    :param name: Name
    :type name: str, optional


    '''
    _standalone_in_db = True

    _eq_is_data_eq = False
    def __init__(self, planet_structures: List[PlanetsStructure], input_speeds: List[Tuple[float, float]], name: str = ''):


        self.planet_structures = planet_structures
        self.number_input = len(input_speeds)

        DessiaObject.__init__(self, name=name)

    def planetaries_possibilities_recursive_function(self, n, planetaries, possibilitie, list_possibilities, meshing_chains,
                                                     meshing_chains_planetary_type, meshing_chains_occupation,
                                                     meshing_chains_planet_index):


        possibilitie_2 = copy.copy(possibilitie)
        list_planetary_type = [['Ring', 'GI', -1], ['Sun', 'GE', 1]]
        meshing_chains_2 = copy.copy(meshing_chains)
        number = 0



        if n == len(planetaries):

            if not 0 in meshing_chains_occupation:
                list_possibilities.append(possibilitie_2)

            return list_possibilities

        planetary = planetaries[n]

        for i, meshing_chain in enumerate(meshing_chains_2):
            number = copy.copy(i)

            if meshing_chains_occupation[i] < 2:

                meshing_chains_occupation_2 = copy.copy(meshing_chains_occupation)
                meshing_chains_occupation_2[i] += 1



                if meshing_chains_planet_index[i] == 0:

                    meshing_chains_planet_index_2 = copy.copy(meshing_chains_planet_index)
                    meshing_chains_planet_index_2[i] = 1

                    for planetary_type in list_planetary_type:
                        if not meshing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            meshing_chain_2 = copy.copy(meshing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            meshing_chains_planetary_type_2 = copy.copy(meshing_chains_planetary_type)
                            meshing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, meshing_chain_2[0], planetary_type[1], number]
                            self.planetaries_possibilities_recursive_function(n+1, planetaries, possibilitie_2,
                                                                              list_possibilities, meshing_chains,
                                                                              meshing_chains_planetary_type_2,
                                                                              meshing_chains_occupation_2,
                                                                              meshing_chains_planet_index_2)


                    meshing_chains_planet_index_2[i] = -1


                    for planetary_type in list_planetary_type:

                        if not meshing_chains_planetary_type[i] == planetary_type[0]:
                            planetary_2 = copy.copy(planetary)
                            meshing_chain_2 = copy.copy(meshing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            meshing_chains_planetary_type_2 = copy.copy(meshing_chains_planetary_type)
                            meshing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, meshing_chain_2[-1], planetary_type[1], number]
                            self.planetaries_possibilities_recursive_function(n+1, planetaries, possibilitie_2,
                                                                              list_possibilities, meshing_chains,
                                                                              meshing_chains_planetary_type_2,
                                                                              meshing_chains_occupation_2,
                                                                              meshing_chains_planet_index_2)



                elif meshing_chains_planet_index[i] == 1:

                    meshing_chains_planet_index_2 = copy.copy(meshing_chains_planet_index)
                    meshing_chains_planet_index_2[i] = -1

                    for planetary_type in list_planetary_type:

                        if not meshing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            meshing_chain_2 = copy.copy(meshing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            meshing_chains_planetary_type_2 = copy.copy(meshing_chains_planetary_type)
                            meshing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, meshing_chain_2[-1], planetary_type[1], number]
                            self.planetaries_possibilities_recursive_function(n+1, planetaries, possibilitie_2,
                                                                              list_possibilities, meshing_chains,
                                                                              meshing_chains_planetary_type_2,
                                                                              meshing_chains_occupation_2,
                                                                              meshing_chains_planet_index_2)


                else:
                    meshing_chains_planet_index_2 = copy.copy(meshing_chains_planet_index)
                    meshing_chains_planet_index_2[i] = -1

                    for planetary_type in list_planetary_type:
                        if not meshing_chains_planetary_type[i] == planetary_type[0]:

                            planetary_2 = copy.copy(planetary)
                            meshing_chain_2 = copy.copy(meshing_chain)

                            planetary_2.planetary_type = planetary_type[0]
                            planetary_2.p = planetary_type[2]

                            meshing_chains_planetary_type_2 = copy.copy(meshing_chains_planetary_type)
                            meshing_chains_planetary_type_2[i] = planetary_type[0]

                            possibilitie_2[n] = [planetary_2, meshing_chain_2[0], planetary_type[1], number]
                            self.planetaries_possibilities_recursive_function(n+1, planetaries, possibilitie_2,
                                                                              list_possibilities, meshing_chains,
                                                                              meshing_chains_planetary_type_2,
                                                                              meshing_chains_occupation_2,
                                                                              meshing_chains_planet_index_2)

        return list_possibilities



    def planetaries_possibilities(self, planetaries, planets_structure, planet_carrier):
        meshing_chains = planets_structure.meshing_chain()



        if len(planetaries) < len(meshing_chains):
            return []

        possibilitie = []
        connection = []
        meshing_chains_planetary_type = []
        meshing_chains_planet_index = []
        meshing_chains_occupation = []

        for i in enumerate(planetaries):
            possibilitie.append(0)

        for i in enumerate(meshing_chains):
            meshing_chains_planetary_type.append(0)
            meshing_chains_occupation.append(0)
            meshing_chains_planet_index.append(0)

        list_possibilities = self.planetaries_possibilities_recursive_function(0, planetaries, possibilitie, [],
                                                                               meshing_chains, meshing_chains_planetary_type,
                                                                               meshing_chains_occupation,
                                                                               meshing_chains_planet_index)


        for double in planets_structure.doubles:

            connection.append(Connection([double.nodes[0], double.nodes[1]], 'D'))

        list_solution = []

        for i, possibilitie_2 in enumerate(list_possibilities):

            connection_2 = copy.copy(connection)
            previous_planetary = []

            for planetary_connection in possibilitie_2:
                connection_2.append(Connection([planetary_connection[:-1][0], planetary_connection[:-1][1]], planetary_connection[:-1][2]))

                if not planetary_connection[3]  in previous_planetary:

                    meshing_chain = meshing_chains[planetary_connection[3]]

                    previous_planetary.append(planetary_connection[3])

                    if len(meshing_chain) > 1:

                        if meshing_chain[0] == planetary_connection[1]:
                            for planet_1, planet_2 in zip(meshing_chain[:-1], meshing_chain[1:]):

                                connection_2.append(Connection([planet_1, planet_2], planetary_connection[2]))


                        elif meshing_chain[-1] == planetary_connection[1]:

                           meshing_chain_2 = meshing_chain[::-1]

                           for planet_1, planet_2 in zip(meshing_chain_2[:-1], meshing_chain_2[1:]):

                                connection_2.append(Connection([planet_1, planet_2], planetary_connection[2]))


            planetary_gear = PlanetaryGear(planetaries, planets_structure.planets, planet_carrier, connection_2, 3, 'PlanetaryGear'+str(i))
            list_path = planetary_gear.path_planetary_to_planetary()
            list_planets = []
            # print(planetary_gear.planets)
            # print(list_path)
            for planet in planetary_gear.planets:
                list_planets.append(planet.name)

            for path in list_path:
                for element in path:

                    if element.name in list_planets:
                        list_planets.remove(element.name)


            if not list_planets:

                list_solution.append(PlanetaryGear(copy.deepcopy(planetaries),
                                                   copy.deepcopy(planets_structure.planets), copy.copy(planet_carrier), connection_2, 'PlanetaryGear'+str(i)))



        return list_solution


    def solution_sort_recursive_function(self, new_first_node, first_node_check, new_graph, graph_check, previous_node,
                                         previous_node_check):

        list_number_false = []
        list_neighbors = list(nx.neighbors(graph_check, first_node_check))

        if len(list(nx.neighbors(graph_check, first_node_check))) != len(list(nx.neighbors(new_graph, new_first_node))):

                    valid = True
                    return valid

        for new_neighbor in nx.neighbors(new_graph, new_first_node):

            valid = False

            if new_neighbor != previous_node:
                number_false = 0


                for neighbor_check in list_neighbors:

                    if (len(previous_node_check) < 2 or neighbor_check != previous_node_check[-2]):

                        if type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor]) == \
                           type(nx.get_node_attributes(graph_check, neighbor_check)[neighbor_check]):

                            if type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor]) == Planetary:


                                if nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor].planetary_type == \
                                   nx.get_node_attributes(graph_check, neighbor_check)[neighbor_check].planetary_type:


                                    valid = self.solution_sort_recursive_function(new_neighbor, neighbor_check, new_graph, graph_check,
                                                                                  new_first_node, previous_node_check + [neighbor_check])

                                    if not valid:
                                        list_neighbors.remove(neighbor_check)

                                    else:
                                        number_false += 1


                                else:

                                    number_false += 1

                            else:
                                valid = self.solution_sort_recursive_function(new_neighbor, neighbor_check, new_graph, graph_check,
                                                                              new_first_node, previous_node_check + [neighbor_check])

                                if not valid:

                                        list_neighbors.remove(neighbor_check)
                                else:
                                    number_false += 1

                        else:
                            number_false += 1



                if number_false >= len(list(nx.neighbors(graph_check, first_node_check)))-(len(previous_node_check) > 1):

                    valid = True
                    return valid


                else:
                    if type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor]) == Planetary:

                            list_number_false.append([number_false,
                                                      nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor].planetary_type])

                    else:
                        list_number_false.append([number_false, type(nx.get_node_attributes(new_graph, new_neighbor)[new_neighbor])])


        return valid



    def solution_sort(self, new_planetary_gear, planetary_gears_check):

        new_graph = new_planetary_gear.graph()
        new_graph.remove_node(str(new_planetary_gear.planet_carrier))

        for i, planet in enumerate(new_planetary_gear.planets):
            new_graph.remove_node('Pv'+str(i))

        for node in nx.nodes(new_graph):

            if len(list(nx.neighbors(new_graph, node))) == 1:
                first_node = node
                break

        list_valid = []
        list_possible_node = []
        for planetary_gear in planetary_gears_check:


            graph_check = planetary_gear.graph()
            graph_check.remove_node(str(planetary_gear.planet_carrier))

            for i, planet in enumerate(planetary_gear.planets):
                graph_check.remove_node('Pv'+str(i))

            possible_first_node_check = []

            for node in nx.nodes(graph_check):

                if len(list(nx.neighbors(graph_check, node))) == 1 and type(nx.get_node_attributes(graph_check, node)[node]) == \
                   type(nx.get_node_attributes(new_graph, first_node)[first_node]):

                    if type(nx.get_node_attributes(graph_check, node)[node]) == Planetary:

                        if nx.get_node_attributes(graph_check, node)[node].planetary_type ==  \
                           nx.get_node_attributes(new_graph, first_node)[first_node].planetary_type:

                           possible_first_node_check.append(node)

                    else:
                        possible_first_node_check.append(node)


            for node in possible_first_node_check:
                valid = self.solution_sort_recursive_function(first_node, node, new_graph, graph_check, first_node, [node])
                list_valid.append(valid)

            if possible_first_node_check:
                list_possible_node.append([possible_first_node_check])


        if  False in list_valid:
            return False
            return False


        return True

    def decision_tree(self) -> List[PlanetaryGear]:
        tree = dt.DecisionTree()


        planet_carrier = PlanetCarrier('PlanetCarrier')
        planetaries = []
        n = 0
        list_solution = []
        list_paths_type = []

        for i in range(self.number_input-1):
            planetaries.append(Planetary(7, 'Sun', 'Planetary_'+str(i)))

        tree.SetCurrentNodeNumberPossibilities(len(self.planet_structures))
        # print(len(self.planet_structures))
        node = tree.NextNode(True)


        while not tree.finished:
            valid = True

            if len(node) == 1:
             planet_architecture = self.planet_structures[node[0]]
             list_planetary_gears = self.planetaries_possibilities(planetaries, planet_architecture, planet_carrier)
             for planetary_gear in list_planetary_gears:
                 planetary_gear.number_group_solution_planet_structure = node[0]

             tree.SetCurrentNodeNumberPossibilities(len(list_planetary_gears))
             list_check = []

            if len(node) == 2:

                planetary_gear = list_planetary_gears[node[1]]


                if self.solution_sort(planetary_gear, list_check):
                    planetary_gear.number_group_solution_architecture = len(list_solution)
                    list_solution.append(planetary_gear)
                    list_check.append(planetary_gear)
                tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(valid)

        return list_solution

class GeneratorPlanetaryGearsZNumber(DessiaObject):
    '''
    A generator of all the number of tooth in a planetary gear

    :param planetary_gear: The planetary gears that we want to generate all the number of tooth possible
    :type planetary_gear: PlanetaryGear

    :param input_speeds: The list of speed range input
    :type input_speeds: List[List[float]]

    :param Z_range_sun: The range of number tooth that can take a normal gear
    :type Z_range_sun: List[int]

    :param Z_range_ring: The range of number tooth that can take a ring
    :type Z_range_ring: List[int]

    :param number_planet: The number of planet which are arround the planetary gear ( exemple: 3,4 or 5)
    :type number_planet: int

    :param name: Name
    :type name: str, optional


    '''
    _standalone_in_db = True

    _eq_is_data_eq = False

    def __init__(self, planetary_gear: PlanetaryGear, input_speeds: List[Tuple[float, float]], input_torques: List[Tuple[float, float]],
                 Z_range_sun: List[int], Z_range_ring: List[int], number_planet: int, speed_max_planet: float = 1000000, name: str = ''):


        self.planetary_gear = planetary_gear
        self.planetary_gear.number_branch_planet = number_planet

        self.input_speeds = input_speeds
        self.input_torques = input_torques
        self.speed_max_planet = speed_max_planet
        self.number_input = len(input_speeds)
        self.Z_range_sun = Z_range_sun
        self.Z_range_ring = Z_range_ring
        self.number_planet = number_planet
        DessiaObject.__init__(self, name=name)

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



    # def test_speed_precision(self,planetary_gear):
    #     for planetary in planetary_gear.planetaries:
    #         range_speed = planetary_gear.speed_range(planetary, planetary_gear.planet_carrier, self.precision)

    #         if not range_speed:
    #             return False
    #     list_planetary=copy.copy(planetary_gear.planetaries)
    #     for planetary in planetary_gear.planetaries:
    #         list_planetary.remove(planetary)
    #         for planetary_2 in list_planetary:
    #             range_speed = planetary_gear.speed_range(planetary, planetary_2, self.precision)

    #             if not range_speed :
    #                 return False

    #     return True



    def test_GCD(self, Z_1, Z_2):

        if m.gcd(Z_1, Z_2) != 1:

                     return False

        return True


    def test_vitesse_and_assembly_condition(self, planetary_gear, begin_meshing_chain, end_meshing_chain,
                                            list_previous_planetary, list_path):

        list_previous_planetary_2 = copy.copy(list_previous_planetary)
        list_previous_planetary_2.remove(begin_meshing_chain)


        range_speed = planetary_gear.speed_range(begin_meshing_chain, planetary_gear.planet_carrier, list_previous_planetary, 1, list_path)

        if range_speed == 'simplex':


        #x=[speed_diff_1,speed_diff_2,speed_1,speed_2]
            c = [-1, 0, 0, 0]
            A = [[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, 1], [0, 1, 0, -1]]
            b = [begin_meshing_chain.speed_input[1], -begin_meshing_chain.speed_input[0],
                 planetary_gear.planet_carrier.speed_input[1], -planetary_gear.planet_carrier.speed_input[0]]
            speed_diff_1_bound = (0, None)
            speed_diff_2_bound = (0, None)
            speed_1_bound = (None, None)
            speed_2_bound = (None, None)

            for i, planetary in enumerate(list_previous_planetary_2):
                speed_input_planetary = planetary.speed_input
                # reason=planetary_gear.reason(list_path[i][0])
                path = list_path[i]
                reason = planetary_gear.reason(path[0])

                if reason < 0:
                    A.extend([[-reason, (1-reason), -reason, -(1-reason)], [-reason, 1-reason, reason, 1-reason]])
                else:

                    if reason < 1:
                        A.extend([[reason, 1-reason, -reason, -(1-reason)], [reason, 1-reason, reason, 1-reason]])
                    else:
                        A.extend([[reason, -(1-reason), -reason, -(1-reason)], [reason, -(1-reason), reason, (1-reason)]])

                b.extend([-speed_input_planetary[0], speed_input_planetary[1]])

            res = op.linprog(c, A_ub=A, b_ub=b, bounds=[speed_diff_1_bound, speed_diff_2_bound, speed_1_bound, speed_2_bound])

            if not res.success:
                return False


        elif not range_speed:
            return False


        elif range_speed[begin_meshing_chain][0] > range_speed[begin_meshing_chain][1]:
            return False

        elif range_speed[planetary_gear.planet_carrier][0] > range_speed[planetary_gear.planet_carrier][1]:
            return False





        if not planetary_gear.test_assembly_condition(self.number_planet, [begin_meshing_chain, end_meshing_chain]):
            return False

        return True

    def test_torque(self, planetary_gear, first_planetary, list_path, list_previous_planetaries):
        num_var = (len(list_previous_planetaries))*2
        num_eq = (len(list_previous_planetaries)+2)*2
        A = np.zeros((num_eq, num_var))
        y = 0
        b = []
        c = [0]*num_var
        c[1] = -1
        second_planetary = list_previous_planetaries[0]

        reason_second_planetary = planetary_gear.reason(list_path[0][0])
        list_previous_planetaries.remove(second_planetary)
        list_path.remove(list_path[0])
        bounds = []

        for i in range(int(num_var/2)):
            A[y][2*i] = 1
            A[y][2*i+1] = 1
            y += 1
            A[y][2*i] = -1
            A[y][2*i+1] = 1
            y += 1

            bounds.extend([(None, None), (0, None)])
            if i != num_var/2-1:
                b.extend([list_previous_planetaries[i].torque_input[1], -list_previous_planetaries[i].torque_input[0]])
            else:
                position_planet_carrier = 2*i
                b.extend([planetary_gear.planet_carrier.torque_input[1], -planetary_gear.planet_carrier.torque_input[0]])






        b.extend([first_planetary.torque_input[1], -first_planetary.torque_input[0], second_planetary.torque_input[1], -second_planetary.torque_input[0]])

        for i, planetary in enumerate(list_previous_planetaries):
            reason_planetary = planetary_gear.reason(list_path[i][0])

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

        if reason_second_planetary == 1:
            coefficient_1 = 10000
        else:

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

        if reason_second_planetary == 1:
            coefficient_2 = 10000
        else:
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

        res = op.linprog(c, A_ub=A, b_ub=b, bounds=bounds)

        if res.success:
            print(1)
            return True
        else:
            return False






    def z_range_mini_max(self, planetary_gear, element, begin_meshing_chain, end_meshing_chain, path, reasons_min_max):


        if not element in path[0]:

            return []
        if not element in reasons_min_max.keys():



            if end_meshing_chain.speed_input[0] != planetary_gear.planet_carrier.speed_input[0]:
                reason_1 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0])/
                               (begin_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0]))

                reason_2 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0])/
                               (begin_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0]))

            else:
                if end_meshing_chain.speed_input[1] == int(end_meshing_chain.speed_input[1]):
                    replace_zero = 0.001
                else:
                    replace_zero = (end_meshing_chain.speed_input[1]-int(end_meshing_chain.speed_input[1]))*0.001

                reason_1 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0])/
                               (replace_zero))
                reason_2 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0])/
                               (replace_zero))




            if end_meshing_chain.speed_input[1] != planetary_gear.planet_carrier.speed_input[0]:
                reason_3 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0])/
                               (begin_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0]))
                reason_4 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0])/
                               (begin_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0]))


            else:
                if end_meshing_chain.speed_input[0] == int(end_meshing_chain.speed_input[0]):
                   replace_zero = 0.001
                else:
                    replace_zero = (end_meshing_chain.speed_input[0]-int(end_meshing_chain.speed_input[0]))*0.001

                reason_3 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[0])/
                               (replace_zero))
                reason_4 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[0])/
                               (replace_zero))



            if end_meshing_chain.speed_input[0] != planetary_gear.planet_carrier.speed_input[1]:
                reason_5 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1])/
                               (begin_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1]))
                reason_6 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1])/
                               (begin_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1]))

            else:
                if end_meshing_chain.speed_input[0] == int(end_meshing_chain.speed_input[0]):
                   replace_zero = 0.001
                else:
                    replace_zero = (end_meshing_chain.speed_input[0]-int(end_meshing_chain.speed_input[0]))*0.001

                reason_5 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1])/
                               (replace_zero))
                reason_6 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1])/
                               (replace_zero))


            if end_meshing_chain.speed_input[1] != planetary_gear.planet_carrier.speed_input[1]:
                reason_7 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1])/
                               (begin_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1]))

                reason_8 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1])/
                               (begin_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1]))

            else:
                if end_meshing_chain.speed_input[1] == int(end_meshing_chain.speed_input[1]):
                   replace_zero = 0.001
                else:
                    replace_zero = (end_meshing_chain.speed_input[1]-int(end_meshing_chain.speed_input[1]))*0.001

                reason_7 = abs((end_meshing_chain.speed_input[0]-planetary_gear.planet_carrier.speed_input[1])/
                               (replace_zero))
                reason_8 = abs((end_meshing_chain.speed_input[1]-planetary_gear.planet_carrier.speed_input[1])/
                               (replace_zero))


            reason_min = min(reason_1, reason_2, reason_3, reason_4, reason_5, reason_6, reason_7, reason_8)
            reason_max = max(reason_1, reason_2, reason_3, reason_4, reason_5, reason_6, reason_7, reason_8)
            reasons_min_max[element] = [reason_min, reason_max]


        reason = planetary_gear.reason_abs(path[0])
        reason_min = reasons_min_max[element][0]
        reason_max = reasons_min_max[element][1]

        if isinstance(element, Planetary):
            if reason_min and reason_max:
                Z_min = int(reason*element.Z/reason_max)-1

            else:
                Z_min = 0

            if reason_max:
                Z_max = int(reason*element.Z/reason_min)+1

            else:
                Z_max = self.Z_range_sun[1]
        else:
            if reason_min and reason_max:
                Z_min = int((element.Z*reason_min)/reason)-1

            else:
                Z_min = 0

            if reason_max:
                Z_max = int((element.Z*reason_max)/reason)+1

            else:
                Z_max = self.Z_range_sun[1]


        Z_range_mini_maxi = [Z_min, Z_max]

        return Z_range_mini_maxi

    def decision_tree_speed_possibilities(self):
        tree = dt.DecisionTree()
        tree.SetCurrentNodeNumberPossibilities(self.number_input)
        node = tree.NextNode(True)
        list_solution = []
        while not tree.finished:
            planetary_gear = copy.deepcopy(self.planetary_gear)

            if len(node) == 1:


                    input_speeds_2 = copy.deepcopy(self.input_speeds)
                    input_speeds_2.remove(self.input_speeds[node[0]])

                    list_possibility_speed = []
                    self.multiplication_possibility_speed(input_speeds_2, 0, [], list_possibility_speed)

                    tree.SetCurrentNodeNumberPossibilities(len(list_possibility_speed))


            if len(node) == 2:

                    possibility_speed = list_possibility_speed[node[1]]
                    for i, planetary in enumerate(planetary_gear.planetaries):
                        planetary.speed_input = possibility_speed[i]
                        planetary.torque_input = self.input_torques[self.input_speeds.index(possibility_speed[i])]
                    planetary_gear.planet_carrier.speed_input = self.input_speeds[node[0]]
                    planetary_gear.planet_carrier.torque_input = self.input_torques[node[0]]
                    list_solution.append(planetary_gear)

                    tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(True)

        return list_solution

    def decision_tree(self) -> List[PlanetaryGear]:
        list_planetary_gears_speed = self.decision_tree_speed_possibilities()
        #print(list_planetary_gears_speed[0].planet_carrier.speed_input)
        list_solution = []
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
            meshing_chains_modif = planetary_gear.meshing_chain()
            meshing_chains = copy.copy(meshing_chains_modif)
            number_element_meshing_chain = []
            numbers_meshing_chain = []
            number_meshing_chain = 0
            totals_element_previous_meshing_chain = []
            total_element_previous_meshing_chain = 0
            flags_meshing_change = []

            flag_gcd = []
            numbers_planetaries_by_meshing_chain = []


            for i, meshing_chain in enumerate(meshing_chains_modif):
                if isinstance(meshing_chain[-1], Planetary):

                    if meshing_chain[-1].planetary_type == 'Ring':

                        meshing_chains_modif[i] = meshing_chain[::-1]
                        meshing_chains[i] = meshing_chain[::-1]

                    if  not isinstance(meshing_chain[0], Planetary) and meshing_chain[-1].planetary_type == 'Sun':

                        meshing_chains_modif[i] = meshing_chain[::-1]
                        meshing_chains[i] = meshing_chain[::-1]
                meshing_chain_2 = copy.copy(meshing_chain)
                number_planetaries = 0
                for element in meshing_chain_2:
                    if isinstance(element, Planetary):
                        number_planetaries += 1
                    if element in list_planet_remove:

                        meshing_chains_modif[i].remove(element)
                numbers_planetaries_by_meshing_chain.append(number_planetaries)


            print(numbers_planetaries_by_meshing_chain)
            if numbers_planetaries_by_meshing_chain[0] == 1:
                if 2 in numbers_planetaries_by_meshing_chain:
                    meshing_chain_1 = meshing_chains_modif[numbers_planetaries_by_meshing_chain.index(2)]
                    meshing_chains_modif[numbers_planetaries_by_meshing_chain.index(2)] = meshing_chains_modif[0]
                    meshing_chains_modif[0] = meshing_chain_1



            for meshing_chain in meshing_chains_modif:

                number_element_meshing_chain.append(len(meshing_chain))
                flags_meshing_change.append(1)

                for i, element in enumerate(meshing_chain):
                    flag_gcd.append(2)

                    if i != 0:
                        flags_meshing_change.append(0)

                    totals_element_previous_meshing_chain.append(total_element_previous_meshing_chain)
                    numbers_meshing_chain.append(number_meshing_chain)

                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':

                        list_tree.append(self.Z_range_ring[1]-self.Z_range_ring[0])
                        list_node_range_data.append(self.Z_range_ring[0])

                    else:
                        list_tree.append(self.Z_range_sun[1]-self.Z_range_sun[0])
                        list_node_range_data.append(self.Z_range_sun[0])

                number_meshing_chain += 1
                total_element_previous_meshing_chain += len(meshing_chain)

            list_planet_remove_neighbour = []


            for i, planet in enumerate(list_planet_remove):
                planet.Z = 1
                list_planet_remove_neighbour.append([planet])

                for meshing in planetary_gear.meshings:

                    if meshing.nodes[0] == planet:
                        list_planet_remove_neighbour[i].append(meshing.nodes[1])

                    if meshing.nodes[1] == planet:
                        list_planet_remove_neighbour[i].append(meshing.nodes[0])


            tree = dt.RegularDecisionTree(list_tree)

            Z_range_mini_maxi = []
            Z_range_mini_maxi_2 = []
            flag_meshing_change = 0
            flag_Z_range_mini_maxi = 0
            flag_Z_range_mini_maxi_2 = 0
            number_max_z_planet = [self.Z_range_sun[1]]*len(meshing_chains_modif)
            list_planetaries_Z_range_mini_maxi = []
            list_path = []
            reason_min_max = {}

            while not tree.finished:

                valid = True
                node = tree.current_node

                number_meshing_chain = numbers_meshing_chain[len(node)-1]

                flag_meshing_change = flags_meshing_change[len(node)-1]
                total_element_previous_meshing_chain = totals_element_previous_meshing_chain[len(node)-1]

                element = meshing_chains_modif[number_meshing_chain][len(node)-total_element_previous_meshing_chain-1]
                element.Z = list_node_range_data[len(node)-1]+ node[len(node)-1]


                if len(node) == 1:


                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':
                        number_max_z_planet[number_meshing_chain] = element.Z





                elif not flag_meshing_change:


                    previous_element = meshing_chains_modif[number_meshing_chain][len(node)-total_element_previous_meshing_chain-2]

                    if len(meshing_chains_modif[number_meshing_chain]) > 2 \
                    and meshing_chains_modif[number_meshing_chain][0].planetary_type == 'Ring':

                        number_tot_z_previous_planets = 0
                        number_max_z_previous_planets = 0
                        number_planet = 0
                        for i in range(len(node)-total_element_previous_meshing_chain-2):
                            previous_planet = meshing_chains_modif[number_meshing_chain][i+1]

                            number_planet += 1
                            number_tot_z_previous_planets += previous_planet.Z
                            if previous_planet.Z > number_max_z_previous_planets:
                                number_max_z_previous_planets = previous_planet.Z



                        if isinstance(element, Planetary):

                            flag_impose_z_number = True

                            for planet in list_planet_remove:
                              if planet in meshing_chains[number_meshing_chain]:
                                  flag_impose_z_number = False
                                  break

                            if flag_impose_z_number and number_planet == 1:
                                    if element.Z != (number_max_z_planet[number_meshing_chain]-number_max_z_previous_planets*2):

                                        valid = False

                            else:
                                if element.Z >= (number_max_z_planet[number_meshing_chain]-number_max_z_previous_planets*2) \
                                or element.Z < (number_max_z_planet[number_meshing_chain]-number_tot_z_previous_planets*2):

                                        valid = False
                        else:

                            if element.Z >= (number_max_z_planet[number_meshing_chain])/2:

                                valid = False
                    else:
                        if element.Z >= (number_max_z_planet[number_meshing_chain])/2:

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
                        number_max_z_planet[number_meshing_chain] = element.Z


                    else:
                        number_max_z_planet[number_meshing_chain] = (self.Z_range_sun[1])*2


                if len(node) == number_element_meshing_chain[number_meshing_chain]+total_element_previous_meshing_chain and valid:

                    begin_meshing_chain = meshing_chains_modif[number_meshing_chain][0]
                    end_meshing_chain = meshing_chains_modif[number_meshing_chain][-1]

                    #planetary_gear = PlanetaryGear(planetary_gear.planetaries, planetary_gear.planets,
                                                   #planetary_gear.planet_carrier, planetary_gear.connections, planetary_gear.name)

                    if Z_range_mini_maxi:

                        if element.Z < Z_range_mini_maxi[0] or element.Z > Z_range_mini_maxi[1]:
                            valid = False

                    if Z_range_mini_maxi_2:

                        if element.Z < Z_range_mini_maxi_2[0]  or element.Z > Z_range_mini_maxi_2[1]:
                            valid = False

                    if valid:

                        if numbers_meshing_chain[len(node)-1] == 0:

                                first_planetary = begin_meshing_chain

                                if not first_planetary in list_planetaries_Z_range_mini_maxi:
                                    list_planetaries_Z_range_mini_maxi.append(first_planetary)


                                if isinstance(begin_meshing_chain, Planetary) and isinstance(end_meshing_chain, Planetary):

                                    if not end_meshing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([begin_meshing_chain, end_meshing_chain]))

                                    if not flag_Z_range_mini_maxi:


                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element, begin_meshing_chain,
                                                                                  end_meshing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)-1],
                                                                                  reason_min_max)

                                        flag_Z_range_mini_maxi = 1



                                        if element.Z > Z_range_mini_maxi[0] and element.Z < Z_range_mini_maxi[1]:
                                            list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                            list_path_speed_test = []
                                            for previous_planetary in list_previous_planetary[1:]:
                                                list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                            valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                             end_meshing_chain,
                                                                                             list_previous_planetary,
                                                                                             list_path_speed_test)

                                        else:
                                            valid = False
                                    else:
                                        list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                        list_path_speed_test = []
                                        for previous_planetary in list_previous_planetary[1:]:
                                            list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                        valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                         end_meshing_chain,
                                                                                         list_previous_planetary,
                                                                                         list_path_speed_test)



                        else:

                            if not begin_meshing_chain in list_planetaries_Z_range_mini_maxi:


                                        list_planetaries_Z_range_mini_maxi.append(begin_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, begin_meshing_chain]))

                            if not flag_Z_range_mini_maxi:

                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element,
                                                                                  first_planetary, begin_meshing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)-1],
                                                                                  reason_min_max)

                                        flag_Z_range_mini_maxi = 1

                                        if Z_range_mini_maxi:
                                            if element.Z > Z_range_mini_maxi[0] and element.Z < Z_range_mini_maxi[1]:
                                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)+1]
                                                list_path_speed_test = []
                                                for previous_planetary in list_previous_planetary[1:]:
                                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                                 begin_meshing_chain,
                                                                                                 list_previous_planetary,
                                                                                                 list_path_speed_test)

                                            else:
                                                valid = False

                            else:
                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)+1]
                                list_path_speed_test = []
                                for previous_planetary in list_previous_planetary[1:]:
                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                 begin_meshing_chain,
                                                                                 list_previous_planetary,
                                                                                 list_path_speed_test)




                            if isinstance(end_meshing_chain, Planetary) and valid:


                                if not end_meshing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, end_meshing_chain]))


                                if not flag_Z_range_mini_maxi_2:

                                        Z_range_mini_maxi_2 = self.z_range_mini_max(planetary_gear, element, first_planetary, end_meshing_chain,
                                                                                    list_path[list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)-1], reason_min_max)
                                        flag_Z_range_mini_maxi_2 = 1
                                        if Z_range_mini_maxi_2:

                                            if element.Z > Z_range_mini_maxi[0] and element.Z < Z_range_mini_maxi[1]:

                                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)]
                                                list_path_speed_test = []

                                                for previous_planetary in list_previous_planetary[1:]:
                                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                                 end_meshing_chain,
                                                                                                 list_previous_planetary,
                                                                                                 list_path_speed_test)

                                            else:
                                                valid = False
                                else:
                                    list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                    list_path_speed_test = []
                                    for previous_planetary in list_previous_planetary[1:]:
                                        list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                    valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                     end_meshing_chain,
                                                                                     list_previous_planetary,
                                                                                     list_path_speed_test)






                        if number_meshing_chain == len(meshing_chains_modif)-1 and valid:
                            list_previous_planetary = list_planetaries_Z_range_mini_maxi
                            list_path_torque_test = []

                            for previous_planetary in list_previous_planetary[1:]:
                                        list_path_torque_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                            if not self.test_torque(planetary_gear, first_planetary, list_path_torque_test, list_previous_planetary[1:]):
                                valid = False

                            list_tree_planetary = []





                            if list_planet_remove and valid:

                                for i in range(len(list_planet_remove)):
                                    list_planet_remove_2 = copy.copy(list_planet_remove)
                                    list_planet_remove_neighbour_2 = copy.copy(list_planet_remove_neighbour)
                                    list_planet_remove[i].Z = 1
                                    valid_planet = True

                                    for meshing_chain in meshing_chains:
                                        if list_planet_remove[i] in meshing_chain and valid_planet:

                                            if meshing_chain[0].planetary_type == 'Ring':
                                                number_planet = 0
                                                for element_2 in meshing_chain:
                                                    if  isinstance(element_2, Planet):
                                                        number_planet += 1

                                                if isinstance(meshing_chain[-1], Planetary) and number_planet == 1:
                                                    number_z_max = (meshing_chain[0].Z - meshing_chain[-1].Z)/2

                                                    if number_z_max < self.Z_range_sun[0]:
                                                        list_planet_remove[i].Z = self.Z_range_sun[0]
                                                        valid_planet = False
                                                    neighbour = list_planet_remove_neighbour[i]
                                                    if not self.test_GCD(number_z_max, neighbour[1].Z):
                                                        valid_planet = False

                                                        break
                                                    if not self.test_GCD(number_z_max, neighbour[2].Z):
                                                        valid_planet = False
                                                        break
                                                    list_planet_remove[i].Z = number_z_max
                                                    list_planet_remove_2.remove(list_planet_remove[i])
                                                    list_planet_remove_neighbour_2.remove(neighbour)

                                                else:
                                                    if isinstance(meshing_chain[-1], Planetary):
                                                        number_z_max = int((meshing_chain[0].Z- meshing_chain[-1].Z)/2)+1
                                                    else:
                                                        number_z_max = meshing_chain[0].Z/2


                                                    if number_z_max - self.Z_range_sun[0] > 0:
                                                        list_tree_planetary.append(number_z_max - self.Z_range_sun[0])

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
                                            planetary_gear.speed_max_planet = planetary_gear.speed_max_planets()
                                            if planetary_gear.speed_max_planet < self.speed_max_planet:

                                                print(planetary_gear)
                                                list_solution.append(copy.deepcopy(planetary_gear))

                                            # print(planetary_gear)
                                            # if list_solution:
                                            #     return list_solution

                                        tree_planet.NextNode(valid_planet)

                            else:
                                planetary_gear.speed_max_planet = planetary_gear.speed_max_planets()
                                if valid and planetary_gear.speed_max_planet < self.speed_max_planet:

                                    list_solution.append(copy.deepcopy(planetary_gear))
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


        # a=0
        # for i in range(len(list_solution)):
        #     if not self.test_speed_precision(list_solution[i+a]):
        #         list_solution.pop(i+a)
        #         a-=1
        return list_solution


class GeneratorPlanetaryGearsZNumberReason(DessiaObject):
    '''
    A generator of all the number of tooth in a planetary gear

    :param planetary_gear: The planetary gears that we want to generate all the number of tooth possible
    :type planetary_gear: PlanetaryGear

    :param input_speeds: The list of speed range input
    :type input_speeds: List[List[float]]

    :param Z_range_sun: The range of number tooth that can take a normal gear
    :type Z_range_sun: List[int]

    :param Z_range_ring: The range of number tooth that can take a ring
    :type Z_range_ring: List[int]

    :param number_planet: The number of planet which are arround the planetary gear ( exemple: 3,4 or 5)
    :type number_planet: int

    :param name: Name
    :type name: str, optional


    '''
    _standalone_in_db = True

    _eq_is_data_eq = False

    def __init__(self, planetary_gear: PlanetaryGear, input_reason: List[Tuple[float, float]],
                 input_speeds: List[Tuple[float, float]], input_torques: List[Tuple[float, float]],
                 speed_planet_carrer: Tuple[float, float], torque_planet_carrer: Tuple[float, float],
                 Z_range_sun: List[int], Z_range_ring: List[int], number_planet: int, speed_max_planet: float = 1000000, name: str = '',
                 number_solution: int = 0, different_solution: bool = False):

        self.input_reason = input_reason
        self.planetary_gear = planetary_gear
        self.planetary_gear.number_branch_planet = number_planet
        self.input_speeds = input_speeds
        self.input_torques = input_torques
        self.speed_max_planet = speed_max_planet
        self.number_planetary_reason_input = 0
        self.number_solution = number_solution
        self.number_input = len(input_speeds)
        self.Z_range_sun = Z_range_sun
        self.Z_range_ring = Z_range_ring
        self.number_planet = number_planet
        self.speed_planet_carrer = speed_planet_carrer
        self.torque_planet_carrer = torque_planet_carrer
        self.different_solution = different_solution
        DessiaObject.__init__(self, name=name)

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



    # def test_speed_precision(self,planetary_gear):
    #     for planetary in planetary_gear.planetaries:
    #         range_speed = planetary_gear.speed_range(planetary, planetary_gear.planet_carrier, self.precision)

    #         if not range_speed:
    #             return False
    #     list_planetary=copy.copy(planetary_gear.planetaries)
    #     for planetary in planetary_gear.planetaries:
    #         list_planetary.remove(planetary)
    #         for planetary_2 in list_planetary:
    #             range_speed = planetary_gear.speed_range(planetary, planetary_2, self.precision)

    #             if not range_speed :
    #                 return False

    #     return True



    def test_GCD(self, Z_1, Z_2):

        if m.gcd(Z_1, Z_2) != 1:

                     return False

        return True

    def test_reason_min_max(self, planetary_gear):
        list_planetaries = copy.copy(planetary_gear.planetaries)
        planetary_input = planetary_gear.planetaries[self.number_planetary_reason_input]
        list_planetaries.remove(planetary_input)

        for planetary in list_planetaries:
            path = planetary_gear.path_planetary_to_planetary([planetary_input, planetary])

            reason = planetary_gear.reason(path[0])
            for i, speed_input in enumerate(self.input_speeds):
                if planetary.speed_input == speed_input:
                    if reason < self.input_reason[i-1][0] or reason > self.input_reason[i-1][1]:
                        # print(reason)
                        # print(self.input_reason)

                        return False

        return True





    def test_vitesse_and_assembly_condition(self, planetary_gear, begin_meshing_chain, end_meshing_chain,
                                            list_previous_planetary, list_path):

        # list_previous_planetary_2 = copy.copy(list_previous_planetary)
        # list_previous_planetary_2.remove(begin_meshing_chain)


        # range_speed = planetary_gear.speed_range(begin_meshing_chain, planetary_gear.planet_carrier, list_previous_planetary, 1, list_path)

        # if range_speed == 'simplex':


        # #x=[speed_diff_1,speed_diff_2,speed_1,speed_2]
        #     c = [-1, 0, 0, 0]
        #     A = [[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, 1], [0, 1, 0, -1]]
        #     b = [begin_meshing_chain.speed_input[1], -begin_meshing_chain.speed_input[0],
        #          planetary_gear.planet_carrier.speed_input[1], -planetary_gear.planet_carrier.speed_input[0]]
        #     speed_diff_1_bound = (0, None)
        #     speed_diff_2_bound = (0, None)
        #     speed_1_bound = (None, None)
        #     speed_2_bound = (None, None)

        #     for i, planetary in enumerate(list_previous_planetary_2):
        #         speed_input_planetary = planetary.speed_input
        #         # reason=planetary_gear.reason(list_path[i][0])
        #         path = list_path[i]
        #         reason = planetary_gear.reason(path[0])

        #         if reason < 0:
        #             A.extend([[-reason, (1-reason), -reason, -(1-reason)], [-reason, 1-reason, reason, 1-reason]])
        #         else:

        #             if reason < 1:
        #                 A.extend([[reason, 1-reason, -reason, -(1-reason)], [reason, 1-reason, reason, 1-reason]])
        #             else:
        #                 A.extend([[reason, -(1-reason), -reason, -(1-reason)], [reason, -(1-reason), reason, (1-reason)]])

        #         b.extend([-speed_input_planetary[0], speed_input_planetary[1]])

        #     res = op.linprog(c, A_ub=A, b_ub=b, bounds=[speed_diff_1_bound, speed_diff_2_bound, speed_1_bound, speed_2_bound])

        #     if not res.success:
        #         return False


        # elif not range_speed:
        #     return False


        # elif range_speed[begin_meshing_chain][0] > range_speed[begin_meshing_chain][1]:
        #     return False

        # elif range_speed[planetary_gear.planet_carrier][0] > range_speed[planetary_gear.planet_carrier][1]:
        #     return False





        if not planetary_gear.test_assembly_condition(self.number_planet, [begin_meshing_chain, end_meshing_chain]):

            return False

        return True








    def z_range_mini_max(self, planetary_gear, element, begin_meshing_chain, end_meshing_chain, path, reasons_min_max):



        if not element in path[0]:

            return []
        if not element in reasons_min_max.keys():
            index_input = planetary_gear.planetaries.index(begin_meshing_chain)
            index_output = planetary_gear.planetaries.index(end_meshing_chain)
            if index_input == 0:
                reason_min_max = self.input_reason[index_output-1]
            elif index_output == 0:
                reason_min_max = [1/self.input_reason[index_input-1][1], 1/self.input_reason[index_input-1][0]]
            else:
                reason_min_max = [self.input_reason[index_output-1][0]/self.input_reason[index_input-1][1], self.input_reason[index_output-1][1]/self.input_reason[index_input-1][0]]

            reasons_min_max[element] = reason_min_max

        reason = planetary_gear.reason(path[0])
        reason_min = reasons_min_max[element][0]
        reason_max = reasons_min_max[element][1]

        if isinstance(element, Planetary):
            # if reason_min and reason_max:
                Z_min = reason*element.Z/reason_max

            # else:
            #     Z_min = 0

            # if reason_max:
                Z_max = reason*element.Z/reason_min

            # else:
            #     Z_max = self.Z_range_sun[1]
        else:
            # if reason_min and reason_max:
                Z_min = (element.Z*reason_min)/reason

            # else:
                # Z_min = 0

            # if reason_max:
                Z_max = (element.Z*reason_max)/reason

            # else:
                # Z_max = self.Z_range_sun[1]


        Z_range_mini_maxi = [min([Z_min, Z_max]), max([Z_max, Z_min])]
        # if Z_max>0 :
        #     print(Z_max)

        return Z_range_mini_maxi

    def decision_tree_speed_possibilities(self):
        tree = dt.DecisionTree()


        list_solution = []
        input_speeds_2 = copy.deepcopy(self.input_speeds)

        self.planetary_gear.planet_carrier.speed_input = self.speed_planet_carrer
        self.planetary_gear.planet_carrier.torque_input = self.torque_planet_carrer
        list_possibility_speed = []
        self.multiplication_possibility_speed(input_speeds_2, 0, [], list_possibility_speed)

        speed_reason_input = input_speeds_2[0]
        tree.SetCurrentNodeNumberPossibilities(len(list_possibility_speed))
        node = tree.NextNode(True)

        while not tree.finished:






            if len(node) == 1:

                    possibility_speed = list_possibility_speed[node[0]]
                    for i, planetary in enumerate(self.planetary_gear.planetaries):
                        planetary.speed_input = possibility_speed[i]

                        planetary.torque_input = self.input_torques[self.input_speeds.index(possibility_speed[i])]
                        if possibility_speed[i] == speed_reason_input:
                            self.number_planetary_reason_input = i
                    planetary_gear = copy.deepcopy(self.planetary_gear)
                    list_solution.append(planetary_gear)

                    tree.SetCurrentNodeNumberPossibilities(0)

            node = tree.NextNode(True)

        return list_solution

    def decision_tree(self) -> List[PlanetaryGear]:
        list_planetary_gears_speed = self.decision_tree_speed_possibilities()
        #print(list_planetary_gears_speed[0].planet_carrier.speed_input)
        list_solution = []
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
            meshing_chains_modif = planetary_gear.meshing_chain()
            meshing_chains = copy.copy(meshing_chains_modif)
            number_element_meshing_chain = []
            numbers_meshing_chain = []
            number_meshing_chain = 0
            totals_element_previous_meshing_chain = []
            total_element_previous_meshing_chain = 0
            flags_meshing_change = []

            flag_gcd = []
            numbers_planetaries_by_meshing_chain = []


            for i, meshing_chain in enumerate(meshing_chains_modif):
                if isinstance(meshing_chain[-1], Planetary):

                    if meshing_chain[-1].planetary_type == 'Ring':

                        meshing_chains_modif[i] = meshing_chain[::-1]
                        meshing_chains[i] = meshing_chain[::-1]

                    if  not isinstance(meshing_chain[0], Planetary) and meshing_chain[-1].planetary_type == 'Sun':

                        meshing_chains_modif[i] = meshing_chain[::-1]
                        meshing_chains[i] = meshing_chain[::-1]
                meshing_chain_2 = copy.copy(meshing_chain)
                number_planetaries = 0
                for element in meshing_chain_2:
                    if isinstance(element, Planetary):
                        number_planetaries += 1
                    if element in list_planet_remove:

                        meshing_chains_modif[i].remove(element)
                numbers_planetaries_by_meshing_chain.append(number_planetaries)


            print(numbers_planetaries_by_meshing_chain)
            if numbers_planetaries_by_meshing_chain[0] == 1:
                if 2 in numbers_planetaries_by_meshing_chain:
                    meshing_chain_1 = meshing_chains_modif[numbers_planetaries_by_meshing_chain.index(2)]
                    meshing_chains_modif[numbers_planetaries_by_meshing_chain.index(2)] = meshing_chains_modif[0]
                    meshing_chains_modif[0] = meshing_chain_1


            for meshing_chain in meshing_chains_modif:

                number_element_meshing_chain.append(len(meshing_chain))
                flags_meshing_change.append(1)

                for i, element in enumerate(meshing_chain):
                    flag_gcd.append(2)

                    if i != 0:
                        flags_meshing_change.append(0)

                    totals_element_previous_meshing_chain.append(total_element_previous_meshing_chain)
                    numbers_meshing_chain.append(number_meshing_chain)

                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':

                        list_tree.append(self.Z_range_ring[1]-self.Z_range_ring[0])
                        list_node_range_data.append(self.Z_range_ring[0])

                    else:
                        list_tree.append(self.Z_range_sun[1]-self.Z_range_sun[0])
                        list_node_range_data.append(self.Z_range_sun[0])

                number_meshing_chain += 1
                total_element_previous_meshing_chain += len(meshing_chain)

            list_planet_remove_neighbour = []


            for i, planet in enumerate(list_planet_remove):
                planet.Z = 1
                list_planet_remove_neighbour.append([planet])

                for meshing in planetary_gear.meshings:

                    if meshing.nodes[0] == planet:
                        list_planet_remove_neighbour[i].append(meshing.nodes[1])

                    if meshing.nodes[1] == planet:
                        list_planet_remove_neighbour[i].append(meshing.nodes[0])


            tree = dt.RegularDecisionTree(list_tree)

            Z_range_mini_maxi = []
            Z_range_mini_maxi_2 = []
            flag_meshing_change = 0
            flag_Z_range_mini_maxi = 0
            flag_Z_range_mini_maxi_2 = 0
            number_max_z_planet = [self.Z_range_sun[1]]*len(meshing_chains_modif)
            list_planetaries_Z_range_mini_maxi = []
            list_path = []
            reason_min_max = {}

            while not tree.finished and len(list_solution) < self.number_solution:

                valid = True
                node = tree.current_node

                number_meshing_chain = numbers_meshing_chain[len(node)-1]

                flag_meshing_change = flags_meshing_change[len(node)-1]
                total_element_previous_meshing_chain = totals_element_previous_meshing_chain[len(node)-1]

                element = meshing_chains_modif[number_meshing_chain][len(node)-total_element_previous_meshing_chain-1]
                element.Z = list_node_range_data[len(node)-1]+ node[len(node)-1]


                if len(node) == 1:


                    if isinstance(element, Planetary) and element.planetary_type == 'Ring':
                        number_max_z_planet[number_meshing_chain] = element.Z





                elif not flag_meshing_change:


                    previous_element = meshing_chains_modif[number_meshing_chain][len(node)-total_element_previous_meshing_chain-2]

                    if len(meshing_chains_modif[number_meshing_chain]) > 2 \
                    and meshing_chains_modif[number_meshing_chain][0].planetary_type == 'Ring':

                        number_tot_z_previous_planets = 0
                        number_max_z_previous_planets = 0
                        number_planet = 0
                        for i in range(len(node)-total_element_previous_meshing_chain-2):
                            previous_planet = meshing_chains_modif[number_meshing_chain][i+1]

                            number_planet += 1
                            number_tot_z_previous_planets += previous_planet.Z
                            if previous_planet.Z > number_max_z_previous_planets:
                                number_max_z_previous_planets = previous_planet.Z



                        if isinstance(element, Planetary):

                            flag_impose_z_number = True

                            for planet in list_planet_remove:
                              if planet in meshing_chains[number_meshing_chain]:
                                  flag_impose_z_number = False
                                  break

                            if flag_impose_z_number and number_planet == 1:
                                    if element.Z != (number_max_z_planet[number_meshing_chain]-number_max_z_previous_planets*2):

                                        valid = False

                            else:
                                if element.Z >= (number_max_z_planet[number_meshing_chain]-number_max_z_previous_planets*2) \
                                or element.Z < (number_max_z_planet[number_meshing_chain]-number_tot_z_previous_planets*2):

                                        valid = False
                        else:

                            if element.Z >= (number_max_z_planet[number_meshing_chain])/2:

                                valid = False
                    else:
                        if element.Z >= (number_max_z_planet[number_meshing_chain])/2:

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
                        number_max_z_planet[number_meshing_chain] = element.Z


                    else:
                        number_max_z_planet[number_meshing_chain] = (self.Z_range_sun[1])*2


                if len(node) == number_element_meshing_chain[number_meshing_chain]+total_element_previous_meshing_chain and valid:

                    begin_meshing_chain = meshing_chains_modif[number_meshing_chain][0]
                    end_meshing_chain = meshing_chains_modif[number_meshing_chain][-1]

                    #planetary_gear = PlanetaryGear(planetary_gear.planetaries, planetary_gear.planets,
                                                   #planetary_gear.planet_carrier, planetary_gear.connections, planetary_gear.name)

                    if Z_range_mini_maxi:

                        if Z_range_mini_maxi[1] < 0:
                            return []
                        if element.Z < Z_range_mini_maxi[0] or element.Z > Z_range_mini_maxi[1]:

                            valid = False

                    if Z_range_mini_maxi_2:

                        if Z_range_mini_maxi_2[1] < 0:
                            return []

                        if element.Z < Z_range_mini_maxi_2[0]  or element.Z > Z_range_mini_maxi_2[1]:
                            valid = False


                    if valid:

                        if numbers_meshing_chain[len(node)-1] == 0:

                                first_planetary = begin_meshing_chain

                                if not first_planetary in list_planetaries_Z_range_mini_maxi:
                                    list_planetaries_Z_range_mini_maxi.append(first_planetary)


                                if isinstance(begin_meshing_chain, Planetary) and isinstance(end_meshing_chain, Planetary):

                                    if not end_meshing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([begin_meshing_chain, end_meshing_chain]))

                                    if not flag_Z_range_mini_maxi:


                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element, begin_meshing_chain,
                                                                                  end_meshing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)-1],
                                                                                  reason_min_max)

                                        flag_Z_range_mini_maxi = 1



                                        if element.Z < Z_range_mini_maxi[0] or element.Z < Z_range_mini_maxi[1]:
                                            list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                            list_path_speed_test = []
                                            for previous_planetary in list_previous_planetary[1:]:
                                                list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                            valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                             end_meshing_chain,
                                                                                             list_previous_planetary,
                                                                                             list_path_speed_test)

                                        else:
                                            valid = False
                                    else:
                                        list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                        list_path_speed_test = []
                                        for previous_planetary in list_previous_planetary[1:]:
                                            list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                        valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                         end_meshing_chain,
                                                                                         list_previous_planetary,
                                                                                         list_path_speed_test)



                        else:

                            if not begin_meshing_chain in list_planetaries_Z_range_mini_maxi:


                                        list_planetaries_Z_range_mini_maxi.append(begin_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, begin_meshing_chain]))

                            if not flag_Z_range_mini_maxi:

                                        Z_range_mini_maxi = self.z_range_mini_max(planetary_gear, element,
                                                                                  first_planetary, begin_meshing_chain,
                                                                                  list_path[list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)-1],
                                                                                  reason_min_max)

                                        flag_Z_range_mini_maxi = 1

                                        if Z_range_mini_maxi:
                                            if element.Z > Z_range_mini_maxi[0] and element.Z < Z_range_mini_maxi[1]:
                                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)+1]
                                                list_path_speed_test = []
                                                for previous_planetary in list_previous_planetary[1:]:
                                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                                 begin_meshing_chain,
                                                                                                 list_previous_planetary,
                                                                                                 list_path_speed_test)

                                            else:
                                                valid = False

                            else:
                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(begin_meshing_chain)+1]
                                list_path_speed_test = []
                                for previous_planetary in list_previous_planetary[1:]:
                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])

                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                 begin_meshing_chain,
                                                                                 list_previous_planetary,
                                                                                 list_path_speed_test)




                            if isinstance(end_meshing_chain, Planetary) and valid:


                                if not end_meshing_chain in list_planetaries_Z_range_mini_maxi:

                                        list_planetaries_Z_range_mini_maxi.append(end_meshing_chain)
                                        list_path.append(planetary_gear.path_planetary_to_planetary([first_planetary, end_meshing_chain]))


                                if not flag_Z_range_mini_maxi_2:

                                        Z_range_mini_maxi_2 = self.z_range_mini_max(planetary_gear, element, first_planetary, end_meshing_chain,
                                                                                    list_path[list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)-1], reason_min_max)
                                        flag_Z_range_mini_maxi_2 = 1
                                        if Z_range_mini_maxi_2:

                                            if element.Z > Z_range_mini_maxi_2[0] and element.Z < Z_range_mini_maxi_2[1]:

                                                list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)]
                                                list_path_speed_test = []

                                                for previous_planetary in list_previous_planetary[1:]:
                                                    list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                                valid = self.test_vitesse_and_assembly_condition(planetary_gear, begin_meshing_chain,
                                                                                                 end_meshing_chain,
                                                                                                 list_previous_planetary,
                                                                                                 list_path_speed_test)

                                            else:
                                                valid = False
                                else:
                                    list_previous_planetary = list_planetaries_Z_range_mini_maxi[0:list_planetaries_Z_range_mini_maxi.index(end_meshing_chain)+1]
                                    list_path_speed_test = []
                                    for previous_planetary in list_previous_planetary[1:]:
                                        list_path_speed_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])
                                    valid = self.test_vitesse_and_assembly_condition(planetary_gear, first_planetary,
                                                                                     end_meshing_chain,
                                                                                     list_previous_planetary,
                                                                                     list_path_speed_test)






                        if number_meshing_chain == len(meshing_chains_modif)-1 and valid:
                            list_previous_planetary = list_planetaries_Z_range_mini_maxi
                            list_path_torque_test = []

                            for previous_planetary in list_previous_planetary[1:]:
                                        list_path_torque_test.append(list_path[list_planetaries_Z_range_mini_maxi.index(previous_planetary)-1])



                            # if not self.test_reason_min_max(planetary_gear):

                            #     valid = False

                            list_tree_planetary = []





                            if list_planet_remove and valid:

                                for i in range(len(list_planet_remove)):
                                    list_planet_remove_2 = copy.copy(list_planet_remove)
                                    list_planet_remove_neighbour_2 = copy.copy(list_planet_remove_neighbour)
                                    list_planet_remove[i].Z = 1
                                    valid_planet = True

                                    for meshing_chain in meshing_chains:
                                        if list_planet_remove[i] in meshing_chain and valid_planet:

                                            if meshing_chain[0].planetary_type == 'Ring':
                                                number_planet = 0
                                                for element_2 in meshing_chain:
                                                    if  isinstance(element_2, Planet):
                                                        number_planet += 1

                                                if isinstance(meshing_chain[-1], Planetary) and number_planet == 1:
                                                    number_z_max = (meshing_chain[0].Z - meshing_chain[-1].Z)/2

                                                    if number_z_max < self.Z_range_sun[0]:
                                                        list_planet_remove[i].Z = self.Z_range_sun[0]
                                                        valid_planet = False
                                                    neighbour = list_planet_remove_neighbour[i]
                                                    if not self.test_GCD(number_z_max, neighbour[1].Z):
                                                        valid_planet = False

                                                        break
                                                    if not self.test_GCD(number_z_max, neighbour[2].Z):
                                                        valid_planet = False
                                                        break
                                                    list_planet_remove[i].Z = number_z_max
                                                    list_planet_remove_2.remove(list_planet_remove[i])
                                                    list_planet_remove_neighbour_2.remove(neighbour)

                                                else:
                                                    if isinstance(meshing_chain[-1], Planetary):
                                                        number_z_max = int((meshing_chain[0].Z- meshing_chain[-1].Z)/2)+1
                                                    else:
                                                        number_z_max = meshing_chain[0].Z/2


                                                    if number_z_max - self.Z_range_sun[0] > 0:
                                                        list_tree_planetary.append(number_z_max - self.Z_range_sun[0])

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
                                            planetary_gear.speed_max_planet = planetary_gear.speed_max_planets()
                                            if planetary_gear.speed_max_planet < self.speed_max_planet:
                                                 if self.different_solution and list_solution:
                                                     valid_different = True
                                                     precedent_planetary_gear = list_solution[-1]
                                                     for precedent_planetary, planetary in zip(precedent_planetary_gear.planetaries, planetary_gear.planetaries):
                                                         if precedent_planetary.Z == planetary.Z:
                                                             valid_different = False
                                                     for precedent_planet, planet in zip(precedent_planetary_gear.planets, planetary_gear.planets):
                                                        if precedent_planet.Z == planet.Z:
                                                             valid_different = False
                                                     if valid_different:
                                                         print(planetary_gear)
                                                         list_solution.append(copy.deepcopy(planetary_gear))
                                                 else:
                                                    print(planetary_gear)
                                                    list_solution.append(copy.deepcopy(planetary_gear))

                                            # print(planetary_gear)
                                            # if list_solution:
                                            #     return list_solution

                                        tree_planet.NextNode(valid_planet)

                            else:
                                planetary_gear.speed_max_planet = planetary_gear.speed_max_planets()
                                if valid and planetary_gear.speed_max_planet < self.speed_max_planet:

                                    if self.different_solution and list_solution:
                                        valid_different = True
                                        precedent_planetary_gear = list_solution[-1]
                                        for precedent_planetary, planetary in zip(precedent_planetary_gear.planetaries, planetary_gear.planetaries):
                                            if precedent_planetary.Z == planetary.Z:
                                                valid_different = False
                                        for precedent_planet, planet in zip(precedent_planetary_gear.planets, planetary_gear.planets):
                                           if precedent_planet.Z == planet.Z:
                                                valid_different = False
                                        if valid_different:
                                            print(planetary_gear)
                                            list_solution.append(copy.deepcopy(planetary_gear))
                                    else:
                                       print(planetary_gear)
                                       list_solution.append(copy.deepcopy(planetary_gear))

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


        # a=0
        # for i in range(len(list_solution)):
        #     if not self.test_speed_precision(list_solution[i+a]):
        #         list_solution.pop(i+a)
        #         a-=1
        return list_solution



class GeneratorPlanetaryGearsGeometry(DessiaObject):
    _standalone_in_db = True

    _eq_is_data_eq = False
    def __init__(self, planetary_gear: PlanetaryGear,
                 number_planet: int, D_min: float,
                 D_max: float,
                 recircle_power_max: float,
                 internal_diameter_min: float = 0,
                 recirculation: bool = False,
                 name: str = '',):
        self.recirculation = recirculation
        self.planetary_gear = planetary_gear
        self.number_planet = number_planet
        self.D_min = D_min
        self.D_max = D_max
        self.position_min_max = PositionMinMaxPlanetaryGear(self.planetary_gear)
        self.d_min = 0
        self.recircle_power_max = recircle_power_max
        self.internal_diameter_min = internal_diameter_min
        self.x0_save = []
        DessiaObject.__init__(self, name=name)



    def function_minimize_equation(self, meshing_chain, X, Y, M):
        planet_double = []
        teta = 2*m.pi/self.number_planet
        f = []
        f2 = []
        index_x = 0
        limit_x = []
        for double in self.planetary_gear.doubles:
            planet_double.append(double.nodes)



        for i, element in enumerate(meshing_chain):
            if isinstance(element, Planet):
                index_planet = self.planetary_gear.planets.index(element)

                if i == 1 and isinstance(meshing_chain[0], Planetary):
                    if meshing_chain[0].planetary_type == 'Ring':
                        f.append(X[index_planet]**2+Y[index_planet]**2-(M*(meshing_chain[0].Z-element.Z)/2)**2)

                    else:
                        f.append(X[index_planet]**2+Y[index_planet]**2-(M*(meshing_chain[0].Z+element.Z)/2)**2)


                    if len(meshing_chain) > 2 and not isinstance(meshing_chain[2], Planetary):

                        index_other_planet = self.planetary_gear.planets.index(meshing_chain[2])
                        f.append((X[index_planet]-X[index_other_planet])**2+(Y[index_planet]-Y[index_other_planet])**2-(M*(meshing_chain[2].Z+element.Z)/2)**2)


                if i == len(meshing_chain)-2 and isinstance(meshing_chain[-1], Planetary):
                    if meshing_chain[-1].planetary_type == 'Ring':
                        f.append(X[index_planet]**2+Y[index_planet]**2-(M*(meshing_chain[-1].Z-element.Z)/2)**2)

                    else:
                        f.append(X[index_planet]**2+Y[index_planet]**2-(M*(meshing_chain[-1].Z+element.Z)/2)**2)

                elif i == 0 or i == len(meshing_chain)-1:
                    if i == 0:
                        index_other_planet = self.planetary_gear.planets.index(meshing_chain[1])
                        f.append((X[index_planet]-X[index_other_planet])**2+(Y[index_planet]-Y[index_other_planet])**2-(M*(meshing_chain[i+1].Z+element.Z)/2)**2)

                else:

                    index_other_planet = self.planetary_gear.planets.index(meshing_chain[i+1])
                    f.append((X[index_planet]-X[index_other_planet])**2+(Y[index_planet]-Y[index_other_planet])**2-(M*(meshing_chain[i+1].Z+element.Z)/2)**2)



                for double in planet_double:
                    if element in double:

                        if element == double[0]:
                            other_double_index = self.planetary_gear.planets.index(double[1])
                        if element == double[1]:
                            other_double_index = self.planetary_gear.planets.index(double[0])

                        f.append(X[other_double_index]-X[index_planet])

                        f.append(Y[other_double_index]-Y[index_planet])

                        planet_double.remove(double)



        # for i in range(len(f)):
        #     f[i]=20*f[i]
        return f

    def function_minimize_inequation_meshing_chain(self, x, meshing_chain, X, Y):
        teta = 2*m.pi/self.number_planet
        X_prime = copy.deepcopy(X)
        Y_prime = copy.deepcopy(Y)
        index_x = 0
        f2 = []

        for i, element in enumerate(meshing_chain):
            if isinstance(element, Planet):
                index_planet = self.planetary_gear.planets.index(element)
                X_prime[index_planet] = X[index_planet]*m.cos(teta)-Y[index_planet]*m.sin(teta)
                Y_prime[index_planet] = X[index_planet]*m.sin(teta)+Y[index_planet]*m.cos(teta)

        for i, element in enumerate(meshing_chain):
            if isinstance(element, Planet):
                index_planet = self.planetary_gear.planets.index(element)
                if i > 1:
                    for element_2 in meshing_chain[:i-1]:

                        if isinstance(element_2, Planetary):
                            f2.append((Y[index_planet])**2+(X[index_planet])**2-x[index_x])
                            index_x += 1

                        else:
                            index_other_planet = self.planetary_gear.planets.index(element_2)
                            f2.append((Y[index_planet]-Y[index_other_planet])**2+(X[index_planet]-X[index_other_planet])**2-x[index_x])
                            index_x += 1



                for element_3 in meshing_chain:
                    if isinstance(element_3, Planet):
                        index_other_planet = self.planetary_gear.planets.index(element_3)
                        f2.append((Y_prime[index_other_planet]-Y[index_planet])**2+(X_prime[index_other_planet]-X[index_planet])**2-x[index_x])
                        # print(X_prime[index_other_planet])
                        # print(Y_prime[index_other_planet])
                        # print(f2)
                        index_x += 1

            elif i == len(meshing_chain)-1:
                    for element_2 in meshing_chain[:i-1]:
                       if isinstance(element_2, Planet):

                            index_other_planet = self.planetary_gear.planets.index(element_2)
                            f2.append((Y[index_other_planet])**2+(X[index_other_planet])**2-x[index_x])

                            index_x += 1

        F = 0

        for f in f2:
            F += f**2

        return F

    def function_minimize_inequation_meshing_chain_min_max(self, meshing_chain, M):
        min_x_max_x = []
        index_x = 0
        x0 = []
        outside_diameter_coefficient = 4*(0.01*M+M*1.2)
        for i, element in enumerate(meshing_chain):
            if isinstance(element, Planet):
                index_planet = self.planetary_gear.planets.index(element)
                if i > 1:
                    for element_2 in meshing_chain[:i-1]:
                        if isinstance(element_2, Planetary):
                            if element_2.planetary_type == 'Ring':
                                min_x_max_x.append([-np.inf, ((((element_2.Z-element.Z))*M-outside_diameter_coefficient)/2)**2])
                                x0.append(((((element_2.Z-element.Z))*M-outside_diameter_coefficient)/2)**2)
                            else:
                                min_x_max_x.append([(((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2, np.inf])
                                x0.append((((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2)

                        else:
                            min_x_max_x.append([(((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2, np.inf])
                            x0.append((((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2)

                for element_3 in meshing_chain:
                    if isinstance(element_3, Planet):
                        min_x_max_x.append([(((element_3.Z+element.Z)*M+outside_diameter_coefficient)/2)**2, np.inf])
                        x0.append((((element_3.Z+element.Z)*M+outside_diameter_coefficient)/2)**2)

            elif i == len(meshing_chain)-1:
                    for element_2 in meshing_chain[:i-1]:
                       if isinstance(element_2, Planet):
                            if element.planetary_type == 'Ring':
                                min_x_max_x.append([-np.inf, (((element.Z-element_2.Z)*M-outside_diameter_coefficient)/2)**2])
                                x0.append((((element.Z-element_2.Z)*M-outside_diameter_coefficient)/2)**2)
                            else:
                                min_x_max_x.append([(((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2, np.inf])
                                x0.append((((element_2.Z+element.Z)*M+outside_diameter_coefficient)/2)**2)



        return min_x_max_x, x0


    def function_inequation_constrain(self, x, planetary_gear):
        meshing_chains = planetary_gear.meshing_chain()
        X = [0]*len(self.planetary_gear.planets)
        Y = [0]*len(self.planetary_gear.planets)
        M = [0]*len(meshing_chains)
        index_x = 0
        f2 = []
        # print(x)
        for i, meshing_chain in enumerate(meshing_chains):
            M[i] = x[index_x]
            index_x += 1

        X[0] = 0
        Y[0] = x[index_x]
        index_x += 1
        for i, planets in enumerate(planetary_gear.planets[1:]):
            X[i+1] = x[index_x]
            index_x += 1
            Y[i+1] = x[index_x]
            index_x += 1


        for i, meshing_chain in enumerate(meshing_chains):
            m2 = M[i]
            if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                     f2.append(self.D_max-m2*meshing_chain[0].Z)
                     f2.append(-self.D_min+m2*meshing_chain[0].Z)




            elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                 f2.append(self.D_max-m2*meshing_chain[-1].Z)
                 f2.append(-self.D_min+m2*meshing_chain[-1].Z)


            else:

                for element in meshing_chain:
                    if isinstance(element, Planet):
                        index_planet = planetary_gear.planets.index(element)
                        f2.append((self.D_max/2)**2-(X[index_planet]**2+Y[index_planet]**2))
                        f2.append(-(self.D_min/2)**2+(X[index_planet]**2+Y[index_planet]**2))
        if self.internal_diameter_min:
            for i, meshing_chain in enumerate(meshing_chains):
                m2 = M[i]
                if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Sun':
                    f2.append(-self.internal_diameter_min+m2*meshing_chain[0].Z)


                elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Sun':
                    f2.append(-self.internal_diameter_min+m2*meshing_chain[-1].Z)

                else:

                    for element in meshing_chain:
                        if isinstance(element, Planet):
                            index_planet = planetary_gear.planets.index(element)

                            f2.append(-(self.internal_diameter_min/2)+((X[index_planet]**2+Y[index_planet]**2)**1/2)-element.Z*m2/2)
        F = 0
        # print(f2)
        for f in f2:
            F += f-abs(f)
        # if F==0:
        #     F=1
        print(F)
        return F


    def verification(self):

        meshing_chains = self.planetary_gear.meshing_chain()


        def function_verification(x, planetary_gear, meshing_chains):

            X = [0]*len(self.planetary_gear.planets)
            Y = [0]*len(self.planetary_gear.planets)
            M = [0]*len(meshing_chains)
            index_x = 0
            f2 = []

            for i, meshing_chain in enumerate(meshing_chains):
                M[i] = x[index_x]
                index_x += 1

            X[0] = 0
            Y[0] = x[index_x]
            index_x += 1

            for i, planets in enumerate(planetary_gear.planets[1:]):
                X[i+1] = x[index_x]
                index_x += 1
                Y[i+1] = x[index_x]
                index_x += 1

            for i, meshing_chain in enumerate(meshing_chains):

                m2 = M[i]
                function_minimize_equation = self.function_minimize_equation(meshing_chain, X, Y, m2)
                for f in function_minimize_equation:
                    f2.append(f*10)

                f2.extend(self.function_minimize_equation(meshing_chain, X, Y, m2))

                min_max_x, x0 = self.function_minimize_inequation_meshing_chain_min_max(meshing_chain, m2)

                res_1 = op.minimize(self.function_minimize_inequation_meshing_chain, x0, bounds=min_max_x, args=(meshing_chain, X, Y))

                if res_1.fun > 0.0000000000001:
                      f2.append((res_1.fun**(1/2))*1000)

                if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                              if m2*meshing_chain[0].Z*1.3 > self.D_max:
                                  f2.append((self.D_max-m2*meshing_chain[0].Z*1.3))
                              elif m2*meshing_chain[0].Z*1.3 < self.D_min:
                                  f2.append((self.D_min-m2*meshing_chain[0].Z*1.3))







                elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                      if m2*meshing_chain[-1].Z*1.3 > self.D_max:
                          f2.append((self.D_max-m2*meshing_chain[-1].Z*1.3))
                      elif m2*meshing_chain[-1].Z*1.3 < self.D_min:
                          f2.append((self.D_min-m2*meshing_chain[-1].Z*1.3))






                else:
                    # print(1)
                    for element in meshing_chain:
                        if isinstance(element, Planet):

                            index_planet = planetary_gear.planets.index(element)
                            d = (X[index_planet]**2+Y[index_planet]**2)**(1/2)+element.Z*m2/2
                            if d > self.D_max/2:
                                f2.append(((self.D_max/2)-(X[index_planet]**2+Y[index_planet]**2)**(1/2)-element.Z*m2/2))








            F = 0

            for f in f2:
                F += (abs(f))
            # print(f2)
            # print(x)
            # print(F)
            return F
        min_max_x_2 = []
        min_max_x_3 = []
        x0 = []
        for meshing_chain in meshing_chains:
            min_max_x_2.append([0, self.D_max])
            min_max_x_3.append([0, self.D_max])
            x0.append(1)

        min_max_x_2.append([-self.D_max/2, self.D_max/2])
        min_max_x_3.append([-self.D_min/2, self.D_min/2])
        x0.append(1)

        for planets in self.planetary_gear.planets[1:]:
            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            min_max_x_3.append([-self.D_min/2, self.D_min/2])
            min_max_x_3.append([-self.D_min/2, self.D_min/2])
            x0.append(1)
            x0.append(1)


        # for i, meshing_chain in enumerate(meshing_chains):
        #         if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
        #               min_max_x_2.append([self.D_min, self.D_max])
        #               min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
        #               x0.append(self.D_min)

        #         elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
        #               min_max_x_2.append([self.D_min, self.D_max])
        #               min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
        #               x0.append(self.D_min)

        #         else:

        #             for element in meshing_chain:
        #                 if isinstance(element, Planet):
        #                     min_max_x_2.append([0, self.D_max])
        #                     min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
        #                     x0.append(self.D_min)
        j = 0
        fx = 0.02
        while (j != 3 and fx > 0.0001):

            x0_2 = pyDOE.lhs(len(x0), 1000)
            for x0_1 in x0_2:

                for i, x in enumerate(x0):
                    x0_1 = random.random()
                    x0[i] = (min_max_x_3[i][1]-min_max_x_3[i][0])*x0_1 + min_max_x_3[i][0]

                # print(function_verification(x0,self.planetary_gear)-3000)
                if function_verification(x0, self.planetary_gear, meshing_chains) < 3:
                    break
            min_x = []
            max_x = []
            for i in range(len(min_max_x_2)):
                min_x.append(min_max_x_2[i][0])
                max_x.append(min_max_x_2[i][1])
            s = 0.1
            if (self.D_max-self.D_min) < 1:
                s = (self.D_max-self.D_min)*0.1

            # constraint={'type':'ineq','fun':self.function_inequation_constrain,'args':(self.planetary_gear,),}
            # res_2=op.minimize(function_verification,x0,bounds=min_max_x_2, method='SLSQP' , args=(self.planetary_gear),constraints=constraint,tol=0.0001, options={'ftol':1e-10,'maxiter':150})
            xra, fx = cma.fmin(function_verification, x0, s, args=(self.planetary_gear, meshing_chains), options={'bounds':[min_x, max_x],
                                                                                                                  'tolfun': 1e-12,
                                                                                                                  'verbose': 3,
                                                                                                                  'ftarget': 1e-24,
                                                                                                                  'maxiter': 1500})[0:2]
            j += 1

            # print(xra,fx)
        # if res_2.success==True and res_2.fun<0.001:
        #     break
        # if fx<0.001:
        #     break

        if fx < 0.0001:




            x = xra
            X = []
            Y = []
            M = []
            index_x = 0
            teta = 2*m.pi/self.number_planet
            for meshing_chain in meshing_chains:
                M.append(x[index_x])
                index_x += 1
            X.append(0)
            Y.append(x[index_x])
            index_x += 1

            for planets in self.planetary_gear.planets[1:]:
                X.append(x[index_x])
                index_x += 1
                Y.append(x[index_x])
                index_x += 1

            X_prime = []
            Y_prime = []
            for y in range(self.number_planet):
                X_prime.append([])
                Y_prime.append([])
                for i, x in enumerate(X):
                        X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                        Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))



            # list_color = ['mediumblue', 'purple', 'green', 'k']
            # plt.figure()
            # ax=plt.subplot(aspect='equal')

            z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

            for r, meshing_chain in enumerate(meshing_chains):
                z2 = z[r]
                m2 = M[r]
                # color=list_color[i]

                for i, planet in enumerate(self.planetary_gear.planets):
                    if planet in meshing_chain:
                        # ax.add_patch(plt.Circle([X[i],Y[i]], m2*planet.Z/2,color=color,fill=False))
                        planet.positions = [(z2, Y[i], X[i])]
                        planet.module = m2
                        d = (((X[i]**2+Y[i]**2)**0.5)*2+m2*planet.Z)
                        if d > self.d_min:
                            self.d_min = d

                        for y in range(self.number_planet):
                            # ax.add_patch(plt.Circle([X_prime[y][i],Y_prime[y][i]], m2*planet.Z/2,color=color,fill=False))
                            planet.positions.append((z2, Y_prime[y][i], X_prime[y][i]))


                for planetary in self.planetary_gear.planetaries:
                    if planetary in meshing_chain:
                       # ax.add_patch(plt.Circle([0,0], m2*planetary.Z/2,color=color,fill=False))
                       planetary.module = m2
                       planetary.position = (z2, 0, 0)
                       d = m2*planetary.Z
                       if d > self.d_min:
                           self.d_min = d

            # ax.relim()
            # ax.autoscale_view()
            # plt.show()
            self.planetary_gear.position = True
            self.x0_save = xra
            return self.planetary_gear


        else:
            self.planetary_gear.sum_Z_planetary = -101
            return self.planetary_gear


    def verification_recirculation(self,):

        meshing_chains = self.planetary_gear.meshing_chain()


        def function_verification(x, planetary_gear, meshing_chains):

            X = [0]*len(self.planetary_gear.planets)
            Y = [0]*len(self.planetary_gear.planets)
            M = [0]*len(meshing_chains)
            index_x = 0
            f2 = []

            for i, meshing_chain in enumerate(meshing_chains):
                M[i] = x[index_x]
                index_x += 1

            X[0] = 0
            Y[0] = x[index_x]
            index_x += 1

            for i, planets in enumerate(planetary_gear.planets[1:]):
                X[i+1] = x[index_x]
                index_x += 1
                Y[i+1] = x[index_x]
                index_x += 1

            for i, meshing_chain in enumerate(meshing_chains):

                m2 = M[i]
                f2.extend(self.function_minimize_equation(meshing_chain, X, Y, m2))

                min_max_x, x0 = self.function_minimize_inequation_meshing_chain_min_max(meshing_chain, m2)

                res_1 = op.minimize(self.function_minimize_inequation_meshing_chain, x0, bounds=min_max_x, args=(meshing_chain, X, Y))

                if res_1.fun > 0.000000000001:
                      f2.append(res_1.fun*100)

                if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                      f2.append(x[index_x]-m2*meshing_chain[0].Z*1.3)
                      index_x += 1

                elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                      f2.append(x[index_x]-m2*meshing_chain[-1].Z*1.3)
                      index_x += 1

                else:

                    for element in meshing_chain:
                        if isinstance(element, Planet):
                            index_planet = planetary_gear.planets.index(element)
                            f2.append(((x[index_x]/2)**2-(X[index_planet]**2+Y[index_planet]**2))/10000000)
                            index_x += 1






            F = 0
            # print(f2)
            for f in f2:
                F += ((abs(f)/self.D_min)**2)
            # print(f2)
            # print(x)
            # print(F)
            return F
        min_max_x_2 = []
        min_max_x_3 = []
        x0 = []
        for meshing_chain in meshing_chains:
            min_max_x_2.append([0, self.D_max/7])
            min_max_x_3.append([0, self.D_min/7])
            x0.append(1)

        min_max_x_2.append([-self.D_max/2, self.D_max/2])
        min_max_x_3.append([-self.D_min/2, self.D_min/2])
        x0.append(1)

        for planets in self.planetary_gear.planets[1:]:
            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            min_max_x_3.append([-self.D_min/2, self.D_min/2])
            min_max_x_3.append([-self.D_min/2, self.D_min/2])
            x0.append(1)
            x0.append(1)


        for i, meshing_chain in enumerate(meshing_chains):
                if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                      min_max_x_2.append([self.D_min, self.D_max])
                      min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
                      x0.append(self.D_min)

                elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                      min_max_x_2.append([self.D_min, self.D_max])
                      min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
                      x0.append(self.D_min)

                else:

                    for element in meshing_chain:
                        if isinstance(element, Planet):
                            min_max_x_2.append([self.D_min, self.D_max])
                            min_max_x_3.append([self.D_min, self.D_min+0.0001*self.D_min])
                            x0.append(self.D_min)

        min_x = []
        max_x = []
        for i in range(len(min_max_x_2)):
            min_x.append(min_max_x_2[i][0])
            max_x.append(min_max_x_2[i][1])
        j = 0
        fx = 0.02
        list_planetary_gears = []
        sum_power_recirculation = []
        while j != 50:

            x0_2 = pyDOE.lhs(len(x0), 100)
            for x0_1 in x0_2:

                for i, x in enumerate(x0):
                    # x0_1 = random.random()

                    x0[i] = (min_max_x_3[i][1]-min_max_x_3[i][0])*x0_1[i] + min_max_x_3[i][0]

                x = x0
                X = []
                Y = []
                M = []
                index_x = 0
                teta = 2*m.pi/self.number_planet
                for meshing_chain in meshing_chains:
                    M.append(x[index_x])
                    index_x += 1
                X.append(0)
                Y.append(x[index_x])
                index_x += 1

                for planets in self.planetary_gear.planets[1:]:
                    X.append(x[index_x])
                    index_x += 1
                    Y.append(x[index_x])
                    index_x += 1

                X_prime = []
                Y_prime = []
                for y in range(self.number_planet):
                    X_prime.append([])
                    Y_prime.append([])
                    for i, x in enumerate(X):
                            X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                            Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




                z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

                for r, meshing_chain in enumerate(meshing_chains):
                    z2 = z[r]
                    m2 = M[r]


                    for i, planet in enumerate(self.planetary_gear.planets):
                        if planet in meshing_chain:

                            planet.positions = [(z2, Y[i], X[i])]
                            planet.module = m2


                            for y in range(self.number_planet):

                                planet.positions.append((z2, Y_prime[y][i], X_prime[y][i]))


                    for planetary in self.planetary_gear.planetaries:
                        if planetary in meshing_chain:

                           planetary.module = m2
                           planetary.position = (z2, 0, 0)

                self.planetary_gear.update_position_mech()


                # print(function_verification(x0,self.planetary_gear)-3000)
                if function_verification(x0, self.planetary_gear, meshing_chains) < 1:
                    power_recirculation = self.planetary_gear.recirculation_power()
                    number_good_loop = 0
                    for power in power_recirculation:
                        if power[1] < self.recircle_power_max:
                            number_good_loop += 1
                    print(power_recirculation)
                    if number_good_loop == len(power_recirculation):
                        s = 0.1
                        if (self.D_max-self.D_min) < 1:
                            s = (self.D_max-self.D_min)*0.1



                        # constraint={'type':'ineq','fun':self.function_inequation_constrain,'args':(self.planetary_gear,),}
                        # res_2=op.minimize(function_verification,x0,bounds=min_max_x_2, method='SLSQP' , args=(self.planetary_gear),constraints=constraint,tol=0.0001, options={'ftol':1e-10,'maxiter':150})
                        xra, fx = cma.fmin(function_verification, x0, s, args=(self.planetary_gear, meshing_chains), options={'bounds':[min_x, max_x],
                                                                                                                              'tolfun': 1e-10,
                                                                                                                              'verbose': 3,
                                                                                                                              'ftarget': 1e-2,
                                                                                                                              'maxiter': 2000})[0:2]
                        x = xra
                        X = []
                        Y = []
                        M = []
                        index_x = 0
                        teta = 2*m.pi/self.number_planet
                        for meshing_chain in meshing_chains:
                            M.append(x[index_x])
                            index_x += 1
                        X.append(0)
                        Y.append(x[index_x])
                        index_x += 1

                        for planets in self.planetary_gear.planets[1:]:
                            X.append(x[index_x])
                            index_x += 1
                            Y.append(x[index_x])
                            index_x += 1

                        X_prime = []
                        Y_prime = []
                        for y in range(self.number_planet):
                            X_prime.append([])
                            Y_prime.append([])
                            for i, x in enumerate(X):
                                    X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                                    Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




                        z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

                        for r, meshing_chain in enumerate(meshing_chains):
                            z2 = z[r]
                            m2 = M[r]


                            for i, planet in enumerate(self.planetary_gear.planets):
                                if planet in meshing_chain:

                                    planet.positions = [(z2, Y[i], X[i])]
                                    planet.module = m2


                                    for y in range(self.number_planet):

                                        planet.positions.append((z2, Y_prime[y][i], X_prime[y][i]))


                            for planetary in self.planetary_gear.planetaries:
                                if planetary in meshing_chain:

                                   planetary.module = m2
                                   planetary.position = (z2, 0, 0)

                        self.planetary_gear.update_position_mech()
                        power_recirculation = self.planetary_gear.recirculation_power()
                        print(power_recirculation)
                        number_good_loop = 0
                        sum_power = 0
                        for power in power_recirculation:
                            if power[1] > sum_power:
                                sum_power = power[1]
                            if power[1] < self.recircle_power_max:
                                number_good_loop += 1
                        if number_good_loop == len(power_recirculation):
                            list_planetary_gears.append(xra)

                            sum_power_recirculation.append(sum_power)
                            if len(list_planetary_gears) > 9:

                                j = 49
                                break

            j += 1

        # sum_power_recirculation=[]
        # list_planetary_gears_2=[]
        # for xra_2 in list_planetary_gears:

        #     xra, fx = cma.fmin(function_verification, xra_2, 0.1, args=(self.planetary_gear, meshing_chains), options={'bounds':[min_x, max_x],
        #                                                                                                                         'tolfun': 1e-10,
        #                                                                                                                         'verbose': 3,
        #                                                                                                                         'ftarget': 1e-10,
        #                                                                                                                         'maxiter': 2000})[0:2]
        #     x = xra
        #     X = []
        #     Y = []
        #     M = []
        #     index_x = 0
        #     teta = 2*m.pi/self.number_planet
        #     for meshing_chain in meshing_chains:
        #         M.append(x[index_x])
        #         index_x += 1
        #     X.append(0)
        #     Y.append(x[index_x])
        #     index_x += 1

        #     for planets in self.planetary_gear.planets[1:]:
        #         X.append(x[index_x])
        #         index_x += 1
        #         Y.append(x[index_x])
        #         index_x += 1

        #     X_prime = []
        #     Y_prime = []
        #     for y in range(self.number_planet):
        #         X_prime.append([])
        #         Y_prime.append([])
        #         for i, x in enumerate(X):
        #                 X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
        #                 Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




        #     z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

        #     for r, meshing_chain in enumerate(meshing_chains):
        #         z2 = z[r]
        #         m2 = M[r]


        #         for i, planet in enumerate(self.planetary_gear.planets):
        #             if planet in meshing_chain:

        #                 planet.positions = [(z2, Y[i], X[i])]
        #                 planet.module = m2


        #                 for y in range(self.number_planet):

        #                     planet.positions.append((z2, Y_prime[y][i], X_prime[y][i]))


        #         for planetary in self.planetary_gear.planetaries:
        #             if planetary in meshing_chain:

        #                planetary.module = m2
        #                planetary.position = (z2, 0, 0)
        #     power_recirculation=self.planetary_gear.recirculation_power()
        #     print(power_recirculation)
        #     number_good_loop=0
        #     sum_power=0
        #     for power in power_recirculation:
        #         if power[1]>sum_power:
        #             sum_power=power[1]
        #         if power[1]<100:
        #             number_good_loop+=1
        #     if number_good_loop==len(power_recirculation):
        #         list_planetary_gears_2.append(xra)
        #         sum_power_recirculation.append(sum_power)


        if sum_power_recirculation:
            min_power = min(sum_power_recirculation)
            xra = list_planetary_gears[sum_power_recirculation.index(min_power)]
        else:
            return self.planetary_gear




            # print(xra,fx)
        # if res_2.success==True and res_2.fun<0.001:
        #     break
        # if fx<0.001:
        #     break

        # if fx < 0.1:




        x = xra
        X = []
        Y = []
        M = []
        index_x = 0
        teta = 2*m.pi/self.number_planet
        for meshing_chain in meshing_chains:
            M.append(x[index_x])
            index_x += 1
        X.append(0)
        Y.append(x[index_x])
        index_x += 1

        for planets in self.planetary_gear.planets[1:]:
            X.append(x[index_x])
            index_x += 1
            Y.append(x[index_x])
            index_x += 1

        X_prime = []
        Y_prime = []
        for y in range(self.number_planet):
            X_prime.append([])
            Y_prime.append([])
            for i, x in enumerate(X):
                    X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                    Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))



        # list_color = ['mediumblue', 'purple', 'green', 'k']
        # plt.figure()
        # ax=plt.subplot(aspect='equal')

        z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

        for r, meshing_chain in enumerate(meshing_chains):
            z2 = z[r]
            m2 = M[r]
            # color=list_color[i]

            for i, planet in enumerate(self.planetary_gear.planets):
                if planet in meshing_chain:
                    # ax.add_patch(plt.Circle([X[i],Y[i]], m2*planet.Z/2,color=color,fill=False))
                    planet.positions = [(z2, Y[i], X[i])]
                    planet.module = m2
                    d = (((X[i]**2+Y[i]**2)**0.5)*2+m2*planet.Z)
                    if d > self.d_min:
                        self.d_min = d

                    for y in range(self.number_planet):
                        # ax.add_patch(plt.Circle([X_prime[y][i],Y_prime[y][i]], m2*planet.Z/2,color=color,fill=False))
                        planet.positions.append((z2, Y_prime[y][i], X_prime[y][i]))


            for planetary in self.planetary_gear.planetaries:
                if planetary in meshing_chain:
                    # ax.add_patch(plt.Circle([0,0], m2*planetary.Z/2,color=color,fill=False))
                    planetary.module = m2
                    planetary.position = (z2, 0, 0)
                    d = m2*planetary.Z
                    if d > self.d_min:
                        self.d_min = d

        # ax.relim()
        # ax.autoscale_view()
        # plt.show()
        self.planetary_gear.position = True
        power_recirculation = self.planetary_gear.recirculation_power()
        print(power_recirculation)

        return self.planetary_gear


        # else:
        #     self.planetary_gear.sum_Z_planetary = -101
        #     return self.planetary_gear




    def optimize_min(self):

        if self.recirculation:
            self.optimize_min_recirculation()

        meshing_chains = self.planetary_gear.meshing_chain()
        position_min_max = self.position_min_max
        i = 0
        dim_max = [0]*len(meshing_chains)


        if self.planetary_gear.position == True:



            def function_verification(x, planetary_gear, meshing_chains, recirculation_power, l, F2, M_temporary, dim_max=dim_max):

                    X = [0]*len(self.planetary_gear.planets)
                    Y = [0]*len(self.planetary_gear.planets)
                    M = [0]*len(meshing_chains)
                    index_x = 0
                    f2 = []
                    for i in range(len(dim_max)):

                        dim_max[i] = 0
                    for i, meshing_chain in enumerate(meshing_chains):
                        M[i] = x[index_x]
                        index_x += 1

                    X[0] = 0
                    Y[0] = x[index_x]
                    index_x += 1
                    for i, planets in enumerate(planetary_gear.planets[1:]):
                        X[i+1] = x[index_x]
                        index_x += 1
                        Y[i+1] = x[index_x]
                        index_x += 1

                    for i, meshing_chain in enumerate(meshing_chains):

                        m2 = M[i]
                        function_minimize_equation = self.function_minimize_equation(meshing_chain, X, Y, m2)
                        for f in function_minimize_equation:
                            f2.append(f*10)
                        min_max_x, x0 = self.function_minimize_inequation_meshing_chain_min_max(meshing_chain, m2)

                        res_1 = op.minimize(self.function_minimize_inequation_meshing_chain, x0, bounds=min_max_x,
                                            args=(meshing_chain, X, Y))

                        if res_1.fun > 0.0000000000001:
                              f2.append((res_1.fun**(1/2))*10000)

                        if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':

                              if m2*meshing_chain[0].Z*1.3 > self.D_max:
                                  f2.append((self.D_max-m2*meshing_chain[0].Z*1.3))
                              elif m2*meshing_chain[0].Z*1.3 < self.D_min:
                                  f2.append((self.D_min-m2*meshing_chain[0].Z*1.3))
                              dim_max[i] = m2*meshing_chain[0].Z*1.3


                              # if f2[-1] > 0:
                              #      f2[-1] = f2[-1]*100



                        elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                              if m2*meshing_chain[-1].Z*1.3 > self.D_max:
                                  f2.append((self.D_max-m2*meshing_chain[-1].Z*1.3))
                              elif m2*meshing_chain[-1].Z*1.3 < self.D_min:
                                  f2.append((self.D_min-m2*meshing_chain[-1].Z*1.3))
                              dim_max[i] = m2*meshing_chain[-1].Z*1.3

                              # if f2[-1] > 0:
                              #      f2[-1] = f2[-1]*100



                        else:
                            # print(1)
                            for element in meshing_chain:
                                if isinstance(element, Planet):

                                    index_planet = planetary_gear.planets.index(element)
                                    d = (X[index_planet]**2+Y[index_planet]**2)**(1/2)+element.Z*m2/2
                                    if d > self.D_max/2:
                                        f2.append(((self.D_max/2)-(X[index_planet]**2+Y[index_planet]**2)**(1/2)-element.Z*m2/2))

                                    if abs(d)*2 > dim_max[i]:
                                        dim_max[i] = abs(d)*2


                                    # if f2[-1] > 0:
                                    #     f2[-1] = f2[-1]*100

                    f3 = []

                    if self.internal_diameter_min:

                        for i, meshing_chain in enumerate(meshing_chains):
                            if l[1] == 100:
                                M_temporary[i] = M[i]

                            m2 = M_temporary[i]
                            if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Sun':
                                f3.append(-self.internal_diameter_min+m2*meshing_chain[0].Z)



                            elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Sun':
                                f3.append(-self.internal_diameter_min+m2*meshing_chain[-1].Z)




                            else:

                                for element in meshing_chain:
                                    if isinstance(element, Planet):

                                        index_planet = planetary_gear.planets.index(element)

                                        # f3.append((-(self.internal_diameter_min/2)+((X[index_planet]**2+Y[index_planet]**2)**(1/2))))
                                        # if f3[-1]>0:
                                        #     print(1)
                                        f3.append((-(self.internal_diameter_min/2)+((X[index_planet]**2+Y[index_planet]**2)**(1/2))-element.Z*m2/2)*100)
                                        # if limit<0:
                                        #    f3.append(-10)

                                        # print(X[index_planet])
                                        # print(Y[index_planet])

                                        # print(X[index_planet]**2+Y[index_planet]**2)


                                        # print(x)
                        if l[1] == 100:
                            l[1] = 0
                        else:
                            l[1] += 1



                    if  recirculation_power and l[0] == 10:

                        X = []
                        Y = []
                        M = []
                        index_x = 0
                        teta = 2*m.pi/self.number_planet
                        for meshing_chain in meshing_chains:
                            M.append(x[index_x])
                            index_x += 1
                        X.append(0)
                        Y.append(x[index_x])
                        index_x += 1


                        for planets in self.planetary_gear.planets[1:]:
                            X.append(x[index_x])
                            index_x += 1
                            Y.append(x[index_x])
                            index_x += 1




                        X_prime = []
                        Y_prime = []
                        teta = 2*m.pi/self.number_planet
                        for y in range(self.number_planet):
                            X_prime.append([])
                            Y_prime.append([])
                            for i, j in enumerate(X):
                                X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                                Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




                        z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

                        for i, meshing_chain in enumerate(meshing_chains):
                            z2 = z[i]
                            m2 = M[i]


                            for i, planet in enumerate(self.planetary_gear.planets):
                                if planet in meshing_chain:

                                    planet.module = m2
                                    list_positions_planet = [(z2, Y[i], X[i])]




                                    for y in range(self.number_planet):

                                        list_positions_planet.append((z2, Y_prime[y][i], X_prime[y][i]))

                                    planet.positions = list_positions_planet


                            for planetary in self.planetary_gear.planetaries:
                                if planetary in meshing_chain:
                                    planetary.module = m2
                                    planetary.position = (z2, 0, 0)
                        self.planetary_gear.update_position_mech()










                    F = 0

                    for f in f2:
                        F += abs(f)
                    # print(f2)
                    # print(f3)
                    # print('aez')
                    for f in f3:

                        F += abs(f-abs(f))


                    if  recirculation_power and l[0] == 20:
                        F2[0] = 0
                        power_recirculation = self.planetary_gear.recirculation_power()
                        power_max = 0
                        for power in power_recirculation:
                            if power_max < power[1]:
                                power_max = power[1]
                        F2[0] += F*((power_max/10000))
                        # if power_max>self.recircle_power_max:

                        #     print(power_max)
                        #     F2[0]+=power_max*1000
                                # self.optimize_min_recirculation()
                        l[0] = 0

                    F += F2[0]

                    l[0] += 1

                    return F

            min_max_x_2 = []
            x0 = []
            M_temporary = []
            for meshing_chain in meshing_chains:
                min_max_x_2.append([0, self.D_max])
                x0.append(meshing_chain[0].module)
                M_temporary.append(0)


            min_max_x_2.append([-self.D_max, self.D_max])
            x0.append(self.planetary_gear.planets[0].positions[0][1])

            for planets in self.planetary_gear.planets[1:]:
                min_max_x_2.append([-self.D_max, self.D_max])
                min_max_x_2.append([-self.D_max, self.D_max])

                x0.append(planets.positions[0][0])
                x0.append(planets.positions[0][1])

            # for i, meshing_chain in enumerate(meshing_chains):
            #     if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
            #           min_max_x_2.append([self.D_min, self.D_max])

            #           x0.append(self.D_min)

            #     elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
            #           min_max_x_2.append([self.D_min, self.D_max])

            #           x0.append(self.D_min)

            #     else:

            #         for element in meshing_chain:
            #             if isinstance(element, Planet):
            #                 min_max_x_2.append([0, self.D_max])

            #                 x0.append(self.D_min)

            min_x = []
            max_x = []
            for i in range(len(min_max_x_2)):
                min_x.append(min_max_x_2[i][0])
                max_x.append(min_max_x_2[i][1])
            if self.recirculation:
                 recirculation_power = [1]
            else:
                recirculation_power = []

            # print(self.planetary_gear.recirculation_power())
            l = [0, 10]
            F2 = [0]
            s = 1
            if self.D_max-self.D_min < 1:
                s = (self.D_max-self.D_min)*0.1 #TODO Link with recirculation power
            # print(x0)

            xra, fx = cma.fmin(function_verification,
                               self.x0_save,
                               s,
                               args=(self.planetary_gear,
                                     meshing_chains,
                                     recirculation_power,
                                     l, F2, M_temporary),
                               options={'bounds':[min_x, max_x],
                                        'tolfun': 1e-15,
                                        'verbose': 3,
                                        'ftarget': 1e-24, #TODO
                                        'maxiter': 4000, 'tolstagnation':20000})[0:2]
            if fx > 0.00001:
                self.planetary_gear.min_Z_planetary = -200
                return None



            if not self.internal_diameter_min:
                dim_max_min = min(dim_max)

                coeff_reduc = self.D_min/dim_max_min

            else:
                coeff_reduc = 1


            x = xra

            X = []
            Y = []
            M = []
            index_x = 0
            teta = 2*m.pi/self.number_planet
            for meshing_chain in meshing_chains:
                print(15)
                print(M)
                M.append(x[index_x]*coeff_reduc)
                index_x += 1
            X.append(0)
            Y.append(x[index_x]*coeff_reduc)
            index_x += 1

            for planets in self.planetary_gear.planets[1:]:
                X.append(x[index_x]*coeff_reduc)
                index_x += 1
                Y.append(x[index_x]*coeff_reduc)
                index_x += 1

            X_prime = []
            Y_prime = []
            for y in range(self.number_planet):
                X_prime.append([])
                Y_prime.append([])
                for i, x in enumerate(X):
                    X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                    Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




            z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

            for i, meshing_chain in enumerate(meshing_chains):
                z2 = z[i]
                m2 = M[i]


                for i, planet in enumerate(self.planetary_gear.planets):
                    if planet in meshing_chain:

                        position_min_max.enter_module(m2, self.planetary_gear, planet, 'Min')
                        list_positions_planet = [(z2, Y[i], X[i])]




                        for y in range(self.number_planet):

                            list_positions_planet.append((z2, Y_prime[y][i], X_prime[y][i]))

                        position_min_max.enter_position(list_positions_planet, self.planetary_gear, planet, 'Min')


                for planetary in self.planetary_gear.planetaries:
                    if planetary in meshing_chain:
                        position_min_max.enter_module(m2, self.planetary_gear, planetary, 'Min')
                        position_min_max.enter_position((z2, 0, 0), self.planetary_gear, planetary, 'Min')









        return position_min_max

    def optimize_min_recirculation(self):

        meshing_chains = self.planetary_gear.meshing_chain()
        position_min_max = self.position_min_max

        if self.planetary_gear.position == True:
            list_previous_value = []
            list_i = []
            def function_verification(x, planetary_gear, meshing_chains, list_previous_value, list_i, l, F2):

                    X = [0]*len(self.planetary_gear.planets)
                    Y = [0]*len(self.planetary_gear.planets)
                    M = [0]*len(meshing_chains)
                    index_x = 0
                    f2 = []

                    for i, meshing_chain in enumerate(meshing_chains):
                        M[i] = x[index_x]
                        index_x += 1
                    X[0] = 0
                    Y[0] = x[index_x]
                    index_x += 1
                    for i, planets in enumerate(planetary_gear.planets[1:]):
                        X[i+1] = x[index_x]
                        index_x += 1
                        Y[i+1] = x[index_x]
                        index_x += 1

                    for i, meshing_chain in enumerate(meshing_chains):

                        m2 = M[i]

                        f2.extend(self.function_minimize_equation(meshing_chain, X, Y, m2)*10)

                        min_max_x, x0 = self.function_minimize_inequation_meshing_chain_min_max(meshing_chain, m2)

                        res_1 = op.minimize(self.function_minimize_inequation_meshing_chain, x0, bounds=min_max_x,
                                            args=(meshing_chain, X, Y))

                        if res_1.fun > 0.000000000001:
                              f2.append(res_1.fun*1000000)

                        if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                              f2.append(x[index_x]-m2*meshing_chain[0].Z*1.3)
                              index_x += 1

                        elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                              f2.append(x[index_x]-m2*meshing_chain[-1].Z*1.3)
                              index_x += 1

                        else:

                            for element in meshing_chain:
                                if isinstance(element, Planet):
                                    index_planet = planetary_gear.planets.index(element)
                                    f2.append(((x[index_x]/2)**2-(X[index_planet]**2+Y[index_planet]**2)))
                                    index_x += 1



                    F = 0

                    for f in f2:
                        F += abs(f)

                    if  l[0] == 20 or F2 == 0:
                        X = []
                        Y = []
                        M = []
                        index_x = 0
                        teta = 2*m.pi/self.number_planet
                        for meshing_chain in meshing_chains:
                            M.append(x[index_x])
                            index_x += 1
                        X.append(0)
                        Y.append(x[index_x])
                        index_x += 1


                        for planets in self.planetary_gear.planets[1:]:
                            X.append(x[index_x])
                            index_x += 1
                            Y.append(x[index_x])
                            index_x += 1




                        X_prime = []
                        Y_prime = []
                        teta = 2*m.pi/self.number_planet
                        for y in range(self.number_planet):
                            X_prime.append([])
                            Y_prime.append([])
                            for i, x in enumerate(X):
                                X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                                Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




                        z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

                        for i, meshing_chain in enumerate(meshing_chains):
                            z2 = z[i]
                            m2 = M[i]


                            for i, planet in enumerate(self.planetary_gear.planets):
                                if planet in meshing_chain:

                                    planet.module = m2
                                    list_positions_planet = [(z2, Y[i], X[i])]




                                    for y in range(self.number_planet):

                                        list_positions_planet.append((z2, Y_prime[y][i], X_prime[y][i]))

                                    planet.positions = list_positions_planet


                            for planetary in self.planetary_gear.planetaries:
                                if planetary in meshing_chain:
                                    planetary.module = m2
                                    planetary.position = (z2, 0, 0)

                        self.planetary_gear.update_position_mech()
                        try:
                            power_recirculation = self.planetary_gear.recirculation_power()
                        except:
                            power_recirculation = [[10000000, 10000000]]

                        index = 0
                        power_max = 0

                        for power in power_recirculation:
                            if power[1] > power_max:
                                power_max = power[1]
                            if power[1] < self.recircle_power_max:
                                index += 1

                        F2[0] = power_max*F
                        l[0] = 0
                    else:
                        l[0] += 1







                    F += F2[0]

                    if int(F) not in list_previous_value:
                        list_previous_value.append(int(F))

                        list_i.append(0)
                    else:
                        index = list_previous_value.index(int(F))
                        list_i[index] += 1

                        if list_i[index] > 400:
                            # if power_max<self.recircle_power_max:
                            #     return 0
                            # else:
                            return 0




                    return F

            min_max_x_2 = []
            x0 = []

            for meshing_chain in meshing_chains:
                min_max_x_2.append([0, self.D_max/7])
                x0.append(meshing_chain[0].module)

            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            x0.append(self.planetary_gear.planets[0].positions[0][1])

            for planets in self.planetary_gear.planets[1:]:
                min_max_x_2.append([-self.D_max/2, self.D_max/2])
                min_max_x_2.append([-self.D_max/2, self.D_max/2])

                x0.append(planets.positions[0][2])
                x0.append(planets.positions[0][1])

            for meshing_chain in meshing_chains:
                min_max_x_2.append([0, self.D_max/7])
                x0.append(meshing_chain[0].module)

            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            x0.append(self.planetary_gear.planets[0].positions[0][1])

            for planets in self.planetary_gear.planets[1:]:
                min_max_x_2.append([-self.D_max/2, self.D_max/2])
                min_max_x_2.append([-self.D_max/2, self.D_max/2])

                x0.append(planets.positions[0][2])
                x0.append(planets.positions[0][1])

            for i, meshing_chain in enumerate(meshing_chains):
                if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                      min_max_x_2.append([self.D_min, self.D_max])

                      x0.append(self.d_min)

                elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                      min_max_x_2.append([self.D_min, self.D_max])

                      x0.append(self.d_min)

                else:

                    for element in meshing_chain:
                        if isinstance(element, Planet):
                            min_max_x_2.append([0, self.D_max])

                            x0.append(self.d_min)

            min_x = []
            max_x = []
            for i in range(len(min_max_x_2)):
                min_x.append(min_max_x_2[i][0])
                max_x.append(min_max_x_2[i][1])

            succes = 0
            x = 1

            while succes == False or x == 1:
                self.verification_recirculation()
                meshing_chains = self.planetary_gear.meshing_chain()
                min_max_x_2 = []
                x0 = []

                for meshing_chain in meshing_chains:
                    min_max_x_2.append([0, self.D_max/7])
                    x0.append(meshing_chain[0].module)

                min_max_x_2.append([-self.D_max/2, self.D_max/2])
                x0.append(self.planetary_gear.planets[0].positions[0][1])

                for planets in self.planetary_gear.planets[1:]:
                    min_max_x_2.append([-self.D_max/2, self.D_max/2])
                    min_max_x_2.append([-self.D_max/2, self.D_max/2])

                    x0.append(planets.positions[0][2])
                    x0.append(planets.positions[0][1])

                for i, meshing_chain in enumerate(meshing_chains):
                    if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                          min_max_x_2.append([self.D_min, self.D_max])

                          x0.append(self.D_min)

                    elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                          min_max_x_2.append([self.D_min, self.D_max])

                          x0.append(self.D_min)

                    else:

                        for element in meshing_chain:
                            if isinstance(element, Planet):
                                min_max_x_2.append([self.D_min, self.D_max])

                                x0.append(self.D_min)

                min_x = []
                max_x = []
                for i in range(len(min_max_x_2)):
                    min_x.append(min_max_x_2[i][0])
                    max_x.append(min_max_x_2[i][1])
                # self.planetary_gear.babylonjs()
                l = [0]
                F2 = [0]
                res = op.minimize(function_verification, x0, args=(self.planetary_gear, meshing_chains, list_previous_value, list_i, l, F2),
                                  bounds=min_max_x_2, options={'maxiter': 100, 'ftol': 1e-12 * np.finfo(float).eps})
                x = res.fun
                succes = res.success
                list_previous_value = []
                list_i = []



            x = res.x

            X = []
            Y = []
            M = []
            index_x = 0
            teta = 2*m.pi/self.number_planet
            for meshing_chain in meshing_chains:
                M.append(x[index_x])
                index_x += 1

            X.append(0)
            Y.append(x[index_x])
            index_x += 1

            for planets in self.planetary_gear.planets[1:]:
                X.append(x[index_x])
                index_x += 1
                Y.append(x[index_x])
                index_x += 1

            X_prime = []
            Y_prime = []
            for y in range(self.number_planet):
                X_prime.append([])
                Y_prime.append([])
                for i, x in enumerate(X):
                    X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                    Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




            z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

            for i, meshing_chain in enumerate(meshing_chains):
                z2 = z[i]
                m2 = M[i]


                for i, planet in enumerate(self.planetary_gear.planets):
                    if planet in meshing_chain:
                        d = (((X[i]**2+Y[i]**2)**0.5)*2+m2*planet.Z)
                        if d > self.d_min:
                            self.d_min = d



                for planetary in self.planetary_gear.planetaries:
                    if planetary in meshing_chain:
                        d = m2*planetary.Z
                        if d > self.d_min:
                           self.d_min = d

















    def optimize_max(self):
        meshing_chains = self.planetary_gear.meshing_chain()
        position_min_max = self.position_min_max
        if self.planetary_gear.position == True:

            def function_verification(x, planetary_gear, meshing_chains):

                    X = [0]*len(self.planetary_gear.planets)
                    Y = [0]*len(self.planetary_gear.planets)
                    M = [0]*len(meshing_chains)
                    index_x = 0
                    f2 = []

                    for i, meshing_chain in enumerate(meshing_chains):
                        M[i] = x[index_x]
                        index_x += 1
                    X[0] = 0
                    Y[0] = x[index_x]
                    index_x += 1
                    for i, planets in enumerate(planetary_gear.planets[1:]):
                        X[i+1] = x[index_x]
                        index_x += 1
                        Y[i+1] = x[index_x]
                        index_x += 1

                    for i, meshing_chain in enumerate(meshing_chains):

                        m2 = M[i]
                        f2.extend(self.function_minimize_equation(meshing_chain, X, Y, m2))

                        min_max_x, x0 = self.function_minimize_inequation_meshing_chain_min_max(meshing_chain, m2)

                        res_1 = op.minimize(self.function_minimize_inequation_meshing_chain, x0, bounds=min_max_x, args=(meshing_chain, X, Y))

                        if res_1.fun > 0.000000000001:
                              f2.append(res_1.fun*100000000)

                        if isinstance(meshing_chain[0], Planetary) and meshing_chain[0].planetary_type == 'Ring':
                              f2.append((self.D_max-m2*meshing_chain[0].Z*1.3)/100)
                              if f2[-1] < 0:
                                  f2[-1] = f2[-1]*100000000
                              index_x += 1

                        elif isinstance(meshing_chain[-1], Planetary) and meshing_chain[-1].planetary_type == 'Ring':
                              f2.append((self.D_max-m2*meshing_chain[-1].Z*1.3)/100)
                              if f2[-1] < 0:
                                  f2[-1] = f2[-1]*100000000
                              index_x += 1

                        else:

                            for element in meshing_chain:
                                if isinstance(element, Planet):
                                    index_planet = planetary_gear.planets.index(element)
                                    f2.append(((self.D_max/2)**2-(X[index_planet]**2+Y[index_planet]**2))/100)
                                    if f2[-1] < 0:
                                        f2[-1] = f2[-1]*100000000
                                    index_x += 1

                    # for planet in self.planetary_gears.planets:
                    #     planet.position







                    F = 0
                    # print(f2)
                    for f in f2:
                        F += abs(f)
                    # print(f2)
                    # print(x)
                    # print(F)
                    return F

            min_max_x_2 = []
            x0 = []
            for meshing_chain in meshing_chains:
                min_max_x_2.append([0, self.D_max/7])
                x0.append(meshing_chain[0].module)

            min_max_x_2.append([-self.D_max/2, self.D_max/2])
            x0.append(self.planetary_gear.planets[0].positions[0][1])

            for planets in self.planetary_gear.planets[1:]:
                min_max_x_2.append([-self.D_max/2, self.D_max/2])
                min_max_x_2.append([-self.D_max/2, self.D_max/2])

                x0.append(planets.positions[0][2])
                x0.append(planets.positions[0][1])

            min_x = []
            max_x = []
            for i in range(len(min_max_x_2)):
                min_x.append(min_max_x_2[i][0])
                max_x.append(min_max_x_2[i][1])
            s = 0.1
            if (self.D_max-self.D_min) < 1:
                s = (self.D_max-self.D_min)*0.1

            xra, fx = cma.fmin(function_verification, x0, s, args=(self.planetary_gear, meshing_chains), options={'bounds':[min_x, max_x],
                                                                                                                  'tolfun': 1e-10,
                                                                                                                  'verbose': 3,
                                                                                                                  'ftarget': 1e-8,
                                                                                                                  'maxiter': 2000})[0:2]
            if fx > 1:
                self.planetary_gear.min_Z_planetary = -200
                return None


            x = xra

            X = []
            Y = []
            M = []
            index_x = 0
            teta = 2*m.pi/self.number_planet
            for meshing_chain in meshing_chains:
                M.append(x[index_x])
                index_x += 1
            X.append(0)
            Y.append(x[index_x])
            index_x += 1
            for planets in self.planetary_gear.planets[1:]:
                X.append(x[index_x])
                index_x += 1
                Y.append(x[index_x])
                index_x += 1

            X_prime = []
            Y_prime = []
            for y in range(self.number_planet):
                X_prime.append([])
                Y_prime.append([])
                for i, x in enumerate(X):
                    X_prime[y].append(X[i]*m.cos(y*teta)-Y[i]*m.sin(y*teta))
                    Y_prime[y].append(X[i]*m.sin(y*teta)+Y[i]*m.cos(y*teta))




            z = self.planetary_gear.meshing_chain_position_z(meshing_chains)

            for i, meshing_chain in enumerate(meshing_chains):
                z2 = z[i]
                m2 = M[i]


                for i, planet in enumerate(self.planetary_gear.planets):
                    if planet in meshing_chain:

                        position_min_max.enter_module(m2, self.planetary_gear, planet, 'Max')
                        list_positions_planet = [(z2, Y[i], X[i])]


                        for y in range(self.number_planet):

                            list_positions_planet.append((z2, Y_prime[y][i], X_prime[y][i]))

                        position_min_max.enter_position(list_positions_planet, self.planetary_gear, planet, 'Max')



                for planetary in self.planetary_gear.planetaries:
                    if planetary in meshing_chain:

                        position_min_max.enter_module(m2, self.planetary_gear, planetary, 'Max')
                        position_min_max.enter_position((z2, 0, 0), self.planetary_gear, planetary, 'Max')





        return position_min_max

class GeneratorPlanetaryGears(DessiaObject):
    _eq_is_data_eq = False
    def __init__(self, number_shafts: int, speed_shafts: List[List[Tuple[float, float]]], torque_shafts: List[List[Tuple[float, float]]],
                 number_max_planet: int, D_min: int, D_max: int, internal_diameter_min: float = 0, reason_min_max: List[Tuple[float, float]] = [],
                 speed_planet_carrer: Tuple[float, float] = [], torque_planet_carrer: Tuple[float, float] = [], number_solution: int = 100000000000):

        self.number_shafts = number_shafts
        self.speed_shafts = speed_shafts
        self.torque_shaft = torque_shafts
        if speed_shafts:
            number_input = len(speed_shafts[0])
        else:
            number_input = len(reason_min_max)+2
        self.number_max_meshing_plan = number_input
        self.number_min_meshing_plan = int(number_input/2)

        self.number_max_planet = number_max_planet
        self.number_junction_max = int(number_input/2)-1
        self.reason_min_max = reason_min_max
        self.D_min = D_min
        self.D_max = D_max
        self.torque_input = []
        for i in range(number_input-1):
            self.torque_input.append((-2, 2))
        self.number_input = number_input
        self.number_solution = number_solution
        # if not number_solution:
        #     self.number_solution = 100000000000
        self.speed_max_planet = 3*40
        self.internal_diameter_min = internal_diameter_min
        # self.speed_planet_carrer=speed_planet_carrer
        # self.torque_planet_carrer=torque_planet_carrer





    def speed_conversion(self):
        speed_possibility = []
        alpha_min_max = []
        for shaft in range(self.number_shafts):
            number_speed_planet_carrer = shaft
            speed_possibility2 = [(10, 20)]
            alpha_min_max2 = []

            number_total = [v for v in range(self.number_shafts)]
            number_total.remove(shaft)

            for element in range(len(number_total)-1):

                number_planetary_input = number_total[0]
                number_planetary_output = number_total[element+1]
                liste_alpha = []
                for i, speed_shaft in enumerate(self.speed_shafts):
                    alpha2 = []
                    for speed_planet_carrer in speed_shaft[number_speed_planet_carrer]:

                        for speed_input in  speed_shaft[number_planetary_input]:

                            for speed_output in speed_shaft[number_planetary_output]:

                                alpha2.append((speed_output-speed_planet_carrer)/(speed_input-speed_planet_carrer))
                                # print(speed_planet_carrer,speed_output,speed_input)
                                # print(alpha2)

                    liste_alpha.append([min(alpha2), max(alpha2)])
                # print(liste_alpha)
                alpha_total_min = -np.inf
                alpha_total_max = np.inf

                for alpha in liste_alpha:
                    if alpha[0] > alpha_total_min:
                        alpha_total_min = alpha[0]
                    if alpha[1] < alpha_total_max:
                        alpha_total_max = alpha[1]
                # print(alpha_total_min,alpha_total_max)
                if alpha_total_min <= alpha_total_max:

                    list_w3 = []
                    for alpha in [alpha_total_min, alpha_total_max]:
                        for w1 in [10, 20]:
                            for w2 in [30, 40]:
                                list_w3.append(alpha*w1+w2*(1-alpha))

                    w3min = min(list_w3)
                    w3max = max(list_w3)

                    speed_possibility2.append((w3min, w3max))
                    alpha_min_max2.append((alpha_total_min, alpha_total_max))

            if len(speed_possibility2) == self.number_shafts:
                speed_possibility.append(speed_possibility2)
                alpha_min_max.append(alpha_min_max2)
        return speed_possibility, alpha_min_max

    def reason_conversion(self):
        speed_possibility = [[10, 20]]
        w_max = 40
        for list_alpha in self.reason_min_max:
            list_w3 = []
            for alpha in list_alpha:
                for w1 in [10, 20]:
                    for w2 in [30, 40]:
                        list_w3.append(alpha*w1+w2*(1-alpha))


            w3min = min(list_w3)
            w3max = max(list_w3)
            if w3max > w_max:
                w_max = w3max
            speed_possibility.append([w3min, w3max])
        self.speed_max_planet = w_max*3
        return speed_possibility







    def generator(self):
        list_solution_planet_strucuture = []
        for m in range(self.number_min_meshing_plan, self.number_max_meshing_plan+1):

            for p in range(self.number_max_planet):

                for j in range(self.number_junction_max+1):

                    generator = GeneratorPlanetsStructure(number_max_planet=p+1, number_max_meshing_plan=m, number_junction=j,
                                                          number_max_junction_by_planet=2, min_planet_branch=1)
                    list_solution_planet_strucuture.extend(generator.decision_tree())

        input_speed = [0]*self.number_input
        generator = GeneratorPlanetaryGearsArchitecture(planet_structures=list_solution_planet_strucuture, input_speeds=input_speed)

        list_solution_planetary_gear_architecture = generator.decision_tree()

        if not self.reason_min_max:
            result = self.speed_conversion()
            speed_input = result[0][0]
            alpha_min_max = result[1][0]
        else:
            speed_input = self.reason_conversion()
            alpha_min_max = self.reason_min_max
        list_planetary_gears = []
        i = 0

        while len(list_planetary_gears) < self.number_solution*100 and i <= len(list_solution_planetary_gear_architecture)-1:
            planetary_gear = list_solution_planetary_gear_architecture[i]
            generator = GeneratorPlanetaryGearsZNumberReason(planetary_gear=planetary_gear, input_reason=alpha_min_max,
                                                             input_speeds=speed_input, input_torques=self.torque_input,
                                                             speed_planet_carrer=(30, 40), torque_planet_carrer=(-2, 2),
                                                             Z_range_sun=[15, 80], Z_range_ring=[40, 100],
                                                             number_planet=3, number_solution=self.number_solution*15,
                                                             speed_max_planet=self.speed_max_planet, different_solution=True)
            solutions = generator.decision_tree()
            if solutions:
                list_planetary_gears.extend(solutions)




            i += 1
        planetary_gear_results = []
        if not list_planetary_gears:
            return None

        else:

            Z_min = [planetary_gear.z_min for planetary_gear in list_planetary_gears]
            print(Z_min)
            sucess = False
            num_planetary_gear = 0
            n = copy.copy(len(list_planetary_gears))
            while not sucess or num_planetary_gear > n-1 or Z_min == []:
                print(n)
                print(Z_min)
                if not Z_min:
                    break
                Z_max = max(Z_min)
                planetary_gear = list_planetary_gears[Z_min.index(Z_max)]
                print(planetary_gear)
                Z_min.remove(Z_max)
                list_planetary_gears.remove(planetary_gear)
                # if len(planetary_gear.planets)>2:
                #         continue
                if self.internal_diameter_min:

                    for planetary in planetary_gear.planetaries:

                        if planetary.planetary_type == 'Sun':
                            print(planetary.Z)
                            if planetary.Z < self.internal_diameter_min/0.005:
                                continue
                print(planetary_gear)
                generator = GeneratorPlanetaryGearsGeometry(planetary_gear=planetary_gear,
                                                            number_planet=3, D_min=self.D_min,
                                                            D_max=self.D_max,
                                                            recircle_power_max=250,
                                                            internal_diameter_min=self.internal_diameter_min)
                generator.verification()
                if generator.planetary_gear.sum_Z_planetary < 0:
                     continue

                position_min_max = generator.optimize_min()


                if generator.planetary_gear.min_Z_planetary > 0:


                    # generator.planetary_gear.babylonjs()
                    planetary_gear_results.append(PlanetaryGearResult(planetary_gear=generator.planetary_gear, position_min_max=position_min_max))

                    planetary_gear_results[-1].update_geometry()
                    planetary_gear_results[-1].planetary_gear.update_length()
                    # planetary_gear_results[-1].babylonjs()
                    # plot_data.plot_canvas(planetary_gear_results[-1].plot_data())
                    if len(planetary_gear_results) == self.number_solution:
                        sucess = True
                num_planetary_gear += 1
            return planetary_gear_results
















class SolutionSort():

    def __init__(self, planetary_gears_list: List[List[PlanetaryGear]]):

        self.planetary_gears_list = planetary_gears_list

    def list_planetary_gears(self):
        planetary_gears = []
        list_solution = []
        Z_planetary = []
        for planetary_gear in self.planetary_gears_list:
            planetary_gears.extend(planetary_gear)
        for planetary_gear in planetary_gears:
            Z = []
            for planetary in planetary_gear.planetaries:
                Z.append(planetary.Z)
            if not Z  in Z_planetary:
                Z_planetary.append(Z)
                list_solution.append(planetary_gear)
                if len(list_solution) > 100:

                    return list_solution
        return list_solution



        # for planetary_gear in self.planetary_gears_list:
        #     if len(planetary_gear)>10:
        #          for planetary in planetary_gear:
        #             planetary.sum_Z_planetary = 0

        #             planetary.max_Z_planetary = 0
        #             planetary.min_Z_planetary = 100000

        #             for planetary2 in planetary.planetaries:
        #                 planetary.sum_Z_planetary += planetary2.Z

        #                 if planetary.max_Z_planetary < planetary2.Z:
        #                     planetary.max_Z_planetary = planetary2.Z

        #                 if planetary.min_Z_planetary > planetary2.Z:
        #                     planetary.min_Z_planetary = planetary2.Z

        #          planetary_gears.extend((planetary_gear[:10]))
        #          for planetary in planetary_gears:
        #              print(planetary)


        #     else:
        #         for planetary in planetary_gear:
        #             planetary.sum_Z_planetary = 0

        #             planetary.max_Z_planetary = 0
        #             planetary.min_Z_planetary = 100000

        #             for planetary2 in planetary.planetaries:
        #                 planetary.sum_Z_planetary += planetary2.Z

        #                 if planetary.max_Z_planetary < planetary2.Z:
        #                     planetary.max_Z_planetary = planetary2.Z

        #                 if planetary.min_Z_planetary > planetary2.Z:
        #                     planetary.min_Z_planetary = planetary2.Z

        #         planetary_gears.extend(copy.deepcopy(planetary_gear))


        #     if len(planetary_gears)>10:
        #         return planetary_gears

        # return planetary_gears

    def solution_sort(self):
        planetary_gears_list = self.planetary_gears_list

        list_solution = []
        Z_planetary = []
        planetary_gears = []
        Z_new_type = 0
        Z_planetary2 = []
        for planetary_gear in planetary_gears_list:
            planetary_gears.extend(planetary_gear)
        for planetary_gear in planetary_gears:
            Z = []

            for planetary in planetary_gear.planetaries:
                # Z.append(planetary.planetary_type)
                Z.append(planetary.planetary_type)

            if not Z  in Z_planetary:
                list_solution2 = []
                Z_new_type = 1
                Z_planetary.append(Z)
            if Z_new_type == 1:
                Z2 = []
                for planetary in planetary_gear.planetaries:
                    Z2.append(planetary.Z)
                if not Z2 in Z_planetary2:
                    Z_planetary2.append(Z2)
                    list_solution2.append(planetary_gear)

                    if len(list_solution2) > 0:
                        print(len(planetary_gears))
                        list_solution.extend(list_solution2)

                        Z_new_type = 0

        return list_solution


# class CatalogClass():

#     def __init__(self, planetary_gears: List[PlanetaryGear]):
#         self.planetary_gears = planetary_gears

#     def catalog(self):
#         variables = ['name', 'sum_Z_planetary', 'd_min', 'sum_speed_planetary', 'speed_planet_carrer', 'min_Z_planetary', 'max_Z_planetary']
#         array = []

#         choice_args = ['sum_Z_planetary', 'd_min', 'sum_speed_planetary', 'speed_planet_carrer', 'min_Z_planetary', 'max_Z_planetary']

#         minimized_attributes = {'sum_Z_planetary':True, 'd_min': True, 'sum_speed_planetary':False, 'speed_planet_carrer':False, 'min_Z_planetary':True, 'max_Z_planetary':True}

#         pareto_settings = ParetoSettings(minimized_attributes=minimized_attributes,
#                                          enabled=True)

#         for i, planetary_gear in enumerate(self.planetary_gears):


#             array.append((planetary_gear.name+str(i), planetary_gear.sum_Z_planetary, planetary_gear.d_min,
#                           planetary_gear.sum_speed_planetary, planetary_gear.speed_planet_carrer, planetary_gear.min_Z_planetary,
#                           planetary_gear.max_Z_planetary))


#         catalog = Catalog(array=array,
#                           variables=variables,
#                           choice_variables=choice_args,
#                           objectives=[],
#                           pareto_settings=pareto_settings,
#                           name='Planetary_gears')

#         return [catalog], self.planetary_gears











