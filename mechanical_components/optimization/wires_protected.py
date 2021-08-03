#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""

"""

from typing import List
import networkx as nx
from .common import RoutingOptimizer
import mechanical_components.wires as wires
import volmdlr.primitives3d as vmp3d
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import volmdlr as vm

# class WiringOptimizer(RoutingOptimizer):

#     def NumberHarnesses(self, paths):
#         G = nx.Graph()
#         G.add_nodes_from(self.line_graph)
#         for path in paths:
#             routes = [(w1, w2) for w1, w2 in zip(path[:-1], path[1:])]
#             for route1, route2 in zip(routes[:-1], routes[1:]):
#                 G.add_edge(route1, route2)
        
#         # Removing nodes with degree == 0
#         for node in self.line_graph.nodes():
#             if G.degree(node) == 0:
#                 G.remove_node(node)
# #        nx.draw_kamada_kawai(G)
#         return nx.number_connected_components(G)
        
#     def route(self, wires_specs: List[wires.RoutingSpec]):
#         shortest_paths = []
#         shortest_paths_lengths = []
#         for wire_spec in wires_specs:
#             # Checking if pat is in cache

#             if (wire_spec.source, wire_spec.destination) in \
#                     self._shortest_paths_cache:
#                 shortest_path = self._shortest_paths_cache[
#                     (wire_spec.source, wire_spec.destination)]
#             elif (wire_spec.source, wire_spec.destination) in \
#                     self._shortest_paths_cache:
#                 shortest_path = self._shortest_paths_cache[
#                     (wire_spec.source, wire_spec.destination)]
#             else:            
#                 shortest_path = nx.shortest_path(self.graph,
#                                                  wire_spec.source,
#                                                  wire_spec.destination,
#                                                  weight='distance')
#                 self._shortest_paths_cache[
#                     (wire_spec.source, wire_spec.destination)] = shortest_path
                
#             shortest_paths.append(shortest_path)
#             length = self.PathLength(shortest_path)
#             shortest_paths_lengths.append(length)
            
#         wires2 = []
#         for ipath, path in enumerate(shortest_paths):
#             wires2.append(wires.Wire(path, wires_specs[ipath].diameter,
#                                      color=wires_specs[ipath].color,
#                                      name=wires_specs[ipath].name))
#         return wires.Wiring(wires2, [])

#     def volmdlr_primitives(self):
#         volumes = []
#         for route in self.routes:
#             point1, point2 = route[0], route[1]
#             if point1 != point2:
#                 cylinder = vmp3d.Cylinder((point1+point2)/2, point2-point1, 0.05,
#                                           (point2-point1).norm())
#                 volumes.append(cylinder)
#         return volumes


class WiringOptimizer(RoutingOptimizer):

    def NumberHarnesses(self, paths):
        G = nx.Graph()
        G.add_nodes_from(self.line_graph)
        for path in paths:
            routes = [(w1, w2) for w1, w2 in zip(path[:-1], path[1:])]
            for route1, route2 in zip(routes[:-1], routes[1:]):
                G.add_edge(route1, route2)
        
        # Removing nodes with degree == 0
        for node in self.line_graph.nodes():
            if G.degree(node) == 0:
                G.remove_node(node)
#        nx.draw_kamada_kawai(G)
        return nx.number_connected_components(G)
        
    def route(self, wires_specs: List[wires.RoutingSpec]):
        shortest_paths = []
        shortest_paths_lengths = []
        for wire_spec in wires_specs:
            # Checking if pat is in cache

            if (wire_spec.source, wire_spec.destination) in \
                    self._shortest_paths_cache:
                shortest_path = self._shortest_paths_cache[
                    (wire_spec.source, wire_spec.destination)]
            elif (wire_spec.source, wire_spec.destination) in \
                    self._shortest_paths_cache:
                shortest_path = self._shortest_paths_cache[
                    (wire_spec.source, wire_spec.destination)]
            else:            
                shortest_path = nx.shortest_path(self.graph,
                                                 wire_spec.source,
                                                 wire_spec.destination,
                                                 weight='distance')
                self._shortest_paths_cache[
                    (wire_spec.source, wire_spec.destination)] = shortest_path
                
            shortest_paths.append(shortest_path)
            length = self.PathLength(shortest_path)
            shortest_paths_lengths.append(length)
            
        wires2 = []
        for ipath, path in enumerate(shortest_paths):
            wires2.append(wires.Wire(path, wires_specs[ipath].diameter,
                                     color=wires_specs[ipath].color,
                                     name=wires_specs[ipath].name))
        return wires.Wiring(wires2, [])


    def route2(self, wires_specs: List[wires.RoutingSpec]):

        counter = 0
        # self.plot()
        for i in range(len(wires_specs)):
            shortest_paths = []
            path_changed = False
            for wire_spec in wires_specs:
                shortest_path = nx.shortest_path(self.graph,
                                         wire_spec.source,
                                         wire_spec.destination,
                                         weight='distance')
                shortest_paths.append(shortest_path)
                if (wire_spec.source, wire_spec.destination) in self._shortest_paths_cache:
                    if self._shortest_paths_cache[(wire_spec.source, wire_spec.destination)] != shortest_path:
                        path_changed = True
                        self._shortest_paths_cache[(wire_spec.source, wire_spec.destination)] = shortest_path
                else:
                    self._shortest_paths_cache[(wire_spec.source, wire_spec.destination)] = shortest_path
                    path_changed = True

            if path_changed:
                print('changed path', i+1, 'times')
                self.plot_routes(shortest_paths, wires_specs)
                over_all_path_cost = self.over_all_path_cost(shortest_paths,self._graph)
                print('over_all_path_cost :', over_all_path_cost)
                update_graph = self.update_graph(shortest_paths)
            else:
                if counter < 3 and i < len(wires_specs) * 1:
                    counter +=1
                    self.update_graph(shortest_paths)
                else:
                    break

        wires2 = []
        for ipath, path in enumerate(shortest_paths):
            wires2.append(wires.Wire(path, wires_specs[ipath].diameter,
                                     color=wires_specs[ipath].color,
                                     name=wires_specs[ipath].name))

        return wires.Wiring(wires2, [])

    def multi_source_multi_destination_routing(self,wires_specs: List[wires.RoutingSpec]):
        
        ''' 
        Calculaes the wiring routing for multi-source-multi-destination problems
        :param wires_specs: list of the wires specifications 
        returns a list of wires for the shortest paths
        '''
        list_shortest_paths = []
        for i in range(len(wires_specs)):
            wires_specs_ = wires_specs[i:] + wires_specs[:i]
            shortest_paths = []
            for wire_spec in wires_specs_:
                shortest_path = nx.shortest_path(self.graph,
                                         wire_spec.source,
                                         wire_spec.destination,
                                         weight='distance')
                self.update_graph([shortest_path])
                shortest_paths.append(shortest_path)

            # self.plot_routes(shortest_paths, wires_specs)
            over_all_path_cost = self.over_all_path_cost(shortest_paths,self._graph)
            list_shortest_paths.append([over_all_path_cost, shortest_paths])
            # print('over_all_path_cost :', over_all_path_cost)
            self.restart_graph()

        list_shortest_paths.sort()
        shortest_paths = list_shortest_paths[0][1]
        # print('final over_all_path_cost :', list_shortest_paths[0][0])

        wires2 = []
        for ipath, path in enumerate(shortest_paths):
            wires2.append(wires.Wire(path, wires_specs[ipath].diameter,
                                     color=wires_specs[ipath].color,
                                     name=wires_specs[ipath].name))

        return wires.Wiring(wires2, [])

    def single_source_multi_destination_routing(self, wires_specs):
        ''' 
        Uses the method EAMDSP - Efficient Algorithm for Multi-Destination Shortest Path. 
        :param wires_specs: list of the wires specifications 
        returns a list of wires for the shortest paths
        '''
        paths = []
        previous_path = []
        shortest_paths =[]
        finished = False
        s = wires_specs[0].source
        destinations = [wire_spec.destination for wire_spec in wires_specs]
        while not finished:
            paths_dict = []
            visiting_nodes = []
            for i, destination in enumerate(destinations):

                shortest_path = nx.shortest_path(self.graph,
                                         s,
                                         destination,
                                         weight='distance')
                path_length = self.PathLength(shortest_path)
                for node in shortest_path:
                    visiting_nodes.append([node, destination])
                paths_dict.append([path_length, destination])

            paths_dict.sort()
            partial_path = []
            e_Dest_Node = paths_dict[0][1]
            for visisting_node, visisting_node_dest in visiting_nodes:
                if visisting_node_dest == e_Dest_Node:
                    if not paths or visisting_node != paths[-1]:
                        paths.append(visisting_node)
                        partial_path.append(visisting_node)
            # print('previous_path :',previous_path)
            # print('partial_path :', partial_path)
            
            path = []
            for node in partial_path:
                if node in previous_path:
                    previous_path.remove(previous_path[-1])
                else:
                    path.append(node)
                    
            shortest_path = previous_path[:] + path[:]
            previous_path = shortest_path
                    
            shortest_paths.append(shortest_path)
            s = e_Dest_Node
            destinations.remove(e_Dest_Node)
            
            if len(destinations) == 0:
                finished = True
        # print('shortest path :', shortest_paths)
        self.plot_routes([shortest_paths[2]], wires_specs)
        over_all_path_cost = self.over_all_path_cost([paths], self.graph)
        print('over_all_path_cost EAMDSP :', over_all_path_cost)
        return paths

    def over_all_path_cost(self, shortest_paths, graph = None):

        '''
        Calculates the over all distance for the shotest paths. The calculation does not take into the same edge more than once

        :param shortest_paths: list of shortest paths
        returns the sum all the paths distances, but considering sharing edges just once
        '''
        if graph == None:
            graph = self.graph
        over_all_path_cost = 0
        dict_shortest_paths_edges = {}
        for shortest_path in shortest_paths:
            for node1, node2 in zip(shortest_path[:-1], shortest_path[1:]):
                if (node1, node2) not in dict_shortest_paths_edges and (node2, node1) not in dict_shortest_paths_edges and (node1, node2) in graph.edges:
                    dict_shortest_paths_edges[(node1, node2)] = graph.edges[node1, node2]['distance']

        over_all_path_cost = sum(list(dict_shortest_paths_edges.values()))

        return over_all_path_cost

    def plot_routes(self, shortest_paths, wires_specs):
        
        '''
        plots the routes resulted from the routing method
        '''
        ax = list(self.graph.nodes)[0].plot()
        for start, end in list(self.graph.edges):
            start.plot(ax=ax, color = 'gray')
            end.plot(ax=ax, color = 'gray')
        for shortest_path in shortest_paths:
            for node1, node2 in zip(shortest_path[:-1], shortest_path[1:]):
                if (node1, node2) in self.graph.edges:
                    line = vm.edges.LineSegment3D(node1, node2)
                    line.plot(ax=ax, color = 'r')
        for wire_spec in wires_specs:
            wire_spec.source.plot(ax=ax, color = 'b')
            wire_spec.destination.plot(ax=ax, color = 'y')
            
    def update_graph(self, shortest_paths):
        '''
        updates the graph
        :param shortest paths: list of shortest paths
        '''

        shortest_paths_edges = []
        shortest_paths_nodes = []
        
        for shortest_path in shortest_paths:
            for node1, node2 in zip(shortest_path[:-1], shortest_path[1:]):
                if self.graph.edges[node1, node2]['distance'] > self._graph.edges[node1, node2]['distance'] * 0.1:
                    self.graph.edges[node1, node2]['distance'] *= 0.8
                if (node1, node2) not in shortest_paths_edges:
                    shortest_paths_edges.append((node1, node2))
                if node1 not in shortest_paths_nodes:
                    shortest_paths_nodes.append(node1)
                if node2 not in shortest_paths_nodes:
                    shortest_paths_nodes.append(node2)
        path_connected_edges_nodes = []
        not_path_connected_edges = []
        for edge in self.graph.edges():
            if edge not in shortest_paths_edges and (edge[1], edge[0]) not in shortest_paths_edges and (edge[0] in shortest_paths_nodes or edge[1] in shortest_paths_nodes):
                self.graph.edges()[edge]['distance'] *= 0.85
                path_connected_edges_nodes.extend(edge)
            else:
                not_path_connected_edges.append(edge)
        
        for edge in not_path_connected_edges:
            if edge[0] in path_connected_edges_nodes or edge[1] in path_connected_edges_nodes:
                self.graph.edges()[edge]['distance'] *= 0.90
     
    def volmdlr_primitives(self):
        volumes = []
        for route in self.routes:
            point1, point2 = route[0], route[1]
            if point1 != point2:
                cylinder = vmp3d.Cylinder((point1+point2)/2, point2-point1, 0.05,
                                          (point2-point1).norm())
                volumes.append(cylinder)
        return volumes



