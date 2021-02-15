#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""
Created on Wed Jul 10 15:21:32 2019

@author: steven
"""

import networkx as nx
from .common import RoutingOptimizer
import mechanical_components.wires as wires

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
        
    
    def Route(self, wires_specs):
        shortest_paths = []
        shortest_paths_lengths = []
        for wire_spec in wires_specs:
            # Checking if pat is in cache
            if (wire_spec['source'], wire_spec['destination']) in self._shortest_paths_cache:
                shortest_path = self._shortest_paths_cache[(wire_spec['source'], wire_spec['destination'])]
            elif (wire_spec['source'], wire_spec['destination']) in self._shortest_paths_cache:
                shortest_path = self._shortest_paths_cache[(wire_spec['source'], wire_spec['destination'])]
            else:            
                shortest_path = nx.shortest_path(self.graph,
                                                 wire_spec['source'],
                                                 wire_spec['destination'],
                                                 weight = 'distance')
                self._shortest_paths_cache[(wire_spec['source'], wire_spec['destination'])] = shortest_path
                
            shortest_paths.append(shortest_path)
            length = self.PathLength(shortest_path)
            shortest_paths_lengths.append(length)
            
#        if self.NumberHarnesses(shortest_paths) == 1:
#            harness_wires = []
#            for ipath, path in enumerate(shortest_paths):
#                if 'name' in wires_specs[ipath]:
#                    name = wires_specs[ipath]['name']
#                else:
#                    name = ''
#                harness_wires.append(wires.Wire(path, wires_specs[ipath]['diameter'], name=name))
#            return wires.Wiring([], [wires.WireHarness(harness_wires)])
#        else:
        wires2 = []
        for ipath, path in enumerate(shortest_paths):
            if 'name' in wires_specs[ipath]:
#                    print()
                name = wires_specs[ipath]['name']
            else:
                name = ''
            wires2.append(wires.Wire(path, wires_specs[ipath]['diameter'], name=name))
        return wires.Wiring(wires2, [])