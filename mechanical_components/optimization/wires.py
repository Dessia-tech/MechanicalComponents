#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:21:16 2018

"""

import networkx as nx
import mechanical_components.wires as wires
#import matplotlib.pyplot as plt
from .common import RoutingOptimizer

class WiringOptimizer(RoutingOptimizer):
    def __init__(self, waypoints, routes):
        RoutingOptimizer.__init__(self, waypoints, routes)
    
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
            
#            print(wire_spec['source'] in self.waypoints)
#            print(wire_spec['destination'] in self.waypoints)
            shortest_path = nx.shortest_path(self.graph,
                                             wire_spec['source'],
                                             wire_spec['destination'],
                                             weight = 'distance')
            shortest_paths.append(shortest_path)
            length = self.PathLength(shortest_path)
            shortest_paths_lengths.append(length)
            
        if self.NumberHarnesses(shortest_paths) == 1:
            harness_wires = []
            for ipath, path in enumerate(shortest_paths):
                if 'name' in wires_specs[ipath]:
                    name = wires_specs[ipath]['name']
                else:
                    name = ''
                harness_wires.append(wires.Wire(path, wires_specs[ipath]['diameter'], name=name))
            return wires.Wiring([], [wires.WireHarness(harness_wires)])
        else:
            wires2 = []
            for ipath, path in enumerate(shortest_paths):
                if 'name' in wires_specs[ipath]:
#                    print()
                    name = wires_specs[ipath]['name']
                else:
                    name = ''
                wires2.append(wires.Wire(path, wires_specs[ipath]['diameter'], name=name))
            return wires.Wiring(wires2, [])
            
            
            
        