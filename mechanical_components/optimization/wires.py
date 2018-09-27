#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:21:16 2018

"""

import networkx as nx
import mechanical_components.wires as wires

class WiringOptimizer:
    def __init__(self, waypoints, routes, wires_specs):
        self.waypoints = waypoints
        self.routes = routes
        self.wires_specs = wires_specs
        
        # Creating graph
        self.graph = nx.Graph()
        self.graph.add_nodes_from(waypoints)
#        self.graph.add_edges_from(routes)
        for waypoint1, waypoint2 in routes:
            self.graph.add_edge(waypoint1, waypoint2, distance=(waypoint2-waypoint1).Norm())
        self.line_graph = nx.line_graph(self.graph)
#        nx.draw_kamada_kawai(self.graph)
        
    def PathLength(self, path):
        length = 0.
        for waypoint1, waypoint2 in zip(path[:-1], path[1:]):
            length += self.graph[waypoint1][waypoint2]['distance']
        return length
    
    def NumberHarnesses(self, paths):
        G = nx.Graph()
        G.add_nodes_from(self.line_graph)
        for path in paths:
            routes = [(w1, w2) for w1, w2 in zip(path[:-1], path[1:])]
            for route1, route2 in zip(routes[:-1], routes[1:]):
                G.add_edge(route1, route2)
        
        # Removing nodes with degree == 0
#        nodes = G.nodes().copy()
        for node in self.line_graph.nodes():
            if G.degree(node) == 0:
                G.remove_node(node)
#        nx.draw_kamada_kawai(G)
        return nx.number_connected_components(G)
        
    
    def Optimize(self, n_harnesses_max = 1):
        shortest_paths = []
        shortest_paths_lengths = []
        for wire_spec in self.wires_specs:
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
                harness_wires.append(wires.Wire(path, self.wires_specs[ipath]['diameter']))
            return wires.Wiring([], [wires.WireHarness(harness_wires)])
        else:
            alternatives_paths = []
            alternatives_paths_length = []
            # trying to Diverting cables in order to create harness
            for i_wire_spec, wire_spec in enumerate(self.wires_specs):
                cutoff = 3* len(shortest_paths[i_wire_spec])
                wire_spec = self.wires_specs[i_wire_spec]
                source = wire_spec['source']
                destination = wire_spec['destination']
#                alternatives_paths_wire = list(nx.all_simple_paths(self.graph, source, destination, cutoff = cutoff))
#                alternatives_paths.append(alternatives_paths_wire)
#                alternatives_paths_length.append(self.PathLength(alternatives_paths_wire))    
            
            
            wires2 = []
            for ipath, path in enumerate(shortest_paths):
                wires2.append(wires.Wire(path, self.wires_specs[ipath]['diameter']))
            return wires.Wiring(wires2, [])
            
            
            
        