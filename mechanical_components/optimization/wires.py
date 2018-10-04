#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:21:16 2018

"""

import networkx as nx
import mechanical_components.wires as wires
import matplotlib.pyplot as plt

class WiringOptimizer:
    def __init__(self, waypoints, routes):
        self.waypoints = waypoints
        self.routes = routes
#        self.wires_specs = wires_specs
        
        # Creating graph
        self.graph = nx.Graph()
        self.graph.add_nodes_from(waypoints)
#        self.graph.add_edges_from(routes)
        for waypoint1, waypoint2 in routes:
            self.graph.add_edge(waypoint1, waypoint2, distance=(waypoint2-waypoint1).Norm())
        self.line_graph = nx.line_graph(self.graph)
        
    def DrawNetwork(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for wpt1, wpt2 in self.routes:
            ax.plot([wpt1[0], wpt2[0]], [wpt1[1], wpt2[1]], 'o-k')
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
        
    
    def RouteAndCreateMaximumHarnesses(self,  wires_specs):
        shortest_paths = []
        shortest_paths_lengths = []
        for wire_spec in wires_specs:
            print(wire_spec['source'] in self.waypoints)
            print(wire_spec['destination'] in self.waypoints)
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
                harness_wires.append(wires.Wire(path, wires_specs[ipath]['diameter']))
            return wires.Wiring([], [wires.WireHarness(harness_wires)])
        else:
            wires2 = []
            for ipath, path in enumerate(shortest_paths):
                wires2.append(wires.Wire(path, wires_specs[ipath]['diameter']))
            return wires.Wiring(wires2, [])
            
            
            
        