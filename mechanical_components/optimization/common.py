#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Common stuff between mc optimization algorithms

"""
import networkx as nx
import matplotlib.pyplot as plt

class RoutingOptimizer:
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
    
