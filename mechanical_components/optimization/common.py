#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Common stuff between mc optimization algorithms

"""
import networkx as nx
import matplotlib.pyplot as plt
import volmdlr as vm


class RoutingOptimizer:
    def __init__(self, waypoints, routes):
        self.waypoints = waypoints
        self.routes = routes
        
        # Creating graph
        self.graph = nx.Graph()
        self.graph.add_nodes_from(waypoints)
        for waypoint1, waypoint2 in routes:
            self.graph.add_edge(waypoint1, waypoint2, distance=(waypoint2-waypoint1).Norm())
        self.line_graph = nx.line_graph(self.graph)
        
    def DrawNetwork(self, x=vm.x3D, y=vm.y3D):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for wpt1, wpt2 in self.routes:
            vm.LineSegment3D(wpt1, wpt2).MPLPlot2D(x, y, ax)
            
    def PathLength(self, path):
        length = 0.
        for waypoint1, waypoint2 in zip(path[:-1], path[1:]):
            length += self.graph[waypoint1][waypoint2]['distance']
        return length
    