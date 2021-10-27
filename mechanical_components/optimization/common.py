#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""

Common stuff between mc optimization algorithms

"""
import dessia_common as dc
import networkx as nx
import matplotlib.pyplot as plt
import volmdlr as vm
import volmdlr.edges as vme
import copy



class RoutingOptimizer(dc.DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['graph', 'line_graph']
    
    def __init__(self, routes):
        # self.waypoints = waypoints
        self.routes = routes
        
        # Setting Cache
        self._shortest_paths_cache = {}
        
        # Creating graph
        self.graph = nx.Graph()
        self._graph = nx.Graph()
        # self.graph.add_nodes_from(waypoints)
        # waypoints = set()
        # for waypoint1, waypoint2 in routes:
        #     waypoints.add(waypoint1)
        #     waypoints.add(waypoint2)
        
        edge_and_distance = []
        for waypoint1, waypoint2 in routes:
            # self.graph.add_nodes_from([waypoint1, waypoint2])
            # self.graph.add_edge(waypoint1, waypoint2, distance=(waypoint2-waypoint1).norm())
            edge_and_distance.append((waypoint1, waypoint2, (waypoint2-waypoint1).norm()))
        self.graph.add_weighted_edges_from(edge_and_distance, weight='distance')
        self._graph.add_weighted_edges_from(edge_and_distance, weight='distance')
            
        self.line_graph = None
        
    def set_line_graph(self):
        self.line_graph = nx.line_graph(self.graph)

    def plot_graph(self):
        pos = nx.kamada_kawai_layout(self.graph)
        nx.draw_networkx_nodes(self.graph, pos)#, node_color=str, node_size=int, nodelist=list)0
        nx.draw_networkx_edges(self.graph, pos)

    def plot(self):
        ax = list(self.graph.edges)[0][0].plot()
        for start, end in list(self.graph.edges):
            start.plot(ax=ax)
            end.plot(ax=ax)
            line = vm.edges.LineSegment3D(start, end)
            line.plot(ax=ax, color = 'g')


    def Draw(self, x=vm.X3D, y=vm.Y3D):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for wpt1, wpt2 in self.routes:
            vme.LineSegment3D(wpt1, wpt2).plot2d(x, y, ax=ax)
            
    def PathLength(self, path):
        length = 0.
        for waypoint1, waypoint2 in zip(path[:-1], path[1:]):
            length += self.graph[waypoint1][waypoint2]['distance']
        return length
        
    def restart_graph(self):
        self.graph = copy.deepcopy(self._graph)
