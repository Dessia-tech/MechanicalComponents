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

    def volmdlr_primitives(self):
        volumes = []
        for route in self.routes:
            point1, point2 = route[0], route[1]
            if point1 != point2:
                cylinder = vmp3d.Cylinder((point1+point2)/2, point2-point1, 0.05,
                                          (point2-point1).norm())
                volumes.append(cylinder)
        return volumes

class WireRouting():
    def __init__(self, edges):
        self.edges = edges
    
    def weigthed_graph(self):
        G = nx.Graph()
        for edge in self.edges:
            distance = edge[0].point_distance(edge[1])
            G.add_edges_from([(edge[0], edge[1],{'distance': distance})])
            #  {'distance' : 2}
        return G


    def minimim_spannig_tree(self):
        G = self.weigthed_graph()
        mst = nx.minimum_spanning_tree(G, 'distance')
        return mst


