#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:11:29 2018

"""

from typing import List
from dessia_common.core import DessiaObject
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d as primitives3d
import matplotlib.pyplot as plt
import math
import networkx as nx

class Wire(DessiaObject):
    """
    :param waypoints: a list of volmdlr.Point3D waypoints
    """
    _standalone_in_db = True
    
    def __init__(self, waypoints:List[vm.Point3D], diameter:float, name:str=''):
        self.waypoints = waypoints
        self.diameter = diameter
        self.name = name
                
        self._utd_path = False
        
    def _Path(self):
        radii = {}
        for i in range(len(self.waypoints)-2):
            # if lines are not colinear
            if vme.Line2D(self.waypoints[i],
                          self.waypoints[i+1]).unit_direction_vector()\
                    .dot(vme.Line2D(self.waypoints[i+1],
                                    self.waypoints[i+2]).unit_direction_vector())!=1:
                radii[i+1] = 4*self.diameter
        return  primitives3d.OpenRoundedLineSegments3D(self.waypoints, radii, adapt_radius = True)
#        return  primitives3d.RoundedLineSegments3D(self.waypoints, {}, adapt_radius = True)        
    
    def _get_path(self):
        if not self._utd_path:
            self._path = self._Path()
            self._utd_path = True
        return self._path
    
    path = property(_get_path)    
    
    def length(self, estimate = False):
        if estimate:
            length_estimate = 0.
            for wpt1, wpt2 in zip(self.waypoints[:-1], self.waypoints[1:]):
                length_estimate += wpt2.point_distance(wpt1)
            return length_estimate
        else:
            return self.path.length()
    
    def Draw(self, ax):
        x = []
        y = []
        for waypoint in self.waypoints:
            x.append(waypoint[0])
            y.append(waypoint[1])
        ax.plot(x, y, '-k')
        ax.plot([x[0], x[-1]], [y[0], y[-1]], 'ok')
    
    def volmdlr_primitives(self):
        section = vmw.Circle2D(vm.O2D, 0.5 * self.diameter)
        return [primitives3d.Sweep(section, self.path, name=self.name)]

class AWGWire(Wire):
    def __init__(self, waypoints, n, name=''):
        diameter = 0.001 * math.exp(2.1104 - 0.11594*n)
        self.n = n
        Wire.__init__(self, waypoints, diameter, name)
        
iec_sections = [0.5e-6,	0.75e-6, 1e-6, 1.5e-6, 2.5e-6, 4e-6, 6e-6, 10e-6, 16e-6,
                25e-6, 35e-6, 50e-6, 70e-6, 95e-6, 120e-6, 150e-6, 185e-6,
                240e-6, 300e-6, 400e-6, 500e-6, 630e-6, 800e-6, 1000e-6,
                1200e-6, 1400e-6, 1600e-6, 1800e-6, 	2000e-6, 2500e-6]

class IECWire(Wire):
    def __init__(self, waypoints, section, name=''):
        self.section = section
        diameter = 2 * math.sqrt(section/math.pi)
        Wire.__init__(self, waypoints, diameter, name)

class WireHarness(DessiaObject):
    _standalone_in_db = True
    
    def __init__(self, wires:List[Wire], name:str=''):
        self.wires = wires
        self.name = name
        
    def length(self):
        length = 0.
        for wire in self.wires:
            length += wire.length()
        return length
        
    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        
        for wire in self.wires:
            wire.Draw(ax)
        
        return ax

class Wiring(DessiaObject):
    """
    Defines a combination of single wires and wire harnesses.
    
    """
    _standalone_in_db = True
    def __init__(self, single_wires:List[Wire],
                 wire_harnesses:List[WireHarness],
                 name:str=''):
        self.single_wires = single_wires
        self.wire_harnesses = wire_harnesses
        self.name = name


        wires = single_wires[:]
        for harness in wire_harnesses:
            wires.extend(harness.wires)
        self.wires = wires
        
        self.wires_from_waypoints = self.WiresFromWaypoints()

    def __getitem__(self, key):
        key = frozenset(key)
        if key in self.wires_from_waypoints:
            return self.wires_from_waypoints[key]
        else:
            return []
        
    def length(self, estimate=False):
        """
        Gives the cumulative length of wires
        
        :param estimate: If set to True, compute the length without the raddi of wires
        
        """
        
        length = 0.
        for wire in self.wires:
            length += wire.length(estimate)
        return length
    
    def Draw(self, x3D=vm.X3D, y3D=vm.Y3D, ax=None):
        wire_sep = 0.005
#        lines = []
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = None
        
        G = self.Graph()
        Gr, wires_from_waypoints = self.CommonRoutes()
        wire_lines = {wire:{} for wire in self.wires}# nested dicts first key wire, second frozenset(waypoint1, waypoint2)
        for w1, w2 in Gr.edges():
   
            # Getting wires in the section
            wires = wires_from_waypoints[frozenset((w1, w2))]
            nwires = len(wires)
            # getting intermediate waypoints
            waypoints = nx.shortest_path(G, source = w1, target = w2)
            for waypoint1, waypoint2 in zip(waypoints[:-1], waypoints[1:]):
                l3D = vme.LineSegment3D(waypoint1, waypoint2)
                l2D = l3D.plane_projection2d(vm.O3D, x3D, y3D)
                if l2D.length() > 0.:                    
                    v2D = l2D.normal_vector()
                    for iwire, wire in enumerate(wires):
                        delta_wire = (iwire - 0.5 * (nwires-1)) * wire_sep * v2D
                        lwire = l2D.translation(delta_wire, True)
                        wire_lines[wire][frozenset((waypoint1, waypoint2))] = lwire
                else:
                    for iwire, wire in enumerate(wires):
                        wire_lines[wire][frozenset((waypoint1, waypoint2))] = l2D
                                        

                        
        for wire in self.wires:
            waypoint0_2D = wire.waypoints[0].plane_projection2d(vm.O3D, x3D, y3D)
            line = wire_lines[wire][frozenset((wire.waypoints[0], wire.waypoints[1]))]
            if line.points[0].point_distance(waypoint0_2D) < line.points[1].point_distance(waypoint0_2D):
                waypoints_draw = [line.points[0]]
            else:
                waypoints_draw = [line.points[1]]

            nwaypoints = len(wire.waypoints)
            for iwaypoint in range(nwaypoints-2):
                waypoint1_2D = wire.waypoints[iwaypoint].plane_projection2d(vm.O3D, x3D, y3D)
                waypoint2_2D = wire.waypoints[iwaypoint+1].plane_projection2d(vm.O3D, x3D, y3D)
                waypoint3_2D = wire.waypoints[iwaypoint+2].plane_projection2d(vm.O3D, x3D, y3D)
                line1 = vme.LineSegment2D(waypoint1_2D, waypoint2_2D)
                line2 = vme.LineSegment2D(waypoint2_2D, waypoint3_2D)
                
                line1_draw = wire_lines[wire][frozenset((wire.waypoints[iwaypoint], wire.waypoints[iwaypoint+1]))]
                line2_draw = wire_lines[wire][frozenset((wire.waypoints[iwaypoint+1], wire.waypoints[iwaypoint+2]))]
                
                if (line1.length() == 0) or (line2.length() == 0):
                    waypoints_draw.append(wire.waypoints[iwaypoint+1])
                else:
                    u1 = line1.unit_direction_vector()
                    u2 = line2.unit_direction_vector()
                    if abs(u1.dot(u2)) != 1:
                        bv = u2 - u1# bissector vector towards inner of corner
                        bl = vme.Line2D(waypoint2_2D, waypoint2_2D+bv)
                        i1 = vm.Point2D.line_intersection(bl, line1_draw)
                        i2 = vm.Point2D.line_intersection(bl, line2_draw)
                        if waypoint2_2D.point_distance(i1) < waypoint2_2D.point_distance(i2):
                            waypoints_draw.append(i2)
                        else:
                            waypoints_draw.append(i1)

                    else:
                        waypoints_draw.append(line2.points[0])

            waypointn_2D = wire.waypoints[-1].plane_projection2d(vm.O3D, x3D, y3D)
            line = wire_lines[wire][frozenset((wire.waypoints[-2], wire.waypoints[-1]))]
            if line.points[0].point_distance(waypointn_2D) < line.points[1].point_distance(waypointn_2D):
                waypoints_draw.append(line.points[0])
            else:
                waypoints_draw.append(line.points[1])
            
            x = [w[0] for w in waypoints_draw]
            y = [w[1] for w in waypoints_draw]
            ax.plot(x, y, 'o-k')

        ax.set_aspect('equal')
        return fig, ax
    
    def CommonRoutes(self):
        wires_from_waypoints = self.WiresFromWaypoints()
        # Computing reduced graph
        Gr = self.Graph()# This needs to be a copy of the graph!
        node_delete = True
        while node_delete:
            node_delete = False
            for waypoint, degree in nx.degree(Gr):
                if degree == 2:
                    # Seeing whats connected
                    waypoint1, waypoint2 = Gr[waypoint]
                    # If there is the same wires on each side
                    wires1 = wires_from_waypoints[frozenset((waypoint1, waypoint))]
                    wires2 = wires_from_waypoints[frozenset((waypoint2, waypoint))]
                    if set(wires1) == set(wires2):
                        # Contracting node from graph
                        Gr.remove_node(waypoint)
                        Gr.add_edge(waypoint1, waypoint2)
                        del wires_from_waypoints[frozenset((waypoint1, waypoint))]
                        del wires_from_waypoints[frozenset((waypoint2, waypoint))]
                        wires_from_waypoints[frozenset((waypoint1, waypoint2))] = wires1
                        node_delete = True
                        break
        
        
        return Gr, wires_from_waypoints
                    
                
    # TODO: Performance caching this and graph
    def WiresFromWaypoints(self):
        wires = {}

        for wire in self.wires:
            for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
                key = frozenset((waypoint1, waypoint2))
                if not key in wires:
                    wires[key] = [wire]
                else:
                    wires[key].append(wire)
                    
#        for wire_harness in self.wire_harnesses:
#            for wire in wire_harness.wires:
#                for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
#                    key = frozenset((waypoint1, waypoint2))
#                    if not key in wires:
#                        wires[key] = [wire]
#                    else:
#                        wires[key].append(wire)

        return wires
    
    def Graph(self):
        G = nx.Graph()
        # Adding nodes
        for wire in self.wires:
            G.add_nodes_from(wire.waypoints)
        for wire_harness in self.wire_harnesses:
            for wire in wire_harness.wires:
                G.add_nodes_from(wire.waypoints)

        # Adding edges
        for wire in self.wires:
            for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
                G.add_edge(waypoint1, waypoint2)
#                G.edges[waypoint1, waypoint2]['wires'].append(wire)

#        for wire_harness in self.wire_harnesses:
#            for wire in wire_harness.wires:
#                for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
#                    G.add_edge(waypoint1, waypoint2)
                
#        nx.draw_kamada_kawai(G)
        return G
    
    def volmdlr_primitives(self):
        wire_volumes = []
        for wire in self.wires:
            wire_volumes.extend(wire.volmdlr_primitives())
        return wire_volumes