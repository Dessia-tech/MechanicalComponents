#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:11:29 2018

"""

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import matplotlib.pyplot as plt
import math
import networkx as nx

class Wire:
    """
    :param waypoints: a list of volmdlr.Point3D waypoints
    """
    def __init__(self, waypoints, diameter):
        self.waypoints = waypoints
        self.diameter = diameter
        
        radii = {}
        for i in range(len(self.waypoints)-2):
            # if lines are not colinear
            if vm.Line2D(waypoints[i],waypoints[i+1]).DirectionVector(unit=True).Dot(vm.Line2D(waypoints[i+1], waypoints[i+2]).DirectionVector(unit=True))!=1:                
                radii[i+1] = 4*self.diameter
        self.path = primitives3D.RoundedLineSegments3D(self.waypoints, radii, adapt_radius = True)        
        
    def Length(self):
        return self.path.Length()
    
    def Draw(self, ax):
        x = []
        y = []
        for waypoint in self.waypoints:
            x.append(waypoint[0])
            y.append(waypoint[1])
        ax.plot(x, y, '-k')
        ax.plot([x[0], x[-1]], [y[0], y[-1]], 'ok')
    
    def CADVolume(self):
        first_dir = self.waypoints[1] - self.waypoints[0]
        section = vm.Contour3D([vm.Circle3D(self.waypoints[0], 0.5 * self.diameter, first_dir)])
        return primitives3D.Sweep(section, self.path)

class AWGWire(Wire):
    def __init__(self, waypoints, n):
        diameter = 0.001 * math.exp(2.1104 - 0.11594*n)
        self.n = n
        Wire.__init__(self, waypoints, diameter)
        
iec_sections = [0.5e-6,	0.75e-6, 1e-6, 1.5e-6, 2.5e-6, 4e-6, 6e-6, 10e-6, 16e-6,
                25e-6, 35e-6, 50e-6, 70e-6, 95e-6, 120e-6, 150e-6, 185e-6,
                240e-6, 300e-6, 400e-6, 500e-6, 630e-6, 800e-6, 1000e-6,
                1200e-6, 1400e-6, 1600e-6, 1800e-6, 	2000e-6, 2500e-6]

class IECWire(Wire):
    def __init__(self, waypoints, section):
        self.section = section
        diameter = 2 * math.sqrt(section/math.pi)
        Wire.__init__(self, waypoints, diameter)

class WireHarness:
    def __init__(self, wires):
        self.wires = wires
        
    def Length(self):
        length = 0.
        for wire in self.wires:
            length += wire.Length()
        return length
        
    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        
        for wire in self.wires:
            wire.Draw(ax)
        
        return ax
        
    def CADVolumes(self):
        volumes = []
        for wire in self.wires:
            volumes.append(wire.CADVolume())
        return volumes

class Wiring:
    def __init__(self, wires, wire_harnesses):
        self.wires = wires
        self.wire_harnesses = wire_harnesses
        self.wires_from_waypoints = self.WiresFromWaypoints()

        
#    def __setitem__(self, key, item):
#        key = frozenset(key)
#        self.wires_from_waypoints[key] = item

    def __getitem__(self, key):
        key = frozenset(key)
        if key in self.wires_from_waypoints:
            return self.wires_from_waypoints[key]
        else:
            return []
        
    def Length(self):
        length = 0.
        for wire in self.wires:
            length += wire.Length()
        for harness in self.wire_harnesses:
            length += harness.Length()
        return length
        
#    def Draw(self, ax=None):
#        if ax is None:
#            fig = plt.figure()
#            ax = fig.add_subplot(111)
#        else:
#            fig = None
#        for wire in self.wires:
#            wire.Draw(ax)
#        for harness in self.wire_harnesses:
#            harness.Draw(ax)
#        return fig, ax
        
    
    
    def Draw(self, x3D=vm.x3D, y3D=vm.y3D, ax=None):
        wire_sep = 0.008
#        lines = []
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = None
        
        G = self.Graph()
        Gr, wires_from_waypoints = self.CommonRoutes()
        for w1, w2 in Gr.edges():
            # Getting wires in the section
            wires = wires_from_waypoints[frozenset((w1, w2))]
            nwires = len(wires)
            # getting intermediate waypoints
            waypoints = nx.shortest_path(G, source = w1, target = w2)
            for waypoint1, waypoint2 in zip(waypoints[:-1], waypoints[1:]):
                l3D = vm.LineSegment3D(waypoint1, waypoint2)
                l2D = l3D.PlaneProjection2D(x3D, y3D)
                v2D = l2D.NormalVector()
                v2D.Normalize()
                for iwire, wire in enumerate(wires):
                    delta_wire = (iwire - 0.5 * (nwires-1)) * wire_sep * v2D
#                    print(nwires)
                    lwire = l2D.Translation(delta_wire, True)
#                    lines.append(l2D)     
                    lwire.MPLPlot(ax)
                
        
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
                    
        for wire_harness in self.wire_harnesses:
            for wire in wire_harness.wires:
                for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
                    key = frozenset((waypoint1, waypoint2))
                    if not key in wires:
                        wires[key] = [wire]
                    else:
                        wires[key].append(wire)

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

        for wire_harness in self.wire_harnesses:
            for wire in wire_harness.wires:
                for waypoint1, waypoint2 in zip(wire.waypoints[:-1], wire.waypoints[1:]):
                    G.add_edge(waypoint1, waypoint2)
#                    G.edges[waypoint1, waypoint2]['wires'].append(wire)
                
#        nx.draw_kamada_kawai(G)
        return G
        
    def CADExport(self, name='An_unnamed_wiring',
                  python_path='python',
                  path_lib_freecad='/usr/lib/freecad/lib/',
                  export_types=['fcstd']):
        groups = []
        wire_volumes = []
        for wire in self.wires:
            wire_volumes.append(wire.CADVolume())
        groups.append(('Wires', wire_volumes))
        
#        harnesses_volumes = []
        for harness in self.wire_harnesses:
            groups.append(('harness', harness.CADVolumes()))
            
        m = vm.VolumeModel(groups)
        m.FreeCADExport(name, python_path, path_lib_freecad, export_types)