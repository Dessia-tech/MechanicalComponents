#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:11:29 2018

"""

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import matplotlib.pyplot as plt

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
        self.path = primitives3D.RoundedLineSegments3D(self.waypoints, radii)        
        
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
        
    def Length(self):
        length = 0.
        for wire in self.wires:
            length += wire.Length()
        for harness in self.wire_harnesses:
            length += harness.Length()
        return length
        
    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        x = []
        y = []
        for wire in self.wires:
            wire.Draw(ax)
        for harness in self.wire_harnesses:
            harness.Draw(ax)
        return ax
        
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