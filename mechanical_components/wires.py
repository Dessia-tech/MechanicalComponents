#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:11:29 2018

"""

from typing import List, Tuple
from dessia_common.core import DessiaObject
from scipy.optimize import bisect, minimize
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d as primitives3d
import matplotlib.pyplot as plt
import math
import networkx as nx


import volmdlr.faces as vmf

class Wire(DessiaObject):
    """
    :param waypoints: a list of volmdlr.Point3D waypoints
    """
    _standalone_in_db = True
    
    def __init__(self, waypoints: List[vm.Point3D], diameter: float,
                 color: Tuple[float, float, float] = None, name: str = ''):
        self.waypoints = waypoints
        self.diameter = diameter
        self.color = color
        self.name = name
                
        self._utd_path = False
        
    def _path(self):
        radii = {}
        for i in range(len(self.waypoints)-2):
            # if lines are not colinear
            if vme.Line2D(self.waypoints[i],
                          self.waypoints[i+1]).unit_direction_vector()\
                    .dot(vme.Line2D(self.waypoints[i+1],
                                    self.waypoints[i+2]).unit_direction_vector())!=1:
                radii[i+1] = 4*self.diameter
        return primitives3d.OpenRoundedLineSegments3D(
            self.waypoints, radii, adapt_radius = True)
#        return  primitives3d.RoundedLineSegments3D(self.waypoints, {}, adapt_radius = True)        
    
    # def _bspline_path(self):
    #     for p1, p2 in 


    def _get_path(self):
        if not self._utd_path:
            self._path = self._path()
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
    
    def Draw(self, ax=None):
        x = []
        y = []
        for waypoint in self.waypoints:
            x.append(waypoint[0])
            y.append(waypoint[1])
        ax.plot(x, y, '-k')
        ax.plot([x[0], x[-1]], [y[0], y[-1]], 'ok')
    
    def volmdlr_primitives(self):
        section = vmw.Circle2D(vm.O2D, 0.5 * self.diameter)
        if self.color is not None:
            return [primitives3d.Sweep(section, self.path,
                                       color=self.color, name=self.name)]
        else:
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

class JunctionWire(Wire):
    def __init__(self,
                 point1: vm.Point3D, tangeancy1:vm.Vector3D,
                 point2: vm.Point3D, tangeancy2:vm.Vector3D,
                 targeted_length:float, diameter:float,
                 minimum_radius_curv:float = 0, name:str=''):
        # diameter = 2 * math.sqrt(section/math.pi)
        self.tangeancy1 = tangeancy1
        self.tangeancy2 = tangeancy2
        self.targeted_length = targeted_length
        
        self.minimum_radius_curv = minimum_radius_curv
        Wire.__init__(self, [point1, point2], diameter, name)

    # def _create_path_from_force(self, force):
    #     point1 = self.waypoints[0]
    #     point2 = self.waypoints[1]
        
    #     point1_t = point1 + force * self.tangeancy1
    #     point2_t = point2 + force * self.tangeancy2
        
    #     points = [point1, point1_t, point2_t, point2]
        
    #     bezier_curve = vme.BezierCurve3D(degree=3,
    #                                       control_points=points,
    #                                       name='bezier curve 1')
    #     return vmw.Wire3D([bezier_curve])

    # def _path(self):
        
    #     def find_force(force):
    #         bezier_curve = self._create_path_from_force(force)
            
    #         return bezier_curve.length() - self.targeted_length
        
    #     # print(self.diameter, self.length())
    #     res = bisect(find_force, self.diameter, self.targeted_length)
    #     # print(res)
    #     bezier_curve = self._create_path_from_force(res)
    #     l = bezier_curve.length()
    #     n = 20
    #     points = [bezier_curve.point_at_abscissa(i * l / n) for i in range (n+1)]
    #     return primitives3d.OpenRoundedLineSegments3D(
    #         points, {}, adapt_radius=True)
    
    def _create_path_from_forces(self, force1, force2):
        point1 = self.waypoints[0]
        point2 = self.waypoints[1]
        
        point1_t = point1 + force1 * self.tangeancy1
        point2_t = point2 + force2 * self.tangeancy2
        
        points = [point1, point1_t, point2_t, point2]
        
        bezier_curve = vme.BezierCurve3D(degree=3,
                                          control_points=points,
                                          name='bezier curve 1')
        return vmw.Wire3D([bezier_curve])
    
    # def _create_path_from_forcesv2(self, force1, force2):
    #     point1 = self.waypoints[0]
    #     point2 = self.waypoints[1]
        
    #     point1_t = point1 + force1 * self.tangeancy1
    #     point2_t = point2 + force2 * self.tangeancy2
        
    #     middle_point = (point1_t + point2_t)/2 + force2*self.tangeancy2/100 + \
    #         force1 * self.tangeancy1/100
        
    #     points = [point1, point1_t, middle_point, point2_t, point2]
        
    #     bezier_curve = vme.BezierCurve3D(degree=3,
    #                                      control_points=points,
    #                                      name='bezier curve 1')
    #     return vmw.Wire3D([bezier_curve])
    
    def _path(self):
        
        def find_forces(forces):
            # if self.minimum_radius_curv == 0:
            #     bezier_curve = self._create_path_from_forces(abs(forces[0]), forces[1])
            
            #     return abs(bezier_curve.length() - self.targeted_length)
            
            # else :
            #     bezier_curve = self._create_path_from_forces(abs(forces[0]), forces[1])
                
            #     points = bezier_curve.primitives[0].points
            #     start = points[0]
                
            #     radius = []
            #     for middle, end in zip(points[1:-1],points[2:]):
            #         circle = vmw.Circle3D.from_3_points(start, middle, end)
            #         radius.append(circle.radius)
            #         start = middle
                  
            #     return abs(bezier_curve.length() - self.targeted_length) + \
            #         100*(self.minimum_radius_curv - min(radius))
            bezier_curve = self._create_path_from_forces(abs(forces[0]), forces[1])
            
            return abs(bezier_curve.length() - self.targeted_length)
                
        
        # print(self.diameter, self.length())
        # res = bisect(find_force, self.diameter, self.targeted_length)
        res = minimize(find_forces, (self.targeted_length, self.targeted_length),
                            options={'eps': 1e-6})
        # print()
        # print(res.x, self.targeted_length)
        # print(res.fun)
        
        # print(res)
        bezier_curve = self._create_path_from_forces(abs(res.x[0]), res.x[1])
        l = bezier_curve.length()
        n = 20
        points = [bezier_curve.point_at_abscissa(i * l / n) for i in range (n+1)]
        return primitives3d.OpenRoundedLineSegments3D(
            points, {}, adapt_radius=True)
    
    def minimum_curvature_radius(self):
        points = self.path.points
        start = points[0]
        
        radius = []
        for middle, end in zip(points[1:-1],points[2:]):
            circle = vmw.Circle3D.from_3_points(start, middle, end)
            radius.append(circle.radius)
            start = middle         
            
        # print('>>>>>>>',radius)
        
        # ax = self.path.plot()
        # points[0].plot(ax=ax, color='g')
        # points[-1].plot(ax=ax, color='m')
        # for pt, rad in zip(points[1:-1], radius) :
        #     if rad < self.minimum_radius_curv :
        #         pt.plot(ax=ax, color='r')
        #     else :
        #         pt.plot(ax=ax)
        return min(radius)
    
    def to_bezier_curv(self):
        
        points = self.path.points
        start, second = points[0], points[1]
         
        ctl_points = [start, second]
        
        # ax = self.path.plot()
        
        for third in points[2:]:
            ax = self.path.plot()
            start.plot(ax=ax, color='r')
            second.plot(ax=ax, color='b')
            # third.plot(ax=ax, color='g')
            
            circle = vmw.Circle3D.from_3_points(start, second, third)
            
            # circle.plot(ax=ax)
            if circle.radius < self.minimum_radius_curv :
                s_to_c = circle.center - second
                s_to_c.normalize()
                new_center = second + s_to_c*self.minimum_radius_curv
                new_center.plot(ax=ax, color='m')
                
                # line_to_cut = vme.LineSegment3D(new_center, third)
                # line_to_cut.plot(ax=ax, color='g')
                new_circle = vmw.Circle3D(vm.Frame3D(new_center, circle.frame.u, 
                                                     circle.frame.v,circle.frame.w), 
                                          self.minimum_radius_curv)
                new_circle.plot(ax=ax, color='m')
                
                new_circle_2d = new_circle.to_2d(new_circle.frame.origin,
                                                 new_circle.frame.u,
                                                 new_circle.frame.v)
                
                third2d = third.to_2d(new_circle.frame.origin,
                                      new_circle.frame.u,
                                      new_circle.frame.v)
                
                c_to_t = third2d - new_circle_2d.center
                
                # line_to_cut2d = line_to_cut.to_2d(new_circle.frame.origin,
                #                                   new_circle.frame.u,
                #                                   new_circle.frame.v)
                line_to_cut2d = vme.LineSegment2D(new_circle_2d.center, third2d + c_to_t)
                
                ax2d = new_circle_2d.plot(color='r')
                line_to_cut2d.plot(ax=ax2d, color='g')
                
                tessel_pts = new_circle_2d.tessellation_points()
                
                #initialization
                start2d = start.to_2d(new_circle.frame.origin,
                                      new_circle.frame.u,
                                      new_circle.frame.v)
                
                second2d = second.to_2d(new_circle.frame.origin,
                                        new_circle.frame.u,
                                        new_circle.frame.v)
                
                start2d.plot(ax=ax2d, color='r')
                second2d.plot(ax=ax2d, color='b')
                
                dist_min = None
                for pt1, pt2 in zip(tessel_pts, tessel_pts[1:]+[tessel_pts[0]]):
                    line_circle = vme.LineSegment2D(pt1, pt2)
                    inter = line_circle.linesegment_intersections(line_to_cut2d)
                    if inter :
                        new_dist = second2d.point_distance(inter[0])
                        if dist_min is None or new_dist < dist_min:
                            dist_min = new_dist
                            soluce = inter[0]
                            new_third = inter[0].to_3d(new_circle.frame.origin,
                                                       new_circle.frame.u,
                                                       new_circle.frame.v)
                
                # inter = vme.Line2D(tessel_pts[0], tessel_pts[1]).line_intersections(line_to_cut2d)
                # dist_min = second2d.point_distance(inter[0])
                # soluce = inter[0]
                
                # new_third = inter[0].to_3d(new_circle.frame.origin,
                #                            new_circle.frame.u,
                #                            new_circle.frame.v)
                
                # for pt1, pt2 in zip(tessel_pts[1:], tessel_pts[2:]+[tessel_pts[0]]):
                #     line_circle = vme.Line2D(pt1, pt2)
                #     inter = line_circle.line_intersections(line_to_cut2d)
                #     new_dist = second2d.point_distance(inter[0])
                #     if new_dist < dist_min :
                #         dist_min = new_dist
                #         soluce = inter[0]
                #         new_third = inter[0].to_3d(new_circle.frame.origin,
                #                                    new_circle.frame.u,
                #                                    new_circle.frame.v)
                # new_third.plot(ax=ax2d)
                ctl_points.append(new_third)
                soluce.plot(ax=ax2d, color='g')
                
            else :
                ctl_points.append(third)
             
                
            ctl_points[-1].plot(ax=ax, color='g')
            
            
            start = second
            second = ctl_points[-1]
            
        # bezier_curve = vme.BezierCurve3D(degree=3, control_points=ctl_points,
        #                                   name='bezier curvature')
        # l = bezier_curve.length()
        # print('>>>>', l)
        # n = 20
        # points = [bezier_curve.point_at_abscissa(i * l / n) for i in range (n+1)]
        return primitives3d.OpenRoundedLineSegments3D(
            ctl_points, {}, adapt_radius=True)
            
            
                
        
                    

class WireHarness(DessiaObject):
    _standalone_in_db = True
    
    def __init__(self, wires: List[Wire], name: str = ''):
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


class RoutingSpec(DessiaObject):
    _standalone_in_db = True

    def __init__(self, source: vm.Point3D, destination: vm.Point3D,
                 diameter: float, color: Tuple[float, float, float] = None,
                 name: str = ''):
        self.source = source
        self.destination = destination
        self.diameter = diameter
        self.color = color
        self.name = name


class Wiring(DessiaObject):
    """
    Defines a combination of single wires and wire harnesses.
    
    """
    _standalone_in_db = True
    _non_serializable_attributes = ['wires_from_waypoints', 'wires']

    def __init__(self, single_wires: List[Wire],
                 wire_harnesses: List[WireHarness],
                 name: str = ''):
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
            length += wire.length(estimate=estimate)
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
                        bv = u2 - u1  # bissector vector towards inner of corner
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
        Gr = self.Graph()  # This needs to be a copy of the graph!
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
                if key not in wires:
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

    def spaced_wires(self):
        spaced_wires = {}
        G, common_routes = self.CommonRoutes()
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_edges(G, pos)

    def volmdlr_primitives(self):
        wire_volumes = []
        for wire in self.wires:
            wire_volumes.extend(wire.volmdlr_primitives())
        return wire_volumes
