#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:20:08 2019

@author: ringhausen
"""


import volmdlr as vm
import volmdlr.primitives2D
import volmdlr.primitives3D as vm3d
import math
import numpy as npy
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy import interpolate
import json
import pkg_resources
from itertools import chain
from itertools import permutations


###############################################################################

def AngleClockwise(u, v):
    v0 = vm.Vector2D((0,0))
    if u == v0 or v == v0:
        return 0
    else:
        dot = u.Dot(v)
        norm_u = u.Norm()
        norm_v = v.Norm()
        cross = u.Cross(v)
        inner_angle = math.acos(dot/(norm_u*norm_v))
#        if cross < 0:
#            return 2*math.pi-inner_angle
#        else:
#            return inner_angle
        if cross < 0:
            return inner_angle
        else:
            return 2*math.pi-inner_angle
        
        
def SommeAngle(points):
        somme = 0
        p = points.copy()
        if p[0] != p[-1]:
            p.append(points[0])
            p.append(points[1])
        else: 
            p.append(points[1])
        nb = 0
#        print('p:', p)
        for i in range(len(p)-2):
            v1 = vm.Vector2D((p[i][0]-p[i+1][0], p[i][1]-p[i+1][1]))
            v2 = vm.Vector2D((p[i+2][0]-p[i+1][0], p[i+2][1]-p[i+1][1]))
#            print('v1', v1)
#            print('v2', v2)
            
        
            if AngleClockwise(v2, v1) < math.pi:
                somme += AngleClockwise(v2, v1)
#                print(AngleClockwise(v2, v1)*180/math.pi)
                nb += 1
            if AngleClockwise(v2, v1) >= math.pi:
#                somme += -AngleClockwise(v2, v1) + math.pi
                somme += AngleClockwise(v2, v1) - 2*math.pi
#                print((AngleClockwise(v2, v1) - 2*math.pi)*180/math.pi)
                nb -= 1
#            print('somme', somme*180/math.pi)
        somme_totale = (nb-2)*180
#        print(somme*180/math.pi, somme_totale)
        return somme*180/math.pi, somme_totale
    

    
def OffsetHatch(points, offset):
    n = len(points)
    vmpoints = []
    for point in points:
        vmpoints.append(vm.Point2D(point))
    line = vm.primitives2D.RoundedLineSegments2D(vmpoints, {})
    offset_line = line.Offset(-offset)
    m = len(offset_line.points)
    xy = npy.zeros((n+m,2))
    for i in range(n):
        xy[i][0] = line.points[i][0]
        xy[i][1] = line.points[i][1]
    for i in range(m):
        xy[-i-1][0] = offset_line.points[i][0]
        xy[-i-1][1] = offset_line.points[i][1]
    polygone = patches.Polygon(xy, hatch='/')
    return polygone


def OffsetSurface(surface, offset):
    line = vm.primitives2D.RoundedLineSegments2D(surface, {})
    offset_line = line.Offset(-offset)
    return offset_line.points

def PrintLines(lines):
    print('{')
    for i in range(len(lines)):
        print('(')
        print(i,':')
        print(lines[i][0])
        print(lines[i][1].points)
        print(lines[i][2])
        print(lines[i][3])
        print(')')
    print('}')
#points = []
##pts = [(0, 0), (1, 0), (1, 1), (0, 2), (0, 1), (-0.75, 0.314)]
##pts = [(0,0), (1,0), (0.5,0.5)]
##pts = [(0,0), (1,0), (0,1), (1,1)]
#pts = [[0.05, 0.02], [0.02, 0.02], [0.,   0.03], [0.,    0.035], [0.03,  0.035], [0.03, 0.04], [0.08, 0.04], [0.08,  0.035], [0.085, 0.035], [0.085, 0.04 ], [0.086, 0.04 ], [0.086, 0.034], [0.05, 0.02]]
##pts  = [[0.05, 0.02], [0.02, 0.02], [0.,   0.03], [0.,    0.035], [0.086, 0.035], [0.086, 0.034], [0.05, 0.02]]
#for i in pts:
#    points.append(vm.Point2D(i))
#SommeAngle(points)
#points.reverse()
#SommeAngle(points)

###############################################################################

class ShaftMaterial:
    def __init__(self, Re, t_yield, C, name='',):
        self.Re = Re            # Ultimate tensile strength (Pa)
        self.t_yield = t_yield  # Yield stength (Pa)
        self.C = C              # Hardness (Brinell Scale) (Bhn)
        self.name = name

###############################################################################

class Accessory:
    def __init__(self, montage='axial'):#, functional_surfaces):
        """
        The montage should be either radial, axial, inbuilt or flat
        """
        self.montage = montage

###############################################################################

class Surface:
    def __init__(self, points, part_left, part_right):
        self.points = points            # A list of coordinates
        self.part_left = part_left      # Is an Accessory
        self.part_right = part_right
    
    @classmethod
    def Instantiate(cls, coordinate, part_left, part_right):
        points = [vm.Point2D(p) for p in coordinate]
        return cls(points, part_left, part_right)
        # functional_surface = [(point1, point2,...), part_left, part_right]
        # left or right according to the oriented path between the points

###############################################################################

class FunctionalSurfaces:
    def __init__(self, surfaces):
        self.surfaces = surfaces    # A list of Surface.Instantiate objects
        self.offset = 0.001
        
    def RightPartSurface(self, part):
        part_functional_surfaces = []
        for surface in self.surfaces:
            if surface.part_right == part:
                part_functional_surfaces.append(surface)
            elif surface.part_left == part:
                right_surface = Surface.Instantiate(surface.points[::-1], 
                                                    surface.part_right, 
                                                    surface.part_left)
                part_functional_surfaces.append(right_surface)
        return part_functional_surfaces
        # The part is on the right-hand side when the points are read
            
    def PartContour(self, part):
        part_functional_surfaces = self.RightPartSurface(part)
        
        points = []
        inf_points = []
        sup_points = []
        left_points = []
        right_points = []
        for surface in part_functional_surfaces:
            first_point = surface.points[0]
            last_point = surface.points[-1]
            if first_point[0] > last_point[0]:
                inf_points.append(surface.points)
            elif first_point[0] < last_point[0]:
                sup_points.append(surface.points)
            else:
                if first_point[1] < last_point[1]:
                    left_points.append(surface.points)
                else:
                    right_points.append(surface.points)
            
        inf_points.sort(key=lambda p: p[0])
        inf_points.reverse()
        sup_points.sort(key=lambda p: p[0])
        left_points.sort(key=lambda p: p[1])
        right_points.sort(key=lambda p: p[1])
        right_points.reverse()
        
        points.append(inf_points)
        points.append(left_points)
        points.append(sup_points)
        points.append(right_points)
        
        points = list(chain.from_iterable(points))
        points = list(chain.from_iterable(points))
        
        if len(part_functional_surfaces) > 1:
            points.append(points[0])
        
        polygone_hatches = []
        for surface in part_functional_surfaces:
            polygone_hatches.append(OffsetHatch(surface.points, self.offset))
        
        line = vm.primitives2D.RoundedLineSegments2D(points, {})
        contour = vm.Contour2D([line])
        return contour, polygone_hatches
    
    def PartPlot(self, part, axe=None):
        part_contour, polygone_hatches = self.PartContour(part)
        f, a = part_contour.MPLPlot(ax=axe, style='')
        
        for polygone in polygone_hatches:
            if polygone_hatches != None:
                a.add_patch(polygone)
        return a
    
    def Parts(self):
        parts = []
        for surface in self.surfaces:
            parts.append(surface.part_left)
            parts.append(surface.part_right)
        single_parts = list(dict.fromkeys(parts))
        return single_parts
    
    def Plot(self):
        parts = self.Parts()
        axe = self.PartPlot(parts[0])
        for part in parts:
            self.PartPlot(part, axe)

    def SurfaceOrder(self, part):
        
        epsilon = 1e-6 
        
        part_functional_surfaces = self.RightPartSurface(part)
        surfaces_points = []
        for surface in part_functional_surfaces:
            surfaces_points.append(surface.points)
        perm = list(permutations(surfaces_points))
        
        # Because the problem is circular, we can remove some of the permutations
        remove_index = []
        for index, elem in enumerate(perm):
            if elem[0] != surfaces_points[0]:
                remove_index.append(index)
        for i in remove_index[::-1]:
            perm.pop(i) 
            
        # Checks the self intersection and the angles of the polygon
        lines = {}
        for i, perm_i in enumerate(perm):
            points = list(chain.from_iterable(perm_i))
            
            polygon_points = points.copy()
            if polygon_points[0] == polygon_points[-1]:
                del polygon_points[-1] 
            polygon = vm.Polygon2D(polygon_points)
            
            if points[0] != points[-1]:         
                points.append(points[0])
            line = vm.primitives2D.RoundedLineSegments2D(points, {}, 
                                                         closed=True)
            angle = SommeAngle(points)            

            lines[i] = perm_i, line, polygon.SelfIntersect(), angle[1]-angle[0]

            ax = self.PartPlot(part)
            line.MPLPlot(ax, style='r')
        
        
        
        # Repairs the polygons that are inside out
        imax = len(lines)
        for i in range(imax):
            if lines[i][3] < epsilon: #and lines[i][3] < epsilon:
                offset_points = []
                repair_points = []
                for surface in lines[i][0]:
                    line_surface = vm.primitives2D.RoundedLineSegments2D(surface, {})
                    repair_points.extend(surface)
                    offset_points.extend(line_surface.Offset(-self.offset).points)
                offset_points.reverse()
                repair_points.extend(offset_points)
                
                print('repair_points', repair_points)
                if repair_points[0] == repair_points[-1]:
                    del repair_points[-1]
                repair_polygon = vm.Polygon2D(repair_points)
                
                repair_points_line = repair_points.copy()
                repair_points_line.append(repair_points_line[0])
                repair_line = vm.primitives2D.RoundedLineSegments2D(repair_points_line, {})
#                intersect = repair_polygon.SelfIntersect()
                repair_line.MPLPlot(style='b')
                angle = SommeAngle(repair_points)
                
                lines[len(lines)] = lines[i][0], repair_line, repair_polygon.SelfIntersect(), angle[1]-angle[0]

                
        # Repair the polygons that are lightly self intersecting, meaning that
        # they can be repaired by extending the contour to some of the offset's
        # points.
        for i in range(imax):
            print('========')

            if lines[i][2][0]: #and lines[i][3] > 360-epsilon and lines[i][3] < 360+epsilon:
                intersect_line1 = lines[i][2][1]
                intersect_line2 = lines[i][2][2]
                p11, p12 = intersect_line1.points
                p21, p22 = intersect_line2.points
                print('p11', p11)
                print('p12', p12)
                print('p21', p21)
                print('p22', p22)                
                    
                for j, surface in enumerate(lines[i][0]):
                    for p1_surface, p2_surface in zip(surface[::], surface[1::]):
                        print('p1_surface', p1_surface)
                        print('p2_surface', p2_surface)
                        
                        # Looking for the functional surface that is intersected
                        if [p11, p12] == [p1_surface, p2_surface]:

                            print('==> 1 <==')
                            
                            troncated_surfaces = lines[i][0][:j+1]
                            
                            offset_points = OffsetSurface(troncated_surfaces[-1], self.offset)
                            
                            i_point = 1
                            intersect0 = True
                            while intersect0 and i_point <= len(offset_points):
                                # New_surface is composed of the points just before the intersection,
                                # of the offset points in order to pass by the intersection, and finally 
                                # of the remaining points that have been troncated.
                                new_surface = list(troncated_surfaces)
                                new_surface.append(list(offset_points[-ii] for ii in range(1,i_point+1)))
                                if len(lines[i][0][j+1:]) != 0:    
                                    new_surface.append(lines[i][0][j+1:][0])
                                    
                                polygon_points = list(chain.from_iterable(new_surface))
                                
                                if polygon_points[0] == polygon_points[-1]:
                                    del polygon_points[-1]
                                
                                rls_points = polygon_points.copy()
                                rls_points.append(rls_points[0])
                                rls = vm.primitives2D.RoundedLineSegments2D(rls_points, {})
                                
                                polygon = vm.Polygon2D(polygon_points)
                                intersect = polygon.SelfIntersect()
                                intersect0 = intersect[0]
                                
                                i_point += 1
                                
                            rls.MPLPlot(style='g')
                            angle = SommeAngle(rls_points)
                            
                            lines[len(lines)] = lines[i][0], rls, intersect, angle[1]-angle[0]

                                
                        if [p12, p11] == [p1_surface, p2_surface]:

                            print('==> 2 <==')
                            
                            troncated_surfaces = lines[i][0][:j+1]
                            
                            offset_points = OffsetSurface(troncated_surfaces[-1], self.offset)
                            
                            i_point = 1
                            intersect0 = True
                            while intersect0 and i_point <= len(offset_points):
                                # New_surface is composed of the points just before the intersection,
                                # of the offset points in order to pass by the intersection, and finally 
                                # of the remaining points that have been troncated.
                                new_surface = list(troncated_surfaces)
                                new_surface.append(list(offset_points[-ii] for ii in range(1,i_point+1)))
                                if len(lines[i][0][j+1:]) != 0:    
                                    new_surface.append(lines[i][0][j+1:][0])
                                    
                                polygon_points = list(chain.from_iterable(new_surface))
                                rls_points = polygon_points.copy()
                                rls_points.append(rls_points[0])
                                
                                rls = vm.primitives2D.RoundedLineSegments2D(rls_points, {})
                                
                                polygon = vm.Polygon2D(polygon_points)
                                intersect = polygon.SelfIntersect()
                                intersect0 = intersect[0]
                                
                                i_point += 1
                                
                            rls.MPLPlot(style='g')
                            angle = SommeAngle(rls_points)
                            
                            lines[len(lines)] = lines[i][0], rls, intersect, angle[1]-angle[0]

                            
                        if [p21, p22] == [p1_surface, p2_surface]:
                            
                            print('==> 3 <==')
                            
                            troncated_surfaces = lines[i][0][:j+1]
                            
                            offset_points = OffsetSurface(troncated_surfaces[-1], self.offset)
                            
                            i_point = 1
                            intersect0 = True
                            while intersect0 and i_point <= len(offset_points):
                                # New_surface is composed of the points just before the intersection,
                                # of the offset points in order to pass by the intersection, and finally 
                                # of the remaining points that have been troncated.
                                new_surface = list(troncated_surfaces)
                                new_surface.append(list(offset_points[-ii] for ii in range(1,i_point+1)))
                                if len(lines[i][0][j+1:]) != 0:    
                                    new_surface.append(lines[i][0][j+1:][0])
                                    
                                polygon_points = list(chain.from_iterable(new_surface))
                                rls_points = polygon_points.copy()
                                rls_points.append(rls_points[0])
                                
                                rls = vm.primitives2D.RoundedLineSegments2D(rls_points, {})
                                
                                polygon = vm.Polygon2D(polygon_points)
                                intersect = polygon.SelfIntersect()
                                intersect0 = intersect[0]
                                
                                i_point += 1
                                
                            rls.MPLPlot(style='g')
                            angle = SommeAngle(rls_points)
                            
                            lines[len(lines)] = lines[i][0], rls, intersect, angle[1]-angle[0]

                            
                        if [p22, p21] == [p1_surface, p2_surface]:
                            
                            print('==> 4 <==')
                            
                            troncated_surfaces = lines[i][0][:j+1]
                            
                            offset_points = OffsetSurface(troncated_surfaces[-1], self.offset)
                            
                            i_point = 1
                            intersect0 = True
                            while intersect0 and i_point <= len(offset_points):
                                # New_surface is composed of the points just before the intersection,
                                # of the offset points in order to pass by the intersection, and finally 
                                # of the remaining points that have been troncated.
                                new_surface = list(troncated_surfaces)
                                new_surface.append(list(offset_points[-ii] for ii in range(1,i_point+1)))
                                if len(lines[i][0][j+1:]) != 0:    
                                    new_surface.append(lines[i][0][j+1:][0])
                                    
                                polygon_points = list(chain.from_iterable(new_surface))
                                rls_points = polygon_points.copy()
                                rls_points.append(rls_points[0])
                                
                                rls = vm.primitives2D.RoundedLineSegments2D(rls_points, {})
                                
                                polygon = vm.Polygon2D(polygon_points)
                                intersect = polygon.SelfIntersect()
                                intersect0 = intersect[0]
                                
                                i_point += 1
                                
                            rls.MPLPlot(style='g')
                            angle = SommeAngle(rls_points)
                            
                            lines[len(lines)] = lines[i][0], rls, intersect, angle[1]-angle[0]
    
                            
                            
            print('############')
            
                    
        
        print()
        PrintLines(lines)
        print()
                
        return lines
    
    
    
    
    def PolygonLink(self, surface1, surface2):
        epsilon = 1e-6 
        # We suppose that the part is on the right-hand side
        # There is always 2 way of linking 2 polygones, but this method only
        # gives you 1 way of linking them, the more direct way.
        n1 = len(surface1)
        n2 = len(surface2)
        polygone1 = OffsetHatch(surface1.points, self.offset)
        polygone2 = OffsetHatch(surface2.points, self.offset)
        
        p1_1 = surface1.points[0]
        p2_1 = surface1.points[-1]
        p3_1 = polygone1.points[n1+1]
        p4_1 = polygone1.points[-1]
        p_1 = [p1_1, p2_1, p3_1, p4_1]
        
        p1_2 = surface2.points[0]
        p2_2 = surface2.points[-1]
        p3_2 = polygone2.points[n2+1]
        p4_2 = polygone2.points[-1]
        p_2 = [p1_2, p2_2, p3_2, p4_2]
        
        intersect = True
        intersect1 = []
        # p1_1 cannot be linked to p1_2 because of their direction
        for i in range(1,4):
            polygon_points = p_2
            polygon_points.insert(i+1, p1_1)
            polygon_points.insert(i+2, p_1[i])
            polygone1 = vm.Polygon2D(polygon_points) 
            intersect = polygone1.SelfIntersect()
            intersect1.append(intersect)
        
        intersect2 = []
        intersect = True
        for i in range(1,4):
            polygon_points = p_1
            polygon_points.insert(i+1, p1_2)
            polygon_points.insert(i+2, p_2[i])
            polygone2 = vm.Polygon2D(polygon_points) 
            intersect = polygone2.SelfIntersect()
            intersect1.append(intersect)
            
        possible_link = []
        for i, boo1 in enumerate(intersect1):
            if not boo1:
                line1 = vm.LineSegment2D(p1_1, p_1[i+1])
                for j, boo2 in enumerate(intersect2):
                    if not boo2:
                        line2 = vm.LineSegment2D(p1_2, p_2[j+1])
                        p, a, b = vm.Point2D.LinesIntersection(line1, line2, True)
                        if not (a > 0+epsilon and a < 1-epsilon and b > 0+epsilon and b < 1-epsilon):
                            possible_link.append((i,j))
                            
                            
                            
                            
        return 
            

###############################################################################


        
        
###############################################################################
###############################################################################    
    
class Shaft:
    def __init__(self, functional_surfaces, shaft_material=None):
        self.functional_surfaces = functional_surfaces
        self.shaft_material = shaft_material

    def ShaftDrawingPoints(self):
        
        points_inf = []
        for accessoire in self.sorted_functional_accessories:
            if accessoire.external_vector[1] > 0:
                points_inf.extend([tuple(vm.Point2D(p) for p in accessoire.functional_points)])
            else:
                if self.offset != 0:
                    points_offset = []
                    for point in accessoire.functional_points:
                        pt = vm.Point2D(point)
                        points_offset.append(pt.Translation(self.offset * vm.Vector2D(accessoire.external_vector)))
                    points_inf.append(tuple(points_offset))
                    
        points_sup = []
        for accessoire in self.sorted_functional_accessories:
            if accessoire.external_vector[1] < 0:
                points_sup.extend([tuple(vm.Point2D(p) for p in accessoire.functional_points)])
            else:
                if self.offset != 0:
                    points_offset = []
                    for point in accessoire.functional_points:
                        pt = vm.Point2D(point)
                        points_offset.append(pt.Translation(self.offset * vm.Vector2D(accessoire.external_vector)))
                    points_sup.append(tuple(points_offset))

        return points_inf, points_sup 
    
    
    
    def ShaftDrawingContour(self):
        points_inf, points_sup = self.ShaftDrawingPoints()
        print('points inf: ', points_inf, '\n')
        print('points sup: ', points_sup)
        
        points_inf = list(chain.from_iterable(points_inf))
        points_sup = list(chain.from_iterable(points_sup))
        
#        points_sup.reverse()
        points = points_inf + points_sup + [points_inf[0]]
        
        line = vm.primitives2D.RoundedLineSegments2D(points, {})
        contour = vm.Contour2D([line])
        return contour

            
    
    def Plot2(self):
        contour = self.ShaftDrawingContour()
        f, a = contour.MPLPlot()
        contour_accessories, contour_inbuilts, linkage_contours = self.functional_accessories.FunctionalContour()
        for contour_accessory in contour_accessories:
            contour_accessory.MPLPlot(ax=a, style='r')
        for contour_inbuilt in contour_inbuilts:
            contour_inbuilt.MPLPlot(ax=a, style='b')
        for contour_linkage in linkage_contours:
            contour_linkage.MPLPlot(ax=a, style='g')
            
        
    def ShaftCheck(self):
        points_inf, points_sup = self.ShaftDrawingPoints()
        
        for i, accessoire in enumerate(self.sorted_functional_accessories[1:-1], 1):
            if accessoire.montage == 'axial':

                if accessoire.external_vector[1] > 0:
                    points_inf_left = points_inf[: i]
                    points_inf_right = points_inf[i+1 :]

                    if points_inf_left != [] and points_inf_right != []:
                        min_inf_left = min([point[1] for point in list(chain.from_iterable(points_inf_left))])
                        min_inf_right = min([point[1] for point in list(chain.from_iterable(points_inf_right))])
                        r_accessoire = min([point[1] for point in accessoire.functional_points])
                        
                        if r_accessoire > min_inf_left and r_accessoire > min_inf_right:
                            return False
                    
                else:
                    points_sup_left = points_sup[: i]
                    points_sup_right = points_sup[i+1 :]

                    if points_sup_left != [] and points_sup_right != []:
                        max_sup_left = max([point[1] for point in list(chain.from_iterable(points_sup_left))])
                        max_sup_right = max([point[1] for point in list(chain.from_iterable(points_sup_right))])
                        R_accessoire = max([point[1] for point in accessoire.functional_points])
                        
                        if R_accessoire < max_sup_left and R_accessoire < max_sup_right:
                            return False
        return True
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    