#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""
Created on Mon Apr  1 17:20:08 2019

@author: ringhausen
"""

from itertools import chain
from itertools import permutations
from itertools import product
import math
import numpy as npy
import scipy as spy
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import dectree
from dessia_common.core import DessiaObject
import volmdlr as vm
from volmdlr import primitives2d
import volmdlr.primitives3d as vm3d
#from scipy import interpolate
#import json
#import pkg_resources

###############################################################################

def AngleClockwise(vector1, vector2):
    """
    Return the clockwise angle in radians between vector1 and vector2.
    """
    vector0 = vm.Vector2D((0, 0))
    if vector0 in (vector1, vector2):
        return 0

    dot = vector1.Dot(vector2)
    norm_vec_1 = vector1.Norm()
    norm_vec_2 = vector2.Norm()
    cross = vector1.Cross(vector2)
    inner_angle = math.acos(dot/(norm_vec_1*norm_vec_2))

    if cross < 0:
        return inner_angle

    return 2*math.pi-inner_angle


def SommeAngle(points):
    """
    The points define a polygon and this function returns the summ of the
    polygon's angles, taken clockwise.
    Return a tuple of the acutal angle summ and of the total summ that the
    polygon should theoretically have, according to the folowing formula :
        theoretical summ = (nb_vertices - 2) * 180          (in degrees)
    All the angles returned by this function are in degrees.
    """
    somme = 0
    pts = points.copy()
    if pts[0] != pts[-1]:
        pts.append(points[0])
        pts.append(points[1])
    else:
        pts.append(points[1])
    nb_vertices = 0
    for i in range(len(pts)-2):
        vector1 = vm.Vector2D((pts[i][0]-pts[i+1][0], pts[i][1]-pts[i+1][1]))
        vector2 = vm.Vector2D((pts[i+2][0]-pts[i+1][0], pts[i+2][1]-pts[i+1][1]))
        angle_clockwise = AngleClockwise(vector2, vector1)
        if angle_clockwise == math.pi:
            pass
        elif angle_clockwise == 2*math.pi:
            pass
        elif angle_clockwise < math.pi:
            somme += AngleClockwise(vector2, vector1)
            nb_vertices += 1
        elif angle_clockwise > math.pi:
            somme += AngleClockwise(vector2, vector1) - 2*math.pi
            nb_vertices -= 1
    somme_totale = (nb_vertices-2)*180
#    somme_totale = (len(pts)-4)*180
    return somme*180/math.pi, somme_totale



def OffsetHatch(points, offset):
    """
    Returns a hatched patch, corresponding to the polygon created by the points
    connected together and their offset in a direction.
    """
    len_points = len(points)
    line = primitives2D.RoundedLineSegments2D(points, {})
    offset_line = line.Offset(-offset)
    len_offset_points = len(offset_line.points)
    xy_coordinate = npy.zeros((len_points+len_offset_points, 2))
    for i in range(len_points):
        xy_coordinate[i][0] = line.points[i][0]
        xy_coordinate[i][1] = line.points[i][1]
    for i in range(len_offset_points):
        xy_coordinate[-i-1][0] = offset_line.points[i][0]
        xy_coordinate[-i-1][1] = offset_line.points[i][1]
    polygone = patches.Polygon(xy_coordinate, hatch='/')
    return polygone


def OffsetSurface(surface, offset):
    """
    Returns the offset points of the surface.
    """
    line = primitives2D.RoundedLineSegments2D(surface, {})
    offset_line = line.Offset(-offset)
    return offset_line.points


def PolygonFragments(polygon1, polygon2):
    """
    Returns the list of points of polygon1, fractioned with the points of
    polygon2, every time an edge of polygon1 intersects a corner of polygon2
    """
    polygon_points = []

    for i, point1 in enumerate(polygon1.points):
        merge_points = []
        line1 = vm.LineSegment2D(polygon1.points[i-1], point1)
        for j, point2 in enumerate(polygon2.points):
            distance = line1.PointDistance(point2)
            if math.isclose(distance, 0):
                # One point of intersection
                curvilinea_ratio = vm.LineSegment2D(polygon1.points[i-1], point2).Length()\
                                 / vm.LineSegment2D(polygon1.points[i-1], point1).Length()
                merge_points.append((j, point2, curvilinea_ratio))

        nb_merge_points = len(merge_points)
        polygon_points.append(polygon1.points[i-1])
        if nb_merge_points != 0:
            # One or more intersection points
            # Sort them according to the nearest to p1
            merge_points = sorted(merge_points, key=lambda p: p[2])
            for index in range(nb_merge_points):
                if  not math.isclose(merge_points[index][2], 0)\
                and not math.isclose(merge_points[index][2], 1):
                    polygon_points.append(merge_points[index][1])

    return polygon_points

def PolygonFusion2(polygon1, polygon2):
    """
    Fuses 2 polygons that are just touching each other.
    The 2 contours have common points.
    Returns a new contour.
    """
    polygon1_points = PolygonFragments(polygon1, polygon2)
    polygon2_points = PolygonFragments(polygon2, polygon1)

    polygon1_points = polygon1_points + polygon1_points
    polygon2_points = polygon2_points + polygon2_points

    fusion_polygon_points = []
    i = 1
    current_points = polygon1_points
    other_points = polygon2_points
    cursor = current_points[i]
    cursor0 = current_points[i-1]
    while cursor != cursor0:
        fusion_polygon_points.append(cursor)
        swap_lists = False
        for j, _ in enumerate(other_points):
            if other_points[j] == cursor:
                swap_lists = True
                break
        if swap_lists:
            current_points, other_points = other_points, current_points
            i = j+1
            cursor = current_points[i]
        else:
            i += 1
            cursor = current_points[i]
    fusion_polygon_points.append(cursor)
    fusion_polygon = vm.Polygon2D(fusion_polygon_points)
    return fusion_polygon

def OverlappingVerticalLinesIntersection(line1_y1, line1_y2, line2_y1, line2_y2):
    """
    Consider having 2 vertical lines with respective y1 and y2 as their
    starting and ending ordiante. These two lines are on the same line.
    This function tells you if your two lines are overlapping or not.
    """
    y1_max = max([line1_y1, line1_y2])
    y1_min = min([line1_y1, line1_y2])
    y2_max = max([line2_y1, line2_y2])
    y2_min = min([line2_y1, line2_y2])

    if y1_min >= y2_max or y1_max <= y2_min:
        return False

    return True

def points_inside_circle(points, circle):
    point_inside_the_circle = False
    for point in points:
        if point.PointDistance(circle.center) < circle.radius:
            point_inside_the_circle = True
            break
    return point_inside_the_circle


def rolling_circle_in_polygon(polygon, interpoints_distance=0.001):
    # discrétisation du polygon
    polygon_mesh = []
    # on parcourt les arrêtes

    for (vertice1, vertice2) in zip(polygon.points, polygon.points[1:]+[polygon.points[0]]):
        side_direction = vm.Vector2D((vertice2[0] - vertice1[0], vertice2[1] - vertice1[1]))
        normalized_side_direction = vm.Vector2D((vertice2[0] - vertice1[0], vertice2[1] - vertice1[1]))
        normalized_side_direction.Normalize()
        pt_number = 0
        # on ajoute les points un par un sans dépasser la longueur du côté
        segment_mesh = []
        while interpoints_distance * pt_number < side_direction.Norm():
            side_point = vertice1 + interpoints_distance * pt_number * normalized_side_direction
            segment_mesh.append(side_point)
            pt_number += 1
        polygon_mesh.append(segment_mesh)

    # prendre un point quelconque
    # construire le plus grand cercle possible
    min_radius = 1e+10
    for index1, segment1 in enumerate(polygon_mesh):
        for index2, segment2 in enumerate(polygon_mesh):
            if index2-index1 >= 2 and index2-index1 < len(polygon_mesh)-1:

                for point1 in segment1[1:]:
                    seg1 = vm.LineSegment2D(segment1[0], segment1[-1])
                    seg2 = vm.LineSegment2D(segment2[0], segment2[-1])

                    circle1, circle2  = seg2.CreateTangentCircle(point1, seg1)

                    polygon_mesh_modified = polygon_mesh[:]
                    polygon_mesh_modified.pop(index1)
                    polygon_mesh_modified = [p for seg in polygon_mesh_modified for p in seg]
                    if circle1 is not None and not points_inside_circle(polygon_mesh_modified, circle1) and circle1.radius < min_radius and polygon.PointBelongs(circle1.center):
                        min_radius = circle1.radius
                        min_segment1_index = index2
                        min_segment2_index = index1

                for point2 in segment2[1:]:
                    seg1 = vm.LineSegment2D(segment1[0], segment1[-1])
                    seg2 = vm.LineSegment2D(segment2[0], segment2[-1])

                    circle1, circle2  = seg1.CreateTangentCircle(point2, seg2)

                    polygon_mesh_modified = polygon_mesh[:]
                    polygon_mesh_modified.pop(index2)
                    polygon_mesh_modified = [p for seg in polygon_mesh_modified for p in seg]
                    if circle1 is not None and not points_inside_circle(polygon_mesh_modified, circle1) and circle1.radius < min_radius and polygon.PointBelongs(circle1.center):
                        min_radius = circle1.radius
                        min_segment1_index = index1
                        min_segment2_index = index2

    # returns the diamter of the smallest circle
    return min_radius*2, min_segment1_index, min_segment2_index

def is_wide_enought(polygon, minimal_width):
        actual_width, edge1_index, edge2_index = rolling_circle_in_polygon(polygon)
        if actual_width > minimal_width or math.isclose(actual_width, minimal_width, abs_tol=1e-08):
            return True, actual_width, None, None
        return False, actual_width, edge1_index, edge2_index


#def WidthFunction(width, closed_rounded_line):
#    offset_contour = closed_rounded_line.Offset(width/2)
##    offset_contour.MPLPlot(style='g')
#    polygon = vm.Polygon2D(offset_contour.points)
#    if polygon.SelfIntersect()[0]:
#        return -1
#    else:
#        return 1
#
#def ClosedRoundedLineWidth(closed_rounded_line):
#    rounded_line_x_length = max(p[0] for p in closed_rounded_line.points)-min(p[0] for p in closed_rounded_line.points)
#    rounded_line_y_length = max(p[1] for p in closed_rounded_line.points)-min(p[1] for p in closed_rounded_line.points)
#    rounded_line_min_dim = min(rounded_line_x_length, rounded_line_y_length)
#    actual_width = spy.optimize.brenth(WidthFunction, 0, rounded_line_min_dim, args=(closed_rounded_line,))
#    return actual_width

# TEST
#plt.close('all')
#points = []
#pts = [(0,0),(1,1),(1,0),(2,2),(0,1)]
#for pt in pts:
#    points.append(vm.Point2D(pt))
#linesegment2D = vm.primitives2D.RoundedLineSegments2D(points, {})
#line = (points[2], points[3])
#offset = 0.7
#fig, axe = linesegment2D.MPLPlot()
#LineOffset(linesegment2D, line, offset).MPLPlot(ax=axe, style='r')

#points = []
#pts = [(0,0),(0.1,0),(0.1,0.1),(0.05,0.1),(0.05,0.2),(0.1,0.2),(0.1,0.3),(0,0.3)]
#for pt in pts:
#    points.append(vm.Point2D(pt))
#linesegment2D = vm.primitives2D.RoundedLineSegments2D(points, {}, closed=True)
#linesegment2D.MPLPlot()
#actual_width = ClosedRoundedLineWidth(linesegment2D)
#print(actual_width)

###############################################################################

class ShaftMaterial(DessiaObject):
    """
    A material for a Part

    :param Re: The utilmate tensil stength.
    :type Re: Pa
    :param t_yield: The yield strength.
    :type y_yield: Pa
    :param C: The hardness on the Brinell Scale.
    :type C: Bhn
    :param name: The name of the ShaftMaterial.
    """
    def __init__(self, Re, t_yield, C, name='',):
        self.Re = Re            # Ultimate tensile strength (Pa)
        self.t_yield = t_yield  # Yield stength (Pa)
        self.C = C              # Hardness (Brinell Scale) (Bhn)
        self.name = name

###############################################################################

class Part(DessiaObject):
    """
    :param montage: How the Part should be mounted, either 'radial', 'axial' \
    or 'inbuilt'
    :param material: A ShaftMaterial Object
    :param name: The Part name
    """
    def __init__(self, montage='axial', material='', name=''):
        self.montage = montage
        self.material = material    # A ShaftMaterial
        self.name = name
        self.surfaces = []          # A list of PartSurface



    def PossibleShafts(self, yz=(0, 0)):
        """
        Computes a list of Shafts created by the different permutations of the
        surfaces discarding the circular permutations.

        :returns: **shafts** A list of Shaft Objects
        """
        perm = list(permutations(self.surfaces))

        # Because the problem is circular, we can remove some of the permutations
        remove_index = []
        for index, perm_i in enumerate(perm):
            if perm_i[0] != self.surfaces[0]:
                remove_index.append(index)
        for i in remove_index[::-1]:
            perm.pop(i)

        print('PossibleShafts : len(perm)', len(perm))

        # Checks the self intersection and the angles of the polygon
        shafts = []
        for i, perm_i in enumerate(list(perm)):
            perm_i = list(perm_i)
#            print('PossibleShafts : len(perm_i)', len(perm_i))
            points = []
            for j, partsurface in enumerate(perm_i):
                if partsurface.same_shaftline or min([p[1] for p in partsurface.points]) >= yz[1]:
                    points = points + partsurface.points
                else:
                    symm_partsurface = PartSurface([vm.Point2D((p[0], 2*yz[1]-p[1])) for p in partsurface.points], partsurface.otherpart, partsurface.part, partsurface.same_shaftline)
                    perm_i[j] = symm_partsurface
                    points = points + symm_partsurface.points

            if points[0] != points[-1]:
                points.append(points[0])
            line = primitives2D.RoundedLineSegments2D(points, {},
                                                         closed=True)
            shaft = Shaft(perm_i, line, yz=yz, montage=self.montage, name='{}_{}'.format(self.name, i))
            shafts.append(shaft)

        return shafts

    def ViableShafts(self, yz=(0, 0)):
        """
        Computes a list of Shafts based on the reparation of the PossibleShafts.

        :returns: **repaired_shafts** A list of viable Shaft Objects.
        """
        shafts = self.PossibleShafts(yz)
        repaired_shafts = []
        print('ViableShafts nb of shafts to repair', len(shafts))
        for shaft in shafts:
            repaired_shaft = shaft.GoodRepair()
            if repaired_shaft:
                repaired_shafts.append(repaired_shaft)

        repaired_shafts = list(chain.from_iterable(repaired_shafts))

        unique_repaired_shafts = []
        for shaft in repaired_shafts:
            boo = True
            for unique_shaft in unique_repaired_shafts:
                if shaft == unique_shaft:
                    boo = False
            if boo:
                unique_repaired_shafts.append(shaft)

        return unique_repaired_shafts

    def QuickViableShaft(self, yz=(0, 0)):

        print(self.name)
        print('nombre de partsurface', len(self.surfaces))

        inscreasing_x_surfaces = sorted(self.surfaces, key=lambda ps: ps.points[0][0])

        points = []
        for j, partsurface in enumerate(inscreasing_x_surfaces):
            if partsurface.same_shaftline or min([p[1] for p in partsurface.points]) >= yz[1]:
                points = points + partsurface.points
            else:
                symm_partsurface = PartSurface([vm.Point2D((p[0], 2*yz[1]-p[1])) for p in partsurface.points], partsurface.otherpart, partsurface.part, partsurface.same_shaftline)
                inscreasing_x_surfaces[j] = symm_partsurface
                points = points + symm_partsurface.points

        if points[0] != points[-1]:
            points.append(points[0])
        line = primitives2D.RoundedLineSegments2D(points, {}, closed=True)

        quick_shaft = Shaft(inscreasing_x_surfaces, line, yz=yz, montage=self.montage, name='{}_{}'.format(self.name, 'quick'))
        print('QuickViableShaft : QuickRepair')
        quick_repair_shaft = quick_shaft.QuickRepair()
        print('------------------------------')

        if quick_repair_shaft is None or not quick_repair_shaft.Viability():
            print('QuickViableShaft : RepairInsideOut')
            quick_repair_shafts = quick_shaft.RepairInsideOut()
            print('----------------------------------')
            if quick_repair_shafts is not None:
                quick_repair_shafts.sort(key= lambda a: a.polygon.Area())
                quick_repair_shaft = quick_repair_shafts[0]
            if quick_repair_shafts is None or not quick_repair_shaft.Viability():
#                print('pts', quick_shaft.polygon.points)
                print('QuickViableShaft : ViableShafts')
                quick_repair_shafts = self.ViableShafts(yz)
                print('-------------------------------')
                quick_repair_shafts.sort(key= lambda a: a.polygon.Area())
                quick_repair_shaft = quick_repair_shafts[0]

        return quick_repair_shaft

###############################################################################

class Surface(DessiaObject):
    """
    :param points: A list of coordinates, that when connected together draws a surface.
    :param part_left: A Part Object, located on the left hand side of the surface.
    :param part_right: A Part Object, located on the right hand side of the surface.
    """
    def __init__(self, points, part_left, part_right, same_shaftline=True):
        self.points = points            # A list of coordinates
        self.part_left = part_left      # Is a Part
        self.part_right = part_right    # Is a Part
        self.same_shaftline = same_shaftline

        if self.part_right is not None:
            self.partsurface_right = PartSurface(self.points, self.part_right, self.part_left, self.same_shaftline)
        if self.part_left is not None:
            self.partsurface_left = PartSurface(self.points[::-1], self.part_left, self.part_right, self.same_shaftline)

    @classmethod
    def Instantiate(cls, coordinate, part_left, part_right, same_shaftline=True):
        """
        A Class Method that convert the coordinate list into a list of volmdlr
        Points2D.
        """
        points = [vm.Point2D(p) for p in coordinate]

        return cls(points, part_left, part_right, same_shaftline)
        # functional_surface = [(point1, point2,...), part_left, part_right]
        # left or right according to the oriented path between the points

###############################################################################

class PartSurface(DessiaObject):
    """
    A rearrangement of Surface Object.
    """
    def __init__(self, points, part, otherpart=None, same_shaftline=True):
        self.points = points        # Is a list of coordiantes
        self.part = part            # Is a Part, it is on the right hand side
        self.otherpart = otherpart  # Is a Part
        self.same_shaftline = same_shaftline

        self.part.surfaces.append(self)

    @classmethod
    def Instantiate(cls, coordinate, part, otherpart=None, same_shaftline=True):
        """
        A Class Method that convert the coordinate list into a list of volmdlr
        Points2D.
        """
        points = [vm.Point2D(p) for p in coordinate]

        return cls(points, part, otherpart, same_shaftline)

###############################################################################

class Shaft(Part):
    """
    For a unique Part Object, different Shaft Objects can be created, because the
    arrangement of the PartSurface Objects influences the drawing of the Shaft's
    contour.

    :param ordered_partsurfaces: A list of PartSurface Objects. The list's order \
    determines the order followed to connect the PartSurfaces. It is only useful \
    if the Shaft Object needs to be repaired.

    :param contour: A closed RoundedLineSegment2D drawing the Shaft's contour. It \
    connects all the surfaces of ordered_partsurfaces.

    :param name: The name of the Shaft Object.
    """
    def __init__(self, ordered_partsurfaces, contour=None, yz=(0, 0),  montage='axial', name=''):
        Part.__init__(self, montage=montage, name=name)
        self.ordered_partsurfaces = ordered_partsurfaces    # A list of PartSurface
        self.contour = contour                              # A RoundedLineSegment
        self.yz = yz

        self.offset = 0.005
        self.EPSILON = 1e-6

        contour_points = []
        if self.contour is None or not self.contour:
            for partsurface in self.ordered_partsurfaces:
                for point in partsurface.points:
                    contour_points.append(point)
            if contour_points[0] != contour_points[-1]:
                contour_points.append(contour_points[0])
            self.contour = vm.primitives2D.RoundedLineSegments2D(contour_points, {}, closed=True)

        if self.ordered_partsurfaces is not None:
            self.ordered_surfaces = []
            for partsurface in self.ordered_partsurfaces:
                self.ordered_surfaces.append(partsurface.points)

        if self.contour is not None:
            polygon_points = self.contour.points
            index_to_del = []
            for i in range(len(polygon_points)):
                if polygon_points[i-1] == polygon_points[i]:
                    index_to_del.append(i-1)
            for index in index_to_del[::-1]:
                del polygon_points[index]
            self.polygon = vm.Polygon2D(polygon_points)

        if self.ordered_partsurfaces is not None:
            self.ordered_surfaces_index_in_contour = []
            for surface_points in self.ordered_surfaces:
                for (surface_point1, surface_point2) in zip(surface_points[:-1], surface_points[1:]):
                    for index, (contour_point1, contour_point2) in enumerate(zip(self.contour.points[:-1], self.contour.points[1:])):
                        if (surface_point1, surface_point2) == (contour_point1, contour_point2):
                            self.ordered_surfaces_index_in_contour.append(index)


    def __eq__(self, other_shaft):
        points1 = self.contour.points[::]
        points2 = other_shaft.contour.points[::]
        for i, pt2 in enumerate(points2):
            if pt2 == points1[0]:
                boo = True
                new_points2 = points2[i:]+points2[:i+1]
                for (pt1, new_pt2) in zip(points1, new_points2):
                    if pt1 != new_pt2:
                        boo = False
                return boo
        return False

    def __hash__(self):
        return len(self.contour.points)

    def Hatch(self):
        """
        Returns a hatched patch of each PartSurface of the Shaft Object.
        """
        polygon_hatches = []
        for surface in self.ordered_surfaces:
            polygon_hatches.append(OffsetHatch(surface, self.offset))
        return polygon_hatches

    def Plot(self, axe=None, hatch=True):
        """
        Plot the Shaft Object.
        :parma hatch: If True, adds to the plot a hatched patch towards the \
        interior of the Shaft for each PartSurfaces.
        """
        if axe is None:
            _, axe = self.contour.MPLPlot()
        else:
            self.contour.MPLPlot(ax=axe)
        axe.set_aspect('equal')
        if hatch and self.ordered_partsurfaces is not None:
            polygon_hatches = self.Hatch()
            for polygon in polygon_hatches:
                if polygon_hatches is not None:
                    axe.add_patch(polygon)

        for surface in self.ordered_surfaces:
            print('=>', surface)
            red_line = vm.primitives2D.RoundedLineSegments2D(surface, {})
            red_line.MPLPlot(ax=axe, style='r')

        return axe

    def AngleDifference(self):
        """
        Computes the angle difference for a Shaft Object between the actual angle
        summ and the theorical angle summ.
        """
        angle = SommeAngle(self.contour.points)
        return angle[1] - angle[0]


    def Viability(self):
        """
        Returns True if the Shaft is viable, meaning if its contour is not
        self intersecting and if the summ of its contour's angles is equal
        to the formula :
            (nb_vertecies - 2) * 180
        """
        angle = self.AngleDifference()
        if not abs(angle) < self.EPSILON:
            print('Viability : inside out', angle)
        if self.polygon.SelfIntersect()[0]:
            print('Viability : intersection')

        viable = abs(angle) < self.EPSILON\
                 and not self.polygon.SelfIntersect()[0]
        return viable


    def RepairInsideOut(self):
        """
        Repairs the polygons that are inside out.
        """
        if self.ordered_partsurfaces is not None:

            repaired_line = None
            repaired_shafts = []
            if not self.Viability():
                if self.AngleDifference() < self.EPSILON:

                    perm = list(permutations(self.ordered_surfaces))
#                    print(self.ordered_surfaces)


#                    # !!! : PAS SUR SI LE PROBLEME EST VRAIMENT CIRCULAIRE
#                    # Because the problem is circular, we can remove some of the permutations
#                    remove_index = []
#                    for index, perm_i in enumerate(perm):
#                        if perm_i[0] != self.ordered_surfaces[0]:
#                            remove_index.append(index)
#                    for i in remove_index[::-1]:
#                        perm.pop(i)

                    print('RepairInsideOut : len(perm)', self.name, len(perm))

                    perm_repaired_line = []
                    for print_i, perm_surfaces in enumerate(perm):
                        offset_points = []
                        repair_points = []
                        for surface in perm_surfaces:
                            line_surface = primitives2D.RoundedLineSegments2D(surface, {})
                            repair_points.extend(surface)
                            offset_points.extend(line_surface.Offset(-self.offset).points)
                        offset_points.reverse()
                        repair_points.extend(offset_points)

                        repair_points_line = repair_points.copy()
                        if repair_points[0] != repair_points[-1]:
                            repair_points_line.append(repair_points_line[0])

                        repaired_line = primitives2D.RoundedLineSegments2D(repair_points_line, {}, closed=True)
                        perm_repaired_line.append(repaired_line)
                        repaired_shaft = Shaft(self.ordered_partsurfaces, repaired_line, self.yz, self.montage, self.name)
                        repaired_shafts.append(repaired_shaft)

            if repaired_line is None:
                print("Repair failed")
                return None

            return repaired_shafts

    def RepairIntersection(self):
        """
        Repair the polygons that are lightly self intersecting, meaning that
        they can be repaired by extending the contour to some of the offset's
        points.
        """
        if self.ordered_partsurfaces is not None:

            repaired_line = None
            if not self.Viability():
                intersection = self.polygon.SelfIntersect()
                # If there is an intersection
                if intersection[0]:
                    intersect_line1 = intersection[1]
                    intersect_line2 = intersection[2]

                    for j, surface in enumerate(self.ordered_surfaces):
                        for p1_surface, p2_surface in zip(surface[::], surface[1::]):

                            # Looking for the functional surface that is intersected
                            # Possibly 4 different cases :
                            if intersect_line1.points == [p1_surface, p2_surface]\
                            or intersect_line2.points == [p2_surface, p1_surface]\
                            or intersect_line1.points == [p2_surface, p1_surface]\
                            or intersect_line2.points == [p1_surface, p2_surface]:

                                troncated_surfaces = self.ordered_surfaces[:j+1]

                                offset_points = OffsetSurface(troncated_surfaces[-1], self.offset)
                                i_point = 1
                                intersect0 = True
                                while intersect0 and i_point <= len(offset_points):
                                # New_surface is composed of the points just before the intersection,
                                # of the offset points in order to pass by the intersection, and
                                # finally of the remaining points that have been troncated.
                                    repaired_surface = list(troncated_surfaces)
                                    repaired_surface.append(list(offset_points[-ii] for ii in range(1, i_point+1)))

                                    if self.ordered_surfaces[j+1:]:
                                        repaired_surface.extend(self.ordered_surfaces[j+1:])

                                    repaired_points = list(chain.from_iterable(repaired_surface))

                                    if repaired_points[0] != repaired_points[-1]:
                                        repaired_points.append(repaired_points[0])
                                    intersect0 = vm.Polygon2D(repaired_points[:-1]).SelfIntersect()[0]
                                    i_point += 1

                                    repaired_line = primitives2D.RoundedLineSegments2D(repaired_points, {}, closed=True)

            if repaired_line is None:
                print("Repair failed")
                return None

            return Shaft(self.ordered_partsurfaces, repaired_line, self.yz, self.montage, self.name)

    def QuickRepair(self):

        if self.ordered_partsurfaces is not None:

            repaired_line = None
            if not self.Viability():
                if self.AngleDifference() < self.EPSILON:

                    perm_repaired_line = []
                    offset_points = []
                    repair_points = []
                    for surface in self.ordered_surfaces:
                        line_surface = primitives2D.RoundedLineSegments2D(surface, {})
                        repair_points.extend(surface)
                        offset_points.extend(line_surface.Offset(-self.offset).points)
                    offset_points.reverse()
                    repair_points.extend(offset_points)

                    repair_points_line = repair_points.copy()
                    if repair_points[0] != repair_points[-1]:
                        repair_points_line.append(repair_points_line[0])

                    repaired_line = primitives2D.RoundedLineSegments2D(repair_points_line, {}, closed=True)
                    perm_repaired_line.append(repaired_line)
                    repaired_shaft = Shaft(self.ordered_partsurfaces, repaired_line, self.yz, self.montage, self.name)
#                    repaired_shafts.append(repaired_shaft)

            if repaired_line is None:
                print("Repair failed")
                return None

            return repaired_shaft

    def Repair(self):
        """
        Repair the Shaft if it is not viable, using either the
        RepairIntersection or RepairInsideOut method according to the issue
        detected.
        """
        shafts = []
        if not self.Viability():
            if self.polygon.SelfIntersect()[0]:
                shaft = self.RepairIntersection()
                if shaft is not None:
                    shafts.append(shaft)
            if self.AngleDifference() < self.EPSILON:
                shaft = self.RepairInsideOut()
                if shaft is not None:
                    shafts.extend(shaft)
        return shafts

#    def Repair(self):
#        """
#        Repair the Shaft if it is not viable, using either the
#        RepairIntersection or RepairInsideOut method according to the issue
#        detected.
#        """
#        shafts = []
#        if not self.Viability():
#            if self.polygon.SelfIntersect()[0]:
#                shaft = self.RepairIntersection()
#                if shaft is not None:
#                    shafts.append(shaft)
#            if self.AngleDifference() < self.EPSILON:
#                shaft = self.RepairInsideOut()
#                if shaft is not None:
#                    shafts.extend(shaft)
#        return shafts


    def GoodRepair(self):
        """
        Repair the Shaft and returns only the viable reparations.
        """
        if self.Viability():
            return [self]
        shafts = self.Repair()
        good_shafts = []
        for shaft in shafts:
            if shaft.Viability():
                good_shafts.append(shaft)
        return good_shafts


    def Wider(self, minimal_width):
        """
        Faire attention à ce que le sens de parcours de self.contour fasse en
        sorte de laisser la matière à droite.
        """
        new_contour = self.contour
        new_polygon = self.polygon

        wide_enough, actual_width, edge1_index, edge2_index = is_wide_enought(new_polygon, minimal_width)

        not_offsetable_edge_indexes = []
        for index in self.ordered_surfaces_index_in_contour:
            not_offsetable_edge_indexes.append(index)
            if index != 0:
                not_offsetable_edge_indexes.append(index-1)
            else:
                not_offsetable_edge_indexes.append(len(self.contour.points))
            if index != len(self.contour.points):
                not_offsetable_edge_indexes.append(index+1)
            else:
                not_offsetable_edge_indexes.append(0)
        not_offsetable_edge_indexes = list(set(not_offsetable_edge_indexes))

        while not wide_enough:
            print('.')
            offset = minimal_width - actual_width
#            offset = - offset
            if edge1_index not in not_offsetable_edge_indexes and edge2_index not in not_offsetable_edge_indexes:
                new_contour = new_contour.OffsetLines([edge1_index], offset/2)
                new_contour = new_contour.OffsetLines([edge2_index], offset/2)
            elif edge1_index in not_offsetable_edge_indexes and edge2_index not in not_offsetable_edge_indexes:
                new_contour = new_contour.OffsetLines([edge2_index], offset)
            elif edge2_index in not_offsetable_edge_indexes and edge1_index not in not_offsetable_edge_indexes:
                new_contour = new_contour.OffsetLines([edge1_index], offset)
            else:
                print("The thinest area is inbetween two functionnal surfaces, so it can't be made wider")
                break

#            new_contour.MPLPlot()

            new_polygon = vm.Polygon2D([p for p in new_contour.points])

            wide_enough, actual_width, edge1_index, edge2_index = is_wide_enought(new_polygon, minimal_width)

        print('final actual_width', actual_width)

        new_shaft = Shaft(self.ordered_partsurfaces, new_contour, self.yz, self.montage, self.name)
#        new_shaft.Plot()
        return new_shaft

    def Contour(self):
        closed_contour = vm.primitives2D.RoundedLineSegments2D(self.contour.points, {}, closed=True)
        contour2D = vm.Contour2D([closed_contour])
        return contour2D

    @classmethod
    def DictToObject(cls, dict_):
        ordered_partsurfaces = dict_['ordered_partsurfaces']
        contour = dict_['contour']
        montage = dict_['montage']
        name = dict_['name']
        shaft = cls(ordered_partsurfaces, contour, montage, name)
        return shaft


    def Dict(self):
        d = {'ordered_partsurfaces': self.ordered_partsurfaces,
             'contour': self.contour,
             'montage': self.montage,
             'name': self.name}
        return d

    def CADVolumes(self, plane_origin=vm.O3D, center=vm.O3D, x=vm.X3D, y=vm.Y3D):
        z = x.Cross(y)
        z.Normalize()

        profile = vm3d.RevolvedProfile(plane_origin, x, y, [self.Contour()], center, vm.X3D)
#        profile.MPLPlot(axe)
        return profile

    def CADExport(self, fcstd_filepath, python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = vm.VolumeModel([('shaft', self.CADVolumes())])
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

###############################################################################

class SurfaceRepertory(DessiaObject):
    """
    Automation of the creation of Viable Shafts.
    self.shafts_product contains the product of all the Viable Shafts, if there
    is more than one for each Shaft.
    """
    def __init__(self, surfaces):
        self.surfaces = surfaces    # A list of Surface

        shafts = []
        parts = list(set([s.part_left  for s in self.surfaces]\
                        +[s.part_right for s in self.surfaces]))
        for part in parts:
            if part is not None:
                shafts.append(part.ViableShafts())

        self.shafts_product = [p for p in product(*shafts)]

###############################################################################

class ShaftLine(DessiaObject):
    """
    An assembly of Shaft objects.
    """
    def __init__(self, shafts, yz=(0, 0), name=''):
        self.shafts = shafts    # A list of Shaft Objects
        self.yz = yz # The y position of the revolved axis for the Shaft Assembly
        self.name = name

        print()
        print('shaftline => ', self.name, self.yz)
        for shaft in self.shafts:
            print(shaft.name, shaft.yz)
        print()

        self.viable_assembly_orders = None # An ordered list of Shaft Objects


    def ShaftConntections(self, shaft):
        """
        Gives you a list of Shaft that are connected to the entry Shaft
        """
        linked_shafts = []
        for shafti in self.shafts:
            if shafti is not shaft:
                for partsurface in shafti.ordered_partsurfaces:
                    if partsurface.otherpart == shaft.ordered_partsurfaces[0].part:
                        linked_shafts.append(shafti)
        return linked_shafts

    def ShaftMountability(self, shaft_to_mount_contour, contour):
        """
        Returns a dictionary that tells you according to the entry contour,
        whether the part is mountable from the left side, the right side, or
        not mountable at all.
        """
        mountable = {'left': True, 'right': True}

        if contour is None:
            return mountable

        part_functional_points = shaft_to_mount_contour.points
        y_max = max(p[1] for p in part_functional_points)
        y_min = min(p[1] for p in part_functional_points)
        x_max = max(p[0] for p in part_functional_points)
        x_min = min(p[0] for p in part_functional_points)

#        print('y_max', y_max)
#        print('y_min', y_min)
#        print('x_max', x_max)
#        print('x_min', x_min)

        # Check if the points of the Shaft contour are not in the axial
        # assembly zone
        for i in range(len(contour.points)-1):
            contour_point1 = contour.points[i]
            contour_point2 = contour.points[i+1]
            if contour_point1[0] <= x_min and contour_point2[0] <= x_min:
                overlap = OverlappingVerticalLinesIntersection(y_max, y_min,
                                                               contour_point1[1],
                                                               contour_point2[1])
                if overlap:
                    mountable['left'] = False
            if contour_point1[0] >= x_max and contour_point2[0] >= x_max:
                overlap = OverlappingVerticalLinesIntersection(y_max, y_min,
                                                               contour_point1[1],
                                                               contour_point2[1])
                if overlap:
                    mountable['right'] = False

        return mountable

    def MountShaftContour(self, shaft_to_mount, contour):
        """
        Returns the new contour of shaft_to_mount when mounted on contour. If
        the shaft is not mountable, returns the entry contour.
        """
        if contour is None:
            return shaft_to_mount.contour

        # In order to mount the parts together, one should be mountable on the
        # other and vis versa.
        mountable1 = self.ShaftMountability(shaft_to_mount.contour, contour)
        mountable2 = self.ShaftMountability(contour, shaft_to_mount.contour)

        mountable = {'left': mountable1['left'] and mountable2['right'],
                     'right': mountable1['right'] and mountable2['left']}

        if  (not mountable['left'])\
        and (not mountable['right'])\
        and shaft_to_mount.montage == 'axial':
            # The Shaft is not mountable
            return contour

        # Create the contour of the Part that will be assembled to Shaft
        if contour.points[0] == contour.points[-1]:
            polygon1_points = contour.points[:-1]
        else:
            polygon1_points = contour.points
        if shaft_to_mount.contour.points[0] == shaft_to_mount.contour.points[-1]:
            polygon2_points = shaft_to_mount.contour.points[:-1]
        else:
            polygon2_points = shaft_to_mount.contour.points

        polygon1 = vm.Polygon2D(polygon1_points)
        polygon2 = vm.Polygon2D(polygon2_points)
        new_polygon = PolygonFusion2(polygon1, polygon2)
        new_contour = primitives2D.RoundedLineSegments2D(new_polygon.points+[new_polygon.points[0]],
                                                            {})

        return new_contour


    def CreateGraph(self, draw=True):
        """
        Returns and plot a networkx graph, where the nodes are the shafts of
        the ShaftLine and the edges represent the contact surfaces between
        the shafts.
        """

        G = nx.Graph()

        G.add_nodes_from(self.shafts)

        edges = []
        for shaft in self.shafts:
            if shaft.ordered_partsurfaces is not None:
                linked_shafts = self.ShaftConntections(shaft)
            else:
                linked_shafts = 'dunno'

            for link_shaft in linked_shafts:
                edges.append((shaft, link_shaft))

        non_duplicated_edges = []
        for node1, node2 in edges:
            if (node1, node2) not in non_duplicated_edges and (node2, node1) not in non_duplicated_edges:
                non_duplicated_edges.append((node1, node2))

        G.add_edges_from(non_duplicated_edges)

#        pos = nx.kamada_kawai_layout(G)
        
        pos = nx.spectral_layout(G)
        if draw:
            plt.figure()
            nx.draw_networkx_nodes(G, pos)
            nx.draw_networkx_edges(G, pos)
            labels = {s:s.name for s in self.shafts}
            nx.draw_networkx_labels(G, pos, labels)

        return G


    def NextShaftsToMount(self, mounted_shafts, viability=False):
        """
        Returns the list of the possible next shafts to mount, taking the
        networkx grpah into account, so that only touching part can be mounted
        together.
        If viability is True, the method also checks that the shafts can
        actually be axially mounted.
        """
        G = self.CreateGraph(draw=False)
        nodes = list(G.nodes())
        potential_next_shafts = [s for s in nodes if s not in mounted_shafts]
        next_shafts = []
        edges = list(G.edges())
        edges_count = {e:0 for e in edges}
        for edge in edges:
            for shaft in mounted_shafts:
                if shaft in edge:
                    edges_count[edge] += 1

        for edge, count in edges_count.items():
            if edge[0] in potential_next_shafts and count == 1:
                next_shafts.append(edge[0])
            if edge[1] in potential_next_shafts and count == 1:
                next_shafts.append(edge[1])
        next_shafts = list(set(next_shafts))

        if viability:
            viable_next_shafts = []
            for next_shaft in next_shafts:
                test_assembly = mounted_shafts + [next_shaft]
                if self.AssemblyViability(test_assembly):
                    viable_next_shafts.append(next_shaft)
            return viable_next_shafts

        return next_shafts


    def OrderedViableMountage(self):
        """
        Returns one or many ordered lists of shafts, showing the order
        according to which the ShaftLine can be mounted.
        """
        self.CreateGraph(draw=False)
        ordered_assemblies = []

        dt = dectree.DecisionTree()

        len_shafts = len(self.shafts)

        while not dt.finished:
            valid = True

            if dt.current_depth == 0:
                dt.SetCurrentNodeDataPossibilities([s for s in self.shafts])

            elif dt.current_depth < len_shafts:

                current_assembly = []
                for index in range(len(dt.current_node)):
                    current_assembly.append(dt.data[tuple(dt.current_node[:index+1])])

                viable_next_shafts = self.NextShaftsToMount(current_assembly, viability=True)
                if not viable_next_shafts:
                    valid = False
                    dt.SetCurrentNodeNumberPossibilities(0)
                else:
                    dt.SetCurrentNodeDataPossibilities(viable_next_shafts)


            elif dt.current_depth == len_shafts:

                current_assembly = []
                for index in range(len(dt.current_node)):
                    current_assembly.append(dt.data[tuple(dt.current_node[:index+1])])
                ordered_assemblies.append(current_assembly)
                dt.SetCurrentNodeNumberPossibilities(0)
                valid = False

            dt.NextNode(valid)

        self.viable_assembly_orders = ordered_assemblies

        return ordered_assemblies


    def CreateContour(self, ordered_shafts=None):
        """
        Returns the contour of the ShaftLine mounted according to the
        order of the entry Shaft list. If the ShaftLine is not mountable
        in this specific order, the method returns None.
        """
        
        if ordered_shafts is None:
            ordered_shafts = self.OrderedViableMountage()[0]
        
        contour = None

        for shaft in ordered_shafts:
            new_contour = self.MountShaftContour(shaft, contour)
            if contour == new_contour:
                return None
            contour = new_contour
        return contour

    def AssemblyViability(self, ordered_shafts):
        """
        We suppose that shafts is ordered according to the assembly patern.
        """
        contour = self.CreateContour(ordered_shafts)
        if contour is None:
            return False
        return True


    def GraphicMountability(self, ordered_assembly=None):
        """
        Plots a simple sequenced visualisation of the mountage following the
        entry order.
        """
        if ordered_assembly is None:
            ordered_assembly = self.OrderedViableMountage()[0]

        # TODO: It's a detail but the subplot grid can be improved
        len_ordered_assembly = len(ordered_assembly)
#        i = int(npy.floor(len_ordered_assembly/2))
#        j = int(npy.ceil(len_ordered_assembly/2))

        i = int(npy.ceil(npy.sqrt(len_ordered_assembly)))

        plt.figure()
        for index in range(len_ordered_assembly):

            contour = self.CreateContour(ordered_assembly[:index+1])
            axe = plt.subplot(i, i, index+1)
            axe.set_aspect('equal')
            contour.MPLPlot(axe)
            ordered_assembly[index].contour.MPLPlot(axe, style='r')

    def Plot(self, axe=None, symm_y_axis=None):
        """
        Plots the ShaftLine Object with delimitation for every Shaft.
        """
#        contour = self.CreateContour(self.viable_assembly_orders[0])
#        _, axe = contour.MPLPlot()
        if axe is None:
            fig, axe = plt.subplots()

        if symm_y_axis is None:
            symm_y_axis = self.yz[1]
        print(symm_y_axis)

        translated_contours = []
        for shaft in self.shafts:
            translated_contour = shaft.contour.Translation(vm.Vector2D((0, symm_y_axis)))
            _, axe = translated_contour.MPLPlot(axe)
            translated_contours.append(translated_contour)

        axe.set_aspect('equal')
        return axe, translated_contours

    def SymmetricPlot(self, symm_y_axis=None, axe=None):
        """
        Plots the whole ShaftLine with its symmetric.
        """
        if symm_y_axis is None:
            symm_y_axis = self.yz[1]

        if axe is None:
            axe, contours = self.Plot(None, symm_y_axis)
        else:
            _, contours = self.Plot(axe, symm_y_axis)
        axe.set_aspect('equal')

        for shaft in self.shafts:
            symm_points_shaft = []
            shaft_contour = shaft.contour.Translation(vm.Vector2D((0, symm_y_axis)))
            for point in shaft_contour.points:
                symm_points_shaft.append(point.Translation(vm.Point2D((0,-2*(point[1]-symm_y_axis)))))
            symm_shaft_contour = vm.primitives2D.RoundedLineSegments2D(symm_points_shaft, {}, closed=True)
            symm_shaft_contour.MPLPlot(axe)

        symm_line = vm.LineSegment2D(vm.Point2D((min([p[0] for shaft in self.shafts for p in shaft.contour.points]), symm_y_axis)),
                                     vm.Point2D((max([p[0] for shaft in self.shafts for p in shaft.contour.points]), symm_y_axis)))
        symm_line.MPLPlot(axe, style='-.k')

        return axe

    def CADVolumes(self, center=None, x=vm.X3D, y=vm.Y3D):
        if center is None:
            center = vm.Point3D((0, self.yz[0], self.yz[1]))
        z = x.Cross(y)
        z.Normalize()
        profiles = []
        for shaft in self.shafts:
            profiles.append(shaft.CADVolumes(center=center))
#            shaft.CADVolumes()[0].MPLPlot()
        return profiles

    def CADExport(self, fcstd_filepath, python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = vm.VolumeModel([('shaftline', self.CADVolumes())])
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

###############################################################################

class ShaftLinesAssembly(DessiaObject):
    def __init__(self, shaftlines, shaftlines_to_position, origin=(0, 0), name=''):
        self.shaftlines = shaftlines
        self.shaftlines_to_position = shaftlines_to_position  # A dictionary
        self.origin = origin
        self.name = name
        
#    def UnfoldPlot(self):
#        print(self.shaftlines_to_position)
#        heights = []
#        heights.append(0)
#        for shaftline1, shaftline2 in zip(self.shaftlines[:-1], self.shaftlines[1:]):
#            d = ((self.shaftlines_to_position[shaftline2][0]-self.shaftlines_to_position[shaftline1][0])**2 + (self.shaftlines_to_position[shaftline2][1]-self.shaftlines_to_position[shaftline1][1])**2)**0.5
#            print(d)
#            heights.append(-d + heights[-1])
#        fig, axe = plt.subplots()
#        for sl_index, shaftline in enumerate(self.shaftlines):
#            shaftline.SymmetricPlot(symm_y_axis=heights[sl_index], axe=axe)
#        
#        return heights

    def Plot(self):
        fig, axe = plt.subplots()
        for shaftline in self.shaftlines:
            shaftline.SymmetricPlot(symm_y_axis=self.shaftlines_to_position[shaftline][1], axe=axe)
            
    def ShaftLineEnvelope(self, shaftline):
        position = self.shaftlines_to_position[shaftline]
        
        restricted_zones = []
        for other_shaftline in self.shaftlines:
            if other_shaftline is not shaftline:
                print('ca arrive')
                contour = other_shaftline.CreateContour()
                print('c parti')
                symm_contour_points = [p.Translation(vm.Point2D((0, -2*p[1]))) for p in contour.points]
                symm_contour = vm.primitives2D.RoundedLineSegments2D(symm_contour_points, {})
                other_position = self.shaftlines_to_position[other_shaftline]
                d = ((position[0]-other_position[0])**2 + (position[1]-other_position[1])**2)**0.5
                translation_vector = vm.Vector2D((position[0]-other_position[0], position[1]-other_position[0]+d))
                translated_contour = vm.primitives2D.RoundedLineSegments2D([p.Translation(vm.Point2D(translation_vector)) for p in contour.points], {})
                translated_symm_contour = vm.primitives2D.RoundedLineSegments2D([p.Translation(vm.Point2D(translation_vector)) for p in symm_contour.points], {})
                restricted_zones.append(translated_symm_contour)
                
        fig, axe = plt.subplots()
        shaftline.Plot(axe=axe, symm_y_axis=position[0])
        for contour in restricted_zones:
            contour.MPLPlot(ax=axe, style='r')
            
    def CADVolumes(self, center=None, x=vm.X3D, y=vm.Y3D):
        if center is None:
            center = self.origin
        z = x.Cross(y)
        z.Normalize()
        profiles = []
        for shaftline in self.shaftlines:
            profiles = profiles + shaftline.CADVolumes(self.shaftlines_to_position[shaftline])
        return profiles

    def CADExport(self, fcstd_filepath, python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = vm.VolumeModel([('shaftlineassembly', self.CADVolumes())])
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

###############################################################################