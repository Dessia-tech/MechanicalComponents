#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:20:08 2019

@author: ringhausen
"""

from itertools import chain
from itertools import permutations
from itertools import product
import math
import numpy as npy
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import dectree
import volmdlr as vm
from volmdlr import primitives2D
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

        if AngleClockwise(vector2, vector1) < math.pi:
            somme += AngleClockwise(vector2, vector1)
            nb_vertices += 1
        if AngleClockwise(vector2, vector1) >= math.pi:
            somme += AngleClockwise(vector2, vector1) - 2*math.pi
            nb_vertices -= 1
    somme_totale = (nb_vertices-2)*180
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

###############################################################################

class ShaftMaterial:
    """
    A material for a Part

    :param Re: The utilmate tensil stength.
    :type Re: Pa
    :param t_yield: The yield strength.
    :type y_yield: Pa
    :param C: The hardness on the Brinell Scale
    :type C: Bhn
    :param name: The name of the ShaftMaterial
    """
    def __init__(self, Re, t_yield, C, name='',):
        self.Re = Re            # Ultimate tensile strength (Pa)
        self.t_yield = t_yield  # Yield stength (Pa)
        self.C = C              # Hardness (Brinell Scale) (Bhn)
        self.name = name

###############################################################################

class Part:
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



    def PossibleShafts(self):
        """
        Computes a list of Shafts created by the different permutations of the
        surfaces discarding the circular permutations.

        :returns: **shafts** A list of Shaft Objects
        """
        perm = list(permutations(self.surfaces))

        # Because the problem is circular, we can remove some of the permutations
        remove_index = []
        for index, elem in enumerate(perm):
            if elem[0] != self.surfaces[0]:
                remove_index.append(index)
        for i in remove_index[::-1]:
            perm.pop(i)

        # Checks the self intersection and the angles of the polygon
        shafts = []
        for i, perm_i in enumerate(perm):
            points = list(chain.from_iterable(j.points for j in perm_i))

            if points[0] != points[-1]:
                points.append(points[0])
            line = primitives2D.RoundedLineSegments2D(points, {},
                                                         closed=True)
            shaft = Shaft(perm_i, line, name='{}_{}'.format(self.name, i))
            shafts.append(shaft)

        return shafts

    def ViableShafts(self):
        """
        Computes a list of Shafts based on the reparation of the PossibleShafts.

        :returns: **repaired_shafts** A list of viable Shaft Objects.
        """
        shafts = self.PossibleShafts()
        repaired_shafts = []
        for shaft in shafts:
            repaired_shaft = shaft.GoodRepair()
            if repaired_shaft:
                repaired_shafts.append(repaired_shaft)

        repaired_shafts = list(chain.from_iterable(repaired_shafts))

        return repaired_shafts

###############################################################################

class Surface:
    """
    :param points: A list of coordinates, that when connected together draws a surface.
    :param part_left: A Part Object, located on the left hand side of the surface.
    :param part_right: A Part Object, located on the right hand side of the surface.
    """
    def __init__(self, points, part_left, part_right):
        self.points = points            # A list of coordinates
        self.part_left = part_left      # Is a Part
        self.part_right = part_right    # Is a Part

        PartSurface(self.points, self.part_right, self.part_left)
        PartSurface(self.points[::-1], self.part_left, self.part_right)

    @classmethod
    def Instantiate(cls, coordinate, part_left, part_right):
        """
        A Class Method that convert the coordinate list into a list of volmdlr
        Points2D.
        """
        points = [vm.Point2D(p) for p in coordinate]
        return cls(points, part_left, part_right)
        # functional_surface = [(point1, point2,...), part_left, part_right]
        # left or right according to the oriented path between the points

###############################################################################

class PartSurface:
    """
    A rearrangement of Surface Object.
    """
    def __init__(self, points, part, otherpart):
        self.points = points        # Is a list of coordiantes
        self.part = part            # Is a Part
        self.otherpart = otherpart  # Is a Part

        self.part.surfaces.append(self)

###############################################################################

class Shaft(Part):
    """
    For a unique Part Object, different Shaft Objects can be created, because the
    arrangement of the PartSurface Objects influences the drawing of the Shafts's
    contour.

    :param ordered_partsurfaces: A list of PartSurface Objects. The list's order \
    determines the order followed to connect the PartSurfaces.

    :param contour: A closed RoundedLineSegment2D drawing the Shaft's contour. It \
    connects all the surfaces of ordered_partsurfaces.

    :param name: The name of the Shaft Object.
    """
    def __init__(self, ordered_partsurfaces, contour, name=''):
        Part.__init__(self, name=name)
        self.ordered_partsurfaces = ordered_partsurfaces    # A list of PartSurface
        self.contour = contour                              # A RoundedLineSegment

        self.offset = 0.001
        self.EPSILON = 1e-6

        self.ordered_surfaces = []
        for partsurface in self.ordered_partsurfaces:
            self.ordered_surfaces.append(partsurface.points)

        if self.contour is not None:
            polygon_points = contour.points
            index_to_del = []
            for i in range(len(polygon_points)):
                if polygon_points[i-1] == polygon_points[i]:
                    index_to_del.append(i-1)
            for index in index_to_del[::-1]:
                del polygon_points[index]
            self.polygon = vm.Polygon2D(polygon_points)



    def Hatch(self):
        """
        Returns a hatched patch of each PartSurface of the Shaft Object.
        """
        polygon_hatches = []
        for surface in self.ordered_surfaces:
            polygon_hatches.append(OffsetHatch(surface, self.offset))
        return polygon_hatches

    def Plot(self, hatch=True):
        """
        Plot the Shaft Object.
        :parma hatch: If True, adds to the plot a hatched patch towards the \
        interior of the Shaft for each PartSurfaces.
        """
        _, axe = self.contour.MPLPlot()
        if hatch:
            polygon_hatches = self.Hatch()
            for polygon in polygon_hatches:
                if polygon_hatches is not None:
                    axe.add_patch(polygon)

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
        viable = abs(angle) < self.EPSILON\
                 and not self.polygon.SelfIntersect()[0]
        return viable


    def RepairInsideOut(self):
        """
        Repairs the polygons that are inside out.
        """
        repaired_line = None
        repaired_shafts = []
        if not self.Viability():
            if self.AngleDifference() < self.EPSILON:

                perm = list(permutations(self.ordered_surfaces))
                perm_repaired_line = []
                for perm_surfaces in perm:
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

                    repaired_line = primitives2D.RoundedLineSegments2D(repair_points_line, {})
                    perm_repaired_line.append(repaired_line)
                    repaired_shaft = Shaft(self.ordered_partsurfaces, repaired_line, self.name)
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

                                repaired_line = primitives2D.RoundedLineSegments2D(repaired_points, {})

        if repaired_line is None:
            print("Repair failed")
            return None

        return Shaft(self.ordered_partsurfaces, repaired_line, self.name)

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

###############################################################################

class SurfaceRepertory:
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
            shafts.append(part.ViableShafts())

        self.shafts_product = [p for p in product(*shafts)]

###############################################################################

class ShaftAssembly:
    """
    An assembly of Shaft objects.
    """
    def __init__(self, shafts):
        self.shafts = shafts    # A list of Shaft

        self.viable_assembly_orders = self.OrderedViableMountage()


    def ShaftConntections(self, shaft):
        """
        Gives you a list of Shaft that are connected to the entry Shaft
        """
        linked_shafts = []
        for shafti in self.shafts:
            if shafti != shaft:
                for partsurface in shafti.ordered_partsurfaces:
                    if partsurface.otherpart == shaft.ordered_partsurfaces[0].part:
                        linked_shafts.append(shafti)
        return linked_shafts



    def ShaftMountability(self, shaft_to_mount, contour):
        """
        Returns a dictionary that tells you according to the entry contour,
        whether the part is mountable from the left side, the right side, or
        not mountable at all.
        """

        mountable = {'left': True, 'right': True}

        if contour is None:
            return mountable

        part_functional_points = shaft_to_mount.contour.points
        y_max = max(p[1] for p in part_functional_points)
        y_min = min(p[1] for p in part_functional_points)
        x_max = max(p[0] for p in part_functional_points)
        x_min = min(p[0] for p in part_functional_points)

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

        mountable = self.ShaftMountability(shaft_to_mount, contour)

        if  (not mountable['left'])\
        and (not mountable['right'])\
        and shaft_to_mount.ordered_partsurfaces[0].part.montage == 'axial':
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
        the ShaftAssembly and the edges represent the contact surfaces between
        the shafts.
        """

        G = nx.Graph()
        G.add_nodes_from(self.shafts)
#        for shaft in self.shafts:
#            G.add_node(shaft.name, attribute=shaft)

        edges = []
        for shaft in self.shafts:
            linked_shafts = self.ShaftConntections(shaft)
            for link_shaft in linked_shafts:
                edges.append((shaft, link_shaft))

        G.add_edges_from(edges)
        pos = nx.kamada_kawai_layout(G)
        if draw:
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
        according to which the ShaftAssembly can be mounted.
        """

        self.CreateGraph()
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

        return ordered_assemblies


    def CreateContour(self, ordered_shafts):
        """
        Returns the contour of the ShaftAssembly mounted according to the
        order of the entry Shaft list. If the ShaftAssembly is not mountable
        in this specific order, the method returns None.
        """
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


    def GraphicMountability(self, ordered_assembly):
        """
        Plots a simple sequenced visualisation of the mountage following the
        entry order.
        """
        len_ordered_assembly = len(ordered_assembly)
        i = int(npy.floor(len_ordered_assembly/2))
        j = int(npy.ceil(len_ordered_assembly/2))

        plt.figure()
        for index in range(len_ordered_assembly):

            contour = self.CreateContour(ordered_assembly[:index+1])
            axe = plt.subplot(i, j, index+1)
            contour.MPLPlot(axe)
            ordered_assembly[index].contour.MPLPlot(axe, style='r')


###############################################################################
###############################################################################