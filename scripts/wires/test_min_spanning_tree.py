import mechanical_components.wires as wires
import mechanical_components.optimization.wires as wires_opt
import numpy as npy
import random
import mechanical_components.optimization.wires_protected as protected_module
import volmdlr as vm
import networkx as nx
import matplotlib.pyplot as plt

point1 = vm.Point3D(1,0,0)
point2 = vm.Point3D(1,0.5,0)
point3 = vm.Point3D(1,1,0)
point4 = vm.Point3D(0,0,0)
point5 = vm.Point3D(-0.5,0.5,0)
point6 = vm.Point3D(-1,1,0)
point7 = vm.Point3D(-1,-0.5,0)
point8 = vm.Point3D(-0.5,-1,0)
point9 = vm.Point3D(-0.5,-0.5,0)
point10 = vm.Point3D(0.5,-0.5,0)
point11 = vm.Point3D(0,1,0)



edges = [(point1, point2), (point2, point3), (point3, point4), (point4, point5), (point5, point7), (point7, point8), (point8, point9), (point9, point10), (point4, point9), (point3, point11), (point11, point5), 
         (point8, point10), (point10, point1)]

graph = nx.Graph()
graph.add_edges_from(edges)
print(list(graph.edges))
ax = list(graph.edges)[0][0].plot()
for start, end in list(graph.edges):
    start.plot(ax=ax)
    end.plot(ax=ax)
    line = vm.edges.LineSegment3D(start, end)
    line.plot(ax=ax, color = 'g')
    

wire_routing = protected_module.WireRouting(edges)
mst = wire_routing.minimim_spannig_tree()
ax = list(mst.edges)[0][0].plot()
for start, end in list(mst.edges):
    start.plot(ax=ax)
    end.plot(ax=ax)
    line = vm.edges.LineSegment3D(start, end)
    line.plot(ax=ax, color = 'r')

