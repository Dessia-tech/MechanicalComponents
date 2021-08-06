import mechanical_components.wires as wires
import mechanical_components.optimization.wires as wires_opt
import numpy as npy
import random

import volmdlr as vm


n_wpts = (15, 5, 4)# Length, width, heightvm.Line2D(waypoints[i),waypoints[i+1]).DirectionVector(unit=True).Dot(vm.Line2D(waypoints[i+1), waypoints[i+2]).DirectionVector(unit=True))==1
grid_size = (0.16, 0.15, 0.08)
min_length_paths = n_wpts[0] + 2
n_wires = 20
connection_probability = 1

waypoints = []
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]):
            grid_point = vm.Point3D(i*grid_size[0], j*grid_size[1], k*grid_size[2])
            waypoints.append(grid_point)

routes = []            
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                                waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k+1]))

for i in range(n_wpts[0]):
    for k in range(n_wpts[2]):
        for j in range(n_wpts[1]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                                waypoints[i*n_wpts[1]*n_wpts[2] + (j+1)*n_wpts[2] +k]))
                
for j in range(n_wpts[1]):
    for k in range(n_wpts[2]):
        for i in range(n_wpts[0]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                                waypoints[(i+1)*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k]))


wo = wires_opt.WiringOptimizer(routes)


wires_specs = []
connected_sources = []
for i in range(n_wires):
    source = random.choice(waypoints[:2*n_wpts[0]])
    destination = random.choice(waypoints[-2*n_wpts[0]:])
    
    wires_specs.append(wires.RoutingSpec(source=source,
                                          destination=destination,
                                          diameter=0.005 + 0.005*random.random()))

wiring = wo.multi_source_multi_destination_routing(wires_specs)
wiring.babylonjs()
  
# source_destination = [(vm.Point3D( 0.0, 0.15, 0.24), vm.Point3D( 2.24, 0.3, 0.0)),
#   (vm.Point3D( 0.0, 0.0, 0.0), vm.Point3D( 2.24, 0.15, 0.08)),
#   (vm.Point3D( 0.0, 0.15, 0.16), vm.Point3D( 2.24, 0.15, 0.16)),
#   (vm.Point3D( 0.0, 0.3, 0.0), vm.Point3D( 2.24, 0.44999999999999996, 0.0)),
#   (vm.Point3D( 0.0, 0.15, 0.24), vm.Point3D( 2.24, 0.6, 0.16))]

# wires_specs =[wires.RoutingSpec(source=source,
#                                             destination=destination,
#                                             diameter=0.005 + 0.005*random.random()) for source, destination in source_destination]

# paths = wo.multi_source_multi_destination_routing(wires_specs)
