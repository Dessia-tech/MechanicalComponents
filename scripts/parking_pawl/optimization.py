import math
import dessia_common.optimization as dc_opt
import mechanical_components.parking_pawl as mcpp

w = w = -3/3.6/(0.73/2)*12
travel = 0.011
mcpp.ParkingPawlOptimizer(wheel_locking_speed=w, locking_mechanism_travel=travel,
                          optimization_bounds=[dc_opt.BoundedAttributeValue('')])