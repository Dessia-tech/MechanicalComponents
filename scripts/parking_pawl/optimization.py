import math
import dessia_common.optimization as dc_opt
import mechanical_components.parking_pawl as mcpp

w = w = -3/3.6/(0.73/2)*12
travel = 0.011

optimization_bounds = [dc_opt.BoundedAttributeValue('wheel_lower_tooth_diameter',
                                                    0.02, 0.05),
                       dc_opt.BoundedAttributeValue('relative_basis_diameter',
                                                    0, 0.02),
                       dc_opt.BoundedAttributeValue('relative_contact_diameter',
                                                    0, 0.02),
                       dc_opt.BoundedAttributeValue('relative_wheel_outer_diameter',
                                                    0, 0.02),
                       dc_opt.BoundedAttributeValue('pawl_offset',
                                                    0.03, 0.09),   
                       dc_opt.BoundedAttributeValue('slope_angle',
                                                    math.radians(20), math.radians(80)),   
                       dc_opt.BoundedAttributeValue('finger_width_ratio',
                                                    0.05, 0.3),   
                       ]

optimizer = mcpp.ParkingPawlOptimizer(wheel_locking_speed=w, locking_mechanism_travel=travel,
                          optimization_bounds=optimization_bounds)

result = optimizer.optimize_gradient()

result.babylonjs()