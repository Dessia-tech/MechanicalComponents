import math
import dessia_common.optimization as dc_opt
import mechanical_components.parking_pawl as mcpp

w = w = -3/3.6/(0.73/2)*12
travel = 0.011

optimization_bounds = [dc_opt.BoundedAttributeValue('wheel_lower_tooth_diameter',
                                                    0.02, 0.05),
                       dc_opt.BoundedAttributeValue('relative_basis_diameter',
                                                    0.0005, 0.02),
                       dc_opt.BoundedAttributeValue('relative_contact_diameter',
                                                    0.0005, 0.02),
                       dc_opt.BoundedAttributeValue('relative_wheel_outer_diameter',
                                                    0.0005, 0.02),
                       dc_opt.BoundedAttributeValue('pawl_offset',
                                                    0.03, 0.09),   
                       dc_opt.BoundedAttributeValue('slope_angle',
                                                    math.radians(20), math.radians(80)),   
                       dc_opt.BoundedAttributeValue('finger_width',
                                                    0.004, 0.015),
                       dc_opt.BoundedAttributeValue('roller_diameter',
                                                    0.010, 0.030),
                       dc_opt.BoundedAttributeValue('axis_outer_diameter',
                                                    0.015, 0.035),
                       dc_opt.BoundedAttributeValue('slope_start_finger_overheight',
                                                    0.005, 0.025),
                       
                       ]

optimizer = mcpp.ParkingPawlOptimizer(wheel_locking_speed=w, locking_mechanism_travel=travel,
                          optimization_bounds=optimization_bounds)

result = optimizer.optimize_gradient()
result.pawl.size_torsion_spring(10 * 9.81)
result_simulation = result.locking_simulation(w)
result_simulation.babylonjs()
result_simulation.mpl_plot()
# result.babylonjs()