import math
import dessia_common.optimization as dc_opt
import mechanical_components.parking_pawl as mcpp

wheels_speed = -3/3.6/(0.73/2)*12

travel = 0.011
wheels_torque = 3764
max_wheels_slack = 0.25
ratios = [1, 4, 12]


optimization_bounds = [dc_opt.BoundedAttributeValue('wheel_lower_tooth_diameter',
                                                    0.035, 0.05),
                       dc_opt.BoundedAttributeValue('pressure_angle',
                                                    math.radians(5), math.radians(40)),
                       dc_opt.BoundedAttributeValue('relative_contact_diameter',
                                                    0.002, 0.02),
                       dc_opt.BoundedAttributeValue('relative_wheel_outer_diameter',
                                                    0.002, 0.02),
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
                                                    0.004, 0.020),
                       dc_opt.BoundedAttributeValue('lower_tooth_ratio',
                                                    0.4, 0.7),
                       ]

parking_pawl_simulations = []
for tn in range(7, 9):
    fixed_parameters = [dc_opt.FixedAttributeValue('teeth_number', tn)]

    for ratio in ratios:
        wheel_speed = ratio*wheels_speed
        wheel_torque = wheels_torque/ratio
        max_engaged_slack = max_wheels_slack/ratio
        optimizer = mcpp.ParkingPawlOptimizer(wheel_locking_speed=wheel_speed,
                                              wheel_torque=wheel_torque,
                                              max_engaged_slack=max_engaged_slack,
                                              locking_mechanism_travel=travel,
                                              fixed_parameters=fixed_parameters,
                                              optimization_bounds=optimization_bounds)
        
        result, objective = optimizer.optimize_cma()
        if result.check():
            # result.babylonjs()
            result.pawl.size_torsion_spring(10 * 9.81)
            result_simulation = result.locking_simulation(wheel_speed)
            
            parking_pawl_simulations.append(result_simulation)
            
simulation_list = mcpp.ParkingPawSimulationlList(parking_pawl_simulations)
simulation_list.plot()