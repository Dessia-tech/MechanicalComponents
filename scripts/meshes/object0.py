#import sys
#del sys.modules['mechanical_components']
import mechanical_components.meshes as meshes_opt
import numpy as npy

# Definition of a Rack
rack1=meshes_opt.Rack(transverse_pressure_angle = 20/180.*npy.pi)
rack1.Update(module = 1)
check,_=rack1.CheckRackViable()
if check:
    rack1.Plot(number_pattern = 10)

# Update of a no ISO rack
rack1.Update(module = 1.21, transverse_pressure_angle = 21/180.*npy.pi, coeff_gear_addendum = 0.9,
               coeff_gear_dedendum = 1.2, coeff_root_radius = 0.2, coeff_circular_tooth_thickness = 0.5)
check,_=rack1.CheckRackViable()
if check:
    RackElem=rack1.Plot(number_pattern = 10)
