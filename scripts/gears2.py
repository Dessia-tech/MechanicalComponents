import mechanical_components.gears_assembly as gears
import numpy as npy

#definition data
erreur=0.02
list_cd=[[0.1,0.101]]
list_gear_set=[(0,1)]
list_speed={0:[1,1000],1:[1,1000]}
list_Z={0:[57,57],1:[106,106]}
list_rack={0:{'name':'Catalogue_A','module':[1*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0:0,1:0}
list_helix_angle={0:[0,0]}
list_material={0:gears.hardened_alloy_steel,1:gears.hardened_alloy_steel}
list_torque={0:200}
list_cycle={0:1e6}

GA_wizard=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,
                                            center_distance=list_cd,rack_list=list_rack,
                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,
                                            Z=list_Z,material=list_material,torque=list_torque,cycle=list_cycle)
GA_wizard.Optimize()
sol=GA_wizard.solutions[-1]

#sol.SVGExport('name.txt',{0:[0,0]})
#v=sol.VolumeModel()
sol.FreeCADExport('Gear')