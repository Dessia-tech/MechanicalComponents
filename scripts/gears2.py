import mechanical_components.gears_assembly as gears
import numpy as npy

#definition data
erreur=0.02
list_cd=[[0.1,0.11]]
list_gear_set=[(0,1)]
list_speed={0:[1,1000],1:[1,1000]}
list_Z={0:[57,57],1:[106,106]}
list_rack={0:{'name':'Catalogue_A','module':[1*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0:0,1:0}
list_helix_angle={0:[0,0]}



GA_wizard=gears.GearAssemblyOptimizerWizard(gear_set=list_gear_set,gear_speed=list_speed,
                                            center_distance=list_cd,rack_list=list_rack,
                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,Z=list_Z)
GA_wizard.Optimize()
sol=GA_wizard.solutions[-1]

sol.SVGGearSet('name.txt',{0:[0,0]})
#v=sol.VolumeModel('3D',{7:[0,0],2:[0.2,0],6:[0.6,0]})
#v.FreeCADExport('python','Gears','/usr/lib/freecad/lib')