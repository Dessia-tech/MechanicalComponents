import mechanical_components.gears_assembly as gears
import numpy as npy

##########################
# script pour AGB Safran
##########################

#definition data
erreur=0.02
list_cd=[[148.6*1e-3,148.6*1e-3+0.01],[130.1*1e-3,130.1*1e-3+0.01],[137*1e-3,137*1e-3+0.01],
         [135*1e-3,135*1e-3+0.01],[142*1e-3,142*1e-3+0.01],[145.4*1e-3,145.4*1e-3+0.01]]
list_gear_set=[(2,4),(4,6),(6,7),(7,0),(0,3),(3,5)]
list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[1000*npy.pi/30,30000*npy.pi/30]}
list_rack={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={2:0,4:0,6:0,7:0,0:0,3:0,5:0}
list_helix_angle={2:[0,0]}
list_material={2:gears.hardened_alloy_steel,4:gears.hardened_alloy_steel}
list_torque={2:1000,6:500,5:300}
list_cycle={2:1e6}

GA_wizard=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,
                                            center_distance=list_cd,rack_list=list_rack,
                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,
                                            material=list_material,torque=list_torque,cycle=list_cycle)
GA_wizard.SearchCenterLine(nb_sol=1)
sol=GA_wizard.solutions_search[-1]

#sol.SVGExport('name.txt',{7:[0,0],2:[0.2,0],6:[0.6,0]})
#v=sol.VolumeModel(name='3D')
sol.FreeCADExport('Gears')