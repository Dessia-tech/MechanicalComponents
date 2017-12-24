import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *
from scipy.interpolate import splprep, splev
from sympy import *
import itertools

from LibSvg import *
from Gear import *

dico_ref={'ratio':{'nom':0.23,'min':None,'max':None,'err':0.05,'prog':None,'prem':None,'pgcd':None},
                     'Z1':{'nom':None,'min':13,'max':20},
                     'Z2':{'nom':None,'min':45,'max':70},
                     'center_distance':{'nom':None,'min':60,'max':70,'err':None,'prog':None},
                     'transverse_pressure_angle':{'nom':None,'min':0.1,'max':0.5,'err':None,'prog':None},
                     'helix_angle':{'nom':None,'min':0.3,'max':0.3,'err':None,'prog':None},
                     'coefficient_profile_shift1':{'nom':None,'min':None,'max':None,'err':None,'prog':None},
                     'coefficient_profile_shift2':{'nom':None,'min':None,'max':None,'err':None,'prog':None},
                     'gear_width':{'nom':None,'min':20,'max':20,'err':None,'prog':None},
                     'maximum_torque':{'nom':None,'min':200,'max':200,'err':None,'prog':None}
                     }

M1=MasterAssemblyGear(**dico_ref)
M1.Optimize()

### Developpement perso
#list_bounds=[(13,13),(54,54),(50,80),(0.01,0.5),(0.3,0.3),(-1,1),(-1,1),(20,20),(200,200)]
#A1=AssemblyGearOptimize1(list_bounds)
#O1=Optimization(A1)
#O1.Optimize()
##AG1=AssemblyGear(*list(O1.solution[-1]))
##xsol=list(npy.array([30.,50.,115.41878235,0.38523225,-0.94639435,0.26922009,0.,20.]))
##xsol=list(npy.array([3.00000000e+01,5.10000000e+01,9.74954025e+01,4.25119842e-01,3.00000000e-01,9.95153575e-01,7.35600926e-02,2.00000000e+01,2.00000000e+02]))
##AG1=AssemblyGear(*xsol)

### Sorties graphiques SVG
for i,j in enumerate(M1.solution):
    AG1=AssemblyGear(*j)
    T1=Trace()
    T1.Rack(AG1.Rack1,5,'rack1-s{}.html'.format(str(i)))
    T1.Rack(AG1.Rack2,5,'rack2-s{}.html'.format(str(i)))
    T1.RackGenere(AG1.Gear1,AG1.Z1,'Cremaillere_Z1-s{}.html'.format(str(i)))
    T1.RackGenere(AG1.Gear2,AG1.Z2,'Cremaillere_Z2-s{}.html'.format(str(i)))
    T1.DetailDent(AG1,AG1.Gear1,AG1.DB1,AG1.DF1,'Creation_Dent1-s{}.html'.format(str(i)))
    T1.DetailDent(AG1,AG1.Gear2,AG1.DB2,AG1.DF2,'Creation_Dent2-s{}.html'.format(str(i)))
    T1.GearAssembly(AG1.Gear1,AG1.Gear2,AG1,AG1.DF1,AG1.DF2,(0,0),(AG1.center_distance,0),'Gear-s{}.html'.format(str(i)))

### Sorties data
dico_export={}
for i,j in enumerate(M1.solution):
    dico_export['s'+str(i)]={}
    AG1=AssemblyGear(*j)
    dico_export['s'+str(i)]['axial_contact_ratio']=AG1.axial_contact_ratio
    dico_export['s'+str(i)]['axial_pitch']=AG1.axial_pitch
    dico_export['s'+str(i)]['center_distance']=AG1.center_distance
    dico_export['s'+str(i)]['circular_tooth_thickness1']=AG1.circular_tooth_thickness1
    dico_export['s'+str(i)]['circular_tooth_thickness2']=AG1.circular_tooth_thickness2
    dico_export['s'+str(i)]['coefficient_profile_shift1']=AG1.coefficient_profile_shift1
    dico_export['s'+str(i)]['coefficient_profile_shift2']=AG1.coefficient_profile_shift2
    dico_export['s'+str(i)]['DB1']=AG1.DB1
    dico_export['s'+str(i)]['DB2']=AG1.DB2
    dico_export['s'+str(i)]['DF1']=AG1.DF1
    dico_export['s'+str(i)]['DF2']=AG1.DF2
    dico_export['s'+str(i)]['gear_width']=AG1.gear_width
    dico_export['s'+str(i)]['helix_angle']=AG1.helix_angle
    dico_export['s'+str(i)]['maximum_torque']=AG1.maximum_torque
    dico_export['s'+str(i)]['normal_circular_pitch']=AG1.normal_circular_pitch
    dico_export['s'+str(i)]['normal_load']=AG1.normal_load
    dico_export['s'+str(i)]['normal_pressure_angle']=AG1.normal_pressure_angle
    dico_export['s'+str(i)]['pitch_diameter_factory1']=AG1.pitch_diameter_factory1
    dico_export['s'+str(i)]['pitch_diameter_factory2']=AG1.pitch_diameter_factory2
    dico_export['s'+str(i)]['radial_contact_ratio']=AG1.radial_contact_ratio
    dico_export['s'+str(i)]['radial_load']=AG1.radial_load
    AG1.SigmaLewis()
    dico_export['s'+str(i)]['sigma_lewis_maximum1']=AG1.sigma_lewis_maximum1
    dico_export['s'+str(i)]['sigma_lewis_maximum2']=AG1.sigma_lewis_maximum2
    dico_export['s'+str(i)]['space_width1']=AG1.space_width1
    dico_export['s'+str(i)]['space_width2']=AG1.space_width2
    dico_export['s'+str(i)]['tangential_load']=AG1.tangential_load
    dico_export['s'+str(i)]['total_contact_ratio']=AG1.total_contact_ratio
    dico_export['s'+str(i)]['transverse_base_pitch']=AG1.transverse_base_pitch
    dico_export['s'+str(i)]['transverse_pressure_angle']=AG1.transverse_pressure_angle
    dico_export['s'+str(i)]['transverse_pressure_angle_rack_T1']=AG1.transverse_pressure_angle_rack_T1
    dico_export['s'+str(i)]['transverse_pressure_angle_rack_T2']=AG1.transverse_pressure_angle_rack_T2
    dico_export['s'+str(i)]['transverse_radial_pitch']=AG1.transverse_radial_pitch
    dico_export['s'+str(i)]['transverse_radial_pitch_rack1']=AG1.transverse_radial_pitch_rack1
    dico_export['s'+str(i)]['transverse_radial_pitch_rack2']=AG1.transverse_radial_pitch_rack2
    dico_export['s'+str(i)]['Z1']=AG1.Z1
    dico_export['s'+str(i)]['Z2']=AG1.Z2
    dico_export['s'+str(i)]['outside_diameter1']=AG1.Gear1.outside_diameter
    dico_export['s'+str(i)]['outside_diameter2']=AG1.Gear2.outside_diameter
    dico_export['s'+str(i)]['root_diameter_active1']=AG1.Gear1.root_diameter_active
    dico_export['s'+str(i)]['root_diameter_active2']=AG1.Gear2.root_diameter_active
    dico_export['s'+str(i)]['root_diameter1']=AG1.Gear1.root_diameter
    dico_export['s'+str(i)]['root_diameter2']=AG1.Gear2.root_diameter
    
    
    
    