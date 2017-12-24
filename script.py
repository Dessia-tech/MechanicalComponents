#import numpy as npy
#import volmdlr as vm
#import volmdlr.primitives3D as primitives3D
#import volmdlr.primitives2D as primitives2D
#import math
#from scipy.linalg import norm,solve,LinAlgError
#from scipy.optimize import *
#from scipy.interpolate import splprep, splev
#from sympy import *
#import itertools

#import LibSvg
import Gears

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

M1=Gears.MasterAssemblyGear(**dico_ref)
M1.Optimize()

####Initialisation
#list_bounds=[(13,13),(54,54),(50,80),(0.01,0.5),(0.3,0.3),(-1,1),(-1,1),(20,20),(200,200)]
#A1=AssemblyGearOptimize1(list_bounds)
#O1=Optimization(A1)
#O1.Optimize()
##AG1=AssemblyGear(*list(O1.solution[-1]))
#
##xsol=list(npy.array([30.,50.,115.41878235,0.38523225,-0.94639435,0.26922009,0.,20.]))
##xsol=list(npy.array([3.00000000e+01,5.10000000e+01,9.74954025e+01,4.25119842e-01,3.00000000e-01,9.95153575e-01,7.35600926e-02,2.00000000e+01,2.00000000e+02]))
##AG1=AssemblyGear(*xsol)
#
#T1=Trace()
#T1.Rack(AG1.Rack1,5,'rack1.html')
#T1.Rack(AG1.Rack2,5,'rack2.html')
#T1.RackGenere(AG1.Gear1,AG1.Z1,'Cremaillere_Z1.html')
#T1.RackGenere(AG1.Gear2,AG1.Z2,'Cremaillere_Z2.html')
#T1.DetailDent(AG1,AG1.Gear1,AG1.DB1,AG1.DF1,'Creation_Dent1.html')
#T1.DetailDent(AG1,AG1.Gear2,AG1.DB2,AG1.DF2,'Creation_Dent2.html')
#T1.GearAssembly(AG1.Gear1,AG1.Gear2,AG1,AG1.DF1,AG1.DF2,(0,0),(AG1.center_distance,0),'Gear.html')

