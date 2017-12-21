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



##Initialisation
list_bounds=[(30,30),(56,56),(60,70),(0.01,0.5),(0.3,0.3),(-1,1),(-1,1),(20,20),(200,200)]
A1=AssemblyGearOptimize1(list_bounds)
O1=Optimization(A1)
O1.Optimize()
AG1=AssemblyGear(*list(O1.solution[-1]))

#xsol=list(npy.array([30.,50.,115.41878235,0.38523225,-0.94639435,0.26922009,0.,20.]))
#xsol=list(npy.array([3.00000000e+01,5.10000000e+01,9.74954025e+01,4.25119842e-01,3.00000000e-01,9.95153575e-01,7.35600926e-02,2.00000000e+01,2.00000000e+02]))
#AG1=AssemblyGear(*xsol)

T1=Trace()
T1.Rack(AG1.Rack1,5,'rack1.html')
T1.Rack(AG1.Rack2,5,'rack2.html')
T1.RackGenere(AG1.Gear1,AG1.Z1,'Cremaillere_Z1.html')
T1.RackGenere(AG1.Gear2,AG1.Z2,'Cremaillere_Z2.html')
T1.DetailDent(AG1,AG1.Gear1,AG1.DB1,AG1.DF1,'Creation_Dent1.html')
T1.DetailDent(AG1,AG1.Gear2,AG1.DB2,AG1.DF2,'Creation_Dent2.html')
T1.GearAssembly(AG1.Gear1,AG1.Gear2,AG1,AG1.DF1,AG1.DF2,(0,0),(AG1.center_distance,0),'Gear.html')

