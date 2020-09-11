import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
#from interval import interval

import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

###############################################################
#input for a specified gear mesh (modification of the modulus)
list_cd=[[0.117,0.16]]
list_gear_set=[[(0,1)]]
list_speed={0:[1000*npy.pi/30,1500*npy.pi/30],1:[4100*npy.pi/30,
               4300*npy.pi/30]}
list_rack={
           0:{'name':'Catalogue_A','module':[1.4*1e-3,1.4*1e-3]},
           1:{'name':'Catalogue_A','module':[1.5*1e-3,1.5*1e-3]}, 
           2:{'name':'Catalogue_A','module':[1.6*1e-3,1.6*1e-3]},
#           3:{'name':'Catalogue_A','module':[2.25*1e-3,1*2.25e-3]},
#           4:{'name':'Catalogue_A','module':[2.5*1e-3,2.5*1e-3]}
           }
list_cycle={0:1e12}

dict_pandas = {'cd': [], 'B': [], 'Z': [], 'df': [], 'ratio': [], 'tq': [], 
               'modulus': [], 'rca': []}

delta_speed = 10
compt = 1
for cd_min in npy.arange(0.07, 0.25, 0.01):
    list_cd=[[cd_min, cd_min + 0.05]]
    for min_speed in npy.arange(100, 310, 10):
        list_speed={0:[min_speed, min_speed + delta_speed],1:[200, 200 + delta_speed]}
        for num_rack, rack in list_rack.items():
            list_rack_choice={0:[num_rack],1:[num_rack]}
            for tq in npy.arange(100, 410, 10):
                list_torque={0: tq, 1: 'output'}
                GA=meshes_opt.MeshAssemblyOptimizer(Z={},
                                   connections=list_gear_set,
                                   gear_speeds=list_speed,
                                   center_distances=list_cd,
                                   rack_list=list_rack,
                                   rack_choice=list_rack_choice,
                                   torques=list_torque,
                                   cycles=list_cycle)
                GA.OptimizeCD(nb_sol=1,verbose=False)
    #            print('Nombre de solutions convergés:',len(GA.solutions))
                if len(GA.solutions) > 0:
                    solution=GA.solutions[-1]
                    dict_pandas['cd'].append(solution.center_distance[0])
                    dict_pandas['B'].append(solution.mesh_assembly[0].gear_width[0])
                    dict_pandas['Z'].append(solution.mesh_assembly[0].Z[0])
                    dict_pandas['df'].append(solution.mesh_assembly[0].DF[0][0])
                    dict_pandas['rca'].append(solution.mesh_assembly[0].radial_contact_ratio[0])
                    ratio = solution.mesh_assembly[0].Z[0]/solution.mesh_assembly[0].Z[1]
                    dict_pandas['ratio'].append(ratio)
                    tq = list(solution.mesh_assembly[0].torque.values())[0]
                    dict_pandas['tq'].append(tq)
                    dict_pandas['modulus'].append(solution.mesh_assembly[0].meshes[0].rack.module)
                    
                    dict_pandas['cd'].append(solution.center_distance[0])
                    dict_pandas['B'].append(solution.mesh_assembly[0].gear_width[1])
                    dict_pandas['Z'].append(solution.mesh_assembly[0].Z[1])
                    dict_pandas['df'].append(solution.mesh_assembly[0].DF[0][1])
                    dict_pandas['rca'].append(solution.mesh_assembly[0].radial_contact_ratio[0])
                    dict_pandas['ratio'].append(1/ratio)
                    dict_pandas['tq'].append(tq/ratio)
                    dict_pandas['modulus'].append(solution.mesh_assembly[0].meshes[1].rack.module)
                    
                    compt += 1
print('Number of solutions: {}'.format(compt))
            
df = pd.DataFrame(dict_pandas)
scatter_matrix(df)
plt.show()

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_trisurf(df.tq, df.cd, df.B, cmap=cm.jet, linewidth=0.2)
plt.show()

sol = df.loc[:,['tq', 'cd', 'ratio', 'modulus', 'rca']].values
poly = PolynomialFeatures(degree=4)
x=poly.fit_transform(sol)

y = df.loc[:,['B']].values
reg = LinearRegression().fit(x, y)
print(reg.score(x, y))
print(reg.coef_)
print(reg.intercept_ )

array2D=poly.fit_transform(sol)
y_interp = reg.predict(array2D)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_trisurf(sol[:,0], sol[:,1], y_interp[:,0], cmap=cm.jet, linewidth=0.1)
#ax.scatter3D(df.tq, df.cd, df.B, cmap='Greens', linewidth=0.2)
plt.show()

plt.figure()
plt.plot(df.B, 1-df.B/(y_interp[:,0]), 'r')

max(df.B-y_interp[:,0])
min(df.B-y_interp[:,0])