#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:13:01 2018

@author: steven
"""

#import itertools as it
from mechanical_components.meshes import MeshAssembly, hardened_alloy_steel
import numpy as npy
import itertools
import networkx as nx
import powertransmission.tools as tools
from scipy.optimize import minimize
import dectree

class ContinuousMeshesAssemblyOptimizer:
    """
    Gear mesh assembly optimizer
    
    :param Z: Dictionary define the tooth number of each mesh {node1: Z1, node2: Z2 ...}
    :param center_distance: List of two elements define gear mesh center-distance [[mesh1_centerdistance_min, mesh1_centerdistance_max], [mesh2_centerdistance_min, mesh2_centerdistance_max] ...] with mesh1 the mesh between node1 and node2 ...    
    :param connections: List of tuple define gear mesh connection [(node1,node2), (node2,node3), (node2,node4)]
    :param transverse_pressure_angle: List of two elements define the transversal pressure angle interval for each mesh [[mesh1_transversepressure_min, mesh1_transversepressure_max], [mesh2_transversepressure_min, mesh2_transversepressure_max] ...]
    :param coefficient_profile_shift: Dictionary defining the minimum and maximum coefficient profile shift for each node {node1: [node1_coeffshift_min,node1_coeffshift_max],node2: [node2_coeffshift_min,node2_coeffshift_max] ...}
    :param gear_graph: NetwokX gear graph connection
    :param rack_list: Dictionary define all admissible rack {rack1:mechanical_components.meshes.Rack, rack2:mechanical_components.meshes.Rack ...}
    :param rack_choice: Dictionary assign for each gear mesh a list of acceptable rack {node1:rack1, node2:rack1, node3:rack4 ...}
    :param list_gear: List of gear mesh in connections [node1, node2, node3, node4]
    :param material: Dictionary defining material for each gear mesh {node1:mechanical_components.meshes.Material, node2:mechanical_components.meshes.Material ...}
    :param torque: Dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
    :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}
    :param safety_factor: Safety factor used for the ISO design
    
    >>> input = {'Z': {0: 53, 1: 180},
                 'center_distance': [[0.117, 0.117]],
                 'coefficient_profile_shift': {0: [-0.8, 0.8], 1: [-0.8, 0.8]},
                 'connections': [(1, 0)],
                 'cycle': {1: 1000000.0},
                 'list_gear': [1, 0],
                 'rack_choice': {0: 0, 1: 0},
                 'rack_list': {0: {'coeff_circular_tooth_thickness': [0.5, 0.5],
                                   'coeff_gear_addendum': [1, 1],
                                   'coeff_gear_dedendum': [1.25, 1.25],
                                   'coeff_root_radius': [0.38, 0.38],
                                   'module': [0.001, 0.003],
                                   'name': 'Optim_Module',
                                   'transverse_pressure_angle_rack': [0.3490658503988659,
                                                                      0.3490658503988659]}},
                 'safety_factor': 1,
                 'torque': {0: 'output', 1: 100},
                 'transverse_pressure_angle': [[0.2617993877991494, 0.5235987755982988]]}
    >>> cmao1 = ContinuousMeshesAssemblyOptimizer(**input)
    """
    def __init__(self, Z, center_distance, connections, transverse_pressure_angle,
                 coefficient_profile_shift, gear_graph,
                 rack_list, rack_choice, list_gear, material,torque,
                 cycle, safety_factor):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear
        Bounds=list(center_distance)
        Bounds.extend(transverse_pressure_angle)
        # TODO: Check if this is a bug
        tp=len(coefficient_profile_shift.keys())
        for i in self.list_gear:
            Bounds.append(coefficient_profile_shift[i])
        self.solutions=[]
        
        self.xi={'Z': Z, 'connections': connections,'gear_graph':gear_graph,'list_gear':list_gear,
                 'material':material,'torque':torque,'cycle':cycle,'safety_factor':safety_factor}
        self.xj,self.dict_xu=self._init()
        self.xt=dict(list(self.xi.items())+list(self.xj.items()))
        
        for k,v in self.dict_xu:
            Bounds.append(self.rack_list[v][k])
        self.Bounds=npy.array(Bounds)
        
        #xk partie non optimisée du vecteur x
        #xu partie optimisée du vecteur x
        
        self.MeshAssembly = MeshAssembly(**self.xt)
        
    # TODO: Pourquoi un deuxième Init ???
    def _init(self):
        xj={'center_distance':self._init_list(self.center_distance),
                 'transverse_pressure_angle':self._init_list(self.transverse_pressure_angle),
                 'coefficient_profile_shift':self._init_item(self.coefficient_profile_shift),
                 'transverse_pressure_angle_rack':[],'coeff_gear_addendum':[],
                 'coeff_gear_dedendum':[],'coeff_root_radius':[],'coeff_circular_tooth_thickness':[]}
        dict_xu=[]
        x_list={}
        for i in list(set(list(self.rack_choice.values()))):
            x_list[i]={'transverse_pressure_angle_rack':[],
                  'coeff_gear_addendum':[],
                  'coeff_gear_dedendum':[],
                  'coeff_root_radius':[],
                  'coeff_circular_tooth_thickness':[]}
            for k,v in self.rack_list[i].items():
                if k not in ['type','name','module']:
                    if v[0]==v[1]:
                        x_list[i][k]=v[0]
                    else:
                        x_list[i][k]=(v[1]-v[0])*float(npy.random.random(1))+v[0]
                        dict_xu.append((k,i))
#        for k in sorted(list(self.rack_choice.keys())):
        for k in self.list_gear:
            v=self.rack_choice[k]
            for k2,v2 in self.rack_list[v].items():
                if k2 not in ['type','name','module']:
                    xj[k2].append(x_list[v][k2])
        return xj,dict_xu
    
    def _convert_xj2Xu(self,xj):
        sol=xj['center_distance'].copy()
        sol.extend(xj['transverse_pressure_angle'])
        # TODO: Check if this is a bug
        tp=len(xj['coefficient_profile_shift'].keys())
        for key in self.list_gear:
            sol.append(xj['coefficient_profile_shift'][key])
        for k,v in self.dict_xu:
            for k2 in sorted(list(self.rack_choice.keys())):
                v2=self.rack_choice[k2]
                if v2==v:
                    ind=k2
            sol.append(xj[k][ind])
        return npy.array(sol)

    def _convert_Xu2xj(self,X):
        X=list(X)
        tp1=len(self.center_distance)
        tp2=len(self.transverse_pressure_angle)
        tp3=len(self.coefficient_profile_shift.keys())
        xj=self.xj
        xj['center_distance']=X[0:tp1]
        xj['transverse_pressure_angle']=X[tp1:tp1+tp2]
        xj['coefficient_profile_shift']={}
        for i in self.list_gear:
            j=self.list_gear.index(i)
            xj['coefficient_profile_shift'][i]=X[tp1+tp2+j]
#        xj['transverse_pressure_angle_rack']=[0]*len(self.rack_choice.keys())
#        xj['coeff_gear_addendum']=[0]*len(self.rack_choice.keys())
#        xj['coeff_gear_dedendum']=[0]*len(self.rack_choice.keys())
#        xj['coeff_root_radius']=[0]*len(self.rack_choice.keys())
#        xj['coeff_circular_tooth_thickness']=[0]*len(self.rack_choice.keys())
        for i,(k,v) in enumerate(self.dict_xu):
            for k2 in sorted(list(self.rack_choice.keys())):
                v2=self.rack_choice[k2]
                if v2==v:
                    ind=k2
                    xj[k][ind]=X[tp1+tp2+tp3+i]
        return xj

    def _init_list(self,list):
        list1=[]
        for li in list:
            list1.append((li[1]-li[0])*float(npy.random.random(1))+li[0])
        return list1
    
    def _init_item(self,dico):
        dict1={}
        for k1,v1 in dico.items():
            dict1[k1]=(v1[1]-v1[0])*float(npy.random.random(1))+v1[0]
        return dict1
    
    def Update(self,xj):
        self.xt=dict(list(self.xi.items())+list(xj.items()))
        self.MeshAssembly.Update(**self.xt)
        return xj
    
    # TODO: implémenter les critères comme des calculs dans le MeshAssembly
    # Les appeller ensuite. Pour un exemple regarder InternalInterference du model3D PWT.
    def Fineq(self,X):
        xj=self._convert_Xu2xj(X)
        xj=self.Update(xj)
        ineq=[]
        #jeu inter-denture
        for lb in self.MeshAssembly.linear_backlash:
            ineq.append(lb)
            ineq.append((4e-4)-lb)
        #contrainte géométrique
        for ne,gs in enumerate(self.MeshAssembly.connections):
            dia1=self.MeshAssembly.meshes[ne][gs[0]].root_diameter_active
            dia2=self.MeshAssembly.meshes[ne][gs[1]].root_diameter_active
            de1=self.MeshAssembly.meshes[ne][gs[0]].outside_diameter
            de2=self.MeshAssembly.meshes[ne][gs[1]].outside_diameter
            cd=self.MeshAssembly.center_distance[ne]
            ineq.append(cd-(de1/2+dia2/2))
            ineq.append(cd-(de2/2+dia1/2))
            oaa1=self.MeshAssembly.meshes[ne][gs[0]].outside_active_angle
            oaa2=self.MeshAssembly.meshes[ne][gs[1]].outside_active_angle
            ineq.append(oaa1)
            ineq.append(oaa2)
            df1=self.MeshAssembly.DF[ne][gs[0]]
            df2=self.MeshAssembly.DF[ne][gs[1]]
            db1=self.MeshAssembly.meshes[ne][gs[0]].DB
            db2=self.MeshAssembly.meshes[ne][gs[1]].DB
            ineq.append(df1-db1)
            ineq.append(df2-db2)
        #contrainte sur le RCA
        for ne,gs in enumerate(self.MeshAssembly.connections):
            rca=self.MeshAssembly.radial_contact_ratio[ne]
            ineq.append(rca-1)
        #contrainte sur le module
        for ne,gs in enumerate(self.MeshAssembly.connections):
            for g in gs:
                mo=self.MeshAssembly.meshes[ne][g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]<lm[1]:
                ineq.append(mo-list_module[0])
                ineq.append(list_module[1]-mo)
        return ineq
    
    def Feq(self,X):
        x=self._convert_Xu2xj(X)
        x=self.Update(x)
        eq=[]
#        for ng in range(self.MeshAssembly.gear_graph.number_of_nodes()):
        for ng in self.list_gear:
            nel=list(self.MeshAssembly.gear_graph.edges(ng))
            if len(nel)>1:
                list_db=[]
                for ne in nel:
                    if (ne[0],ne[1]) in self.MeshAssembly.connections:
                        nes=self.MeshAssembly.connections.index((ne[0],ne[1]))
                    elif (ne[1],ne[0]) in self.MeshAssembly.connections:
                        nes=self.MeshAssembly.connections.index((ne[1],ne[0]))
                    list_db.append(self.MeshAssembly.meshes[nes][ng].DB)
#                    print(self.MeshAssembly.meshes[nes][ng].DB,nes,ng,list_db)
                list1=itertools.combinations(list_db,2)
                for n1,n2 in list1:
                    eq.append(n1-n2)
#        for ne,gs in enumerate(self.MeshAssembly.connections):
#            for g in gs:
#                mo=self.MeshAssembly.meshes[ne][g].rack.module
#                lm=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]==lm[1]:
#                    eq.append(mo-lm[0])
        if eq==[]:
            eq=[0]
        return eq
        
    def Objective(self,X):
        x=self._convert_Xu2xj(X)
        x=self.Update(x)
        fineq=self.Fineq(X)
        feq=self.Feq(X)
        obj=0
        #Maximisation du module pour avoir des pignons avec un faible gear_width
        for ne,gs in enumerate(self.MeshAssembly.connections):
            for g in gs:
                mo=self.MeshAssembly.meshes[ne][g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
                if list_module[0]<list_module[1]:
                    obj+=100*(list_module[1]-mo)**2
                
        #Minimisation des entraxes sur la borne inf
        for num_engr,list_cd in enumerate(self.center_distance):
            obj+=100*(list_cd[0]-self.MeshAssembly.center_distance[num_engr])**2
            
        for lb in self.MeshAssembly.linear_backlash:
            obj+=100*(lb)**2
            
        for i in fineq:
            if i < 0:
                obj+=1000*i**2
            else:
                obj+=0.000001*i
        for i in feq:
            if i<0:
                obj+=1000*i**2
            else:
                obj+=1000*i**2
        return obj
    
    def Optimize(self, verbose = False):
        """ Optimizer function
        
        >>> cmao1.Optimize(verbose = True)
        Iteration n°0 with status 0, min(fineq):-4.7435658082073395e-08, max(eq):0
        Mesh sections: {1: 0.0039733070033111731, 0: 0.0039733070033111731}
        Numbers of teeth: {1: 59, 0: 19}
        Center distances: [0.11700000000000001]
        """
        max_iter=5
        i=0
        arret=0 
        while i<max_iter and arret==0:
            xj0,dict_xu=self._init()
            xj0=self.Update(xj0)
            X0=self._convert_xj2Xu(xj0)
            if self.Feq(X0)==[0]:
                cons = {'type': 'ineq','fun' : self.Fineq}
            else:
                cons = ({'type': 'eq','fun' : self.Feq},{'type': 'ineq','fun' : self.Fineq})
            cx = minimize(self.Objective, X0, bounds=self.Bounds,constraints=cons)
            Xsol=cx.x
            xsol=self._convert_Xu2xj(Xsol)
            xsol=self.Update(xsol)
            if verbose:
                print('Iteration n°{} with status {}, min(fineq):{}, max(eq):{}'.format(i,
                      cx.status,min(self.Fineq(Xsol)),max(npy.abs(self.Feq(Xsol)))))
            if min(self.Fineq(Xsol))>-1e-5 and max(npy.abs(self.Feq(Xsol)))<1e-5:
                self.solutions.append(xsol)
                arret=1
            i=i+1


class MeshAssemblyOptimizer:
    """
    Gear mesh assembly optimizer supervisor
    
    :param connections: List of tuple define gear mesh connection [(node1,node2), (node2,node3)...]
    :param gear_speed: Dictionary defining minimum and maximum speed for each gear mesh {node1: [speed_min,speed_max], node2: [speed_min,speed_max]...}
    :param center_distance: List of two elements define gear mesh connection [[mesh1_centerdistance_min, mesh1_centerdistance_max], [mesh2_centerdistance_min, mesh2_centerdistance_max] ...] with mesh1 the mesh between node1 and node2 ...
    :param transverse_pressure_angle: List of two elements define the transversal pressure angle interval for each mesh [[mesh1_transversepressure_min, mesh1_transversepressure_max], [mesh2_transversepressure_min, mesh2_transversepressure_max] ...]
    :param helix_angle: Dictionary to define for one mesh the minimum and maximum helix angle {node2: [mesh1_helixangle_min, mesh1_helixangle_max]}
    :param gear_width: Dictionary to define for each gear mesh the minimum and maximum gear width {node1: [node1_gearwidth_min,node1_gearwidth_max] ,node2: [node2_gearwidth_min,node2_gearwidth_max] ...}
    :param frequency: List of two elements defining unacceptable frequency interval [[freq1_min,freq1_max], [freq2_min,freq2_max]...]
    :param coefficient_profile_shift: Dictionary defining the minimum and maximum coefficient profile shift for each node {node1: [node1_coeffshift_min,node1_coeffshift_max],node2: [node2_coeffshift_min,node2_coeffshift_max] ...}
    :param rack_list: Dictionary define all admissible rack {rack1:mechanical_components.meshes.Rack, rack2:mechanical_components.meshes.Rack ...}
    :param rack_choice: Dictionary assign for each gear mesh a list of acceptable rack {node1:[rack1,rack2], node2:[rack1],node3:[rack4,rack5] ...}
    :param material: Dictionary defining material for each gear mesh {node1:mechanical_components.meshes.Material, node2:mechanical_components.meshes.Material ...}
    :param torque: Dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
    :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}
    :param safety_factor: Safety factor used for the ISO design
    
    >>> list_cd=[[0.117,0.117],[0.12,0.13]]
    >>> list_gear_set=[(1,0),(0,2)]
    >>> list_speed={1:[1000*npy.pi/30.,1500*npy.pi/30.],0:[4100*npy.pi/30.,
                   4300*npy.pi/30.],2:[200*npy.pi/30.,400*npy.pi/30.]}
    >>> GA = meshes_opt.MeshAssemblyOptimizer(connections = list_gear_set,
                                gear_speed = list_speed,
                                center_distance = list_cd)
    """
    def __init__(self, connections, gear_speed, center_distance, Z=None,
                 transverse_pressure_angle=None, helix_angle=None,
                 gear_width=None, frequency=[[0,0]],
                 coefficient_profile_shift=None, rack_list=None,
                 rack_choice=None, material=None, torque=None, cycle=None,
                 safety_factor=1,verbose=False):

        list_gear=[]
        for gs in connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        nb_set=len(connections)
                    
        if transverse_pressure_angle==None:
            transverse_pressure_angle=[]
            for i in range(nb_set):
                transverse_pressure_angle.append([15/180.*npy.pi,30/180.*npy.pi])
            
        if helix_angle==None:
            helix_angle={list_gear[0]:[15/180.*npy.pi,25/180.*npy.pi]}
        
        if gear_width==None:
            gear_width={list_gear[0]:[15*1e-3,25*1e-3]}
        gw_min=npy.inf
        gw_max=-npy.inf
        for ne in gear_width.keys():
            if gear_width[ne][0]<gw_min:
                gw_min=gear_width[ne][0]
            if gear_width[ne][1]>gw_max:
                gw_max=gear_width[ne][1]
        for ne in list_gear:
            if ne not in gear_width.keys():
                gear_width[ne]=[gw_min,gw_max]
                
        if coefficient_profile_shift==None:
            coefficient_profile_shift={list_gear[0]:[-0.8,0.8]}
        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne]=[-0.8,0.8]
                
        if rack_list==None:
            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,3*1e-3],
                          'transverse_pressure_angle_rack':[20*npy.pi/180.,20*npy.pi/180.],
                          'coeff_gear_addendum':[1,1],
                          'coeff_gear_dedendum':[1.25,1.25],
                          'coeff_root_radius':[0.38,0.38],
                          'coeff_circular_tooth_thickness':[0.5,0.5]}}
                        
            
        if rack_choice==None:
            rack_choice={list_gear[0]:[list(rack_list.keys())[0]]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne]=[list(rack_list.keys())[0]]
                
        if material==None:
            material={list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne]=hardened_alloy_steel
        
        if torque==None:
            torque={list_gear[0]:100,list_gear[1]:'output'}
            
        if cycle==None:
            cycle={list_gear[0]:1e6}
        
        if Z==None:
            Z={}
        self.Z=Z
        self.connections=connections
        self.gear_speed=gear_speed
        self.frequency=frequency
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear
        self.material=material
        self.torque=torque
        self.cycle=cycle
        self.safety_factor=safety_factor

        self.nb_gear=len(list_gear)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.connections)
        self.gear_graph=gear_graph
        
        self.nb_rack=len(self.rack_list.keys())

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.connections_dfs=list(nx.dfs_edges(gear_graph,self.node_init))
        
        if self.Z=={}:
            var_Z=self.AnalyseZ()
            self.Z=var_Z

        self.plex_calcul=self.AnalyzeCombination()
        
        for i,plex in enumerate(self.plex_calcul):
            plex['gear_graph']=self.gear_graph
            plex['rack_list']=self.rack_list
            plex['list_gear']=self.list_gear
            plex['material']=self.material
            plex['torque']=self.torque
            plex['cycle']=self.cycle
            plex['connections']=self.connections
#            plex['center_distance']=self.center_distance
            plex['transverse_pressure_angle']=self.transverse_pressure_angle
            plex['coefficient_profile_shift']=self.coefficient_profile_shift
            plex['safety_factor']=safety_factor
            self.plex_calcul[i]=plex
            
        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]
        
    def AnalyseZ(self):
        """ Analyse of the minimum and maximum admissible teeth number for each gear mesh
        
        :results: dictionary define the minimum and maximum teeth number for each gear mesh {node1 : [Z1_min, Z1_max], node2 : [Z2_min, Z2_max] ...}
        
        >>> list_z = GA.AnalyseZ()
        >>> print(list_z)
        {1:[13,45], 2: [37,56]}
        """
        Z=self.Z
        for i,(engr1,engr2) in enumerate(self.connections):
            cd_min=self.center_distance[i][0]
            cd_max=self.center_distance[i][1]
            module1_min,module1_max=(npy.inf,0)
            for rack_num in self.rack_choice[engr1]:
                mod_min,mod_max=self.rack_list[rack_num]['module']
                module1_min,module1_max=(min(module1_min,mod_min),max(module1_max,mod_max))
            module2_min,module2_max=(npy.inf,0)
            for rack_num in self.rack_choice[engr2]:
                mod_min,mod_max=self.rack_list[rack_num]['module']
                module2_min,module2_max=(min(module2_min,mod_min),max(module2_max,mod_max))
            demul_min=self.gear_speed[engr1][0]/self.gear_speed[engr2][1]
            demul_max=self.gear_speed[engr1][1]/self.gear_speed[engr2][0]
            DF1_max=2*cd_max/(1+demul_min)
            Z1_max=int(DF1_max/module1_min)+1
            DF2_max=2*cd_max*demul_max/(1+demul_max)
            Z2_max=int(DF2_max/module2_min)+1
            DF1_min=2*cd_min/(1+demul_max)
            Z1_min=int(DF1_min/module1_max)-1
            DF2_min=2*cd_min*demul_min/(1+demul_min)
            Z2_min=int(DF2_min/module2_max)-1
            
            if engr1 not in Z.keys():
                Z[engr1]=[Z1_min,Z1_max]
            else:
                Z[engr1]=[max(Z1_min,Z[engr1][0]),min(Z1_max,Z[engr1][1])]
            if engr2 not in Z.keys():
                Z[engr2]=[Z2_min,Z2_max]
            else:
                Z[engr2]=[max(Z2_min,Z[engr2][0]),min(Z2_max,Z[engr2][1])]
        return Z

    def AnalyzeCombination(self, verbose = False):
        """ Analyse with decision tree all admissible configuration
        
        :results: list of all admissible solutions 
        
        >>> list_plex = GA.AnalyzeCombination()
        >>> print(list_plex[0])
        {'Z': {0: 53, 1: 180},
         'center_distance': [[0.117, 0.117]],
         'coefficient_profile_shift': {0: [-0.8, 0.8], 1: [-0.8, 0.8]},
         'connections': [(1, 0)],
         'cycle': {1: 1000000.0},
         'list_gear': [1, 0],
         'rack_choice': {0: 0, 1: 0},
         'rack_list': {0: {'coeff_circular_tooth_thickness': [0.5, 0.5],
                           'coeff_gear_addendum': [1, 1],
                           'coeff_gear_dedendum': [1.25, 1.25],
                           'coeff_root_radius': [0.38, 0.38],
                           'module': [0.001, 0.003],
                           'name': 'Optim_Module',
                           'transverse_pressure_angle_rack': [0.3490658503988659,
                                                              0.3490658503988659]}},
         'safety_factor': 1,
         'torque': {0: 'output', 1: 100},
         'transverse_pressure_angle': [[0.2617993877991494, 0.5235987755982988]]}
        """
        n1=self.node_init
        list_node=[n1]
        for (n1,n2) in self.connections_dfs:
            if n2 not in list_node:
                list_node.append(n2)
                
        np=[]
        list_gear=[]
        for engr_num in list_node:
            np.append(self.Z[engr_num][1]-self.Z[engr_num][0]+1)
            list_gear.append(npy.arange(self.Z[engr_num][0],self.Z[engr_num][1]+1))
        np.extend([self.nb_rack]*self.nb_gear)

        list_rack=list(self.rack_list.keys())
#        print(list_gear)
#        print(list_node)
        
        demul_int_min=1/9.
        demul_int_max=9
#        print(np)
        dt = dectree.RegularDecisionTree(np)
        
        incr=0
        plex_calcul=[]
        self.fonctionnel=[]
        self.fonctionnel_module=[]

        def pgcd(a,b) :
            while a%b != 0 :
                a, b = b, a%b
            return b
        while not dt.finished:
            valid=True
#            print(dt.current_node)
            Z_node = [list_gear[i][node_value] for i,node_value in enumerate(dt.current_node[:self.nb_gear])]
#            print('Znode: ', Z_node)
#            print('dt.current_depth', dt.current_depth, dt.current_node, len(dt.current_node))
            if (dt.current_depth<=(self.nb_gear-1)) and (dt.current_depth>0):
                z1=list_gear[dt.current_depth-1][dt.current_node[dt.current_depth-1]]
                z2=list_gear[dt.current_depth][dt.current_node[dt.current_depth]]
#                print(z1,z2)
                #analyse ACV engrenage 2 à 2
                if (pgcd(z1,z2)!=1):
                    valid=False
                #analyse demul interne
                if valid:
                    demul=z1/z2
#                    print(demul)
                    if (demul > demul_int_max) or (demul < demul_int_min):
                        valid=False
                #analyse ACV de l'ensemble des engrenages entre eux
#                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
#                    for n in liste_node[0:dt.current_depth]:
#                        i=liste_node.index(n)
#                        z=liste_gear[dt.current_node[i]]
#                        if pgcd(z,z2)!=1:
#                            valid=False
                #analyse des vitesses du CDC
                if valid:
                    v0_min,v0_max=self.gear_speed[list_node[0]]
                    z0=list_gear[0][dt.current_node[0]]
                    for engr_index,engr_num in enumerate(list_node[1:dt.current_depth+1]):
                        engr_index +=1
                        if engr_num in self.gear_speed.keys():
                            z=list_gear[engr_index][dt.current_node[engr_index]]
                            demul=z0/z
                            vp_min=self.gear_speed[engr_num][0]/demul
                            vp_max=self.gear_speed[engr_num][1]/demul
                            v0_min=max(v0_min,vp_min)
                            v0_max=min(v0_max,vp_max)
                            if (v0_min>v0_max):
                                valid=False
                                break
                #analyse frequence
                if valid:
                    for freq in self.frequency:
                        zm=list_gear[0][dt.current_node[0]]
                        for engr_num,engr_ind in enumerate(dt.current_node):
                            z=list_gear[engr_num][dt.current_node[engr_num]]
                            f_min=(2*npy.pi*v0_min*zm/z)/z
                            f_max=(2*npy.pi*v0_max*zm/z)/z
                            if (max(f_min,f_max)>freq[0]) and (min(f_min,f_max)<freq[1]):
                                valid=False
                                break

            elif (dt.current_depth>(self.nb_gear-1)):
                # Analyse de la faisabilité des cremailleres
                rack_num=list_rack[dt.current_node[-1]]
                rack_pos=list_node[dt.current_depth-self.nb_gear]
                if rack_num not in self.rack_choice[rack_pos]:
                    valid=False
                    
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)):
                 
                """
                Analyse de la viabilité module/entraxe et construction 
                d'une fonctionnelle mesurant l'écart entre 
                les entraxes estimées et les entraxes mini spécifiées
                """
                liste_DF_min={}
                module_minmax={}
                module_inf,module_sup=(0,npy.inf)
                for tree_pos,tree_val in enumerate(dt.current_node[0:self.nb_gear]):
                    engr_num=list_node[tree_pos]
                    rack_num=list_rack[dt.current_node[tree_pos+self.nb_gear]]
                    module_minmax[engr_num]=self.rack_list[rack_num]['module']
                    module_inf,module_sup=(max(module_inf,module_minmax[engr_num][0]),min(module_sup,module_minmax[engr_num][1]))
                for tree_pos,tree_val in enumerate(dt.current_node[0:self.nb_gear]):
                    z=list_gear[tree_pos][tree_val]
                    engr_num=list_node[tree_pos]
                    liste_DF_min[engr_num]=z*module_inf
                liste_pente_cd_module=[]
                for set_num,(eng1,eng2) in enumerate(self.connections):
                    cd_min=(liste_DF_min[eng1]+liste_DF_min[eng2])/2.
                    liste_pente_cd_module.append(cd_min/module_inf)
                cd_minmax_nv=[]
                module_optimal=0
                for set_num,cd in enumerate(self.center_distance):
                    module_optimal=max(module_optimal,cd[0]/liste_pente_cd_module[set_num])
                if module_optimal>module_sup:
                    valid=False
                module_optimal=max(module_optimal,module_inf)
                fonctionnel=0
                for set_num,cd in enumerate(self.center_distance):
                    cd_optimal=liste_pente_cd_module[set_num]*module_optimal
                    if (cd_optimal)>(cd[1]):
                        valid=False
                        break
                    else:
                        fonctionnel+=(cd_optimal-cd[0])**2
                        cd_minmax_nv.append([cd_optimal,min(cd[1],cd_optimal*1.2)])
#            
#                # analyse coherence DB, demul et angle de pression de la cremaillere
#                for set_num,(engr1,engr2) in enumerate(self.connections):
#                    engr1_pos=liste_node.index(engr1)
#                    rack1_num=liste_rack[dt.current_node[self.nb_gear+engr1_pos]]
#                    transverse_pressure_angle=self.rack_list[rack1_num]['transverse_pressure_angle_rack']
#                    Z1=liste_gear[engr1_pos][dt.current_node[engr1_pos]]
#                    DB1_min=npy.cos(transverse_pressure_angle[0])*Z1*module_optimal
#                    DB1_max=npy.cos(transverse_pressure_angle[1])*Z1*module_sup
#                    engr2_pos=liste_node.index(engr2)
#                    rack2_num=liste_rack[dt.current_node[self.nb_gear+engr2_pos]]
#                    transverse_pressure_angle=self.rack_list[rack2_num]['transverse_pressure_angle_rack']
#                    Z2=liste_gear[engr2_pos][dt.current_node[engr2_pos]]
#                    DB2_min=npy.cos(transverse_pressure_angle[0])*Z2*module_optimal
#                    DB2_max=npy.cos(transverse_pressure_angle[1])*Z2*module_sup
#                    if ((Z1/Z2)<(DB1_min/DB2_max)) or ((Z1/Z2)>(DB1_max/DB2_min)):
#                        valid=False
#                        break
        
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)) & (valid==True):
                gear={}
                rack={}
                for n in list_node:
                    i=list_node.index(n)
                    gear[n]=list_gear[i][dt.current_node[i]]
                    rack[n]=list_rack[dt.current_node[i+self.nb_gear]]
                Export={}
                Export['Z']=gear
                Export['rack_choice']=rack
                Export['center_distance']=cd_minmax_nv
                Export['dw'] = v0_max - v0_min
                self.fonctionnel.append(fonctionnel)
                self.fonctionnel_module.append(module_optimal)
                plex_calcul.append(Export)
                incr+=1
#            print('valid sent to dt', valid)
            dt.NextNode(valid)
        if incr>1:
            if verbose:
                print('Number of combination found: {}'.format(incr))
        else:
            if verbose:
                print('No teeth combination found: increase center distances')
        return plex_calcul


    def Optimize(self,nb_sol=1, list_sol=None, verbose=False):
        """ Gear mesh assembly optimization for given plex configuration
        
        :param nb_sol: number of solution desired, if list_sol = None we take nb_sol = len(plex) the optimize function analyse all the plex possibility
        :param list_sol: list of number define the plex configuration [2,4,5]
        
        >>> GA.Optimize(list_sol=[1,2,3,4], verbose=True)
        """
        compt_nb_sol=0
        if list_sol==None:
            liste_plex=self.plex_calcul
            # Using all combinatoric
            # Sorting data, delta w max first
            dw = [s['dw'] for s in liste_plex]
            liste_plex2 = []
            for i in npy.argsort(dw)[::-1]:
                del liste_plex[i]['dw']
                liste_plex2.append(liste_plex[i])
            liste_plex = liste_plex2
            
        else:
            liste_plex=[]
            for ind_plex in list_sol:
                del self.plex_calcul[ind_plex]['dw']
                liste_plex.append(self.plex_calcul[ind_plex])
            nb_sol=len(liste_plex)
        
        for plex in liste_plex:
            
            ga=ContinuousMeshesAssemblyOptimizer(**plex)
            try:
                ga.Optimize(verbose)
            except ValueError:
                print('Convergence problem')
            if len(ga.solutions)>0:
                xsol=ga.solutions[-1]
                xt=dict(list(ga.xi.items())+list(xsol.items()))
                self.solutions.append(MeshAssembly(**xt))
                compt_nb_sol+=1
                if verbose:
                    print('Mesh sections: {}'.format(self.solutions[-1].gear_width))
                    print('Numbers of teeth: {}'.format(plex['Z']))
                    print('Center distances: {}'.format(self.solutions[-1].center_distance))
                if compt_nb_sol==nb_sol:
                    break

    def SearchOptimumCD(self, nb_sol = 1, verbose = False,
                        progress_callback = lambda x:x):
        """ Gear mesh assembly optimization of the nearest solution with the specifications
        
            * near the minimum (cd1_min, cd2_min ...) of center-distance list [[cd1_min, cd1_max], [cd2_min, cd2_max] ...]
            * maximum module of the specified rack
        
        :param nb_sol: number of solution desired, if -1 we analyse all the possibilities
        
        >>> GA.SearchOptimumCD(nb_sol=1, verbose=True)
        """
        list_fonctionnel=npy.array(self.fonctionnel)
        
        #En cours de construction
        list_fonctionnel_module=npy.array(self.fonctionnel_module)
        sort_fonct_module=npy.argsort(list_fonctionnel_module)
        compt_fonct_mod=0
        fonct_entraxe=[]
        for ind_plex in sort_fonct_module[::-1]:
            fonct_entraxe.append(list_fonctionnel[ind_plex])
            if compt_fonct_mod==5*nb_sol:
                break
            compt_fonct_mod+=1
        list_fonct_entraxe=npy.array(fonct_entraxe)
        sort_list_fonct_entraxe=npy.argsort(list_fonct_entraxe)
        plex_analyse=[]
        for ind_plex in range(nb_sol):
            plex_analyse.append(sort_fonct_module[::-1][sort_list_fonct_entraxe[ind_plex]])
        
        compt_nb_sol=0
        liste_solutions=[]
#        for num_plex in npy.argsort(list_fonctionnel):
        for num_plex in plex_analyse:
            plex=self.plex_calcul[num_plex]
            ga=ContinuousMeshesAssemblyOptimizer(**plex)
            try:
                ga.Optimize(verbose = verbose)
            # BUG: Définir l'erreur
            except ValueError:
                print('Convergence Problem')
            if len(ga.solutions)>0:
                xsol=ga.solutions[-1]
                xt=dict(list(ga.xi.items())+list(xsol.items()))
                solutions=MeshAssembly(**xt)
                valid_cd=True
                for engr_num,cd in enumerate(self.center_distance):
                    if (solutions.center_distance[engr_num])<(cd[0]):
                        valid_cd=False
                    elif (solutions.center_distance[engr_num])>(cd[1]):
                        valid_cd=False
                if valid_cd:
                    liste_solutions.append(solutions)
                    compt_nb_sol+=1
                    if verbose:
                        print('valid solution n°{}'.format(compt_nb_sol))
                        print('Module: {}'.format(solutions.meshes[0][list(solutions.meshes[0].keys())[0]].rack.module))
                    if compt_nb_sol==nb_sol:
                        break
                else:
                    if verbose: 
                        print('unvalid solution')
            else:
                if verbose:
                    print('Solution non convergée')
        # TODO: commentaire en anglais
        #Tri des solutions convergées en fonction de la fonctionnelle actualisée
        list_fonctionnel_nv=[]
        for solutions in liste_solutions:
            fonctionnel=0
            for engr_num,cd in enumerate(self.center_distance):
                fonctionnel+=(solutions.center_distance[engr_num]-cd[0])**2
            list_fonctionnel_nv.append(fonctionnel)
        for indice_plex,num_plex in enumerate(npy.argsort(npy.array(list_fonctionnel_nv))):
            self.solutions.append(liste_solutions[num_plex])
            if verbose and indice_plex==0:
                print('Mesh sections: {}'.format(self.solutions[-1].gear_width))
                print('Center distances: {}'.format(self.solutions[-1].center_distance))
                print('Module: {}'.format(self.solutions[-1].meshes[0][list(self.solutions[-1].meshes[0].keys())[0]].rack.module))

