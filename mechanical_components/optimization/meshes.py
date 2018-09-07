#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:13:01 2018

@author: Pierre-Emmanuel Dumouchel
"""

#import itertools as it
from mechanical_components.meshes import MeshAssembly, hardened_alloy_steel
import numpy as npy
import itertools
import networkx as nx
import powertransmission.tools as tools
from scipy.optimize import minimize
import copy

class ContinuousMeshesAssemblyOptimizer:
    """
    Gear mesh assembly optimizer
    
    :param Z: Dictionary define the tooth number of each mesh {node1: Z1, node2: Z2 ...}
    :param center_distance: List of two elements define gear mesh center-distance [[mesh1_centerdistance_min, mesh1_centerdistance_max], [mesh2_centerdistance_min, mesh2_centerdistance_max] ...] with mesh1 the mesh between node1 and node2 ...    
    :param connections: List of tuple define gear mesh connection [(node1,node2), (node2,node3), (node2,node4)]
    :param transverse_pressure_angle: List of two elements define the transversal pressure angle interval for each mesh [[mesh1_transversepressure_min, mesh1_transversepressure_max], [mesh2_transversepressure_min, mesh2_transversepressure_max] ...]
    :param coefficient_profile_shift: Dictionary defining the minimum and maximum coefficient profile shift for each node {node1: [node1_coeffshift_min,node1_coeffshift_max],node2: [node2_coeffshift_min,node2_coeffshift_max] ...}
    :param rack_list: Dictionary define all admissible rack {rack1:mechanical_components.meshes.Rack, rack2:mechanical_components.meshes.Rack ...}
    :param rack_choice: Dictionary assign for each gear mesh a list of acceptable rack {node1:rack1, node2:rack1, node3:rack4 ...}
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
                 coefficient_profile_shift,
                 rack_list, rack_choice, material,torque,
                 cycle, safety_factor,db,verbose=False):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.Z=Z
        self.connections=connections
        
        # Initailization
        self.solutions=[]
        
        # NetworkX graph construction
        list_gear=[]
        for gs in connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(connections)
        self.gear_graph=gear_graph
        self.list_gear=list_gear
        self.connections_dfs=list(nx.dfs_edges(gear_graph,connections[0][0]))
        
        # Search of the min/max of the base diameter of the first gear mesh
        mesh1,mesh2=connections[0]
        Z1=Z[mesh1]
        Z2=Z[mesh2]
        cd1_min,cd1_max=self.center_distance[0]
        tpa1_min,tpa1_max=self.transverse_pressure_angle[0]
        df1_min=2*cd1_min*Z1/(Z1+Z2)
        db1_min=df1_min*npy.cos(tpa1_max)
        df1_max=2*cd1_max*Z1/(Z1+Z2)
        db1_max=df1_max*npy.cos(tpa1_min)
        self.db=[db1_min,db1_max]
        
        # Search of unknown parameters (borne_min different of borne_max)
        dict_unknown={'db':[],'transverse_pressure_angle':[],
                 'coefficient_profile_shift':[],'transverse_pressure_angle_rack':[],
                 'coeff_gear_addendum':[],'coeff_gear_dedendum':[],
                 'coeff_root_radius':[],'coeff_circular_tooth_thickness':[]}
        dict_global=copy.deepcopy(dict_unknown)
        for i in list(set(list(self.rack_choice.values()))):
            for k,v in self.rack_list[i].items():
                if k not in ['type','name','module']:
                    dict_global[k].append(i)
                    if not v[0]==v[1]:
                        dict_unknown[k].append(i)
        dict_global['db']=[connections[0][0]]
        if not db1_min==db1_max:
            dict_unknown['db']=[connections[0][0]]
        for num_mesh,list_tpa in enumerate(self.transverse_pressure_angle):
            dict_global['transverse_pressure_angle'].append(num_mesh)
            if not list_tpa[0]==list_tpa[1]:
                dict_unknown['transverse_pressure_angle'].append(num_mesh)
        for num_gear in self.list_gear:
            cps=coefficient_profile_shift[num_gear]
            dict_global['coefficient_profile_shift'].append(num_gear)
            if not cps[0]==cps[1]:
                dict_unknown['coefficient_profile_shift'].append(num_gear)
        self.dict_unknown=dict_unknown
        self.dict_global=dict_global
        number_unknown=0
        for key,list_unknown in dict_unknown.items():
            number_unknown+=len(list_unknown)
#        if verbose:
#            print('The total number of unknown for the gear mesh assembly optimization is {}'.format(number_unknown))
        
        # Definition of the Bound matrix for the optimizer
        Bounds=[]
        list_order_unknown=[]
        for key,list_unknown in dict_unknown.items():
            if len(list_unknown)>0:
                if key=='db':
                    Bounds.append(self.db)
                    list_order_unknown.append(key)
                elif key=='transverse_pressure_angle':
                    list_order_unknown.append(key)
                    for num_mesh in list_unknown:
                        Bounds.append(self.transverse_pressure_angle[num_mesh])
                elif key=='coefficient_profile_shift':
                    list_order_unknown.append(key)
                    for num_gear in list_unknown:
                        Bounds.append(self.coefficient_profile_shift[num_gear])
                elif key in ['transverse_pressure_angle_rack','coeff_gear_addendum',
                             'coeff_gear_dedendum','coeff_root_radius',
                             'coeff_circular_tooth_thickness']:
                    list_order_unknown.append(key)
                    for num_rack in list_unknown:
                        Bounds.append(self.rack_list[num_rack][key])
        self.Bounds=npy.array(Bounds)
        self.list_order_unknown=list_order_unknown
        
        # Definition initial condition
        self.X0=self.CondInit()
        self.general_data={'Z': Z, 'connections': connections,
                 'material':material,'torque':torque,'cycle':cycle,
                 'safety_factor':safety_factor,'verbose':verbose}
        self.optimizer_data=self._convert_X2x(self.X0)
        self.xt=dict(list(self.optimizer_data.items())+list(self.general_data.items()))
        self.MeshAssembly = MeshAssembly(**self.xt)
        
        self.save=copy.deepcopy(self.optimizer_data)
        
    def CondInit(self):
        X0=[]
        for interval in self.Bounds:
            X0.append((interval[1]-interval[0])*float(npy.random.random(1))+interval[0])
        return X0
    
    def _convert_X2x(self,X):
        optimizer_data={'center_distance':[],'transverse_pressure_angle':[],
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}} # dictionary of possibily/currently optimization data
        # copy and transfert data to the dictionary optimizer_data
        liste_transverse_pressure_angle=[]
        for key,list_ind in self.dict_global.items():
            for num_elem in list_ind:
                if key in ['transverse_pressure_angle']:
                    if num_elem in self.dict_unknown[key]:
                        value=X[self._position_X(key,num_elem)]
                    else:
                        value=self.transverse_pressure_angle[num_elem][0]
                    liste_transverse_pressure_angle.append(value)
                    if num_elem==0:
                        optimizer_data[key].append(value)
                elif key in ['coefficient_profile_shift']:
                    if num_elem in self.dict_unknown[key]:
                        value=X[self._position_X(key,num_elem)]
                    else:
                        value=self.coefficient_profile_shift[num_elem][0]
                    optimizer_data[key][num_elem]=value
                elif key in ['transverse_pressure_angle_rack',
                             'coeff_gear_addendum','coeff_gear_dedendum',
                             'coeff_root_radius','coeff_circular_tooth_thickness']:
                    if num_elem in self.dict_unknown[key]:
                        value=X[self._position_X(key,num_elem)]
                    else:
                        value=self.rack_list[num_elem][key][0]
                    for num_engr in self.list_gear:
                        if self.rack_choice[num_engr]==num_elem:
                            optimizer_data[key][num_engr]=value
                elif key in ['db']:
                    if num_elem in self.dict_unknown[key]:
                        value=X[self._position_X(key,num_elem)]
                    else:
                        value=self.db[0]
                    db_init=value
        # calculation of the center_distance data
        engr_init=self.dict_global['db'][0]
        db={engr_init:db_init}
        Z=self.Z
        for engr1,engr2 in self.connections_dfs:
            db1=db[engr1]
            db2=Z[engr2]/Z[engr1]*db1
            db[engr2]=db2
        for num_mesh,(engr1,engr2) in enumerate(self.connections):
            tpa=liste_transverse_pressure_angle[num_mesh]
            df1=db[engr1]/npy.cos(tpa)
            df2=db[engr2]/npy.cos(tpa)
            optimizer_data['center_distance'].append((df1+df2)/2)
        return optimizer_data
                    
    def _position_X(self,key,indice):
        # search of the position in the vector X of the data define by (key,indice) in dict_global
        position=0
        for evol_key in self.list_order_unknown:
            val=self.dict_unknown[evol_key]
            if evol_key!=key:
                position+=len(val)
            else:
                break
        for ind in val:
            if ind!=indice:
                position+=1
            else:
                break
        return position
    
    def Update(self,X):
        optimizer_data=self._convert_X2x(X)
        xt=dict(list(optimizer_data.items())+list(self.general_data.items()))
        if self.save!=optimizer_data:
            self.MeshAssembly.Update(**xt)
            self.save=copy.deepcopy(optimizer_data)
    
    def Fineq(self,X):
        
        self.Update(X)
        ineq=[]
        ineq.extend(self.MeshAssembly.ListeIneq())

        #geometric constraint
        for num_mesh,(engr1,engr2) in enumerate(self.MeshAssembly.connections):
            dia1=self.MeshAssembly.meshes[engr1].root_diameter_active
            dia2=self.MeshAssembly.meshes[engr2].root_diameter_active
            de1=self.MeshAssembly.meshes[engr1].outside_diameter
            de2=self.MeshAssembly.meshes[engr2].outside_diameter
            cd=self.MeshAssembly.center_distance[num_mesh]
            ineq.append(cd-(de1/2+dia2/2))
            ineq.append(cd-(de2/2+dia1/2))
            oaa1=self.MeshAssembly.meshes[engr1].outside_active_angle
            oaa2=self.MeshAssembly.meshes[engr2].outside_active_angle
            ineq.append(oaa1)
            ineq.append(oaa2)
            df1=self.MeshAssembly.DF[num_mesh][engr1]
            df2=self.MeshAssembly.DF[num_mesh][engr2]
            db1=self.MeshAssembly.meshes[engr1].db
            db2=self.MeshAssembly.meshes[engr2].db
            ineq.append(df1-db1)
            ineq.append(df2-db2)
        # modulus constraint (not in the open-source part because this parameter is not a ddl)
        for ne,gs in enumerate(self.MeshAssembly.connections):
            for g in gs:
                mo=self.MeshAssembly.meshes[g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
                ineq.append(mo-list_module[0])
                ineq.append(list_module[1]-mo)
        # center-distance constraint (not in the open-source part because this parameter is not a ddl)
        for num_mesh,(engr1,engr2) in enumerate(self.MeshAssembly.connections):
            cd=self.MeshAssembly.center_distance[num_mesh]
            limit_cd=self.center_distance[num_mesh]
            ineq.append(cd-limit_cd[0])
            ineq.append(limit_cd[1]-cd)
        return ineq
        
    def Objective(self,X):
        self.Update(X)
#        fineq=self.Fineq(X)
        obj=0
        obj+=self.MeshAssembly.Functional()
        # maximization of the gear modulus
#        for ne,gs in enumerate(self.MeshAssembly.connections):
#            for g in gs:
#                mo=self.MeshAssembly.meshes[g].rack.module
#                list_module=self.rack_list[self.rack_choice[g]]['module']
#                if list_module[0]<list_module[1]:
#                    obj+=100*(list_module[1]-mo)**2
            
        # Minimisation of center-distance on the minimum bound specified
#        for num_engr,list_cd in enumerate(self.center_distance):
#            obj+=100*(list_cd[0]-self.MeshAssembly.center_distance[num_engr])**2
            
#        for i in fineq:
#            if i < 0:
#                obj+=1000*i**2
#            else:
#                obj+=0.000001*i
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
            X0=self.CondInit()
            self.Update(X0)
            cons = {'type': 'ineq','fun' : self.Fineq}
            cx = minimize(self.Objective, X0, bounds=self.Bounds,constraints=cons)
            Xsol=cx.x
            self.Update(Xsol)
            if verbose:
                print('Iteration n°{} with status {}, min(fineq):{}'.format(i,
                      cx.status,min(self.Fineq(Xsol))))
            if min(self.Fineq(Xsol))>-1e-5:
                optimizer_data=self._convert_X2x(Xsol)
                xt=dict(list(optimizer_data.items())+list(self.general_data.items()))
                self.solutions.append(MeshAssembly(**xt))
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
            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,2.5*1e-3],
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
            torque=[{list_gear[0]:100,list_gear[1]:'output'}]
            
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

        self.plex_calcul=self.AnalyzeCombination(verbose)
        
        for i,plex in enumerate(self.plex_calcul):
            plex['rack_list']=self.rack_list
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
        
        demul_int_min=1/9.
        demul_int_max=9
        dt=tools.RegularDecisionTree(np)
        
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
            Z_node = [list_gear[i][node_value] for i,node_value in enumerate(dt.current_node[:self.nb_gear])]
            if (dt.current_depth<=(self.nb_gear-1)) and (dt.current_depth>0):
                z1=list_gear[dt.current_depth-1][dt.current_node[dt.current_depth-1]]
                z2=list_gear[dt.current_depth][dt.current_node[dt.current_depth]]

                if (pgcd(z1,z2)!=1):
                    valid=False
                # gear ratio analysis
                if valid:
                    demul=z1/z2
                    if (demul > demul_int_max) or (demul < demul_int_min):
                        valid=False
                # NVH analysis of gear mesh each other
#                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
#                    for n in liste_node[0:dt.current_depth]:
#                        i=liste_node.index(n)
#                        z=liste_gear[dt.current_node[i]]
#                        if pgcd(z,z2)!=1:
#                            valid=False
                # speed specification analysis
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
                # rack feasibility analysis
                rack_num=list_rack[dt.current_node[-1]]
                rack_pos=list_node[dt.current_depth-self.nb_gear]
                if rack_num not in self.rack_choice[rack_pos]:
                    valid=False
                    
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)):
                 
                # feasibility analysis of the modulus toward center-distance
                liste_DF_min={}
                module_minmax={}
                module_inf,module_sup=(0,npy.inf)
                # intersection between all module specified by each rack
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
                
                # Analysis of the impact of the new center-distance on the initial base diamter interval
                # Initialisation of the dfs graph with the gear mesh used for the unknown db
                connections_dfs=list(nx.dfs_edges(self.gear_graph,self.connections[0][0]))
                db={}
                for ind,(eng1,eng2) in enumerate(connections_dfs[::-1]):
                    if (eng1,eng2) in self.connections:
                        num_mesh=self.connections.index((eng1,eng2))
                    elif (eng2,eng1) in self.connections:
                        num_mesh=self.connections.index((eng2,eng1))
                    if len(cd_minmax_nv)>num_mesh:
                        cd_min,cd_max=cd_minmax_nv[num_mesh]
                        z1_pos=list_node.index(eng1)
                        z1=list_gear[z1_pos][dt.current_node[z1_pos]]
                        z2_pos=list_node.index(eng2)
                        z2=list_gear[z2_pos][dt.current_node[z2_pos]]
                        tpa1_min,tpa1_max=self.transverse_pressure_angle[num_mesh]
                        df2_min=2*cd_min-2*cd_min*z1/(z1+z2)
                        db2_min=df2_min*npy.cos(tpa1_max)
                        df2_max=2*cd_max-2*cd_max*z1/(z1+z2)
                        db2_max=df2_max*npy.cos(tpa1_min)
                        try:
                            db[eng2]=[max(db2_min,db[eng2][0]),min(db2_max,db[eng2][1])]
                        except KeyError:
                            db[eng2]=[db2_min,db2_max]
                        db[eng1]=[db[eng2][0]*z1/z2,db[eng2][1]*z1/z2]
                    else:
                        break
        
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
                Export['center_distance'] = cd_minmax_nv
                Export['db'] = db
                Export['dw'] = v0_max - v0_min

                self.fonctionnel.append(fonctionnel)
                self.fonctionnel_module.append(module_optimal)
                plex_calcul.append(Export)
                incr+=1
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
#                print('dw: ', liste_plex[i]['dw'])
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
                sol1=ga.solutions[-1]
                self.solutions.append(sol1)
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
        for num_plex in plex_analyse:
            plex=self.plex_calcul[num_plex]
            del plex['dw']
            ga=ContinuousMeshesAssemblyOptimizer(**plex)
            try:
                ga.Optimize(verbose)
            except ValueError:
                print('Convergence Problem')
            if len(ga.solutions)>0:
                solutions=ga.solutions[-1]
                valid_cd=True
                for engr_num,cd in enumerate(self.center_distance):
                    if (solutions.center_distance[engr_num])<(cd[0]*0.999):
                        valid_cd=False
                    elif (solutions.center_distance[engr_num])>(cd[1]*1.001):
                        valid_cd=False
                if valid_cd:
                    liste_solutions.append(solutions)
                    compt_nb_sol+=1
                    if verbose:
                        print('valid solution n°{}'.format(compt_nb_sol))
                        print('Module: {}'.format(solutions.meshes[list(solutions.meshes.keys())[0]].rack.module))
                    if compt_nb_sol==nb_sol:
                        break
                else:
                    if verbose: 
                        print('unvalid solution')
            else:
                if verbose:
                    print('Convergence Problem')
        # admissible solution are sorted toward the updated functional
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
                print('Module: {}'.format(self.solutions[-1].meshes[list(self.solutions[-1].meshes.keys())[0]].rack.module))

