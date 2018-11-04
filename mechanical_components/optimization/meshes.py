#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:13:01 2018

@author: Pierre-Emmanuel Dumouchel
"""

#import itertools as it
from mechanical_components.meshes import MeshAssembly, hardened_alloy_steel
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import itertools
import networkx as nx
#import powertransmission.tools as tools
import dectree

from scipy.optimize import minimize
import copy

class CompositeMeshAssembly:
    def __init__(self,center_distance, connections, transverse_pressure_angle,
                 coefficient_profile_shift,transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness, Z, strong_link=None, material=None,
                 torque=None, cycle=None,
                 safety_factor=1,verbose=False):
        
        self.connections=connections
        self.center_distance=center_distance
        self.mesh_assembly = []
        self.general_data = []
        
        # 
        for num_graph,list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh=0
            general_data={'Z': {}, 'connections': [],
                 'material':{},'torque':{},'cycle':{},
                 'safety_factor':safety_factor,'verbose':verbose}
            input_data={'center_distance':[],'transverse_pressure_angle':[],
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}}
            li_connection=[]
            for num_cd,list_connection in enumerate(connections):
                for num_mesh_iter,gs in enumerate(list_connection):
                    if (gs in list_sub_graph) or (gs[::-1] in list_sub_graph):
                        li_connection.append(gs)
                        for num_gear in gs:
                            if num_gear in coefficient_profile_shift.keys():
                                input_data['coefficient_profile_shift'][num_gear]=coefficient_profile_shift[num_gear]
                            if num_gear in transverse_pressure_angle_rack.keys():
                                input_data['transverse_pressure_angle_rack'][num_gear]=transverse_pressure_angle_rack[num_gear]
                            if num_gear in coeff_gear_addendum.keys():
                                input_data['coeff_gear_addendum'][num_gear]=coeff_gear_addendum[num_gear]
                            if num_gear in coeff_gear_dedendum.keys():
                                input_data['coeff_gear_dedendum'][num_gear]=coeff_gear_dedendum[num_gear]
                            if num_gear in coeff_root_radius.keys():
                                input_data['coeff_root_radius'][num_gear]=coeff_root_radius[num_gear]
                            if num_gear in coeff_circular_tooth_thickness.keys():
                                input_data['coeff_circular_tooth_thickness'][num_gear]=coeff_circular_tooth_thickness[num_gear]                                
                            if num_gear in Z.keys():
                                general_data['Z'][num_gear]=Z[num_gear]
                            if num_gear in material.keys():
                                general_data['material'][num_gear]=material[num_gear]
                        input_data['transverse_pressure_angle'].append(transverse_pressure_angle[num_mesh])    
                    num_mesh+=1
                input_data['center_distance'].append(center_distance[num_cd])
            general_data['connections']=li_connection
            for (eng1,eng2) in list_sub_graph:
                if (eng1,eng2) in torque.keys():
                    general_data['torque'][(eng1,eng2)]=torque[(eng1,eng2)]
                if (eng2,eng1) in torque.keys():
                    general_data['torque'][(eng2,eng1)]=torque[(eng2,eng1)]
                if eng1 not in general_data['cycle'].keys():
                    general_data['cycle'][eng1]=cycle[eng1]
                if eng2 not in general_data['cycle'].keys():
                    general_data['cycle'][eng2]=cycle[eng2]
                    
            self.general_data.append(general_data)
            xt=dict(list(input_data.items())+list(general_data.items()))
            self.mesh_assembly.append(MeshAssembly(**xt))
            
    def _get_graph_dfs(self):
        _graph_dfs,_=gear_graph_simple(self.connections)
        return _graph_dfs
    sub_graph_dfs=property(_get_graph_dfs)
    
    def _get_list_gear(self):
        _,_list_gear=gear_graph_simple(self.connections)
        return _list_gear
    list_gear=property(_get_list_gear)
    
    def SVGExport(self,name,position):
        
        centers=self.PosAxis(position)
        L=[]
        for mesh_assembly_iter in self.mesh_assembly:
            position_svg={}
            for num_gear,pos in centers.items():
                if num_gear in mesh_assembly_iter.Z.keys():
                    position_svg[num_gear]=pos
            L.extend(mesh_assembly_iter.SVGExport('gear',position_svg))
        G1=vm.Contour2D(L)
        G1.MPLPlot()        
        
    def FreeCADExport(self, fcstd_filepath, centers = {}, axis = (1,0,0), export_types=['fcstd'], python_path = 'python',
                      path_lib_freecad = '/usr/lib/freecad/lib'):
        """ Export 3D volume to FreeCAD
        
        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation
        
        :results: export of a FreeCAD file
        """
        for ma in self.mesh_assembly:
            ma.FreeCADExport(fcstd_filepath, centers, axis, python_path, path_lib_freecad, export_types)
        
    def Update(self,optimizer_data):
        output_x=[]
        for num_graph,list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh=0
            input_data={'center_distance':[],'transverse_pressure_angle':[],
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}}
            li_connection=[]
            for num_cd,list_connection in enumerate(self.connections):
                for num_mesh_iter,(eng1,eng2) in enumerate(list_connection):
                    if ((eng1,eng2) in list_sub_graph) or ((eng2,eng1) in list_sub_graph):
                        li_connection.append((eng1,eng2))
                        for key,list_value in optimizer_data.items():
                            if key in ['coefficient_profile_shift',
                                       'transverse_pressure_angle_rack',
                                       'coeff_gear_addendum','coeff_gear_dedendum',
                                       'coeff_root_radius','coeff_circular_tooth_thickness']:
                                input_data[key][eng1]=optimizer_data[key][eng1]
                                input_data[key][eng2]=optimizer_data[key][eng2]
                            elif key in ['center_distance']:
                                input_data[key].append(optimizer_data[key][num_cd])
                            elif key in ['transverse_pressure_angle']:
                                input_data[key].append(optimizer_data[key][num_mesh])
                    num_mesh+=1
            xt=dict(list(input_data.items())+list(self.general_data[num_graph].items()))
            output_x.append(xt)
            
#            if self.save!=optimizer_data:
            self.mesh_assembly[num_graph].Update(**xt)
        return output_x
        
    def PosAxis(self,position):
        # Definition of the initial center for all gear (when not given by the user)
        
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(self.list_gear)
        for num_cd,list_connections in enumerate(self.connections):
            (eng1_m,eng2_m)=list_connections[0]
            if len(list_connections)>1:
                for (eng1,eng2) in list_connections[1:]:
                    gear_graph.add_edges_from([(eng1_m,eng1),(eng2_m,eng2)])        
                    eng1_m=eng1
                    eng2_m=eng2
        list_line=list(nx.connected_component_subgraphs(gear_graph))
        dict_line={}
        for num_line,list_num_eng in enumerate(list_line):
            for num_eng in list_num_eng:
                dict_line[num_eng]=num_line
        def fun(x):
            obj=0
            for num_cd,list_connections in enumerate(self.connections):
                eng1=dict_line[list_connections[0][0]]
                eng2=dict_line[list_connections[0][1]]
                obj+=(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num_cd])**2
            return obj
        def eg(x):
            ine=[]
            for k,val in position.items():
                key=dict_line[k]
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine=[]
            for num_cd,list_connections in enumerate(self.connections):
                eng1=dict_line[list_connections[0][0]]
                eng2=dict_line[list_connections[0][1]]
                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*self.center_distance[num_cd])
                ine.append(1.001*self.center_distance[num_cd]-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
            return ine
        cons = ({'type': 'eq','fun' : eg},{'type': 'ineq','fun' : ineg})
        drap=1
        while drap==1:
            x0=tuple(npy.random.random(2*len(list_line))*1)
            Bound=[[0,1]]*(len(list_line)*2)
            res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
            if (min(ineg(res.x))>0) and (max(eg(res.x))<1e-7):
                drap=0
        x_opt=res.x
        centers={}
        for num_pos,num_eng in enumerate(self.list_gear):
            opt_pos=dict_line[num_eng]
            centers[num_eng]=[x_opt[2*opt_pos],x_opt[2*opt_pos+1]]
        return centers


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
    def __init__(self, Z, center_distance, connections,strong_link, transverse_pressure_angle,
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
        self.strong_link=strong_link
        self.torque=torque
        self.cycle=cycle
        # Initailization
        self.solutions=[]
        
#        # NetworkX graph construction
#        list_gear=[] # list of all gears
#        compt_mesh=0 # number of gear mesh
#        for gs in connections:
#            for (eng1,eng2) in gs:
#                compt_mesh+=1
#                if eng1 not in list_gear:
#                    list_gear.append(eng1)
#                if eng2 not in list_gear:
#                    list_gear.append(eng2)
#        self.list_gear=list_gear
#        # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
#        self.nb_gear=len(list_gear)
#        gear_graph=nx.Graph()
#        gear_graph.add_nodes_from(list_gear)
#        for list_edge in self.connections:
#            gear_graph.add_edges_from(list_edge)                  
#        self.gear_graph=gear_graph
#        sub_graph=list(nx.connected_component_subgraphs(gear_graph))
#        self.sub_graph_dfs=[]
#        for s_graph in sub_graph:
#            node_init=list(s_graph.nodes())[0]
#            self.sub_graph_dfs.append(list(nx.dfs_edges(s_graph,node_init)))
        
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
        # search the unknown list db and transverse_pressure_angle to define in the optimizer
        list_analyze_cd={}
        list_analyze_line_gear=[]
        num_mesh=0
        list_unknown=[]
        for num_cd,list_connections in enumerate(connections):
            for (eng1,eng2) in list_connections:
                for num_line_iter,list_dfs in enumerate(self.sub_graph_dfs): # search the number of the line gear analyze
                    if ((eng1,eng2) in list_dfs) or ((eng2,eng1) in list_dfs):
                        num_line=num_line_iter
                if num_cd not in list_analyze_cd.keys():
                    list_analyze_cd[num_cd]=num_mesh
                    if num_line not in list_analyze_line_gear:
                        list_analyze_line_gear.append(num_line)
                        interval_db=db[eng1]
                        if interval_db[0]!=interval_db[1]:
                            dict_unknown['db'].append(eng1)
                        interval_tpa=self.transverse_pressure_angle[num_mesh]
                        if interval_tpa[0]!=interval_tpa[1]:
                            dict_unknown['transverse_pressure_angle'].append(num_mesh)
                    else:
                        interval_tpa=self.transverse_pressure_angle[num_mesh]
                        if interval_tpa[0]!=interval_tpa[1]:
                            dict_unknown['transverse_pressure_angle'].append(num_mesh)
                else:
                    if num_line not in list_analyze_line_gear:
                        list_analyze_line_gear.append(num_line)
                        interval_tpa=self.transverse_pressure_angle[num_mesh]
                        if interval_tpa[0]!=interval_tpa[1]:
                            dict_unknown['transverse_pressure_angle'].append(num_mesh)
                num_mesh+=1
        for num_gear in self.list_gear:
            dict_global['db'].append(num_gear)
            dict_global['transverse_pressure_angle'].append(num_gear)
        #search another unknown
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
        Bounds = []
        list_order_unknown = []
        for key,list_unknown in dict_unknown.items():
            if len(list_unknown)>0:
                if key=='db':
                    list_order_unknown.append(key)
                    for num_gear in list_unknown:
                        Bounds.append(db[num_gear])
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
        self.Bounds = npy.array(Bounds)
        self.list_order_unknown = list_order_unknown
        
        # Definition initial condition
        self.X0 = self.CondInit()
        
        self.db = db
        optimizer_data = self._convert_X2x(self.X0)
        dic_torque,dic_cycle = self.TorqueCycleMeshAssembly()
        
        self.general_data={'Z': Z, 'connections': connections,
                 'material':material,'torque':dic_torque,'cycle':dic_cycle,
                 'safety_factor':safety_factor,'verbose':verbose}
        input_dat=dict(list(optimizer_data.items())+list(self.general_data.items()))
        self.mesh_assembly=CompositeMeshAssembly(**input_dat)
        
        self.save=copy.deepcopy(optimizer_data)
    
    def _get_graph_dfs(self):
        _graph_dfs,_=gear_graph_simple(self.connections)
        return _graph_dfs
    sub_graph_dfs=property(_get_graph_dfs)
    
    def _get_list_gear(self):
        _,_list_gear=gear_graph_simple(self.connections)
        return _list_gear
    list_gear=property(_get_list_gear)
        
    def TorqueCycleMeshAssembly(self):
        # construction of a graph with gea_mesh and same_speed for analyze with dfs the order of torque calculation
        torque_graph=nx.Graph()
        torque_graph.add_nodes_from(self.list_gear)
        for list_edge in self.connections:
            torque_graph.add_edges_from(list_edge,typ='gear_mesh')
            li_shaft1=[]
            li_shaft2=[]
            for eng1,eng2 in list_edge:
                li_shaft1.append(eng1)
                li_shaft2.append(eng2)
            if len(li_shaft1)>1:
                for pos_gear,num_gear in enumerate(li_shaft1[1:]):
                    valid_strong_ling=False
                    for list_strong_link in self.strong_link:
                        if (num_gear in list_strong_link) and (li_shaft1[pos_gear] in list_strong_link):
                            valid_strong_ling=True
                    if valid_strong_ling:
                        torque_graph.add_edges_from([(num_gear,li_shaft1[pos_gear])],typ='same_speed')
            if len(li_shaft2)>1:
                for pos_gear,num_gear in enumerate(li_shaft2[1:]):
                    valid_strong_ling=False
                    for list_strong_link in self.strong_link:
                        if (num_gear in list_strong_link) and (li_shaft2[pos_gear] in list_strong_link):
                            valid_strong_ling=True
                    if valid_strong_ling:
                        torque_graph.add_edges_from([(num_gear,li_shaft2[pos_gear])],typ='same_speed')
        for num_gear,tq in self.torque.items():
            if tq=='output':
                node_output=num_gear
        torque_graph_dfs=list(nx.dfs_edges(torque_graph,node_output))
        order_torque_calculation=[(eng2,eng1) for (eng1,eng2) in torque_graph_dfs[::-1]]
        # calculation torque distribution
        temp_torque={}
        for eng1 in self.list_gear:
            temp_torque[eng1]=0
        for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
            if eng1 in self.torque.keys():
                temp_torque[eng1]+=self.torque[eng1]
            if torque_graph[eng1][eng2]['typ']=='gear_mesh':
                temp_torque[eng2]+=-temp_torque[eng1]*self.Z[eng2]/self.Z[eng1]
            else:
                temp_torque[eng2]+=temp_torque[eng1]
        dic_torque={}
        for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
            dic_torque[(eng1,eng2)]=temp_torque[eng1]
        # calculation cycle distribution
        dic_cycle=dict(self.cycle)
        node_init=list(self.cycle.keys())[0]
        cycle_graph_dfs=list(nx.dfs_edges(torque_graph,node_init))
        for num_mesh_cy,(eng1,eng2) in enumerate(cycle_graph_dfs):
            if torque_graph[eng1][eng2]['typ']=='gear_mesh':
                dic_cycle[eng2]=dic_cycle[eng1]*self.Z[eng1]/self.Z[eng2]
            else:
                dic_cycle[eng2]=dic_cycle[eng1]
        return dic_torque,dic_cycle
                
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
                if key in ['coefficient_profile_shift']:
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
        
        dict_db={} # temp storage of db dict
        dict_df={} # temp storage of pitch diam
        dict_cd={} # temp storage of center_distance
        dict_tpa={} # temp storage of transversale_pressure_angle
        for num_gear in self.dict_unknown['db']:
            dict_db[num_gear]=X[self._position_X('db',num_gear)]
        for num_mesh in self.dict_unknown['transverse_pressure_angle']:
            dict_tpa[num_mesh]=X[self._position_X('transverse_pressure_angle',num_mesh)]
        num_mesh=0
        for num_cd,list_connection in enumerate(self.connections):
            for num_mesh_iter,(eng1,eng2) in enumerate(list_connection):
                if num_mesh_iter==0: # in this case we must define center_distance
                    if eng1 in dict_db.keys():
                        if num_mesh in dict_tpa.keys():
                            db2=self.Z[eng2]/self.Z[eng1]*dict_db[eng1]
                            dict_db[eng2]=db2
                            dict_df[eng1]=dict_db[eng1]/npy.cos(dict_tpa[num_mesh])
                            dict_df[eng2]=db2/npy.cos(dict_tpa[num_mesh])
                            dict_cd[num_cd]=(dict_df[eng1]+dict_df[eng2])/2
                    elif eng2 in dict_db.keys():
                        if num_mesh in dict_tpa.keys():
                            db1=self.Z[eng1]/self.Z[eng2]*db2
                            dict_db[eng1]=db1
                            dict_df[eng1]=db1/npy.cos(dict_tpa[num_mesh])
                            dict_df[eng2]=dict_db[eng2]/npy.cos(dict_tpa[num_mesh])
                            dict_cd[num_cd]=(dict_df[eng1]+dict_df[eng2])/2
                else: # the center_distance is define
                    if eng1 in dict_db.keys():
                        cd=dict_cd[num_cd]
                        dict_db[eng2]=self.Z[eng2]/self.Z[eng1]*dict_db[eng1]
                        dict_df[eng2]=2*cd*self.Z[eng2]/(self.Z[eng1]+self.Z[eng2])
                        dict_df[eng1]=2*cd-dict_df[eng2]
                        dict_tpa[num_mesh]=npy.arccos(dict_db[eng2]/dict_df[eng2])
                    elif eng2 in dict_db.keys():
                        cd=dict_cd[num_cd]
                        dict_db[eng1]=self.Z[eng1]/self.Z[eng2]*dict_db[eng2]
                        dict_df[eng1]=2*cd*self.Z[eng1]/(self.Z[eng1]+self.Z[eng2])
                        dict_df[eng2]=2*cd-dict_df[eng1]
                        dict_tpa[num_mesh]=npy.arccos(dict_db[eng1]/dict_df[eng1])
                    elif num_mesh in dict_tpa.keys():
                        cd=dict_cd[num_cd]
                        dict_df[eng2]=2*cd*self.Z[eng2]/(self.Z[eng1]+self.Z[eng2])
                        dict_df[eng1]=2*cd-dict_df[eng2]
                        dict_db[eng1]=dict_df[eng1]*npy.cos(dict_tpa[num_mesh])
                        dict_db[eng2]=dict_df[eng2]*npy.cos(dict_tpa[num_mesh])
                    else:
                        cd=dict_cd[num_cd]
                        try:
                            dict_db[eng2]=self.Z[eng2]/self.Z[eng1]*dict_db[eng1]
                        except:
                            dict_db[eng1]=self.Z[eng1]/self.Z[eng2]*dict_db[eng2]
                        dict_df[eng2]=2*cd*self.Z[eng2]/(self.Z[eng1]+self.Z[eng2])
                        dict_df[eng1]=2*cd-dict_df[eng2]
                        dict_tpa[num_mesh]=npy.arccos(dict_db[eng1]/dict_df[eng1])
                num_mesh+=1
        num_mesh=0
        for num_cd,list_connection in enumerate(self.connections):
            optimizer_data['center_distance'].append(dict_cd[num_cd])
            for num_mesh_iter,(eng1,eng2) in enumerate(list_connection):
                optimizer_data['transverse_pressure_angle'].append(dict_tpa[num_mesh])
                num_mesh+=1
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
        
        optimizer_data = self._convert_X2x(X)
        output_x = self.mesh_assembly.Update(optimizer_data)
        return optimizer_data
    
    def Fineq(self,X):
        
        x=self.Update(X)
        ineq=[]
        for mesh_assembly_iter in self.mesh_assembly.mesh_assembly:
            ineq.extend(mesh_assembly_iter.ListeIneq())
            #geometric constraint
            for num_mesh,(engr1,engr2) in enumerate(mesh_assembly_iter.connections):
                dia1=mesh_assembly_iter.meshes[engr1].root_diameter_active
                dia2=mesh_assembly_iter.meshes[engr2].root_diameter_active
                de1=mesh_assembly_iter.meshes[engr1].outside_diameter
                de2=mesh_assembly_iter.meshes[engr2].outside_diameter
                cd=mesh_assembly_iter.center_distance[num_mesh]
                ineq.append(cd-(de1/2+dia2/2))
                ineq.append(cd-(de2/2+dia1/2))
                oaa1=mesh_assembly_iter.meshes[engr1].outside_active_angle
                oaa2=mesh_assembly_iter.meshes[engr2].outside_active_angle
                ineq.append(oaa1)
                ineq.append(oaa2)
                df1=mesh_assembly_iter.DF[num_mesh][engr1]
                df2=mesh_assembly_iter.DF[num_mesh][engr2]
                db1=mesh_assembly_iter.meshes[engr1].db
                db2=mesh_assembly_iter.meshes[engr2].db
                ineq.append(df1-db1)
                ineq.append(df2-db2)
            # modulus constraint (not in the open-source part because this parameter is not a ddl)
            for ne,gs in enumerate(mesh_assembly_iter.connections):
                for g in gs:
                    mo=mesh_assembly_iter.meshes[g].rack.module
                    list_module=self.rack_list[self.rack_choice[g]]['module']
                    ineq.append(mo-list_module[0])
                    ineq.append(list_module[1]-mo)
            # center-distance constraint (not in the open-source part because this parameter is not a ddl)
            for num_mesh,(engr1,engr2) in enumerate(mesh_assembly_iter.connections):
                cd=mesh_assembly_iter.center_distance[num_mesh]
                limit_cd=self.center_distance[num_mesh]
                ineq.append(cd-limit_cd[0])
                ineq.append(limit_cd[1]-cd)
        return ineq
        
    def Objective(self,X):
        x=self.Update(X)
        fineq=self.Fineq(X)
        obj=0
        for mesh_assembly_iter in self.mesh_assembly.mesh_assembly:
            obj+=mesh_assembly_iter.Functional()
        # maximization of the gear modulus
        for ne,mesh_assembly_iter in enumerate(self.mesh_assembly.mesh_assembly):
            for gs in mesh_assembly_iter.connections:
                for g in gs:
                    mo=mesh_assembly_iter.meshes[g].rack.module
                    list_module=self.rack_list[self.rack_choice[g]]['module']
                    if list_module[0]<list_module[1]:
                        obj+=100*(list_module[1]-mo)**2
            
        # Minimisation of center-distance on the minimum bound specified
        for num_engr,list_cd in enumerate(self.center_distance):
            obj+=100*(list_cd[0]-self.mesh_assembly.center_distance[num_engr])**2
            
        for i in fineq:
            if i < 0:
                obj+=1000*i**2
            else:
                obj+=0.000001*i
        return obj
    
    def Optimize(self, verbose = False):
        """ Optimizer function
        
        >>> cmao1.Optimize(verbose = True)
        Iteration n°0 with status 0, min(fineq):-4.7435658082073395e-08, max(eq):0
        Mesh sections: {1: 0.0039733070033111731, 0: 0.0039733070033111731}
        Numbers of teeth: {1: 59, 0: 19}
        Center distances: [0.11700000000000001]
        """
        max_iter=1
        i=0
        arret=0 
        while i<max_iter and arret==0:
            X0=self.CondInit()
            x_temp=self.Update(X0)
            cons = {'type': 'ineq','fun' : self.Fineq}
            cx = minimize(self.Objective, X0, bounds=self.Bounds,constraints=cons)
            Xsol=cx.x
            output_x=self.Update(Xsol)
            if verbose:
                print('Iteration n°{} with status {}, min(fineq):{}'.format(i,
                      cx.status,min(self.Fineq(Xsol))))
            if min(self.Fineq(Xsol))>-1e-5:
                input_dat=dict(list(output_x.items())+list(self.general_data.items()))
                self.solutions.append(CompositeMeshAssembly(**input_dat))
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
    def __init__(self, connections, gear_speed, center_distance, strong_link=[], Z=None,
                 transverse_pressure_angle={}, helix_angle=None,
                 gear_width=None, frequency=[[0,0]],
                 coefficient_profile_shift=None, rack_list=None,
                 rack_choice=None, material=None, torque=None, cycle=None,
                 safety_factor=1,verbose=False):

        self._drap_gear_graph=False
        self._drap_list_gear=False
        self._drap_connections_kinematic_dfs=False
        self._drap_sub_graph_dfs=False
        
        list_gear=[] # list of all gears
        for gs in connections:
            for (eng1,eng2) in gs:
                if eng1 not in list_gear:
                    list_gear.append(eng1)
                if eng2 not in list_gear:
                    list_gear.append(eng2)
        number_mesh=0
        for gs in connections:
            for (eng1,eng2) in gs:
                number_mesh+=1
                  
        # default parameters
        if len(transverse_pressure_angle.keys())<number_mesh:
            for num_mesh in list_gear:
                if num_mesh not in transverse_pressure_angle.keys():
                    transverse_pressure_angle[num_mesh]=[15/180*npy.pi,30/180*npy.pi]

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
                
        speed_min,speed_max=[npy.inf,-npy.inf]
        # definition min/max absolute speed
        for num_engr,(speed_interval_min,speed_interval_max) in gear_speed.items():
            if speed_interval_min<speed_min:
                speed_min=speed_interval_min
            if speed_interval_max>speed_max:
                speed_max=speed_interval_max
            
        if rack_list==None:
#            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,2.5*1e-3],
#                          'transverse_pressure_angle_rack':[20*npy.pi/180.,20*npy.pi/180.],
#                          'coeff_gear_addendum':[1,1],
#                          'coeff_gear_dedendum':[1.25,1.25],
#                          'coeff_root_radius':[0.38,0.38],
#                          'coeff_circular_tooth_thickness':[0.5,0.5]}}
            rack_list={0:{}}
        for num_rack, rack in rack_list.items():
            if 'module' not in rack.keys():
                rack_list[num_rack]['module'] = [1*1e-3,2.5*1e-3]
            if 'transverse_pressure_angle_rack' not in rack.keys():
                rack_list[num_rack]['transverse_pressure_angle_rack'] = [20*npy.pi/180.,20*npy.pi/180.]
            if 'coeff_gear_addendum' not in rack.keys():
                rack_list[num_rack]['coeff_gear_addendum'] = [1,1]
            if 'coeff_gear_dedendum' not in rack.keys():
                rack_list[num_rack]['coeff_gear_dedendum'] = [1.25,1.25]
            if 'coeff_root_radius' not in rack.keys():
                rack_list[num_rack]['coeff_root_radius'] = [0.38,0.38]
            if 'coeff_circular_tooth_thickness' not in rack.keys():
                rack_list[num_rack]['coeff_circular_tooth_thickness'] = [0.5,0.5]
            if 'name' not in rack.keys():
                rack_list[num_rack]['name'] = 'Optim_Module'
                        
            
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
        self.material=material
        self.torque=torque
        self.cycle=cycle
        self.safety_factor=safety_factor
        self.strong_link=strong_link
            
        self.nb_rack=len(self.rack_list.keys())
        
        if self.Z=={}:
            var_Z=self.AnalyseZ()
            self.Z=var_Z
            
        print(self.Z)
        self.plex_calcul = self.AnalyzeCombination(verbose)
        
        for i,plex in enumerate(self.plex_calcul):
            plex['rack_list']=self.rack_list
            plex['material']=self.material
            plex['torque']=self.torque
            plex['cycle']=self.cycle
            plex['connections']=self.connections
            plex['strong_link']=self.strong_link
#            plex['center_distance']=self.center_distance
            plex['transverse_pressure_angle']=self.transverse_pressure_angle
            plex['coefficient_profile_shift']=self.coefficient_profile_shift
            plex['safety_factor']=safety_factor
            self.plex_calcul[i]=plex
            
        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]
        
    def _get_graph_dfs(self):
        if not self._drap_sub_graph_dfs:
            _graph_dfs,_=gear_graph_simple(self.connections)
            self._drap_sub_graph_dfs=True
            self.__cache_drap_sub_graph_dfs=_graph_dfs
            return _graph_dfs
        else:
            return self.__cache_drap_sub_graph_dfs
    sub_graph_dfs=property(_get_graph_dfs)
    
    def _get_list_gear(self):
        if not self._drap_list_gear:
            _,_list_gear=gear_graph_simple(self.connections)
            self._drap_list_gear=True
            self.__cache_drap_list_gear=_list_gear
            return _list_gear
        else:
            return self.__cache_drap_list_gear
    list_gear=property(_get_list_gear)
    
    def _get_gear_graph(self):
        if not self._drap_gear_graph:
            _,_,_gear_graph=gear_graph_complex(self.connections,self.strong_link)
            self._drap_gear_graph=True
            self.__cache_gear_graph=_gear_graph
            return _gear_graph
        else:
            return self.__cache_gear_graph
    gear_graph=property(_get_gear_graph)
    
    def _get_connections_dfs(self):
        _connections_dfs,_,_=gear_graph_complex(self.connections,self.strong_link)
        return _connections_dfs
    connections_dfs=property(_get_connections_dfs)
    
    def _get_connections_kinematic_dfs(self):
        if not self._drap_connections_kinematic_dfs:
            _,_connections_kinematic_dfs,_=gear_graph_complex(self.connections,self.strong_link)
            self._drap_connections_kinematic_dfs=True
            self.__cache_connections_kinematic_dfs=_connections_kinematic_dfs
            return _connections_kinematic_dfs
        else:
            return self.__cache_connections_kinematic_dfs
    connections_kinematic_dfs=property(_get_connections_kinematic_dfs)
        
    def AnalyseZ(self):
        """ Analyse of the minimum and maximum admissible teeth number for each gear mesh
        
        :results: dictionary define the minimum and maximum teeth number for each gear mesh {node1 : [Z1_min, Z1_max], node2 : [Z2_min, Z2_max] ...}
        
        >>> list_z = GA.AnalyseZ()
        >>> print(list_z)
        {1:[13,45], 2: [37,56]}
        """
        Z_max = 130
        
        Z=self.Z
        for i,shaft_mesh in enumerate(self.connections):
            cd_min=self.center_distance[i][0]
            cd_max=self.center_distance[i][1]
            for engr1,engr2 in shaft_mesh:
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
                Z1_max = min(Z1_max, Z_max)
                DF2_max=2*cd_max*demul_max/(1+demul_max)
                Z2_max=int(DF2_max/module2_min)+1
                Z2_max = min(Z2_max, Z_max)
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
        nb_gear=len(self.list_gear)
        list_node=[]
        for eng1,eng2 in self.connections_dfs:
            if eng1 not in list_node:
                list_node.append(eng1)
            if eng2 not in list_node:
                list_node.append(eng2)
                    
        np=[]
        list_Z=[]
        for engr_num in list_node:
            np.append(self.Z[engr_num][1]-self.Z[engr_num][0]+1)
            list_Z.append(npy.arange(self.Z[engr_num][0],self.Z[engr_num][1]+1))
        np.extend([self.nb_rack]*nb_gear)

        list_rack=list(self.rack_list.keys())
        
        demul_int_min=1/9.
        demul_int_max=9
        dt=dectree.RegularDecisionTree(np)

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
            if (dt.current_depth<=(nb_gear-1)) and (dt.current_depth>0):
                current_node=list_node[dt.current_depth]
                list_neighbors=self.gear_graph.neighbors(current_node)
                list_analyze_pos=[dt.current_depth]
                list_analyze_Z=[list_Z[dt.current_depth][dt.current_node[dt.current_depth]]]
                for neighbors in list_neighbors:
                    pos_neighbors=list_node.index(neighbors)
                    if pos_neighbors<dt.current_depth:
                        list_analyze_pos.append(pos_neighbors)
                        z=list_Z[pos_neighbors][dt.current_node[pos_neighbors]]
                        list_analyze_Z.append(z)
                
                # NVH analysis of gear mesh
                for z in list_analyze_Z[1:]:
                    if (pgcd(list_analyze_Z[0],z)!=1):
                        valid=False

                # gear ratio analysis
                if valid:
                    for i,(pos_z,z) in enumerate(zip(list_analyze_pos,list_analyze_Z)):
                        if i>0:
                            if self.gear_graph[list_node[list_analyze_pos[0]]][list_node[pos_z]]['typ']=='gear_mesh':
                                demul=list_analyze_Z[0]/z
                                if (demul > demul_int_max) or (demul < demul_int_min):
                                    valid=False
                # NVH analysis of gear mesh each other
#                if (valid) & (dt.current_depth<=(nb_gear-1)):
#                    for n in liste_node[0:dt.current_depth]:
#                        i=liste_node.index(n)
#                        z=liste_gear[dt.current_node[i]]
#                        if pgcd(z,z2)!=1:
#                            valid=False
                # speed specification analysis
                if valid:
                    v0_min, v0_max = self.gear_speed[list_node[0]]
#                    z0 = list_Z[0][dt.current_node[0]]
                    demul={self.list_gear[0]:1}
                    for (eng1,eng2) in self.connections_kinematic_dfs:
                        pos_eng1=list_node.index(eng1)
                        pos_eng2=list_node.index(eng2)
                        if (pos_eng1>dt.current_depth) or (pos_eng2>dt.current_depth):
                            break
                        z1=list_Z[pos_eng1][dt.current_node[pos_eng1]]
                        z2=list_Z[pos_eng2][dt.current_node[pos_eng2]]
                        if self.gear_graph[eng1][eng2]['typ']=='gear_mesh':
                            demul[eng2]=demul[eng1]*z1/z2
                        else:
                            demul[eng2]=demul[eng1]
                        vsi_min=self.gear_speed[eng2][0] # Specified speed
                        vsi_max=self.gear_speed[eng2][1]
                        vai_min = max(v0_min*demul[eng2],vsi_min)
                        vai_max = min(v0_max*demul[eng2],vsi_max)
                        if vai_min>vai_max:
                            valid = False
                            break
                        v0_min, v0_max = vai_min/demul[eng2],vai_max/demul[eng2]
                    
                # frequncy analysis
#                if valid:
#                    for freq in self.frequency:
#                        zm=list_Z[0][dt.current_node[0]]
#                        for engr_num,engr_ind in enumerate(dt.current_node):
#                            z=list_Z[engr_num][dt.current_node[engr_num]]
#                            f_min=(2*npy.pi*v0_min*zm/z)/z
#                            f_max=(2*npy.pi*v0_max*zm/z)/z
#                            if (max(f_min,f_max)>freq[0]) and (min(f_min,f_max)<freq[1]):
#                                valid=False
#                                break


            elif (dt.current_depth>(nb_gear-1)) and (valid==True):
                # rack feasibility analysis
                rack_num=list_rack[dt.current_node[-1]]
                rack_pos=list_node[dt.current_depth-nb_gear]
                if rack_num not in self.rack_choice[rack_pos]:
                    valid=False
                    
            if (dt.current_depth==(nb_gear+nb_gear-1)) and (valid==True):
                # feasibility analysis of the modulus toward center-distance
                liste_DF_min={}
                module_minmax={}
                module_inf,module_sup=(0,npy.inf)
                # intersection between all module specified by each rack
                for tree_pos,tree_val in enumerate(dt.current_node[0:nb_gear]):
                    engr_num=list_node[tree_pos]
                    rack_num=list_rack[dt.current_node[tree_pos+nb_gear]]
                    module_minmax[engr_num]=self.rack_list[rack_num]['module']
                    module_inf,module_sup=(max(module_inf,module_minmax[engr_num][0]),min(module_sup,module_minmax[engr_num][1]))
                for tree_pos,tree_val in enumerate(dt.current_node[0:nb_gear]):
                    z=list_Z[tree_pos][tree_val]
                    engr_num=list_node[tree_pos]
                    liste_DF_min[engr_num]=z*module_inf
                liste_pente_cd_module=[]
                for set_num,gs in enumerate(self.connections):
                    cd_min=-npy.inf
                    for eng1,eng2 in gs:
                        cd_min=max(cd_min,(liste_DF_min[eng1]+liste_DF_min[eng2])/2)

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
                db={}
                for s_graph in self.sub_graph_dfs:
                    for ind,(eng1,eng2) in enumerate(s_graph):
                        num_mesh=0 #position of the gear mesh (need to define position in transverse_pressure_angle)
                        for num_cd_temp,list_connections in enumerate(self.connections):
                            if (eng1,eng2) in list_connections:
                                num_cd=num_cd_temp #num_cd define number of center-distance to use in the connections list
                                num_mesh+=list_connections.index((eng1,eng2))
                            elif (eng2,eng1) in list_connections:
                                num_cd=num_cd_temp
                                num_mesh+=list_connections.index((eng2,eng1))
                            else:
                                num_mesh+=len(list_connections)
                        if len(cd_minmax_nv)>num_cd:
                            cd_min,cd_max=cd_minmax_nv[num_cd]
                            z1_pos=list_node.index(eng1)
                            z1=list_Z[z1_pos][dt.current_node[z1_pos]]
                            z2_pos=list_node.index(eng2)
                            z2=list_Z[z2_pos][dt.current_node[z2_pos]]
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
        
            if (dt.current_depth==(nb_gear+nb_gear-1)) and (valid==True):
                gear={}
                rack={}
                for n in list_node:
                    i=list_node.index(n)
                    gear[n]=list_Z[i][dt.current_node[i]]
                    rack[n]=list_rack[dt.current_node[i+nb_gear]]
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
#                if verbose:
#                    print('Mesh sections: {}'.format(self.solutions[-1].gear_width))
#                    print('Numbers of teeth: {}'.format(plex['Z']))
#                    print('Center distances: {}'.format(self.solutions[-1].center_distance))
                if compt_nb_sol==nb_sol:
                    break

    def OptimizeCD(self, nb_sol = 1, verbose = False,
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
#            if compt_fonct_mod==5*nb_sol:
#                break
            compt_fonct_mod+=1
        list_fonct_entraxe=npy.array(fonct_entraxe)
        sort_list_fonct_entraxe=npy.argsort(list_fonct_entraxe)
        plex_analyse=[]
        for ind_plex in range(min(5*nb_sol, len(self.fonctionnel))):
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
                
                for num_graph,list_sub_graph in enumerate(ga.sub_graph_dfs):
                    num_cd_sub_graph=0
                    for num_cd,list_connection in enumerate(self.connections):
                        for num_mesh_iter,(eng1,eng2) in enumerate(list_connection):
                            if ((eng1,eng2) in list_sub_graph) or ((eng2,eng1) in list_sub_graph):
                                if (solutions.center_distance[num_cd])<(self.center_distance[num_cd][0]*0.999):
                                    valid_cd=False
                                elif (solutions.center_distance[num_cd])>(self.center_distance[num_cd][1]*1.001):
                                    valid_cd=False
                                num_cd_sub_graph+=1
#                for engr_num,cd in enumerate(self.center_distance):
#                    if (solutions.center_distance[engr_num])<(cd[0]*0.999):
#                        valid_cd=False
#                    elif (solutions.center_distance[engr_num])>(cd[1]*1.001):
#                        valid_cd=False
                if valid_cd:
                    liste_solutions.append(solutions)
                    compt_nb_sol+=1
                    if verbose:
                        print('valid solution n°{}'.format(compt_nb_sol))
                        for num_graph in range(len(solutions.mesh_assembly)):
                            print('Module gear set n°{}: {}'.format(num_graph,solutions.mesh_assembly[num_graph].meshes[list(solutions.mesh_assembly[num_graph].meshes.keys())[0]].rack.module))
                    
                else:
                    if verbose: 
                        print('unvalid solution')
            else:
                if verbose:
                    print('Convergence Problem')
            if compt_nb_sol==nb_sol:
                break
        # admissible solution are sorted toward the updated functional
        list_fonctionnel_nv=[]
        for solutions in liste_solutions:
            fonctionnel=0
            
            for num_cd,cd in enumerate(self.center_distance):
                fonctionnel+=(solutions.center_distance[num_cd]-self.center_distance[num_cd][0])**2

            list_fonctionnel_nv.append(fonctionnel)
        for indice_plex,num_plex in enumerate(npy.argsort(npy.array(list_fonctionnel_nv))):
            self.solutions.append(liste_solutions[num_plex])
            if verbose:
                print('solution sorted n°{}'.format(indice_plex))
                for num_graph in range(len(liste_solutions[num_plex].mesh_assembly)):
                    print('Mesh sections gear set n°{}: {}'.format(num_graph,self.solutions[-1].mesh_assembly[num_graph].gear_width))
                    print('Center distances gear set n°{}: {}'.format(num_graph,self.solutions[-1].mesh_assembly[num_graph].center_distance))
                    print('Module gear set n°{}: {}'.format(num_graph,self.solutions[-1].mesh_assembly[num_graph].meshes[list(self.solutions[-1].mesh_assembly[num_graph].meshes.keys())[0]].rack.module))
      
def gear_graph_simple(connections):
    # NetworkX graph construction
    list_gear=[] # list of all gears
    compt_mesh=0 # number of gear mesh
    for gs in connections:
        for (eng1,eng2) in gs:
            compt_mesh+=1
            if eng1 not in list_gear:
                list_gear.append(eng1)
            if eng2 not in list_gear:
                list_gear.append(eng2)
    # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
    gear_graph=nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge)                  
    sub_graph=list(nx.connected_component_subgraphs(gear_graph))
    sub_graph_dfs=[]
    for s_graph in sub_graph:
        node_init=list(s_graph.nodes())[0]
        sub_graph_dfs.append(list(nx.dfs_edges(s_graph,node_init)))
    return sub_graph_dfs,list_gear

def gear_graph_complex(connections,strong_link):
    # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
    _,list_gear=gear_graph_simple(connections)
    gear_graph=nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge,typ='gear_mesh')
        li_shaft1=[]
        li_shaft2=[]
        for eng1,eng2 in list_edge:
            li_shaft1.append(eng1)
            li_shaft2.append(eng2)
        if len(li_shaft1)>1:
            for pos_gear,num_gear in enumerate(li_shaft1[1:]):
                valid_strong_ling=False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft1[pos_gear] in list_strong_link):
                        valid_strong_ling=True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear,li_shaft1[pos_gear])],typ='same_speed')
                else:
                    gear_graph.add_edges_from([(num_gear,li_shaft1[pos_gear])],typ='same_shaft')
        if len(li_shaft2)>1:
            for pos_gear,num_gear in enumerate(li_shaft2[1:]):
                valid_strong_ling=False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft2[pos_gear] in list_strong_link):
                        valid_strong_ling=True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear,li_shaft2[pos_gear])],typ='same_speed')                   
                else:
                    gear_graph.add_edges_from([(num_gear,li_shaft2[pos_gear])],typ='same_shaft')                   
    connections_dfs=list(nx.dfs_edges(gear_graph,list_gear[0]))
    # construction of a graph without same_shaft attribute
    gear_graph_kinematic=copy.deepcopy(gear_graph)
    for edge,typ in nx.get_edge_attributes(gear_graph_kinematic,'typ').items():
        if typ=='same_shaft':
            gear_graph_kinematic.remove_edges_from([edge])
    connections_kinematic_dfs=list(nx.dfs_edges(gear_graph_kinematic,list_gear[0]))
    return connections_dfs,connections_kinematic_dfs,gear_graph
