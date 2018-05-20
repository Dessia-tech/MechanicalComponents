import numpy as npy
import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
#import cma
#from sympy import *
import itertools
from jinja2 import Environment, PackageLoader, select_autoescape
import networkx as nx
#from interval import interval, inf, imath

import mechanical_components.LibSvgD3 as LibSvg
import powertransmission.tools as tools

import persistent
#from dessia_common import ResultsDBClient
#import pyDOE

class Rack():
    def __init__(self,transverse_pressure_angle=None,transverse_radial_pitch=None):
        if transverse_pressure_angle==None:
            self.transverse_pressure_angle=20*npy.pi/180

    def RackParam(self,transverse_pressure_angle,transverse_radial_pitch):

        self.transverse_radial_pitch=transverse_radial_pitch
        self.module=self.transverse_radial_pitch/math.pi
        self.gear_addendum=self.module
        self.gear_dedendum=1.25*self.module
        self.root_radius=0.38*self.module
        self.circular_tooth_thickness=transverse_radial_pitch/2

        self.tooth_space=self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth=self.gear_addendum+self.gear_dedendum
        self.clearance=self.root_radius-self.root_radius*npy.sin(self.transverse_pressure_angle)

        #paramètre pour la trochoide
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle)-self.root_radius*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle)/(npy.sin(self.transverse_pressure_angle))))
        self.b=self.gear_dedendum-self.root_radius

    def Update(self,transverse_pressure_angle,transverse_radial_pitch):
        self.RackParam(transverse_pressure_angle,transverse_radial_pitch)

    def Dict(self):
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        return d

    def CSVExport(self):
        d=self.__dict__.copy()
        return list(d.keys()),list(d.values())

class Gear():
    def __init__(self,z,db,cp):
        self.Z=z
        self.DB=db
        self.rack=Rack()
        transverse_pressure_angle_rack=self.rack.transverse_pressure_angle
        self.DFF=self.DB/npy.cos(transverse_pressure_angle_rack)
        transverse_radial_pitch_rack=npy.pi*self.DFF/self.Z
        self.rack.Update(transverse_pressure_angle_rack,transverse_radial_pitch_rack)
        self.coefficient_profile_shift=cp
        
        self.outside_diameter=self.DFF+2*(self.rack.gear_addendum+self.rack.module*self.coefficient_profile_shift)
        self.alpha_outside_diameter=npy.arccos(self.DB/self.outside_diameter)
        self.root_diameter=self.DFF-2*(self.rack.gear_dedendum-self.rack.module*self.coefficient_profile_shift)
        self.root_diameter_active=self._RootDiameterActive()
        
        self.alpha_pitch_diameter=npy.arccos(self.DB/self.DFF)
        self.circular_tooth_thickness=self.rack.circular_tooth_thickness+self.rack.module*self.coefficient_profile_shift*npy.tan(self.rack.transverse_pressure_angle)+self.rack.module*self.coefficient_profile_shift*npy.tan(self.rack.transverse_pressure_angle)
        self.outside_active_angle=2*self.circular_tooth_thickness/self.DFF-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter-npy.tan(self.alpha_pitch_diameter)+self.alpha_pitch_diameter)
        
        
    def GearSection(self,diameter):
        #epaisseur de la dent au diameter
        alpha_diameter=npy.arccos(self.DB/diameter)
        theta1=(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(npy.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)
        
    def _RootDiameterActive(self):
        #Analyse diam pied de dent actif
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
#        rho=self.rack.root_radius_T
        phi=-(a+b*npy.tan(npy.pi/2-self.rack.transverse_pressure_angle))/r
        root_diameter_active=2*norm(self._Trochoide(phi))
        return root_diameter_active
        
    def _Trochoide(self,phi,ind='T'):
        if ind=='T':
            drap=1
        else:
            drap=-1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
        rho=self.rack.root_radius
        x2=rho*npy.sin(npy.arctan((a-r*phi)/b)-phi)+a*npy.cos(phi)-b*npy.sin(phi)+r*(npy.sin(phi)-phi*npy.cos(phi))
        y2=-rho*npy.cos(npy.arctan((a-r*phi)/b)-phi)-a*npy.sin(phi)-b*npy.cos(phi)+r*(npy.cos(phi)+phi*npy.sin(phi))
        sol=(y2,x2)
        return sol

class GearAssembly():
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph):
        self.center_distance=center_distance
        self.gear_set=gear_set
        self.transverse_pressure_angle_0=transverse_pressure_angle
        self.gear_graph=gear_graph

        self.DF,DB,self.transverse_pressure_angle,self.node_dfs=self.GearAssemblyParam(Z)

        #instantiation des objets Gears
        self.gears={}
        for g,z in Z.items():
            db=DB[g]
            cp=coefficient_profile_shift[str(g)]
            self.gears[g]=Gear(z,db,cp)
            
        self.linear_backlash=[]
        for g1,g2 in self.node_dfs:
            circular_tooth_thickness1=self.gears[g1].GearSection(self.DF[g1])
            circular_tooth_thickness2=self.gears[g2].GearSection(self.DF[g2])
            transverse_radial_pitch1=npy.pi*self.DF[g1]/self.gears[g1].Z
            space_width1=transverse_radial_pitch1-circular_tooth_thickness1
            space_width2=transverse_radial_pitch1-circular_tooth_thickness2    
            self.linear_backlash.append(min(space_width1-circular_tooth_thickness2,space_width2-circular_tooth_thickness1))

    def GearAssemblyParam(self,Z):
        DF={}
        DB={}
        transverse_pressure_angle=[self.transverse_pressure_angle_0]
        for gs,cd in zip(self.gear_set,self.center_distance):
            Z1=Z[gs[0]]
            Z2=Z[gs[1]]
            #définition DF
            DF1=2*cd*Z1/Z2/(1+Z1/Z2)
            DF2=2*cd-DF1
            DF[gs[0]]=DF1
            DF[gs[1]]=DF2
        DF1=DF[self.gear_set[0][0]]
        DF2=DF[self.gear_set[0][1]]
        DB1=DF1*npy.cos(self.transverse_pressure_angle_0)
        DB2=DF2*npy.cos(self.transverse_pressure_angle_0)
        DB[self.gear_set[0][0]]=DB1
        DB[self.gear_set[0][1]]=DB2
        node_dfs=list(nx.edge_dfs(self.gear_graph, [self.gear_set[0][0],self.gear_set[0][1]]))
        for nd in node_dfs[1:]:
            DB1=DB[nd[0]]
            DF1=DF[nd[0]]
            DF2=DF[nd[1]]
            transverse_pressure_angle.append(npy.arccos(DB1/DF1))
            DB2=DF2*npy.cos(transverse_pressure_angle[-1])
            DB[nd[1]]=DB2
        return DF,DB,transverse_pressure_angle,node_dfs


class ContinuousGearAssemblyOptimizer:
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph):
        self.xk={'Z':Z,'gear_set':gear_set,'gear_graph':gear_graph}
        self.xu={'center_distance':self._init_list(center_distance),'transverse_pressure_angle':self._init_list(transverse_pressure_angle)[0],'coefficient_profile_shift':self._init_item(coefficient_profile_shift)}
        self.xt=dict(list(self.xk.items())+list(self.xu.items()))
        self.GearAssembly=GearAssembly(**self.xt)

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
    


class GearAssemblyOptimizer:
    def __init__(self,gear_width,center_distance,coefficient_profile_shift,gear_speed,Z,gear_set,transverse_pressure_angle,helix_angle,frequency):
        self.Z=Z['data']
        self.gear_set=gear_set['data']
        self.gear_speed=gear_speed['data']
        self.frequency=frequency['data']
        self.center_distance=center_distance['data']
        self.transverse_pressure_angle=transverse_pressure_angle['data']
        self.coefficient_profile_shift=coefficient_profile_shift['data']

        self.nb_gear=npy.array(self.gear_set).max()+1
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(range(self.nb_gear))
        gear_graph.add_edges_from(self.gear_set)
        self.gear_graph=gear_graph

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.node_dfs=list(nx.dfs_edges(gear_graph,self.node_init))

        self.AnalyzeCombination()

    def AnalyzeCombination(self):
        np=[80]*self.nb_gear
        liste_gear=npy.arange(20,20+80)

        demul_int_min=1/1.9
        demul_int_max=1.9

        dt=tools.RegularDecisionTree(np)
        node=dt.current_node
        n1=self.node_init
        liste_node=[n1]

        for (n1,n2) in self.node_dfs:
            if n2 not in liste_node:
                liste_node.append(n2)
        incr=0
        self.plex_calcul=[]

        def pgcd(a,b) :
            while a%b != 0 :
                a, b = b, a%b
            return b

        while not dt.finished:
            valid=True
            #Analyse du Z initial
            if dt.current_depth==0:
                z=liste_gear[dt.current_node[0]]
                if (z<self.Z[str(liste_node[0])][0]) or (z>self.Z[str(liste_node[0])][1]):
                    valid=False
            if dt.current_depth>0:
                (n1,n2)=self.node_dfs[dt.current_depth-1]
                i1=liste_node.index(n1)
                i2=liste_node.index(n2)
                z1=liste_gear[dt.current_node[i1]]
                z2=liste_gear[dt.current_node[i2]]
                #analyse ACV engrenage 2 à 2
                if pgcd(z1,z2)!=1:
                    valid=False
                #analyse ACV de l'ensemble des engrenages entre eux
                if valid:
                    for n in liste_node[0:dt.current_depth]:
                        i=liste_node.index(n)
                        z=liste_gear[dt.current_node[i]]
                        if pgcd(z,z2)!=1:
                            valid=False
                #analyse demul interne
                if valid:
                    demul=liste_gear[dt.current_node[i1]]/liste_gear[dt.current_node[i2]]
                    if (demul > demul_int_max) or (demul < demul_int_min):
                        valid=False
                #Analyse des bornes sur Z
                if valid:
                    if (z2<self.Z[str(n2)][0]) or (z2>self.Z[str(n2)][1]):
                        valid=False
                #analyse des vitesses du CDC
                if valid:
                    v1=self.gear_speed[str(liste_node[0])][0]
                    v2=self.gear_speed[str(liste_node[0])][1]
                    for n in liste_node[1:dt.current_depth+1]:
                        if valid==True:
                            i=liste_node.index(n)
                            if str(n) in self.gear_speed.keys():
                                demul=liste_gear[dt.current_node[0]]/liste_gear[dt.current_node[i]]
                                v1p=self.gear_speed[str(n)][0]/demul
                                v2p=self.gear_speed[str(n)][1]/demul
                                v1=max(v1,v1p)
                                v2=min(v2,v2p)
                                if (v1>v2):
                                    valid=False
                #analyse frequence
                if valid:
                    for freq in self.frequency:
                        zm=liste_gear[dt.current_node[0]]
                        for i in dt.current_node:
                            f1=(60*v1*zm/liste_gear[i])/liste_gear[i]
                            f2=(60*v2*zm/liste_gear[i])/liste_gear[i]
                            zm=liste_gear[i]
                            if (max(f1,f2)>freq[0]) and (min(f1,f2)<freq[1]):
                                valid=False
            if (dt.current_depth==(self.nb_gear-1)) & (valid==True):
                sol={}
                for n in liste_node:
                    i=liste_node.index(n)
                    sol[n]=liste_gear[dt.current_node[i]]
                Temp={}
                Temp['Z']=sol
                Temp['gear_set']=self.gear_set
                Temp['center_distance']=self.center_distance
                Temp['transverse_pressure_angle']=self.transverse_pressure_angle
                Temp['coefficient_profile_shift']=self.coefficient_profile_shift

                self.plex_calcul.append(Temp)
                incr+=1
            dt.NextNode(valid)
        print('Nombre de combinaison trouvées: {}'.format(incr))


    def Optimize(self,callback=lambda x:x):
        lpx=len(self.plex_calcul)
        for ii,i in enumerate(self.plex_calcul):
            # print('{}%'.format(ii/lpx*100))
            callback(ii/lpx)
            plex=i
            plex['gear_graph']=self.gear_graph
            A1=ContinuousGearAssemblyOptimizer(**plex)
            # A1.Optimize()
            # try:
            #     xsol=npy.transpose([A1.solutions[-1]])
            #     A1.DefXU(xsol)
            #     xt=dict(list(A1.xk.items())+list(A1.xu.items())+list(A1.xo.items()))
            #     self.solutions.append(GearAssembly(**xt))
            # except:
            #     pass


class GearAssemblyOptimizerWizard:

    def __init__(self,data):
        dico={}
        for lign in data:
            dico[lign['type']]={}
            for key,val in lign.items():
                if 'type' not in key:
                    dico[lign['type']][key]=val

        #passage en unité SI
        dico['center_distance']['data']=tuple([[0.001*p[0],0.001*p[1]] for p in dico['center_distance']['data']])
        dico['transverse_pressure_angle']['data']=tuple([[npy.pi/180*p[0],npy.pi/180*p[1]] for p in dico['transverse_pressure_angle']['data']])
        dico['helix_angle']['data']={k:[npy.pi/180*dico['helix_angle']['data'][k][0],npy.pi/180*dico['helix_angle']['data'][k][1]] for k in dico['helix_angle']['data'].keys()}
        dico['gear_speed']['data']={k:[2*npy.pi/60*dico['gear_speed']['data'][k][0],2*npy.pi/60*dico['gear_speed']['data'][k][1]] for k in dico['gear_speed']['data'].keys()}

        dico=self.DefaultDataSet(dico)
        self.dico_gear_assembly_optimizer=dico


    def Optimize(self):

        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        M1.Optimize()
#        self.solutions=M1.solutions

    def DefaultDataSet(self,dico):

        #consolidation Z et coefficient_profile_shift
        nb_gear=npy.array(dico['gear_set']['data']).max()+1
        for i in range(nb_gear):
            if str(i) not in dico['Z']['data'].keys():
                dico['Z']['data'][str(i)]=[20,80]
            if str(i) not in dico['coefficient_profile_shift']['data'].keys():
                dico['coefficient_profile_shift']['data'][str(i)]=[-1.2,1.2]

        #consolidation gear_width
        nb_gear_set=len(dico['gear_set']['data'])
        parcourt={}
        for i in range(nb_gear):
            parcourt[str(i)]=[]
        for it1,it2 in dico['gear_set']['data']:
            parcourt[str(it1)].append(str(it2))
            parcourt[str(it2)].append(str(it1))
        know_gw=dico['gear_width']['data'].keys()
        for i in range(nb_gear):
            if str(i) not in know_gw:
                temp=[0,0]
                for j in know_gw:
                    drap=0
                    ptm=j
                    ptmm=j
                    while drap==0:
                        ptp=[]
                        for pt in ptm:
                            ptp=ptp+parcourt[pt]
                            try:
                                ptp.remove(pt)
                            except:
                                pass
                        for pt in ptmm:
                            try:
                                ptp.remove(pt)
                            except:
                                pass
                        if str(i) in ptp:
                            drap=2
                        if ptp==ptmm:
                            drap=1
                        ptmm=ptm
                        ptm=ptp
                    if drap==2:
                        temp[0]=max(temp[0],dico['gear_width']['data'][j][0])
                        temp[1]=max(temp[1],dico['gear_width']['data'][j][1])
                if temp!=[0,0]:
                    dico['gear_width']['data'][str(i)]=temp
                else:
                    dico['gear_width']['data'][str(i)]=[0.02,0.03]
        return dico
