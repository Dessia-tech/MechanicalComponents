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
    def __init__(self,transverse_pressure_angle):
        self.transverse_pressure_angle=transverse_pressure_angle

    def RackParam(self,transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):

        self.transverse_pressure_angle=transverse_pressure_angle
        self.transverse_radial_pitch=self.module*math.pi
        self.gear_addendum=coeff_gear_addendum*self.module
        self.gear_dedendum=coeff_gear_dedendum*self.module
        self.root_radius=coeff_root_radius*self.module
        self.circular_tooth_thickness=coeff_circular_tooth_thickness*self.transverse_radial_pitch

        self.tooth_space=self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth=self.gear_addendum+self.gear_dedendum
        self.clearance=self.root_radius-self.root_radius*npy.sin(self.transverse_pressure_angle)

        #paramètre pour la trochoide
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle)-self.root_radius*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle)/(npy.sin(self.transverse_pressure_angle))))
        self.b=self.gear_dedendum-self.root_radius

    def Update(self,module,transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.module=module
        self.RackParam(transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        

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
    def __init__(self,z,db,cp,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.rack=Rack(transverse_pressure_angle_rack)
        self.GearParam(z,db,cp,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    def Update(self,z,db,cp,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.GearParam(z,db,cp,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    def GearParam(self,z,db,cp,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.Z=z
        self.DB=db
        self.DFF=self.DB/npy.cos(transverse_pressure_angle_rack)
        module_rack=self.DFF/self.Z
        self.rack.Update(module_rack,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        self.coefficient_profile_shift=cp
        
        self.outside_diameter=self.DFF+2*(self.rack.gear_addendum+self.rack.module*self.coefficient_profile_shift)
        self.alpha_outside_diameter=npy.arccos(self.DB/self.outside_diameter)
        self.root_diameter=self.DFF-2*(self.rack.gear_dedendum-self.rack.module*self.coefficient_profile_shift)
        self.root_diameter_active,self.phi_trochoide=self._RootDiameterActive()
        self.alpha_root_diameter_active=npy.arccos(self.DB/self.root_diameter_active)
        
        self.alpha_pitch_diameter=npy.arccos(self.DB/self.DFF)
        self.circular_tooth_thickness=self.rack.circular_tooth_thickness+self.rack.module*self.coefficient_profile_shift*npy.tan(self.rack.transverse_pressure_angle)+self.rack.module*self.coefficient_profile_shift*npy.tan(self.rack.transverse_pressure_angle)
        self.tooth_space=self.rack.transverse_radial_pitch-self.circular_tooth_thickness
        self.outside_active_angle=2*self.circular_tooth_thickness/self.DFF-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter-npy.tan(self.alpha_pitch_diameter)+self.alpha_pitch_diameter)
        self.base_circular_tooth_thickness=self.DB/2*(2*self.circular_tooth_thickness/self.DFF+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter))
        
        self.root_angle=self.tooth_space/(self.DFF/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
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
        return root_diameter_active,phi
    
    def GearContours(self,discret=10,list_number=[None]):
        #Analytical tooth profil
        if list_number==[None]:
            list_number=npy.arange(int(self.Z))
        L=[self._OutsideTrace(0)]
        L.append(self._InvoluteTrace(discret,0,'T'))
        L.append(self._TrochoideTrace(2*discret,0,'T'))
        L.append(self._RootCircleTrace(0))
        L.append(self._TrochoideTrace(2*discret,0,'R'))
        L.append(self._InvoluteTrace(discret,1,'R'))
        for i in list_number[1::]:
            L.append(self._OutsideTrace(i))
            L.append(self._InvoluteTrace(discret,i,'T'))
            L.append(self._TrochoideTrace(2*discret,i,'T'))
            L.append(self._RootCircleTrace(i))
            L.append(self._TrochoideTrace(2*discret,i,'R'))
            L.append(self._InvoluteTrace(discret,i+1,'R'))
        return L
    
    def _InvoluteTrace(self,discret,number,ind='T'):
        
        if ind=='T':
            drap=1
            theta=npy.linspace(npy.tan(self.alpha_outside_diameter),
                               npy.tan(self.alpha_root_diameter_active),discret)
        else:
            drap=-1
            theta=npy.linspace(npy.tan(self.alpha_root_diameter_active),
                               npy.tan(self.alpha_outside_diameter),discret)
        
#        theta=npy.linspace(0,npy.tan(self.alpha_outside_diameter),discret)
        sol=self._Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            self.rac=L.points[-1]
        else:
            L=ref.Rotation(vm.Point2D((0,0)),
                           self.base_circular_tooth_thickness*2/self.DB)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            L.points[0]=self.rac
        return L
        
    def _Involute(self,tan_alpha):
        
        sol=(self.DB/2*npy.cos(tan_alpha)+self.DB/2*tan_alpha*npy.sin(tan_alpha),
             self.DB/2*npy.sin(tan_alpha)-self.DB/2*tan_alpha*npy.cos(tan_alpha))
        return sol
    
    def _TrochoideTrace(self,discret,number,ind='T'):
        if ind=='T':
            drap=1
        else:
            drap=-1
        
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        phi0=a/(self.DFF/2)
        
        ref=[]
        if ind=='R':
            theta=npy.linspace(phi0,drap*self.phi_trochoide,discret)
        else:
            theta=npy.linspace(drap*self.phi_trochoide,phi0,discret)
        for t in theta:
            ref.append(vm.Point2D((self._Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        if ind=='T':
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            L1.points[0]=self.rac
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            self.rac=L1.points[-1]

        return L1
    
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
    
    def _RootCircleTrace(self,number):

        drap=1
        a=drap*self.rack.a
        phi0=a/(self.DFF/2)
        p1=vm.Point2D((self._Trochoide(phi0,'T')))
        p1=p1.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        drap=-1
        a=drap*self.rack.a
        phi0=a/(self.DFF/2)
        p2=vm.Point2D((self._Trochoide(phi0,'R')))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        ref=primitives2D.RoundedLines2D([p1,p2],{},False)
        L2=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
        return L2
    
    def _OutsideTrace(self,number):
        #trace du sommet des dents en arc de cercle
        theta4=npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*npy.cos(theta4),self.outside_diameter/2*npy.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        #ref=vm.Arc2D(p1,p2,p3)
        ref=primitives2D.RoundedLines2D([p3,p2,p1],{},False)
        L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
        return L
    
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
        
        d['rack']=self.rack.Dict()
        return d

class GearAssembly():
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness,list_gear):
        self.center_distance=center_distance
        self.gear_set=gear_set
        self.transverse_pressure_angle=transverse_pressure_angle
        self.gear_graph=gear_graph
        self.list_gear=list_gear

        self.DF,DB,self.node_dfs=self.GearAssemblyParam1(Z)
        #instantiation des objets Gears
        self.gears={}
        for ne,ns in enumerate(self.gear_set):
            self.gears[ne]={}
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[ng]
                ngp=self.list_gear.index(ng)
                tpa=transverse_pressure_angle_rack[ngp]
                cga=coeff_gear_addendum[ngp]
                cgd=coeff_gear_dedendum[ngp]
                crr=coeff_root_radius[ngp]
                cct=coeff_circular_tooth_thickness[ngp]
                self.gears[ne][ng]=Gear(z,db,cp,tpa,cga,cgd,crr,cct)
                
#        for g,z in Z.items():
#            db=DB[g]
#            cp=coefficient_profile_shift[str(g)]
#            self.gears[g]=Gear(z,db,cp)
        self.linear_backlash,self.radial_contact_ratio=self.GearAssemblyParam2(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
            
    def Update(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness,list_gear):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.DF,DB,self.node_dfs=self.GearAssemblyParam1(Z)
        self.linear_backlash,self.radial_contact_ratio=self.GearAssemblyParam2(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    def GearAssemblyParam2(self,Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        for ne,ns in enumerate(self.gear_set):
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[ng]
                ngp=self.list_gear.index(ng)
                tpa=transverse_pressure_angle_rack[ngp]
                cga=coeff_gear_addendum[ngp]
                cgd=coeff_gear_dedendum[ngp]
                crr=coeff_root_radius[ngp]
                cct=coeff_circular_tooth_thickness[ngp]
                self.gears[ne][ng].Update(z,db,cp,tpa,cga,cgd,crr,cct)
            
        linear_backlash=[]
        radial_contact_ratio=[]
        for g1,g2 in self.node_dfs:
            ne=self.gear_set.index((g1,g2))
            circular_tooth_thickness1=self.gears[ne][g1].GearSection(self.DF[ne][g1])
            circular_tooth_thickness2=self.gears[ne][g2].GearSection(self.DF[ne][g2])
            transverse_radial_pitch1=npy.pi*self.DF[ne][g1]/self.gears[ne][g1].Z
            space_width1=transverse_radial_pitch1-circular_tooth_thickness1
            space_width2=transverse_radial_pitch1-circular_tooth_thickness2    
            linear_backlash.append(min(space_width1-circular_tooth_thickness2,space_width2-circular_tooth_thickness1))
            transverse_pressure_angle1=self.transverse_pressure_angle[ne]
            center_distance1=self.center_distance[ne]
            radial_contact_ratio.append(1/2*(npy.sqrt(self.gears[ne][g1].outside_diameter**2-self.gears[ne][g1].DB**2)+npy.sqrt(self.gears[ne][g2].outside_diameter**2-self.gears[ne][g2].DB**2)-2*center_distance1*npy.sin(transverse_pressure_angle1))/(transverse_radial_pitch1*npy.cos(transverse_pressure_angle1)))
        return linear_backlash,radial_contact_ratio

    def GearAssemblyParam1(self,Z):
        DF={}
        DB={}
#        transverse_pressure_angle=[self.transverse_pressure_angle_0]
        for ne,(gs,cd) in enumerate(zip(self.gear_set,self.center_distance)):
            Z1=Z[gs[0]]
            Z2=Z[gs[1]]
            #définition DF
            DF1=2*cd*Z1/Z2/(1+Z1/Z2)
            DF2=2*cd-DF1
            DF[ne]={}
            DF[ne][gs[0]]=DF1
            DF[ne][gs[1]]=DF2
            DB1=float(DF1*npy.cos(self.transverse_pressure_angle[ne]))
            DB2=float(DF2*npy.cos(self.transverse_pressure_angle[ne]))
            DB[ne]={}
            DB[ne][gs[0]]=DB1
            DB[ne][gs[1]]=DB2
        
        node_dfs=list(nx.edge_dfs(self.gear_graph, [self.gear_set[0][0],self.gear_set[0][1]]))

        return DF,DB,node_dfs
    
    def GearAssemblyTrace(self,list_gear,list_center,list_rot):
        TG=[]
        for (i,j,k) in zip(list_gear,list_center,list_rot):
            temp=[]
            for m in i:
                temp1=m.Translation(j)
                temp2=temp1.Rotation(vm.Point2D(j),k)
                temp.append(temp2)
            TG.append(temp)
        return TG
    
    def InitialPosition(self,ne,ens):        
        fun = (lambda tan_alpha : (norm(self.gears[ne][ens[0]]._Involute(tan_alpha))-(self.center_distance[ne]-self.DF[ne][ens[1]]/2))**2)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle1=xsol[0]-npy.arctan(xsol[0])
        fun = (lambda tan_alpha : (norm(self.gears[ne][ens[1]]._Involute(tan_alpha))-(self.DF[ne][ens[1]]/2))**2)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle2=xsol[0]-npy.arctan(xsol[0])        
        Angle1=npy.arccos(self.gears[ne][ens[0]].DB/self.DF[ne][ens[0]])
        Angle2=npy.arccos(self.gears[ne][ens[1]].DB/self.DF[ne][ens[1]])
        Gear1Angle=-(npy.tan(Angle1)-Angle1)
        Gear2Angle=-(npy.tan(Angle2)-Angle2)+npy.pi        
        return [Gear1Angle,Gear2Angle]
    
#    def FreeCADExport(self,name,position,file_path,export_types):
#
#        RIM1=self.Gear1.RimContour()
#        RIM2=self.Gear2.RimContour()
##        RIM1.MPLPlot()
##        RIM2.MPLPlot()
#        gear1=primitives3D.RevolvedProfile(vm.Point3D((position1[0],position1[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
#                                           vm.Vector3D((0,1,0)),[RIM1],vm.Vector3D((position1[0],position1[1],0)),
#                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim1')
#        gear2=primitives3D.RevolvedProfile(vm.Point3D((position2[0],position2[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
#                                           vm.Vector3D((0,1,0)),[RIM2],vm.Vector3D((position2[0],position2[1],0)),
#                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim2')
#        
#        # Teeth
#        TG1=self.Gear1.GearContours(5)
#        TG2=self.Gear2.GearContours(5)
#        list_rot=self.InitialPosition()
#        L1=self.GearAssemblyTrace([TG1,TG2],[position1,position2],list_rot)
#        C1=vm.Contour2D(L1[0])
#        C2=vm.Contour2D(L1[1])
##        C1.MPLPlot()
#        h1=0.5*self.Gear1.alpha_rim*(self.Gear1.outside_diameter-self.Gear1.root_diameter)
#        h2=0.5*self.Gear2.alpha_rim*(self.Gear2.outside_diameter-self.Gear2.root_diameter)
#        r1=self.Gear1.root_diameter/2-h1/2
#        r2=self.Gear2.root_diameter/2-h2/2
##        print(h1,h2,r1,r2)
#        c1=primitives3D.Cylinder((position1[0],position1[1],0.48*self.gear_width),
#                                 (0,0,1),r1,1.05*self.gear_width)
#        c2=primitives3D.Cylinder((position2[0],position2[1],0.48*self.gear_width),
#                                 (0,0,1),r2,1.05*self.gear_width)
#        
##        C1int=vm.Contour2D([vm.Circle2D(vm.Point2D(position1),r1)])
##        C2int=vm.Contour2D([vm.Circle2D(vm.Point2D(position2),r2)])
#        if self.helix_angle==0.:            
#            t1=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
#                                            vm.Vector3D((0,1,0)),[C1],(0,0,self.gear_width))
#            t2=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
#                                            vm.Vector3D((0,1,0)),[C2],(0,0,self.gear_width))
#        else:
#            t1=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
#                                                   vm.Vector3D((0,1,0)),(position1[0],position1[1],0),
#                                                   (0,0,self.gear_width),self.DF1*mt.pi/mt.tan(self.helix_angle),
#                                                   C1)
#            t2=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
#                                                   vm.Vector3D((0,1,0)),(position2[0],position2[1],0),
#                                                   (0,0,self.gear_width),-self.DF2*mt.pi/mt.tan(self.helix_angle),
#                                                   C2)
#
#        # Creating holes in teeths
#        t1=primitives3D.Cut(t1,c1,name='Teeth1')
#        t2=primitives3D.Cut(t2,c2,name='Teeth2')
#
#        model=vm.VolumeModel([gear1,t1,gear2,t2])
##        model=vm.VolumeModel([gear1,t1,gear2])
#        model.FreeCADExport('python',file_path,'/usr/lib/freecad/lib',export_types)
    
    def SVGGearSet(self,name,position):
        #tuple1 et 2 correspondent a la position des centres
#        TG1=self.gears[0].GearContours(5)
#        G1=vm.Contour2D(TG1)
#        G1.MPLPlot()
            
        #optimisation pour le placement des axes des engrenages
        def fun(x):
            obj=0
            for num,it in enumerate(self.gear_set):
                eng1=(self.list_gear).index(it[0])
                eng2=(self.list_gear).index(it[1])
                obj+=(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num])**2
            return obj
        def eg(x):
            ine=[]
            for k,val in position.items():
                key=(self.list_gear).index(k)
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine=[]
            for num,it in enumerate(self.gear_set):
                eng1=(self.list_gear).index(it[0])
                eng2=(self.list_gear).index(it[1])
                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*self.center_distance[num])
                ine.append(1.001*self.center_distance[num]-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
            return ine
        cons = ({'type': 'eq','fun' : eg},{'type': 'ineq','fun' : ineg})
        drap=1
        while drap==1:
            x0=tuple(npy.random.random(2*self.gear_graph.number_of_nodes())*1)
            Bound=[[0,1]]*(self.gear_graph.number_of_nodes()*2)
            res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
            if (min(ineg(res.x))>0) and (max(eg(res.x))<1e-7):
                drap=0
        x_opt=res.x
        
        TG={}
        L1=[]
        Struct=[]
        Rot={}
        for num,en in enumerate(self.node_dfs):
            ens=[self.list_gear.index(en[0]),self.list_gear.index(en[1])]
            position1=(x_opt[2*ens[0]],x_opt[2*ens[0]+1])   
            position2=(x_opt[2*ens[1]],x_opt[2*ens[1]+1])
            #tuple1 et 2 correspondent a la position des centres
            ne=self.gear_set.index(en)
            Rot[ne]={}
            if num==0:
                TG[en[0]]=self.gears[ne][en[0]].GearContours(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2))
            TG[en[1]]=self.gears[ne][en[1]].GearContours(5)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[ne][en[1]]/2))
            #Definition de la position angulaire initiale
            list_rot=self.InitialPosition(ne,en)
            if position2[0]==position1[0]:
                if position2[1]-position1[1]>0:
                    angle=npy.pi/2
                else:
                    angle=-npy.pi/2
            else:
                angle=-npy.arctan((position2[1]-position1[1])/(position2[0]-position1[0]))
            if num==0:
                Rot[ne][en[0]]=list_rot[0]-angle
                Rot[ne][en[1]]=list_rot[1]-angle
            else:
                for k1,v1 in Rot.items():
                    if en[0] in v1.keys():
                        Rot[ne][en[0]]=v1[en[0]]
                        delta_rot=Rot[ne][en[0]]-(list_rot[0]-angle)
                Rot[ne][en[1]]=list_rot[1]-angle-delta_rot*((self.gears[ne][en[0]].Z)/(self.gears[ne][en[1]].Z))
            sol=self.GearAssemblyTrace([TG[en[0]],TG[en[1]]],[position1,position2],list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
            if num==0:
                L1.extend(sol[0])
            L1.extend(sol[1])
        L1.extend(Struct)
        G1=vm.Contour2D(L1)
        G1.MPLPlot()
        
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
        
        for ne,gs in enumerate(self.gear_set):
            d['Gear'+str(gs[0]+1)]=self.gears[ne][gs[0]].Dict()
            d['Gear'+str(gs[1]+1)]=self.gears[ne][gs[1]].Dict()

        return d

class ContinuousGearAssemblyOptimizer:
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,
                 gear_graph,cond_init,rack_list,rack_choice,list_gear):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.cond_init=cond_init
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear
        Bounds=list(center_distance)
        Bounds.extend(transverse_pressure_angle)
        tp=len(coefficient_profile_shift.keys())
        for i in self.list_gear:
            Bounds.append(coefficient_profile_shift[i])
        self.solutions=[]
        
        self.xi={'Z':Z,'gear_set':gear_set,'gear_graph':gear_graph,'list_gear':list_gear}
        self.xj,self.dict_xu=self._init()
        self.xt=dict(list(self.xi.items())+list(self.xj.items()))
        
        for k,v in self.dict_xu:
            Bounds.append(self.rack_list[v][k])
        self.Bounds=npy.array(Bounds)
        
        #xk partie non optimisée du vecteur x
        #xu partie optimisée du vecteur x
        
        self.GearAssembly=GearAssembly(**self.xt)
        
    def _init(self):
        xj={'center_distance':self._init_list(self.center_distance),
                 'transverse_pressure_angle':self._init_list(self.transverse_pressure_angle),
                 'coefficient_profile_shift':self._init_item(self.coefficient_profile_shift),
                 'transverse_pressure_angle_rack':[],'coeff_gear_addendum':[],
                 'coeff_gear_dedendum':[],'coeff_root_radius':[],'coeff_circular_tooth_thickness':[]}
        dict_xu=[]
        x_list={}
        for i in list(set(list(self.rack_choice.values()))):
            x_list[i]={'transverse_pressure_angle_rack':[],'coeff_gear_addendum':[],'coeff_gear_dedendum':[],'coeff_root_radius':[],'coeff_circular_tooth_thickness':[]}
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
        self.GearAssembly.Update(**self.xt)
        return xj
        
    def Fineq(self,X):
        xj=self._convert_Xu2xj(X)
        xj=self.Update(xj)
        ineq=[]
        #jeu inter-denture
        for lb in self.GearAssembly.linear_backlash:
            ineq.append(lb)
            ineq.append((4e-4)-lb)
        #contrainte géométrique
        for ne,gs in enumerate(self.GearAssembly.gear_set):
            dia1=self.GearAssembly.gears[ne][gs[0]].root_diameter_active
            dia2=self.GearAssembly.gears[ne][gs[1]].root_diameter_active
            de1=self.GearAssembly.gears[ne][gs[0]].outside_diameter
            de2=self.GearAssembly.gears[ne][gs[1]].outside_diameter
            cd=self.GearAssembly.center_distance[ne]
            ineq.append(cd-(de1/2+dia2/2))
            ineq.append(cd-(de2/2+dia1/2))
            oaa1=self.GearAssembly.gears[ne][gs[0]].outside_active_angle
            oaa2=self.GearAssembly.gears[ne][gs[1]].outside_active_angle
            ineq.append(oaa1)
            ineq.append(oaa2)
            df1=self.GearAssembly.DF[ne][gs[0]]
            df2=self.GearAssembly.DF[ne][gs[1]]
            db1=self.GearAssembly.gears[ne][gs[0]].DB
            db2=self.GearAssembly.gears[ne][gs[1]].DB
            ineq.append(df1-db1)
            ineq.append(df2-db2)
        #contrainte sur le RCA
        for ne,gs in enumerate(self.GearAssembly.gear_set):
            rca=self.GearAssembly.radial_contact_ratio[ne]
            ineq.append(rca-1)
        for ne,gs in enumerate(self.GearAssembly.gear_set):
            for g in gs:
                mo=self.GearAssembly.gears[ne][g].rack.module
                lm=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]<lm[1]:
                ineq.append(mo-lm[0])
                ineq.append(lm[1]-mo)
        return ineq
    
    def Feq(self,X):
        x=self._convert_Xu2xj(X)
        x=self.Update(x)
        eq=[]
#        for ng in range(self.GearAssembly.gear_graph.number_of_nodes()):
        for ng in self.list_gear:
            nel=list(self.GearAssembly.gear_graph.edges(ng))
            if len(nel)>1:
                list_db=[]
                for ne in nel:
                    if (ne[0],ne[1]) in self.GearAssembly.gear_set:
                        nes=self.GearAssembly.gear_set.index((ne[0],ne[1]))
                    elif (ne[1],ne[0]) in self.GearAssembly.gear_set:
                        nes=self.GearAssembly.gear_set.index((ne[1],ne[0]))
                    list_db.append(self.GearAssembly.gears[nes][ng].DB)
#                    print(self.GearAssembly.gears[nes][ng].DB,nes,ng,list_db)
                list1=itertools.combinations(list_db,2)
                for n1,n2 in list1:
                    eq.append(n1-n2)
#        for ne,gs in enumerate(self.GearAssembly.gear_set):
#            for g in gs:
#                mo=self.GearAssembly.gears[ne][g].rack.module
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
        for i in fineq:
            if i < 0:
                obj+=-1000*i
            else:
                obj+=0.0001*i
        for i in feq:
            if i<0:
                obj+=-1000*i
            else:
                obj+=1000*i
        return obj
    
    def Optimize(self):
        boucle=3
        i=0
        arret=0
        while i<boucle and arret==0:
#            print('Boucle d\'itération locale {}'.format(i))
            xj0,dict_xu=self._init()
            xj0=self.Update(xj0)
            X0=self._convert_xj2Xu(xj0)
            cons = ({'type': 'eq','fun' : self.Feq},{'type': 'ineq','fun' : self.Fineq})
#            cons = {'type': 'ineq','fun' : self.Fineq}
#            print(33,X0)
            cx = minimize(self.Objective, X0, bounds=self.Bounds,constraints=cons)
#            cx = minimize(self.Objective, X0, bounds=self.Bounds)
            Xsol=cx.x
            xsol=self._convert_Xu2xj(Xsol)
            xsol=self.Update(xsol)
#            print(i,cx.status,min(self.Fineq(Xsol)),max(npy.abs(self.Feq(Xsol))),self.Objective(Xsol))
#            print(self.Fineq(Xsol))
            if min(self.Fineq(Xsol))>-1e-5 and max(npy.abs(self.Feq(Xsol)))<1e-5:
#            if min(self.Fineq(Xsol))>-1e-5:
                self.solutions.append(xsol)
                arret=1
#                print('Convergence atteinte avec le statut {}, Valeur de la fonctionnelle {}'.format(cx.status,cx.fun))
            i=i+1


class GearAssemblyOptimizer:
    def __init__(self,gear_width,center_distance,coefficient_profile_shift,gear_speed,Z,gear_set,
                 transverse_pressure_angle,helix_angle,frequency,rack_list,rack_choice,
                 list_gear):
        self.Z=Z
        self.gear_set=gear_set
        self.gear_speed=gear_speed
        self.frequency=frequency
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear

        self.nb_gear=len(list_gear)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.gear_set)
        self.gear_graph=gear_graph
        
        self.nb_rack=len(self.rack_list.keys())

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.node_dfs=list(nx.dfs_edges(gear_graph,self.node_init))
        
        if self.Z=={}:
            self.Z=self.AnalyseZ()
        self.AnalyzeCombination()
        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]
        
    def AnalyseZ(self):
        #nombre de dents adaptatif
        Z=self.Z
        for i,gs in enumerate(self.gear_set):
            cd_min=self.center_distance[i][0]
            cd_max=self.center_distance[i][1]
            module1_min=self.rack_list[self.rack_choice[gs[0]]]['module'][0]
            module1_max=self.rack_list[self.rack_choice[gs[0]]]['module'][1]
            module2_min=self.rack_list[self.rack_choice[gs[1]]]['module'][0]
            module2_max=self.rack_list[self.rack_choice[gs[1]]]['module'][1]
            demul_min=self.gear_speed[gs[0]][0]/self.gear_speed[gs[1]][1]
            demul_max=self.gear_speed[gs[0]][1]/self.gear_speed[gs[1]][0]
            DF1_max=2*cd_max/(1+demul_min)
            Z1_max=int(DF1_max/module1_min)+1
            DF2_max=2*cd_max*demul_max/(1+demul_max)
            Z2_max=int(DF2_max/module2_min)+1
            DF1_min=2*cd_min/(1+demul_max)
            Z1_min=int(DF1_min/module1_max)-1
            DF2_min=2*cd_min*demul_min/(1+demul_min)
            Z2_min=int(DF2_min/module2_max)-1
            if gs[0] not in Z.keys():
                Z[gs[0]]=[Z1_min,Z1_max]
            else:
                Z[gs[0]]=[min(Z1_min,Z[gs[0]][0]),max(Z1_max,Z[gs[0]][1])]
            if gs[1] not in Z.keys():
                Z[gs[1]]=[Z2_min,Z2_max]
            else:
                Z[gs[1]]=[min(Z2_min,Z[gs[1]][0]),max(Z2_max,Z[gs[1]][1])]
        return Z

    def AnalyzeCombination(self):
        #recherche des Zmin et Zmax
        Zmin=npy.inf
        Zmax=0
        for n1,n2 in self.gear_set:
            if min(self.Z[n1][0],self.Z[n2][0])<Zmin:
                Zmin=min(self.Z[n1][0],self.Z[n2][0])
            if max(self.Z[n1][1],self.Z[n2][1])>Zmax:
                Zmax=max(self.Z[n1][1],self.Z[n2][1])
        np=[Zmax+1-Zmin]*self.nb_gear+[self.nb_rack]*self.nb_gear
        liste_gear=npy.arange(Zmin,Zmax+1)
        liste_rack=list(self.rack_list.keys())
        
        demul_int_min=1/1.9
        demul_int_max=1.9
        demul_int_min=1/3.
        demul_int_max=3

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
                if (z<self.Z[liste_node[0]][0]) or (z>self.Z[liste_node[0]][1]):
                    valid=False
            if dt.current_depth>0:
                if dt.current_depth<=(self.nb_gear-1):
#                    print(dt.current_node,self.node_dfs,dt.current_depth)
                    (n1,n2)=self.node_dfs[dt.current_depth-1]
                    i1=liste_node.index(n1)
                    i2=liste_node.index(n2)
                    z1=liste_gear[dt.current_node[i1]]
                    z2=liste_gear[dt.current_node[i2]]
                #analyse ACV engrenage 2 à 2
                if (pgcd(z1,z2)!=1) & (dt.current_depth<=(self.nb_gear-1)):
                    valid=False
                #analyse demul interne
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    demul=liste_gear[dt.current_node[i1]]/liste_gear[dt.current_node[i2]]
                    if (demul > demul_int_max) or (demul < demul_int_min):
                        valid=False
                #Analyse des bornes sur Z
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    if (z2<self.Z[n2][0]) or (z2>self.Z[n2][1]):
                        valid=False
                #analyse ACV de l'ensemble des engrenages entre eux
#                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
#                    for n in liste_node[0:dt.current_depth]:
#                        i=liste_node.index(n)
#                        z=liste_gear[dt.current_node[i]]
#                        if pgcd(z,z2)!=1:
#                            valid=False
                #analyse des vitesses du CDC
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    v1=self.gear_speed[liste_node[0]][0]
                    v2=self.gear_speed[liste_node[0]][1]
                    for n in liste_node[1:dt.current_depth+1]:
                        if valid==True:
                            i=liste_node.index(n)
                            if n in self.gear_speed.keys():
                                demul=liste_gear[dt.current_node[0]]/liste_gear[dt.current_node[i]]
                                v1p=self.gear_speed[n][0]/demul
                                v2p=self.gear_speed[n][1]/demul
                                v1=max(v1,v1p)
                                v2=min(v2,v2p)
                                if (v1>v2):
                                    valid=False
                #analyse frequence
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    for freq in self.frequency:
                        zm=liste_gear[dt.current_node[0]]
                        for i in dt.current_node:
                            f1=(60*v1*zm/liste_gear[i])/liste_gear[i]
                            f2=(60*v2*zm/liste_gear[i])/liste_gear[i]
                            zm=liste_gear[i]
                            if (max(f1,f2)>freq[0]) and (min(f1,f2)<freq[1]):
                                valid=False
#                #analyse faisabilité DF et DB
#                if (valid) & (dt.current_depth==(self.nb_gear-1)):
#                    #optimisation pour le placement des axes des engrenages
#                    def fun(x):
#                        obj=(1/x[0])**2
#                        return obj
#                    def ineg(x):
#                        db={}
#                        ine=[]
#                        for nes,(n1,n2) in enumerate(self.node_dfs):
#                            if (n1 in liste_node[0:dt.current_depth+1]) and (n2 in liste_node[0:dt.current_depth+1]):
#                                i1=liste_node.index(n1)
#                                i2=liste_node.index(n2)
#                                z1=liste_gear[dt.current_node[i1]]
#                                z2=liste_gear[dt.current_node[i2]]
#                                ne=self.gear_set.index((n1,n2))
#                                cd=x[ne+1]
#                                df1=2*cd*z1/z2/(1+z1/z2)
#                                df2=2*cd-df1
#                                if nes==0:
#                                    db[n1]=df1*npy.cos(x[0])
#                                db[n2]=df2/df1*db[n1]
#                                ine.append(df1-db[n1])
#                                ine.append(df2-db[n2])
#                        return ine
#                    cons = ({'type': 'ineq','fun' : ineg})
#                    drap=1
#                    boucle=2
#                    i=0
#                    while drap==1 and i<boucle:
#                        Bound=list([self.transverse_pressure_angle[0]])+list(self.center_distance)
#                        Bound_npy=npy.array(Bound)
#                        x0=(Bound_npy[:,1]-Bound_npy[:,0])*npy.random.random(len(Bound_npy[:,1]))+Bound_npy[:,0]
#                        res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
#                        if (min(ineg(res.x))>0):
#                            drap=0
#                        i+=1
#                    if drap==1:
#                        valid=False
                #Analyse de la faisabilité des cremailleres
#                if (valid) & (dt.current_depth>(self.nb_gear-1)):
#                    r1=liste_rack[dt.current_node[-1]]
#                    e1=dt.current_depth-self.nb_gear
#                    if str(r1) not in self.rack_choice[str(e1)]:
#                        valid=False
                #analyse coherence DB, demul et angle de pression de la cremaillere
#                if (valid) & (dt.current_depth==(self.nb_gear+self.nb_gear-1)):
#                    for ne,ns in enumerate(self.gear_set):
#                        e1=liste_node.index(ns[0])
#                        r1=liste_rack[dt.current_node[self.nb_gear+e1]]
#                        transverse_pressure_angle=self.rack_list[r1]['transverse_pressure_angle_rack']
#                        Z1=liste_gear[dt.current_node[e1]]
#                        module=self.rack_list[r1]['module']
#                        DB1_min=npy.cos(transverse_pressure_angle[0])*Z1*module[0]
#                        DB1_max=npy.cos(transverse_pressure_angle[1])*Z1*module[1]
#                        e2=liste_node.index(ns[1])
#                        r2=liste_rack[dt.current_node[self.nb_gear+e2]]
#                        transverse_pressure_angle=self.rack_list[r2]['transverse_pressure_angle_rack']
#                        Z2=liste_gear[dt.current_node[e2]]
#                        module=self.rack_list[r2]['module']
#                        DB2_min=npy.cos(transverse_pressure_angle[0])*Z2*module[0]
#                        DB2_max=npy.cos(transverse_pressure_angle[1])*Z2*module[1]
#                        if ((Z1/Z2)<(DB1_min/DB2_max)) or ((Z1/Z2)>(DB1_max/DB2_min)):
#                            valid=False
                #Analyse
                
                        
                        
#                        
#                    self.DFF=self.DB/npy.cos(transverse_pressure_angle_rack)
#                    transverse_radial_pitch_rack=npy.pi*self.DFF/self.Z
#                    DB=npy.cos(transverse_pressure_angle_rack)*self.Z*module
        
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)) & (valid==True):
                gear={}
                rack={}
                for n in liste_node:
                    i=liste_node.index(n)
                    gear[n]=liste_gear[dt.current_node[i]]
                    rack[n]=liste_rack[dt.current_node[i+self.nb_gear]]
                    
                Temp={}
                Temp['Z']=gear
                Temp['gear_set']=self.gear_set
                Temp['center_distance']=self.center_distance
                Temp['transverse_pressure_angle']=self.transverse_pressure_angle
                Temp['coefficient_profile_shift']=self.coefficient_profile_shift
#                Temp['cond_init']=res.x
                Temp['cond_init']=0
                Temp['rack_list']=self.rack_list
                Temp['rack_choice']=rack
                Temp['list_gear']=self.list_gear
                self.plex_calcul.append(Temp)
                incr+=1
#                if incr==3:
#                    break
            dt.NextNode(valid)
        print('Nombre de combinaison trouvées: {}'.format(incr))


    def Optimize(self,callback=lambda x:x):
        lpx=len(self.plex_calcul)
        for ii,i in enumerate(self.plex_calcul):
#            print('{}%'.format(ii/lpx*100))
            callback(ii/lpx)
            plex=i
            plex['gear_graph']=self.gear_graph
            A1=ContinuousGearAssemblyOptimizer(**plex)
#            try:
            A1.Optimize()
#            except:
#                print('Problème de convergence')
            if len(A1.solutions)>0:
                xsol=A1.solutions[-1]
#                print(11,xsol)
#                print(22,A1.xj)
                xt=dict(list(A1.xi.items())+list(xsol.items()))
                self.solutions.append(GearAssembly(**xt))
#                break
                
    def SearchCenterLine(self,nb_sol,callback=lambda x:x):
        #recherche de l'ensemble des entraxes
        search1=[]
        for ipl,pl in enumerate(self.plex_calcul):
            pl['DF']={}
            pl['Z_data']={}
            for k,v in pl['Z'].items():
                nr=pl['rack_choice'][k]
                mod=pl['rack_list'][nr]['module']
                module=(mod[0]+mod[1])/2
                pl['DF'][k]=module*v
                pl['Z_data'][k]=[v,v]
            pl['center_distance']=[]
            for ne in self.gear_set:
                t1=(pl['DF'][ne[0]]+pl['DF'][ne[1]])/2
                pl['center_distance'].append([t1,t1])
            del pl['DF']
            search1.append(pl)
        #selection des solutions avec entraxes compatible au CDC
        compt=0
        search2=[]
        for i,el in enumerate(search1):
            Z_data=search1[i]['Z_data'].copy()
            cd_data=search1[i]['center_distance'].copy()
            valid=True
            cd_input=[]
            fonctionnel=0
            for j,cd in enumerate(self.center_distance):
                if (cd_data[j][0]*0.99)<(cd[0]):
                    valid=False
                fonctionnel+=cd_data[j][0]-cd[0]
                cd_input.append([cd_data[j][0]*(0.99),cd_data[j][0]*1.1])
            if valid:
                search2.append([fonctionnel,cd_input,Z_data])
                compt+=1
        print('Nombre de solution avec entraxe supérieur au CDC:',compt)
        #Optimisation des nb_sol premières solutions
        search2_np=npy.array(search2)
        compt=0
        for ind in npy.argsort(search2_np[:,0]):
            cd_input=search2_np[ind][1]
            Z_data=search2_np[ind][2]
            ga=GearAssemblyOptimizerWizard(gear_set=self.gear_set,gear_speed=self.gear_speed,
                                            center_distance=cd_input,Z=Z_data,rack_list=self.rack_list,
                                            rack_choice=self.rack_choice)
            ga.Optimize()
            if len(ga.solutions)>0:
                valid2=True
                for j,v in enumerate(self.center_distance):
                    if (ga.solutions[-1].center_distance[j]*1e3)<(v[0]):
                        valid2=False
                if valid2:
                    self.solutions_search.append(ga.solutions[-1])
                    compt+=1
                    if compt==nb_sol:
                        break


class GearAssemblyOptimizerWizard:

    def __init__(self,gear_set,gear_speed,center_distance,Z={},transverse_pressure_angle=None,
                 helix_angle=None,gear_width=None,frequency=[[0,0]],coefficient_profile_shift=None,
                 rack_list=None,rack_choice=None):

        list_gear=[]
        for gs in gear_set:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        nb_gear=len(list_gear)
        nb_set=len(gear_set)
                    
        if transverse_pressure_angle==None:
            transverse_pressure_angle=[]
            for i in range(nb_set):
                transverse_pressure_angle.append([15/180*npy.pi,30/180*npy.pi])
            
        if helix_angle==None:
            helix_angle={list_gear[0]:[15/180*npy.pi,30/180*npy.pi]}
        
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
            coefficient_profile_shift={list_gear[0]:[-1.2,1.2]}
        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne]=[-1.2,1.2]
                
        if rack_list==None:
            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,2*1e-3],'transverse_pressure_angle_rack':[2*npy.pi0,20*npy.pi],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
            
        if rack_choice==None:
            rack_choice={list_gear[0]:list(rack_list.keys())[0]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne]=list(rack_list.keys())[0]
                
        dico={'gear_set':gear_set,'gear_speed':gear_speed,'center_distance':center_distance,
              'Z':Z,'transverse_pressure_angle':transverse_pressure_angle,'helix_angle':helix_angle,
              'gear_width':gear_width,'frequency':frequency,'coefficient_profile_shift':coefficient_profile_shift,
              'rack_list':rack_list,'rack_choice':rack_choice,'list_gear':list_gear}

        self.dico_gear_assembly_optimizer=dico


    def Optimize(self):
        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        self.analyse=M1
        M1.Optimize()
        self.solutions=M1.solutions
        
    def SearchCenterLine(self,nb_sol=1):
        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        self.analyse=M1
        M1.SearchCenterLine(nb_sol)
        self.solutions_search=M1.solutions_search

