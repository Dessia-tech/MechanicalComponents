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
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.center_distance=center_distance
        self.gear_set=gear_set
        self.transverse_pressure_angle=transverse_pressure_angle
        self.gear_graph=gear_graph

        self.DF,DB,self.node_dfs=self.GearAssemblyParam1(Z)
        #instantiation des objets Gears
        self.gears={}
        for ne,ns in enumerate(self.gear_set):
            self.gears[ne]={}
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[str(ng)]
                tpa=transverse_pressure_angle_rack[ng]
                cga=coeff_gear_addendum[ng]
                cgd=coeff_gear_dedendum[ng]
                crr=coeff_root_radius[ng]
                cct=coeff_circular_tooth_thickness[ng]
                self.gears[ne][ng]=Gear(z,db,cp,tpa,cga,cgd,crr,cct)
                
#        for g,z in Z.items():
#            db=DB[g]
#            cp=coefficient_profile_shift[str(g)]
#            self.gears[g]=Gear(z,db,cp)
        self.linear_backlash,self.radial_contact_ratio=self.GearAssemblyParam2(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
            
    def Update(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.DF,DB,self.node_dfs=self.GearAssemblyParam1(Z)
        self.linear_backlash,self.radial_contact_ratio=self.GearAssemblyParam2(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    def GearAssemblyParam2(self,Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        for ne,ns in enumerate(self.gear_set):
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[str(ng)]
                tpa=transverse_pressure_angle_rack[ng]
                cga=coeff_gear_addendum[ng]
                cgd=coeff_gear_dedendum[ng]
                crr=coeff_root_radius[ng]
                cct=coeff_circular_tooth_thickness[ng]
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
    
    def FreeCADExport(self,name,position,file_path,export_types):

        RIM1=self.Gear1.RimContour()
        RIM2=self.Gear2.RimContour()
#        RIM1.MPLPlot()
#        RIM2.MPLPlot()
        gear1=primitives3D.RevolvedProfile(vm.Point3D((position1[0],position1[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[RIM1],vm.Vector3D((position1[0],position1[1],0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim1')
        gear2=primitives3D.RevolvedProfile(vm.Point3D((position2[0],position2[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[RIM2],vm.Vector3D((position2[0],position2[1],0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim2')
        
        # Teeth
        TG1=self.Gear1.GearContours(5)
        TG2=self.Gear2.GearContours(5)
        list_rot=self.InitialPosition()
        L1=self.GearAssemblyTrace([TG1,TG2],[position1,position2],list_rot)
        C1=vm.Contour2D(L1[0])
        C2=vm.Contour2D(L1[1])
#        C1.MPLPlot()
        h1=0.5*self.Gear1.alpha_rim*(self.Gear1.outside_diameter-self.Gear1.root_diameter)
        h2=0.5*self.Gear2.alpha_rim*(self.Gear2.outside_diameter-self.Gear2.root_diameter)
        r1=self.Gear1.root_diameter/2-h1/2
        r2=self.Gear2.root_diameter/2-h2/2
#        print(h1,h2,r1,r2)
        c1=primitives3D.Cylinder((position1[0],position1[1],0.48*self.gear_width),
                                 (0,0,1),r1,1.05*self.gear_width)
        c2=primitives3D.Cylinder((position2[0],position2[1],0.48*self.gear_width),
                                 (0,0,1),r2,1.05*self.gear_width)
        
#        C1int=vm.Contour2D([vm.Circle2D(vm.Point2D(position1),r1)])
#        C2int=vm.Contour2D([vm.Circle2D(vm.Point2D(position2),r2)])
        if self.helix_angle==0.:            
            t1=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                            vm.Vector3D((0,1,0)),[C1],(0,0,self.gear_width))
            t2=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                            vm.Vector3D((0,1,0)),[C2],(0,0,self.gear_width))
        else:
            t1=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                                   vm.Vector3D((0,1,0)),(position1[0],position1[1],0),
                                                   (0,0,self.gear_width),self.DF1*mt.pi/mt.tan(self.helix_angle),
                                                   C1)
            t2=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                                   vm.Vector3D((0,1,0)),(position2[0],position2[1],0),
                                                   (0,0,self.gear_width),-self.DF2*mt.pi/mt.tan(self.helix_angle),
                                                   C2)

        # Creating holes in teeths
        t1=primitives3D.Cut(t1,c1,name='Teeth1')
        t2=primitives3D.Cut(t2,c2,name='Teeth2')

        model=vm.VolumeModel([gear1,t1,gear2,t2])
#        model=vm.VolumeModel([gear1,t1,gear2])
        model.FreeCADExport('python',file_path,'/usr/lib/freecad/lib',export_types)
    
    def SVGGearSet(self,name,position):
        #tuple1 et 2 correspondent a la position des centres
#        TG1=self.gears[0].GearContours(5)
#        G1=vm.Contour2D(TG1)
#        G1.MPLPlot()
            
        #optimisation pour le placement des axes des engrenages
        def fun(x):
            obj=0
            for num,it in enumerate(self.gear_set):
                eng1=it[0]
                eng2=it[1]
                obj+=(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num])**2
            return obj
        def eg(x):
            ine=[]
            for key,val in position.items():
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine=[]
            for num,it in enumerate(self.gear_set):
                eng1=it[0]
                eng2=it[1]
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
        for num,ens in enumerate(self.node_dfs):
            position1=(x_opt[2*ens[0]],x_opt[2*ens[0]+1])
            position2=(x_opt[2*ens[1]],x_opt[2*ens[1]+1])
            #tuple1 et 2 correspondent a la position des centres
            ne=self.gear_set.index(ens)
            Rot[ne]={}
            if num==0:
                TG[ens[0]]=self.gears[ne][ens[0]].GearContours(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][ens[0]]/2))
            TG[ens[1]]=self.gears[ne][ens[1]].GearContours(5)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[ne][ens[1]]/2))
            #Definition de la position angulaire initiale
            list_rot=self.InitialPosition(ne,ens)
            if position2[0]==position1[0]:
                if position2[1]-position1[1]>0:
                    angle=npy.pi/2
                else:
                    angle=-npy.pi/2
            else:
                angle=-npy.arctan((position2[1]-position1[1])/(position2[0]-position1[0]))
            if num==0:
                Rot[ne][ens[0]]=list_rot[0]-angle
                Rot[ne][ens[1]]=list_rot[1]-angle
            else:
                for k1,v1 in Rot.items():
                    if ens[0] in v1.keys():
                        Rot[ne][ens[0]]=v1[ens[0]]
                        delta_rot=Rot[ne][ens[0]]-(list_rot[0]-angle)
                Rot[ne][ens[1]]=list_rot[1]-angle-delta_rot*((self.gears[ne][ens[0]].Z)/(self.gears[ne][ens[1]].Z))
            sol=self.GearAssemblyTrace([TG[ens[0]],TG[ens[1]]],[position1,position2],list_rot=[Rot[ne][ens[0]],Rot[ne][ens[1]]])
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
    def __init__(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,cond_init,rack_list,rack_choice):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.cond_init=cond_init
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        Bounds=list(center_distance)
        Bounds.extend(transverse_pressure_angle)
        tp=len(coefficient_profile_shift.keys())
        for i in range(tp):
            Bounds.append(coefficient_profile_shift[str(i)])
        self.solutions=[]
        
        self.xi={'Z':Z,'gear_set':gear_set,'gear_graph':gear_graph}
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
        for k in sorted(list(self.rack_choice.keys())):
            v=self.rack_choice[k]
            for k2,v2 in self.rack_list[v].items():
                if k2 not in ['type','name','module']:
                    xj[k2].append(x_list[v][k2])
        return xj,dict_xu
    
    def _convert_xj2Xu(self,xj):
        sol=xj['center_distance'].copy()
        sol.extend(xj['transverse_pressure_angle'])
        tp=len(xj['coefficient_profile_shift'].keys())
        for key in range(tp):
            sol.append(xj['coefficient_profile_shift'][str(key)])
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
        for i in range(tp3):
            xj['coefficient_profile_shift'][str(i)]=X[tp1+tp2+i]
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
        for ng in range(self.GearAssembly.gear_graph.number_of_nodes()):
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
    def __init__(self,gear_width,center_distance,coefficient_profile_shift,gear_speed,Z,gear_set,transverse_pressure_angle,helix_angle,frequency,rack_list,rack_choice):
        self.Z=Z['data']
        self.gear_set=gear_set['data']
        self.gear_speed=gear_speed['data']
        self.frequency=frequency['data']
        self.center_distance=center_distance['data']
        self.transverse_pressure_angle=transverse_pressure_angle['data']
        self.coefficient_profile_shift=coefficient_profile_shift['data']
        self.rack_list=rack_list['data']
        self.rack_choice=rack_choice['data']

        self.nb_gear=npy.array(self.gear_set).max()+1
        self.nb_gear=len(self.gear_speed.keys())
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(range(self.nb_gear))
        gear_graph.add_edges_from(self.gear_set)
        self.gear_graph=gear_graph
        
        self.nb_rack=len(self.rack_list.keys())

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.node_dfs=list(nx.dfs_edges(gear_graph,self.node_init))
        
        self.AnalyzeCombination()
        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]

    def AnalyzeCombination(self):
        #recherche des Zmin et Zmax
        Zmin=npy.inf
        Zmax=0
        for n1,n2 in self.gear_set:
            if min(self.Z[str(n1)][0],self.Z[str(n2)][0])<Zmin:
                Zmin=min(self.Z[str(n1)][0],self.Z[str(n2)][0])
            if max(self.Z[str(n1)][1],self.Z[str(n2)][1])>Zmax:
                Zmax=max(self.Z[str(n1)][1],self.Z[str(n2)][1])
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
                if (z<self.Z[str(liste_node[0])][0]) or (z>self.Z[str(liste_node[0])][1]):
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
                    if (z2<self.Z[str(n2)][0]) or (z2>self.Z[str(n2)][1]):
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
            try:
                A1.Optimize()
            except:
                print('Problème de convergence')
            if len(A1.solutions)>0:
                xsol=A1.solutions[-1]
#                print(11,xsol)
#                print(22,A1.xj)
                xt=dict(list(A1.xi.items())+list(xsol.items()))
                self.solutions.append(GearAssembly(**xt))
#                break
                
    def SearchCenterLine(self,callback=lambda x:x):
        for ipl,pl in enumerate(self.plex_calcul):
            pl['DF']={}
            pl['Z_data']={}
            for k,v in pl['Z'].items():
                nr=pl['rack_choice'][k]
                mod=pl['rack_list'][nr]['module']
                module=(mod[0]+mod[1])/2
                pl['DF'][k]=module*v
                pl['Z_data'][str(k)]=[v,v]
            pl['center_distance']=[]
            for ne in self.gear_set:
                t1=(pl['DF'][ne[0]]+pl['DF'][ne[1]])/2
                pl['center_distance'].append([1e3*t1,1e3*t1])
            del pl['DF']
            self.solutions_search.append(pl)
                
            


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
        dico['rack_list']['data']={k:{'type':dico['rack_list']['data'][k]['type'],'name':dico['rack_list']['data'][k]['name'],'module':[dico['rack_list']['data'][k]['module'][0]*(1e-3),dico['rack_list']['data'][k]['module'][1]*(1e-3)],'transverse_pressure_angle_rack':[dico['rack_list']['data'][k]['transverse_pressure_angle_rack'][0]/180*npy.pi,dico['rack_list']['data'][k]['transverse_pressure_angle_rack'][1]/180*npy.pi],'coeff_gear_addendum':dico['rack_list']['data'][k]['coeff_gear_addendum'],'coeff_gear_dedendum':dico['rack_list']['data'][k]['coeff_gear_dedendum'],'coeff_root_radius':dico['rack_list']['data'][k]['coeff_root_radius'],'coeff_circular_tooth_thickness':dico['rack_list']['data'][k]['coeff_circular_tooth_thickness']} for k in dico['rack_list']['data'].keys()}

        dico=self.DefaultDataSet(dico)
        self.dico_gear_assembly_optimizer=dico


    def Optimize(self):
        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        self.analyse=M1
        M1.Optimize()
        self.solutions=M1.solutions
        
    def SearchCenterLine(self):
        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        self.analyse=M1
        M1.SearchCenterLine()
        self.solutions_search=M1.solutions_search

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
