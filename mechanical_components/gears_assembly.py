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

#from dessia_common import ResultsDBClient
#import pyDOE

#data_coeff_YB_Iso
evol_coeff_yb_iso={'data':[[0.0,1.0029325508401201],[4.701492563229561,0.9310850480431024],[9.104477651442416,0.8782991233021732],[14.5522388104227,0.8240469255458759],[19.328358165913905,0.784457481990179],[23.955224059269884,0.7609970656504502],[28.507462609617956,0.7521994155347822],[34.029850499665855,0.7507331425194141],[40.0,0.7492668574805859]],'x':'Linear','y':'Linear'}
#data_wholer_curve
wholer_hardened_alloy_steel={'data':[[4.296196199237153,1.9797762011105589],[4.824840106199563,1.9413306094362142],[5.3344338175705674,1.908892154601565],[6.115493253679078,1.8632380197445122],[6.596511629990596,1.8560294765618042],[7.144205815889171,1.8536266428508523],[7.691899918442984,1.8524252154829133],[8.010991340520903,1.8524252154829133]],'x':'Log','y':'Log'}
wholer_nitrided_alloy_steel={'data':[[4.104865629699472,1.9252942042661974],[4.568697315952783,1.8521640228225367],[4.887581626297173,1.8046294185503593],[5.438381821440599,1.7900033864666123],[6.402282079596832,1.7918316299646175],[7.264719174616821,1.7918316299646175],[7.989456220850952,1.793659894487549]],'x':'Log','y':'Log'}
wholer_through_hardened_steel={'data':[[4.172369719531124,1.895676495604088],[4.677200861168087,1.7983611100752137],[4.9677168648417585,1.741894170956562],[5.329671247836526,1.6842258044699714],[5.439210101685194,1.672211551815507],[6.091680488353632,1.6734129791834462],[7.139443246155129,1.671010124447568],[8.00146620105282,1.6758158339193243]],'x':'Log','y':'Log'}
wholer_surface_hardened_steel={'data':[[4.281908490035029,1.7611169667937343],[4.701013626493532,1.6998443182033265],[5.015342395492649,1.6553916107142128],[5.358246582896013,1.6109389032250994],[5.620187251510196,1.5857089915731581],[6.020242109032534,1.5748961873115592],[6.567936294931109,1.5748961873115592],[7.263269725861159,1.5772990210225108],[7.996703631318779,1.5772990210225108]],'x':'Log','y':'Log'}
wholer_carbon_steel={'data':[[4.307791955971963,1.6419147590563592],[5.242702822291173,1.535876005424268],[5.938450393343521,1.4700588400224806],[6.518240063668731,1.431665495290182],[7.221234961844144,1.4334937598131132],[7.989456220850952,1.4353220033111185]],'x':'Log','y':'Log'}
wholer_cast_iron={'data':[[4.307791955971963,1.6419147590563592],[5.242702822291173,1.535876005424268],[5.938450393343521,1.4700588400224806],[6.518240063668731,1.431665495290182],[7.221234961844144,1.4334937598131132],[7.989456220850952,1.4353220033111185]],'x':'Log','y':'Log'}
wholer_bronze={'data':[[4.307791955971963,1.6419147590563592],[5.242702822291173,1.535876005424268],[5.938450393343521,1.4700588400224806],[6.518240063668731,1.431665495290182],[7.221234961844144,1.4334937598131132],[7.989456220850952,1.4353220033111185]],'x':'Log','y':'Log'}
wholer_grey_iron={'data':[[4.307791955971963,1.6419147590563592],[5.242702822291173,1.535876005424268],[5.938450393343521,1.4700588400224806],[6.518240063668731,1.431665495290182],[7.221234961844144,1.4334937598131132],[7.989456220850952,1.4353220033111185]],'x':'Log','y':'Log'}
#data_gear_material
sigma_hardened_alloy_steel={'data':[[1.8422104370714443,1.4645831828946267],[1.948612010770208,1.5219116983411152],[2.0605171321606295,1.5810895335609718],[2.141235568740199,1.6254729099758645]],'x':'Log','y':'Log'}
sigma_nitrided_alloy_steel={'data':[[1.8458794622934307,1.4349942652846983],[1.943108482795906,1.488624180937243],[2.0201578941534892,1.5274596179084272],[2.128393990321924,1.5866374531282839]],'x':'Log','y':'Log'}
sigma_through_hardened_steel={'data':[[1.7798371068844516,1.292597616678765],[1.921094370898698,1.3850629693024938],[2.032999472571764,1.4627338829976548],[2.1650841833897223,1.5533499158480155]],'x':'Log','y':'Log'}
sigma_surface_hardened_steel={'data':[[1.8312033811228403,1.115064130895591],[1.932101426847302,1.200132264055036],[2.038503000546066,1.2852003773380847]],'x':'Log','y':'Log'}
sigma_carbon_steel={'data':[[1.677104538690319,1.1002696720906269],[1.7633265032441903,1.1723926463420797],[1.8385414118494579,1.2389677010262203],[1.8844041581135444,1.2796524577707729]],'x':'Log','y':'Log'}
sigma_cast_iron={'data':[[1.4734739247717241,0.922736186307453],[1.5468543306246763,0.9837633214242817],[1.6073931580593532,1.0336946174064863],[1.6404143456225206,1.0688314545837265]],'x':'Log','y':'Log'}
sigma_bronze={'data':[[1.313871566195314,0.7858874572688317],[1.3890864826875238,0.8487638922826322],[1.4294457009773085,0.8802021097895326],[1.4551288380965028,0.9097910273994609]],'x':'Log','y':'Log'}
sigma_grey_iron={'data':[[1.354230792372041,0.7100658633470387],[1.4276111785076375,0.7766409180311793],[1.4936535339166166,0.84691459238566],[1.5431853054026896,0.8986951882648367],[1.5725374677438706,0.933832025442077]],'x':'Log','y':'Log'}

class Material:
    def __init__(self,data_coeff_YB_Iso,data_wholer_curve,data_gear_material):
        self.data_coeff_YB_Iso=data_coeff_YB_Iso
        self.data_wholer_curve=data_wholer_curve
        self.data_gear_material=data_gear_material
        
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
hardened_alloy_steel=Material(evol_coeff_yb_iso,wholer_hardened_alloy_steel,sigma_hardened_alloy_steel)
nitrided_alloy_steel=Material(evol_coeff_yb_iso,wholer_nitrided_alloy_steel,sigma_nitrided_alloy_steel)
through_hardened_steel=Material(evol_coeff_yb_iso,wholer_through_hardened_steel,sigma_through_hardened_steel)
surface_hardened_steel=Material(evol_coeff_yb_iso,wholer_surface_hardened_steel,sigma_surface_hardened_steel)
carbon_steel=Material(evol_coeff_yb_iso,wholer_carbon_steel,sigma_carbon_steel)
cast_iron=Material(evol_coeff_yb_iso,wholer_cast_iron,sigma_cast_iron)
bronze=Material(evol_coeff_yb_iso,wholer_bronze,sigma_bronze)
grey_iron=Material(evol_coeff_yb_iso,wholer_grey_iron,sigma_grey_iron)

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
    def __init__(self, z, db, cp, transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness):
        
        self.rack=Rack(transverse_pressure_angle_rack)
        self.GearParam(z,db,cp,transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum, 
                       coeff_root_radius, coeff_circular_tooth_thickness)
        
    def Update(self,z,db,cp,transverse_pressure_angle_rack,
               coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,
               coeff_circular_tooth_thickness):
        self.GearParam(z,db,cp, transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum,
                       coeff_root_radius, coeff_circular_tooth_thickness)
        
    def GearParam(self,z, db, cp, transverse_pressure_angle_rack,
                  coeff_gear_addendum, coeff_gear_dedendum,
                  coeff_root_radius, coeff_circular_tooth_thickness):
        
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
        self.root_gear_angle=self.circular_tooth_thickness/(self.DFF/2)+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
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
    
    def Contour(self,discret=10,list_number=[None]):
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
    
    ### Stress
    
    def _CoeffYSIso(self,s_thickness_iso):
        #facteur de concentration de contrainte en pied de dent
        rho_f=self.rack.root_radius+self.rack.b**2/(self.DFF/2+self.rack.b)
        coeff_ys_iso=1+0.15*s_thickness_iso/rho_f
        return coeff_ys_iso
    
    def GearSectionISO(self,angle):
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
        theta0=fsolve((lambda theta:a + b*npy.tan(theta) + r*(-angle - self.root_angle/2 - theta + npy.pi/2)) ,0)[0]
        phi0=(a-b*npy.tan(theta0))/r
        pt_iso=self._Trochoide(phi0)
        angle0=npy.arctan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        angle_iso=self.root_gear_angle-2*angle0
        diameter_iso=2*norm(pt_iso)
        s_thickness_iso=diameter_iso*npy.sin(angle_iso/2)
        h_height_iso=(s_thickness_iso/2)/npy.tan(angle)
        return s_thickness_iso,h_height_iso
        
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
    def __init__(self,Z, center_distance, gear_set,transverse_pressure_angle,
                 coefficient_profile_shift,gear_graph, transverse_pressure_angle_rack,
                 coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,
                 coeff_circular_tooth_thickness,list_gear,material,torque,cycle):
        
        self.center_distance=center_distance
        self.gear_set=gear_set
        self.transverse_pressure_angle=transverse_pressure_angle
        self.gear_graph=gear_graph
        self.list_gear=list_gear
        self.material=material
        self.helix_angle=0

        self.DF,DB,self.gear_set_dfs=self.GearGeometryParameter(Z)
        self.cycle=self.CycleParameter(cycle,Z)
        self.torque1,self.normal_load,self.tangential_load,self.radial_load=self.GearTorque(Z,torque,DB)
        
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
                
        self.linear_backlash,self.radial_contact_ratio=self.GearContactRatioParameter(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        self.gear_width,self.sigma_iso,self.sigma_lim=self.GearWidthDefinition()
            
    def Update(self,Z,center_distance,gear_set,transverse_pressure_angle,coefficient_profile_shift,gear_graph,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness,list_gear,material,torque,cycle):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.DF,DB,self.gear_set_dfs=self.GearGeometryParameter(Z)
        self.linear_backlash,self.radial_contact_ratio=self.GearContactRatioParameter(Z,coefficient_profile_shift,DB,transverse_pressure_angle_rack,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    def GearContactRatioParameter(self,Z,coefficient_profile_shift,DB,
                           transverse_pressure_angle_rack,coeff_gear_addendum,
                           coeff_gear_dedendum,coeff_root_radius,
                           coeff_circular_tooth_thickness):
        
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
        for g1,g2 in self.gear_set_dfs:
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

    def GearGeometryParameter(self,Z):
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
        
        gear_set_dfs=list(nx.edge_dfs(self.gear_graph, [self.gear_set[0][0],self.gear_set[0][1]]))

        return DF,DB,gear_set_dfs
    
    def GearTorque(self,Z,torque,DB):
        ggd=self.gear_graph.degree(self.list_gear)
        for ne,nb_connexion in ggd:
            if (ne in torque.keys()) and (nb_connexion==1):
                ne_init=ne
        gs_init=list(self.gear_graph.edges(ne_init))[0]
        gs_torque_dfs=list(nx.edge_dfs(self.gear_graph, gs_init))
        torque1={}
        normal_load=[]
        tangential_load=[]
        radial_load=[]
        torque_input_m=torque[ne_init]
        for (eng1,eng2) in gs_torque_dfs:
            torque2=-torque_input_m*Z[eng2]/Z[eng1]
            if eng2 in torque.keys():
                torque2+=torque[eng2]
            torque1[eng1]=torque_input_m
            torque_input_m=torque2
        for i,(eng1,eng2) in enumerate(gs_torque_dfs):
            try:
                ne=self.gear_set.index((eng1,eng2))
            except:
                ne=self.gear_set.index((eng2,eng1))
            tq1=torque1[eng1]
            normal_load.append(abs(tq1)*2/(DB[ne][eng1]))
            tangential_load.append(abs(tq1)*2/(self.DF[ne][eng1]))
            radial_load.append(npy.tan(self.transverse_pressure_angle[ne])*tangential_load[-1])
        return torque1,normal_load,tangential_load,radial_load
    
    def CycleParameter(self,cycle,Z):
        eng_init=list(cycle.keys())[0]
        for eng in self.list_gear:
            if eng not in cycle.keys():
                cycle[eng]=cycle[eng_init]*Z[eng_init]/Z[eng]
        return cycle
                    
    ### Stress
    
    def SigmaLewis(self):
        
        self.sigma_lewis_maximum1=6*self.tangential_load*self.gear_height_lewis1/(self.gear_width*self.Gear1.root_gear_length**2)
        self.sigma_lewis_maximum2=6*self.tangential_load*self.gear_height_lewis2/(self.gear_width*self.Gear2.root_gear_length**2)
        
    def GearWidthDefinition(self):
        coeff_yf_iso=self._CoeffYFIso()
        coeff_ye_iso=self._CoeffYEIso()
        coeff_yb_iso=self._CoeffYBIso()
        sigma_lim=self.SigmaMaterialISO()
        gear_width={}
        for eng in self.list_gear:
            gear_width[eng]=0
        for ne,(eng1,eng2) in enumerate(self.gear_set):
            gear_width1=abs(self.tangential_load[ne]/(sigma_lim[ne][eng1]*self.gears[ne][eng1].rack.module)*coeff_yf_iso[ne][eng1]*coeff_ye_iso[ne]*coeff_yb_iso[ne][eng1])
            gear_width2=abs(self.tangential_load[ne]/(sigma_lim[ne][eng2]*self.gears[ne][eng2].rack.module)*coeff_yf_iso[ne][eng2]*coeff_ye_iso[ne]*coeff_yb_iso[ne][eng2])
            gear_width_set=max(gear_width1,gear_width2)
            gear_width[eng1]=max(gear_width[eng1],gear_width_set)
            gear_width[eng2]=max(gear_width[eng2],gear_width_set)
        sigma_iso=sigma_lim
        return gear_width,sigma_iso,sigma_lim
        
    def SigmaMaterialISO(self):
        safety_factor=4
        angle=30/180*npy.pi
        sigma_lim={}
        for ne,(eng1,eng2) in enumerate(self.gear_set):
            sigma_lim[ne]={}
            
            matrice_wholer=self.material[eng1].data_wholer_curve
            matrice_material=self.material[eng1].data_gear_material
            sgla=self.material[eng1].FunCoeff(self.cycle[eng1],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl1=self.material[eng1].FunCoeff(sgla,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_1,h_height_iso_1=self.gears[ne][eng1].GearSectionISO(angle)
            coeff_ys_iso=self.gears[ne][eng1]._CoeffYSIso(s_thickness_iso_1)
            sigma_lim[ne][eng1]=(sgl1/(safety_factor*coeff_ys_iso))*10**7
            
            matrice_wholer=self.material[eng2].data_wholer_curve
            matrice_material=self.material[eng2].data_gear_material
            sglb=self.material[eng2].FunCoeff(self.cycle[eng2],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl2=self.material[eng2].FunCoeff(sglb,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_2,h_height_iso_2=self.gears[ne][eng2].GearSectionISO(angle)
            coeff_ys_iso=self.gears[ne][eng2]._CoeffYSIso(s_thickness_iso_2)
            sigma_lim[ne][eng2]=(sgl2/(safety_factor*coeff_ys_iso))*10**7
        return sigma_lim

        
    def _CoeffYFIso(self):
        #facteur de forme pour la contrainte ISO
        angle=30/180*npy.pi
        coeff_yf_iso={}
        for ne,(eng1,eng2) in enumerate(self.gear_set):
            coeff_yf_iso[ne]={}
            s_thickness_iso_1,h_height_iso_1=self.gears[ne][eng1].GearSectionISO(angle)
            s_thickness_iso_2,h_height_iso_2=self.gears[ne][eng2].GearSectionISO(angle)
            coeff_yf_iso[ne][eng1]=(6*(h_height_iso_1/self.gears[ne][eng1].rack.module)*npy.cos(self.transverse_pressure_angle[ne]))/((s_thickness_iso_1/self.gears[ne][eng1].rack.module)**2*npy.cos(self.gears[ne][eng1].rack.transverse_pressure_angle))
            coeff_yf_iso[ne][eng2]=(6*(h_height_iso_2/self.gears[ne][eng2].rack.module)*npy.cos(self.transverse_pressure_angle[ne]))/((s_thickness_iso_2/self.gears[ne][eng2].rack.module)**2*npy.cos(self.gears[ne][eng2].rack.transverse_pressure_angle))
        return coeff_yf_iso
        
    def _CoeffYEIso(self):
        #facteur de conduite pour la contrainte ISO
        coeff_ye_iso=[]
        for ne,eng in enumerate(self.gear_set):
            coeff_ye_iso.append(1/self.radial_contact_ratio[ne])
        return coeff_ye_iso
        
    def _CoeffYBIso(self):
        #facteur de contrefort pour la contrainte ISO
        coeff_yb_iso={}
        for ne,(eng1,eng2) in enumerate(self.gear_set):
            coeff_yb_iso[ne]={}
            matrice_YB=self.material[eng1].data_coeff_YB_Iso
            coeff_yb_iso[ne][eng1]=self.material[eng1].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
            matrice_YB=self.material[eng2].data_coeff_YB_Iso
            coeff_yb_iso[ne][eng2]=self.material[eng2].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
        return coeff_yb_iso
        
    ### Fonction de trace et export
    
    def GearRotate(self,list_gear,list_center,list_rot):
        export=[]
        for (i,center,k) in zip(list_gear,list_center,list_rot):
            model_export=[]
            for m in i:
                model_trans=m.Translation(center)
                model_trans_rot=model_trans.Rotation(vm.Point2D(center),k)
                model_export.append(model_trans_rot)
            export.append(model_export)
        return export
    
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
    
    def VolumeModel(self, centers = [], axis = (1,0,0), name = ''):
        
        x = vm.Vector3D(axis)
        y = x.RandomUnitNormalVector()
        z = vm.Vector3D(npy.cross(x.vector, y.vector))  
        
        if len(centers)==0:
            centers = []
            pos_axis = self.PosAxis({self.list_gear[0]:[0,0]})
            for i in range(int(len(pos_axis)/2)):
                centers.append(tuple(pos_axis[2*i]*y.vector+pos_axis[2*i+1]*z.vector))
        else:
            center_var=[]
            for c in centers:
                center_var.append((0,npy.dot(c,y.vector),npy.dot(c,z.vector)))
            centers=center_var
            
        TG={}#
        Struct=[]
        Rot={}
        primitives=[]
        
#        angles=[0]
        
        for num,en in enumerate(self.gear_set_dfs):
            
            ens=[self.list_gear.index(en[0]),self.list_gear.index(en[1])]
            position1 = centers[ens[0]]
            position2 = centers[ens[1]]
            
            #tuple1 et 2 correspondent a la position des centres
            ne=self.gear_set.index(en)
            Rot[ne]={}
            if num==0:
                TG[en[0]]=self.gears[ne][en[0]].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2))
            TG[en[1]]=self.gears[ne][en[1]].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[ne][en[1]]/2))
            #Definition de la position angulaire initiale
            list_rot=self.InitialPosition(ne,en)
            
            if position2[1]==position1[1]:
                if position2[2]-position1[2]>0:
                    angle=npy.pi/2
                else:
                    angle=-npy.pi/2
            else:
                angle=-npy.arctan((position2[2]-position1[2])/(position2[1]-position1[1]))
            if num==0:
                Rot[ne][en[0]]=list_rot[0]-angle
                Rot[ne][en[1]]=list_rot[1]-angle
            else:
                for k1,v1 in Rot.items():
                    if en[0] in v1.keys():
                        Rot[ne][en[0]]=v1[en[0]]
                        delta_rot=Rot[ne][en[0]]-(list_rot[0]-angle)
                Rot[ne][en[1]]=list_rot[1]-angle-delta_rot*((self.gears[ne][en[0]].Z)/(self.gears[ne][en[1]].Z))
            sol=self.GearRotate([TG[en[0]],TG[en[1]]],[(position1[1::]),(position2[1::])],
                                       list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
        
            C1=vm.Contour2D(sol[0])
            C2=vm.Contour2D(sol[1])
            
            extrusion_vector1 = (self.gear_width[en[0]]*x).vector
            extrusion_vector2 = (self.gear_width[en[1]]*x).vector
            
            if num==0:
                t1=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),y,z,[C1],extrusion_vector1)
                primitives.append(t1)
        
            t2=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),y,z,[C2],extrusion_vector2)
            primitives.append(t2)
#            print(primitives)
        model=vm.VolumeModel(primitives)
        return model

    def FreeCADExport(self, file_path, export_types=['fcstd'], python_path = 'python',
                      freecad_path = '/usr/lib/freecad/lib'):
        
        model = self.VolumeModel()
        model.FreeCADExport(python_path ,file_path, freecad_path, export_types)
        
    def PosAxis(self,position):
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
        return x_opt
    
    def SVGExport(self,name,position):
        x_opt=self.PosAxis(position)
        TG={}
        L1=[]
        Struct=[]
        Rot={}
        for num,en in enumerate(self.gear_set_dfs):
            ens=[self.list_gear.index(en[0]),self.list_gear.index(en[1])]
            position1=(x_opt[2*ens[0]],x_opt[2*ens[0]+1])   
            position2=(x_opt[2*ens[1]],x_opt[2*ens[1]+1])
            #tuple1 et 2 correspondent a la position des centres
            ne=self.gear_set.index(en)
            Rot[ne]={}
            if num==0:
                TG[en[0]]=self.gears[ne][en[0]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2))
            TG[en[1]]=self.gears[ne][en[1]].Contour(5)
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
            sol=self.GearRotate([TG[en[0]],TG[en[1]]],[position1,position2],list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
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
                 gear_graph,cond_init,rack_list,rack_choice,list_gear,material,torque,cycle):
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
        
        self.xi={'Z':Z,'gear_set':gear_set,'gear_graph':gear_graph,'list_gear':list_gear,'material':material,'torque':torque,'cycle':cycle}
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
        #contrainte sur le module
        for ne,gs in enumerate(self.GearAssembly.gear_set):
            for g in gs:
                mo=self.GearAssembly.gears[ne][g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]<lm[1]:
                ineq.append(mo-list_module[0])
                ineq.append(list_module[1]-mo)
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
#        #Maximisation du module pour avoir des pignons avec un faible gear_width
#        for ne,gs in enumerate(self.GearAssembly.gear_set):
#            for g in gs:
#                mo=self.GearAssembly.gears[ne][g].rack.module
#                obj+=100*(1/mo)**2
                
        for lb in self.GearAssembly.linear_backlash:
            obj+=100*((1e-4)-lb)**2
            
        for i in fineq:
            if i < 0:
                obj+=-100*i
            else:
                obj+=0.0001*i
        for i in feq:
            if i<0:
                obj+=-100*i
            else:
                obj+=100*i
        return obj
    
    def Optimize(self):
        boucle=100
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
    def __init__(self,gear_set,gear_speed,center_distance,Z={},transverse_pressure_angle=None,
                 helix_angle=None,gear_width=None,frequency=[[0,0]],coefficient_profile_shift=None,
                 rack_list=None,rack_choice=None,material=None,torque=None,cycle=None):
        
        # Initialisation
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
            coefficient_profile_shift={list_gear[0]:[-1,1]}
        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne]=[-1,1]
                
        if rack_list==None:
            rack_list={0:{'name':'Optim_Module','module':[0.5*1e-3,3*1e-3],'transverse_pressure_angle_rack':[20*npy.pi/180,20*npy.pi/180],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
            
        if rack_choice==None:
            rack_choice={list_gear[0]:list(rack_list.keys())[0]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne]=list(rack_list.keys())[0]
                
        if material==None:
            material={list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne]=hardened_alloy_steel
        
        if torque==None:
            torque={list_gear[0]:100}
            
        if cycle==None:
            cycle={list_gear[0]:1e6}
        
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
        self.material=material
        self.torque=torque
        self.cycle=cycle

        self.nb_gear=len(list_gear)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.gear_set)
        self.gear_graph=gear_graph
        
        self.nb_rack=len(self.rack_list.keys())

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.gear_set_dfs=list(nx.dfs_edges(gear_graph,self.node_init))
        
        if self.Z=={}:
            self.Z=self.AnalyseZ()
        
        self.AnalyzeCombination()
        for i,plex in enumerate(self.plex_calcul):
            plex['gear_graph']=self.gear_graph
            plex['rack_list']=self.rack_list
            plex['list_gear']=self.list_gear
            plex['material']=self.material
            plex['torque']=self.torque
            plex['cycle']=self.cycle
            plex['gear_set']=self.gear_set
            plex['center_distance']=self.center_distance
            plex['transverse_pressure_angle']=self.transverse_pressure_angle
            plex['coefficient_profile_shift']=self.coefficient_profile_shift
            self.plex_calcul[i]=plex
            
            
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
        print('np: ', np)
        dt=tools.RegularDecisionTree(np)
        node=dt.current_node
        n1=self.node_init
        liste_node=[n1]

        for (n1,n2) in self.gear_set_dfs:
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
#                    print(dt.current_node,self.gear_set_dfs,dt.current_depth)
                    (n1,n2)=self.gear_set_dfs[dt.current_depth-1]
                    i1=liste_node.index(n1)
                    i2=liste_node.index(n2)
                    z1=liste_gear[dt.current_node[i1]]
                    z2=liste_gear[dt.current_node[i2]]
                #analyse ACV engrenage 2 à 2
                if (pgcd(z1,z2)!=1) & (dt.current_depth<=(self.nb_gear-1)):
                    valid=False
                #Analyse des bornes sur Z
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    if (z2<self.Z[n2][0]) or (z2>self.Z[n2][1]):
                        valid=False
                #analyse demul interne
                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
                    demul=liste_gear[dt.current_node[i1]]/liste_gear[dt.current_node[i2]]
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
#                        for nes,(n1,n2) in enumerate(self.gear_set_dfs):
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
#                Temp['cond_init']=res.x
                Temp['cond_init']=0
                Temp['rack_choice']=rack
                self.plex_calcul.append(Temp)
                incr+=1
#                if incr==3:
#                    break
            dt.NextNode(valid)
        if incr>1:
            print('Nombre de combinaison trouvées: {}'.format(incr))


    def Optimize(self,callback=lambda x:x):
        lpx=len(self.plex_calcul)
        for ii,i in enumerate(self.plex_calcul):
#            print('{}%'.format(ii/lpx*100))
            callback(ii/lpx)
            plex=i
            
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
                print('Largueur denture des dentures convergées: {}'.format(self.solutions[-1].gear_width))
                print('Nombre de dent des dentures convergées: {}'.format(plex['Z']))
                print('Entraxe des dentures convergées: {}'.format(self.solutions[-1].center_distance))
                break
                
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
            ga=GearAssemblyOptimizer(gear_set=self.gear_set,gear_speed=self.gear_speed,
                                            center_distance=cd_input,Z=Z_data,rack_list=self.rack_list,
                                            rack_choice=self.rack_choice,
                                            torque=self.torque,cycle=self.cycle,material=self.material)
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

