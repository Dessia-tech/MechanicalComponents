import numpy as npy
import volmdlr as vm
#import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
#import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
#from sympy import *
import itertools

import mechanical_components.LibSvg as LibSvg

import persistent
#from dessia_common import ResultsDBClient
#import pyDOE

class Rack(persistent.Persistent):
    def __init__(self,transverse_radial_pitch,transverse_pressure_angle_T=None,circular_tooth_thickness=None,gear_addendum=None,
                 gear_dedendum=None,root_radius_T=None,root_radius_R=None,transverse_pressure_angle_R=None):
        
        self.RackParam(transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                       root_radius_T,root_radius_R,transverse_pressure_angle_R)
        
    def RackParam(self,transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                  root_radius_T,root_radius_R,transverse_pressure_angle_R):
    
        self.transverse_pressure_angle_T=transverse_pressure_angle_T
        self.transverse_pressure_angle_R=transverse_pressure_angle_R
        self.transverse_radial_pitch=transverse_radial_pitch
        self.circular_tooth_thickness=circular_tooth_thickness
        self.gear_addendum=gear_addendum
        self.gear_dedendum=gear_dedendum
        self.root_radius_T=root_radius_T
        self.root_radius_R=root_radius_R
        self.transverse_pressure_angle_R=transverse_pressure_angle_R
        self.module=self.transverse_radial_pitch/npy.pi
        
        #if nothing define, we put the ISO definition for the rack
        if gear_addendum==None:
            self.gear_addendum=self.module
        if gear_dedendum==None:
            self.gear_dedendum=1.25*self.module
        if root_radius_T==None:
            self.root_radius_T=0.38*self.module
        if circular_tooth_thickness==None:
            self.circular_tooth_thickness=transverse_radial_pitch/2
        if transverse_pressure_angle_T==None:
            self.transverse_pressure_angle_T=20/180*npy.pi
        if transverse_pressure_angle_R==None:
            self.transverse_pressure_angle_R=self.transverse_pressure_angle_T
        if root_radius_R==None:
            self.root_radius_R=self.root_radius_T
            
        self.tooth_space=self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth=self.gear_addendum+self.gear_dedendum
        self.clearance_T=self.root_radius_T-self.root_radius_T*npy.sin(self.transverse_pressure_angle_T)
        self.clearance_R=self.root_radius_R-self.root_radius_R*npy.sin(self.transverse_pressure_angle_R)
        self.clearance=npy.max([self.clearance_T,self.clearance_R])
        
        #paramètre pour la trochoide
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T)-self.root_radius_T*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle_T)/(npy.sin(self.transverse_pressure_angle_T))))
        self.b=self.gear_dedendum-self.root_radius_T
        
    def Update(self,transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness=None,gear_addendum=None,
               gear_dedendum=None,root_radius_T=None,root_radius_R=None,transverse_pressure_angle_R=None):
        
        self.RackParam(transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                  root_radius_T,root_radius_R,transverse_pressure_angle_R)
        
    def RackContours(self,list_number):
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*npy.tan(self.transverse_pressure_angle_R),self.gear_addendum))
        p4=p1.Translation((self.circular_tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*npy.tan(self.transverse_pressure_angle_T),self.gear_addendum))
        p5=p4.Translation((self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T),-self.gear_dedendum))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_R),-self.gear_dedendum))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.root_radius_T,5:self.root_radius_R},False)
        Rack_Elem=[]
        for i in list_number:
            Rack_Elem.append(L.Translation(((i)*(p7.vector-p1.vector))))
        return Rack_Elem
    
    def RackComplete(self,list_number):
        
        L=self.RackContours(list_number)
        p1=L[0].points[0]
        p6=L[-1].points[-1]
        p2=p1.Translation((-self.circular_tooth_thickness,0))
        p3=p2.Translation((0,2*self.whole_depth))
        p5=p6.Translation((self.circular_tooth_thickness,0))
        p4=p5.Translation((0,2*self.whole_depth))
        Lbis=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6],{},False)
        L.append(Lbis)
        return L
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def CriteriaIneq(self):
        crit=[self.transverse_radial_pitch-self.circular_tooth_thickness-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T)
        -self.gear_dedendum*npy.tan(self.transverse_pressure_angle_R)-(self.root_radius_T*npy.cos(self.transverse_pressure_angle_T)
        -npy.tan(self.transverse_pressure_angle_T)*self.root_radius_T*(1-npy.sin(self.transverse_pressure_angle_T)))
        -(self.root_radius_R*npy.cos(self.transverse_pressure_angle_R)-npy.tan(self.transverse_pressure_angle_R)
        *self.root_radius_R*(1-npy.sin(self.transverse_pressure_angle_R)))]
        crit.append(self.circular_tooth_thickness-(self.gear_addendum*npy.tan(self.transverse_pressure_angle_T)+self.gear_addendum*npy.tan(self.transverse_pressure_angle_R)))
        return crit
        
    def Mass(self):
        pass
    
    def VolumeModel(self):
        pass
    
    def SVGExport(self,number,name):
        
        R1=self.RackComplete(npy.arange(number))
        Ref=[vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,0)),vm.Point2D(((number+1)*self.transverse_radial_pitch,0)))]
        Ref.append(vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,self.gear_addendum)),vm.Point2D(((number+1)*self.transverse_radial_pitch,self.gear_addendum))))
        Ref.append(vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,-self.gear_dedendum)),vm.Point2D(((number+1)*self.transverse_radial_pitch,-self.gear_dedendum))))
        SVG1=LibSvg.SVGTrace(700,0)
        SVG1.Convert(R1,'R1','black',0.02,0)
        SVG1.Convert(Ref,'Ref','black',0.01,1,'0.01px, 0.08px')
        SVG1.Export(name)
        
    def Dict(self):
        d=self.__dict__.copy()
        return d
    
    def CSVExport(self):
        d=self.__dict__.copy()
        return list(d.keys()),list(d.values())
    
class Gear(persistent.Persistent):
    def __init__(self,tooth_number,rack_data,coefficient_profile_shift=0):
        
        self.rack=rack_data
        
        save=[self.rack.transverse_radial_pitch,tooth_number,coefficient_profile_shift]
        self.save=save[:]
        self.GearParam(tooth_number,rack_data,coefficient_profile_shift)
        self.RootDiameterActive()
        
    def GearParam(self,tooth_number,rack_data,coefficient_profile_shift):
        
        self.tooth_number=tooth_number
        self.coefficient_profile_shift=coefficient_profile_shift
        self.pitch_diameter_factory=rack_data.transverse_radial_pitch*tooth_number/npy.pi
        
        self.outside_diameter=self.pitch_diameter_factory+2*(rack_data.gear_addendum+rack_data.module*self.coefficient_profile_shift)
        self.root_diameter=self.pitch_diameter_factory-2*(rack_data.gear_dedendum-rack_data.module*self.coefficient_profile_shift)
        
        #data sur le cercle primitif ()
        self.circular_tooth_thickness=rack_data.circular_tooth_thickness+rack_data.module*self.coefficient_profile_shift*npy.tan(rack_data.transverse_pressure_angle_T)+rack_data.module*self.coefficient_profile_shift*npy.tan(rack_data.transverse_pressure_angle_R)
        self.tooth_space=rack_data.transverse_radial_pitch-self.circular_tooth_thickness
        
        self.transverse_pressure_angle_T_factory=rack_data.transverse_pressure_angle_T
        circular_tooth_thickness_angle=self.circular_tooth_thickness/(self.pitch_diameter_factory/2)
        tooth_space_angle=self.tooth_space/(self.pitch_diameter_factory/2)
        self.base_diameter=self.pitch_diameter_factory*npy.cos(self.transverse_pressure_angle_T_factory)
        temp=self.pitch_diameter_factory-2*rack_data.gear_dedendum+2*rack_data.root_radius_T-2*rack_data.root_radius_T*npy.sin(rack_data.transverse_pressure_angle_T)-2*rack_data.module*self.coefficient_profile_shift
        
        #prise de la valeur absolue pour gérer les cas de changement de signe lié au déport
        theta2=abs(fsolve((lambda theta:self.base_diameter/2*npy.cos(rack_data.transverse_pressure_angle_T-theta)/npy.cos(theta)-temp/2) ,0)[0])
        
        #print(theta2)
        self.theta0=npy.tan(theta2)-rack_data.transverse_pressure_angle_T
        self.alpha_root_diameter_active=theta2
        self.tan_alpha_root_diameter_active=npy.tan(self.alpha_root_diameter_active)
        self.alpha_outside_diameter=npy.arccos(self.base_diameter/self.outside_diameter)
        self.alpha_pitch_diameter=npy.arccos(self.base_diameter/self.pitch_diameter_factory)
        self.inv_alpha_pitch_diameter=npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter
        self.base_circular_tooth_thickness=self.base_diameter/2*(2*self.circular_tooth_thickness/self.pitch_diameter_factory+2*self.inv_alpha_pitch_diameter)
        self.transverse_base_pitch=rack_data.transverse_radial_pitch*npy.cos(rack_data.transverse_pressure_angle_T)
        
        self.root_active_angle=self.tooth_space/(self.pitch_diameter_factory/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter-npy.tan(self.alpha_root_diameter_active)+self.alpha_root_diameter_active)
        self.root_angle=self.tooth_space/(self.pitch_diameter_factory/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle=self.circular_tooth_thickness/(self.pitch_diameter_factory/2)+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        
        self.outside_active_angle=2*self.circular_tooth_thickness/self.pitch_diameter_factory-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter-npy.tan(self.alpha_pitch_diameter)+self.alpha_pitch_diameter)
        
        
#        self.root_diameter_active=self.base_diameter/npy.cos(theta2)
#        self.alt_internal_contact=self.root_diameter_active/2*npy.sin(rack_data.transverse_pressure_angle_T-theta2)
#        self.RootDiameterActive()
        
        save=[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]
        if not self.save==[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]:
            self.RootDiameterActive()
#            print(save,self.root_diameter_active)
            
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        self.root_diameter_active_bis=2*(rho+b**2/(r+b))
        
        self.save=save[:]
        
        
    def Update(self,tooth_number,rack_data,coefficient_profile_shift):
        
        self.GearParam(tooth_number,rack_data,coefficient_profile_shift)
        save=[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]
        if not self.save==[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]:
            self.RootDiameterActive()
#            print(save,self.root_diameter_active)
        self.save=save[:]
        
    def CriteriaIneq(self):
        crit=[self.root_diameter_active-self.base_diameter-0.3*self.transverse_base_pitch]
        crit.append(self.outside_diameter-self.root_diameter-0.3*self.transverse_base_pitch)
        #Distance top of the rack
        #DistRI=self.DistRackInvolute()
        #crit.append(DistRI)
        #Top of the gear must be positive
        crit.append(self.base_circular_tooth_thickness*2/self.base_diameter-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter))
        return crit
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def GearSection(self,diameter):
        #epaisseur de la dent au diameter
        alpha_diameter=npy.arccos(self.base_diameter/diameter)
        theta1=(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(npy.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)
        

    def GearContours(self,discret=10,list_number=[None]):
        #Analytical involute profil
        self.RootDiameterActive()
        if list_number==[None]:
            list_number=npy.arange(int(self.tooth_number))
        L=[self.InvoluteTrace(discret,0,'T')]
        L.extend(self.TrochoideTrace(2*discret,0,'T'))
        L.append(self.InvoluteTrace(discret,0,'R'))
        L.extend(self.TrochoideTrace(2*discret,0,'R'))
        L.append(self.OutsideTrace(0))
        #print(L)
        for i in list_number:
            L.append(self.InvoluteTrace(discret,i,'T'))
            L.extend(self.TrochoideTrace(2*discret,i,'T'))
            L.append(self.InvoluteTrace(discret,i,'R'))
            L.extend(self.TrochoideTrace(2*discret,i,'R'))
            L.append(self.OutsideTrace(i))
        return L
        
    def InvoluteTrace(self,discret,number,ind='T'):
        
        if ind=='T':
            drap=1
        else:
            drap=-1
        theta=npy.linspace(self.tan_alpha_root_diameter_active,npy.tan(self.alpha_outside_diameter),discret)
#        theta=npy.linspace(0,npy.tan(self.alpha_outside_diameter),discret)
        sol=self.Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        else:
            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        return L
    
    def Trace(self,x,y):
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,len(x)):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        return ref
        
    def Involute(self,tan_alpha):
        
        sol=(self.base_diameter/2*npy.cos(tan_alpha)+self.base_diameter/2*tan_alpha*npy.sin(tan_alpha),
             self.base_diameter/2*npy.sin(tan_alpha)-self.base_diameter/2*tan_alpha*npy.cos(tan_alpha))
        return sol
    
    def Trochoide(self,phi,ind='T'):
        
        if ind=='T':
            drap=1
        else:
            drap=-1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        
        x2=rho*npy.sin(npy.arctan((a-r*phi)/b)-phi)+a*npy.cos(phi)-b*npy.sin(phi)+r*(npy.sin(phi)-phi*npy.cos(phi))
        y2=-rho*npy.cos(npy.arctan((a-r*phi)/b)-phi)-a*npy.sin(phi)-b*npy.cos(phi)+r*(npy.cos(phi)+phi*npy.sin(phi))
        sol=(y2,x2)
        return sol
    
    def RootDiameterActive(self):
        
        fun = (lambda t : (norm(npy.array(self.Involute(t[0]))-npy.dot(npy.array([[npy.cos(-self.root_angle/2),-npy.sin(-self.root_angle/2)],[npy.sin(-self.root_angle/2),npy.cos(-self.root_angle/2)]]),(self.Trochoide(t[1])))))**2)
        bnds = ((0, 1), (-1,1))
        sol=minimize(fun,[1,-1], bounds=bnds, method='SLSQP', tol=1e-10)
        
        self.alpha_root_diameter_active=npy.arctan(sol.x[0])
        self.root_diameter_active=norm(npy.array(self.Involute(sol.x[0])))*2
        self.tan_alpha_root_diameter_active=npy.tan(self.alpha_root_diameter_active)
        self.phi_trochoide=sol.x[1]
#        print(fun(sol.x),sol.x)
        
        #corde de la dent au diam de pied actif
        self.root_gear_length=npy.sin((self.root_gear_angle-2*(npy.tan(self.alpha_root_diameter_active)-self.alpha_root_diameter_active))/2)*(self.root_diameter_active/2)*2
        
    
    def RootTrace(self,discret,number,ind='T'):
        #trace simplifie de pied de dent
        if ind=='T':
            drap=1
        else:
            drap=-1
        theta=npy.linspace(self.tan_alpha_root_diameter_active,npy.tan(self.alpha_root_diameter_active+0.0001),2)
        sol=self.Involute(drap*theta)
        p1=vm.Point2D((sol[0][0],sol[1][0]))
        pt=vm.Point2D((sol[0][1],sol[1][1]))
        p2=p1+p1-pt
        theta3=self.root_active_angle/2-(npy.tan(self.alpha_root_diameter_active)-self.alpha_root_diameter_active)
        p4=vm.Point2D((self.root_diameter/2*npy.cos(theta3),-drap*self.root_diameter/2*npy.sin(theta3)))
        p3=p4.Rotation(vm.Point2D((0,0)),drap*0.0001)
        
        tck,ptx=self.BSpline([p1,p2,p3,p4])   
        a=[]
        for t in npy.linspace(0,1,discret):
            a.append(vm.Point2D((splev(t,tck)[0],splev(t,tck)[1])))
        
        ref=primitives2D.RoundedLines2D(a,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        else:
            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        #L1=[L]+[vm.CompositePrimitive2D([p1,p3,p4])]
        return L
    
    def TrochoideSecondary(self,list_number,ind,angle,number=10):
        
        ref=[]
        for t in npy.linspace(-angle,angle,number):
            ref.append(vm.Point2D((self.Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        L=[]
        for i in list_number:
            L.append(ref.Rotation(vm.Point2D((0,0)),-i*2*npy.pi/self.tooth_number))
        return L
    
    def TrochoideTrace(self,discret,number,ind='T'):
        if ind=='T':
            drap=1
        else:
            drap=-1
        
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        self.phi0=npy.arctan((a)/(self.pitch_diameter_factory-b))
        self.phi0=a/(self.pitch_diameter_factory/2)
        
        ref=[]
        self.RootDiameterActive()
        for t in npy.linspace(self.phi0,drap*self.phi_trochoide,discret):
            ref.append(vm.Point2D((self.Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        if ind=='T':
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            
        theta4=-self.root_angle/2
        p1=vm.Point2D((self.root_diameter/2*npy.cos(theta4),self.root_diameter/2*npy.sin(theta4)))
        p2=vm.Point2D((self.Trochoide(self.phi0,ind)))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        ref=primitives2D.RoundedLines2D([p1,p2],{},False)
        if ind=='T':
            L2=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L2=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        #ref=vm.Arc2D(p1,p2,p3)
        L=[L1,L2]

        return L
    
    def OutsideTrace(self,number):
        #trace du sommet des dents en arc de cercle
        theta4=npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*npy.cos(theta4),self.outside_diameter/2*npy.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        #ref=vm.Arc2D(p1,p2,p3)
        ref=primitives2D.RoundedLines2D([p1,p2,p3],{},False)
        L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        return L
        
    def UpdateProfilRack(self):
        
        #Initialisation rack position
        repere=self.PosInitRack()
        
        L=self.InvoluteOptimize()
        L.append(self.rack_profil0)
        return L
        
    def PosRack(self,alpha,list_number=[1]):
        
        angle=npy.tan(alpha)-self.rack.transverse_pressure_angle_T
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle)
        ptT=vm.Point2D(self.Involute(npy.tan(alpha)))
        DistY=npy.dot(ptT.vector-repere[0].vector,repere[2].vector-repere[0].vector)
        DistX=npy.dot(ptT.vector-repere[0].vector,repere[1].vector-repere[0].vector)
        decal=DistY+self.rack.circular_tooth_thickness-npy.tan(self.rack.transverse_pressure_angle_T)*DistX+npy.tan(self.rack.transverse_pressure_angle_T)*(self.rack.module*self.coefficient_profile_shift)
        temp=self.rack.RackComplete(list_number)
        for (i,j) in enumerate(temp):
            j=j.Rotation(vm.Point2D((0,0)),-npy.pi/2)
            j=j.Translation((self.pitch_diameter_factory/2+self.rack.module*self.coefficient_profile_shift,decal))
            j=j.Rotation(vm.Point2D((0,0)),angle)
            temp[i]=j
        rack_profil=temp
        return repere,rack_profil
        
    def EvolRack(self,alpha,repere,rack_profil):
        
        #Angular evolution of the rack
        delta=alpha*self.pitch_diameter_factory/2
        deltaU=-(repere[2].vector-repere[0].vector)*delta
        for j in range(3):
            repere[j]=repere[j].Rotation(vm.Point2D((0,0)),alpha)
        temp=rack_profil
        temp=temp.Translation((deltaU))
        rack_profil=temp.Rotation(vm.Point2D((0,0)),alpha)
        return repere,rack_profil 

    def BSpline(self,pt):
        x=[]
        y=[]
        for i in pt:
            x.append(i.vector[0])
            y.append(i.vector[1])
        tck, u=splprep([x, y], s=0)
        pas=0.01
        NewPT=splev(npy.arange(0,1+pas,pas), tck)
        pt=[]
        for i,j in zip(NewPT[0],NewPT[1]):
            pt.append(vm.Point2D((i,j)))
        return tck,pt
        
    def Mass(self):
        return self.VolumeModel().Mass()
    
    def VolumeModel(self):
        p=vm.Point3D((0,0,0))
        x=vm.Vector3D((1,0,0))
        y=vm.Vector3D((0,1,0))
        z=vm.Vector3D((0,0,self.width))
        wheel3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        lever3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        return vm.VolumeModel([wheel3D,lever3D])
    
    def GearGenerationSVGExport(self,name):

        L1=self.GearContours(20,[int(self.tooth_number)-1,0,1])
        L2=[]
        for i in npy.linspace(self.alpha_root_diameter_active,self.alpha_outside_diameter,5):
            L=self.PosRack(i,[-1,0,1])
            L2.extend(L[1])
        L3=self.TrochoideSecondary([0,self.tooth_number-1],'T',0.8,1000)
        Temp=self.TrochoideSecondary([0,self.tooth_number-1],'R',0.8,1000)
        L3.extend(Temp)
        #G1=vm.Contour2D(L1)
        #G1.MPLPlot()
        SVG1=LibSvg.SVGTrace(700,0,-npy.pi/2)
        SVG1.Convert(L1,'Rack','black',0.02,0)
        SVG1.Convert(L2,'Rack','black',0.01,0,'0.01px, 0.08px')
        SVG1.Convert(L3,'Rack','blue',0.03,0)
        SVG1.Export(name)
    
    def Dict(self):
        d=self.__dict__.copy()
        d['rack']=self.rack.Dict()
        return d
    
    def CSVExport(self):
        d=self.__dict__.copy()
        (d1,d2)=self.rack.CSVExport()
        del d['rack']
        del d['save']
        return list(d.keys())+d1,list(d.values())+d2
        
class GearAssembly(persistent.Persistent):
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle=0,coefficient_profile_shift1=0,
                 coefficient_profile_shift2=0,gear_width=20,maximum_torque=100,
                 transverse_pressure_angle_rack_T1=None,
                 transverse_pressure_angle_rack_T2=None,circular_tooth_thickness_rack1=None,
                 circular_tooth_thickness_rack2=None,
                 gear_addendum_rack1=None,gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
                 root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
                 root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                   coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                   transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                   gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                   root_radius_T1,root_radius_T2,root_radius_R1,
                   root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        self.Rack1=Rack(self.transverse_radial_pitch_rack1,self.transverse_pressure_angle_rack_T1,circular_tooth_thickness_rack1,
                        gear_addendum_rack1,gear_dedendum_rack1,root_radius_T1,root_radius_R1,transverse_pressure_angle_rack_R1)
        self.Rack2=Rack(self.transverse_radial_pitch_rack2,self.transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack2,
                        gear_addendum_rack2,gear_dedendum_rack2,root_radius_T2,root_radius_R2,transverse_pressure_angle_rack_R2)
        
        self.Gear1=Gear(self.Z1,self.Rack1,self.coefficient_profile_shift1)
        self.Gear2=Gear(self.Z2,self.Rack2,self.coefficient_profile_shift2)
        
        self.GearAssemblyParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
    def GearAssemblyParam(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.Z1=Z1
        self.Z2=Z2
        self.center_distance=center_distance
        self.DF1=2*self.center_distance*self.Z1/self.Z2/(1+self.Z1/self.Z2)
        self.DF2=2*self.center_distance-self.DF1
        self.gear_width=gear_width
        
        self.transverse_pressure_angle=transverse_pressure_angle
        self.transverse_radial_pitch=npy.pi*self.DF1/self.Z1
        self.helix_angle=helix_angle
        if npy.round(helix_angle,1)==0:
            self.normal_circular_pitch=self.transverse_radial_pitch
            self.axial_pitch=None
            self.normal_pressure_angle=self.transverse_pressure_angle
        else:
            self.normal_circular_pitch=self.transverse_radial_pitch*npy.cos(self.helix_angle)
            self.axial_pitch=self.normal_circular_pitch/npy.sin(helix_angle)
            self.normal_pressure_angle=npy.arctan(npy.tan(self.transverse_pressure_angle)*npy.cos(self.helix_angle))
            
        self.DB1=self.DF1*npy.cos(self.transverse_pressure_angle)
        self.DB2=self.DF2*npy.cos(self.transverse_pressure_angle)
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2

        self.transverse_pressure_angle_rack_T1=transverse_pressure_angle_rack_T1
        self.transverse_pressure_angle_rack_T2=transverse_pressure_angle_rack_T2
        if transverse_pressure_angle_rack_T1==None:
            self.transverse_pressure_angle_rack_T1=20/180*npy.pi
        if transverse_pressure_angle_rack_T2==None:
            self.transverse_pressure_angle_rack_T2=20/180*npy.pi
            
        self.pitch_diameter_factory1=self.DB1/npy.cos(self.transverse_pressure_angle_rack_T1)
        self.pitch_diameter_factory2=self.DB2/npy.cos(self.transverse_pressure_angle_rack_T2)
        self.transverse_radial_pitch_rack1=npy.pi*self.pitch_diameter_factory1/self.Z1
        self.transverse_radial_pitch_rack2=npy.pi*self.pitch_diameter_factory2/self.Z2
        
        self.maximum_torque=maximum_torque
        self.normal_load=self.maximum_torque*2/(self.DB1)
        self.tangential_load=self.maximum_torque*2/(self.DF1)
        self.radial_load=npy.tan(self.transverse_pressure_angle)*self.tangential_load
        
    def GearAssemblyParam2(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.transverse_base_pitch=self.transverse_radial_pitch*npy.cos(self.transverse_pressure_angle)
        
        #self.transverse_radial_pitch_rack1=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T1)
        #self.transverse_radial_pitch_rack2=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T2)
        self.circular_tooth_thickness1=self.Rack1.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T1)/npy.cos(self.transverse_pressure_angle)
        self.circular_tooth_thickness1=self.Gear1.GearSection(self.DF1)
        self.circular_tooth_thickness2=self.Rack2.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T2)/npy.cos(self.transverse_pressure_angle)
        self.circular_tooth_thickness2=self.Gear2.GearSection(self.DF2)
        self.space_width1=self.transverse_radial_pitch-self.circular_tooth_thickness1
        self.space_width2=self.transverse_radial_pitch-self.circular_tooth_thickness2
        
        self.linear_backlash=self.space_width1-self.circular_tooth_thickness2
        self.linear_backlash2=self.space_width2-self.circular_tooth_thickness1
        
        #formule spur gear
        self.radial_contact_ratio=1/2*(npy.sqrt(self.Gear1.outside_diameter**2-self.DB1**2)+npy.sqrt(self.Gear2.outside_diameter**2-self.DB2**2)-2*self.center_distance*npy.sin(self.transverse_pressure_angle))/(self.transverse_radial_pitch*npy.cos(self.transverse_pressure_angle))
        self.axial_contact_ratio=self.gear_width*npy.sin(self.helix_angle)/self.normal_circular_pitch
        self.total_contact_ratio=self.radial_contact_ratio+self.axial_contact_ratio
        
        #calcul de la hauteur pour la contrainte de Lewis
        self.gear_height_lewis1=self.Gear1.outside_diameter/2-self.Gear1.root_diameter_active/2-(self.Gear1.outside_active_angle/2*self.Gear1.outside_diameter/2*npy.tan(self.transverse_pressure_angle))
        self.gear_height_lewis2=self.Gear2.outside_diameter/2-self.Gear2.root_diameter_active/2-(self.Gear2.outside_active_angle/2*self.Gear2.outside_diameter/2*npy.tan(self.transverse_pressure_angle))
        
    def Update(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        self.Rack1.Update(self.transverse_radial_pitch_rack1,transverse_pressure_angle_rack_T1,circular_tooth_thickness_rack1,
                        gear_addendum_rack1,gear_dedendum_rack1,root_radius_T1,root_radius_R1,transverse_pressure_angle_rack_R1)
        self.Rack2.Update(self.transverse_radial_pitch_rack2,transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack2,
                        gear_addendum_rack2,gear_dedendum_rack2,root_radius_T2,root_radius_R2,transverse_pressure_angle_rack_R2)
        self.Gear1.Update(self.Z1,self.Rack1,coefficient_profile_shift1)
        self.Gear2.Update(self.Z2,self.Rack2,coefficient_profile_shift2)
        
        self.GearAssemblyParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
    def CriteriaIneq(self):
#        print(self.radial_contact_ratio,self.transverse_radial_pitch,self.transverse_pressure_angle)
        crit=[
              2*self.center_distance-self.Gear1.root_diameter_active-self.Gear2.outside_diameter,
              2*self.center_distance-self.Gear2.root_diameter_active-self.Gear1.outside_diameter,
              self.space_width1-self.circular_tooth_thickness2,
              self.space_width2-self.circular_tooth_thickness1,
              self.Gear1.outside_diameter-self.DF1,
              self.Gear2.outside_diameter-self.DF2,
              self.DF1-self.DB1,
              self.DF2-self.DB2,
              self.radial_contact_ratio-1,
              self.linear_backlash-0.05,
              0.3-self.linear_backlash
              ]
        return crit
        
    def CriteriaEq(self):
        
        crit=[
              self.Gear1.tooth_number-int(self.Gear1.tooth_number),
              self.Gear2.tooth_number-int(self.Gear2.tooth_number),
              ]
        return crit
        
    def InitialPosition(self):
        
        fun = (lambda tan_alpha : (norm(self.Gear1.Involute(tan_alpha))-(self.center_distance-self.DF2/2))**2)
#        bnds = (0,1)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle1=xsol[0]-npy.arctan(xsol[0])
        fun = (lambda tan_alpha : (norm(self.Gear2.Involute(tan_alpha))-(self.DF2/2))**2)
#        bnds = (0,1)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle2=xsol[0]-npy.arctan(xsol[0])
        
        Angle1=npy.arccos(self.Gear1.base_diameter/self.DF1)
        Angle2=npy.arccos(self.Gear2.base_diameter/self.DF2)
        
        #Gear2Angle=-Angle2+npy.pi-Angle1*self.Z1/self.Z2
        Gear1Angle=-(npy.tan(Angle1)-Angle1)
        Gear2Angle=-(npy.tan(Angle2)-Angle2)+npy.pi
        #Gear1Angle=0
        
        return [Gear1Angle,Gear2Angle]
        
    def WheelAssembly(self):
        
        Gear2Angle=self.InitialPosition()
        L1=self.Gear1.WheelContours(0,[-2,-1,0,1,2])
        L2=self.Gear2.WheelContours(Gear2Angle,[-2,-1,0,1,2])
        TG=[]
        for i in L2:
            TG.append(i.Translation((self.center_distance,0)))
        TG.extend(L1)
        return TG
        
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
    
    def SigmaLewis(self):
        
        self.sigma_lewis_maximum1=6*self.tangential_load*self.gear_height_lewis1/(self.gear_width*self.Gear1.root_gear_length**2)
        self.sigma_lewis_maximum2=6*self.tangential_load*self.gear_height_lewis2/(self.gear_width*self.Gear2.root_gear_length**2)
        
    def Dict(self):
        self.SigmaLewis()

        d=self.__dict__.copy()
        del d['Rack1']
        del d['Rack2']
        
        d['Gear1']=self.Gear1.Dict()
        d['Gear2']=self.Gear2.Dict()

        return d
    
    def CSVExport(self):
        self.SigmaLewis()
        d=self.__dict__.copy()
        (d1,d2)=self.Gear1.CSVExport()
        d1=['gear1_'+i for i in d1]
        (d3,d4)=self.Gear2.CSVExport()
        d3=['gear2_'+i for i in d3]
        return list(d.keys())+d1+d3,list(d.values())+d2+d4
    
    def SVGExport(self,name,position1,position2):
        #tuple1 et 2 correspondent a la position des centres
        TG1=self.Gear1.GearContours(10)
        TG2=self.Gear2.GearContours(10)
        list_rot=self.InitialPosition()
        L1=self.GearAssemblyTrace([TG1,TG2],[(0,0),(0,0)],list_rot)
        L2=[]
        L2.append(vm.Circle2D(vm.Point2D(position1),self.DF1/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.DF2/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.Gear2.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.pitch_diameter_factory1/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.pitch_diameter_factory2/2))
        L3=[]
        L3.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.root_diameter_active/2))
        L3.append(vm.Circle2D(vm.Point2D(position2),self.Gear2.root_diameter_active/2))
        #G1=vm.Contour2D(LR)
        #G1.MPLPlot()
        SVG1=LibSvg.SVGTrace(700)
        SVG1.Convert(L1[0],'G1','black',0.04,0)
        SVG1.Convert(L1[1],'G2','red',0.04,0)
        SVG1.Convert(L2,'Construction','blue',0.06,0,'0.1px, 0.3px')
        SVG1.Convert(L3,'Construction','red',0.03,0,'0.1px, 0.4px')
        SVG1.Export(name,{'G1':{'R':[2*npy.pi/self.Gear1.tooth_number,0,0]},'G2':{'R':[-2*npy.pi/self.Gear2.tooth_number,self.center_distance,0]}})
        
    def MeshingSVGExport(self,name,gear):
        if gear=='Z1':
            dent=self.Gear1
            diam_fonct=self.DF1
        elif gear=='Z2':
            dent=self.Gear2
            diam_fonct=self.DF2
        Ldev=[dent.InvoluteTrace(20,0,'T')]
        Ltroc=[dent.TrochoideTrace(40,0,'T')[0]]
        Lpied=[dent.TrochoideTrace(40,0,'T')[1]]
        Lout=[dent.OutsideTrace(0)]
        Lcomplet=dent.GearContours(10,[-1,0,1])
        sol=dent.Involute(npy.linspace(0,npy.tan(dent.alpha_outside_diameter),20))
        Lconst=[dent.Trace(sol[0],sol[1])]
        Lconst.extend(dent.TrochoideSecondary([0],'T',0.8,1000))
        L2=[vm.Circle2D(vm.Point2D((0,0)),dent.base_diameter/2)]
        L2.append(vm.Circle2D(vm.Point2D((0,0)),diam_fonct/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter_active/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.outside_diameter/2))
        L3=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2,0)))]
        L4=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2),dent.root_diameter/2*npy.sin(-dent.root_angle/2))))]
        L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2-dent.phi0),dent.root_diameter/2*npy.sin(-dent.root_angle/2-dent.phi0)))))
        L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2*npy.cos(-dent.root_angle/2-dent.phi_trochoide),dent.outside_diameter/2*npy.sin(-dent.root_angle/2-dent.phi_trochoide)))))
        SVG1=LibSvg.SVGTrace(700,1,-npy.pi/2)
        SVG1.Convert(Ldev,'Ldev','black',0.02,0)
        SVG1.Convert(Ltroc,'Ltroc','blue',0.02,0)
        SVG1.Convert(Lpied,'Lpied','red',0.02,0)
        SVG1.Convert(Lout,'Lout','green',0.02,0)
        SVG1.Convert(Lcomplet,'Lcomplet','black',0.01,1,'0.01px, 0.08px')
        SVG1.Convert(Lconst,'Gc','blue',0.01,1,'0.01px, 0.08px')
        SVG1.Convert(L2,'L2','black',0.01,1,'0.01px, 0.08px')
        SVG1.Convert(L3,'L3','black',0.03,1)
        SVG1.Convert(L4,'L4','black',0.01,1)
        SVG1.Export(name)



class ContinuousGearAssemblyOptimizer:
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           transverse_pressure_angle_rack_T1=None,
           transverse_pressure_angle_rack_T2=None,
           circular_tooth_thickness_rack1=None,circular_tooth_thickness_rack2=None,
           gear_addendum_rack1=None,
           gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
           root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
           root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,
           circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,
           gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)

        self.bounds={'center_distance':center_distance,
                 'transverse_pressure_angle':transverse_pressure_angle,
                 'helix_angle':helix_angle,
                 'coefficient_profile_shift1':coefficient_profile_shift1,
                 'coefficient_profile_shift2':coefficient_profile_shift2,
                 'gear_width':gear_width}
        self.Bounds=npy.array([center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,coefficient_profile_shift2,gear_width])
        
        self.xk={'Z1':self.Z1,'Z2':self.Z2,'maximum_torque':self.maximum_torque}
        self.xo={'transverse_pressure_angle_rack_T1':transverse_pressure_angle_rack_T1,
                 'transverse_pressure_angle_rack_T2':transverse_pressure_angle_rack_T2,
                 'circular_tooth_thickness_rack1':circular_tooth_thickness_rack1,
                 'circular_tooth_thickness_rack2':circular_tooth_thickness_rack2,
                 'gear_addendum_rack1':gear_addendum_rack1,
                 'gear_addendum_rack2':gear_addendum_rack2,
                 'gear_dedendum_rack1':gear_dedendum_rack1,
                 'gear_dedendum_rack2':gear_dedendum_rack2,
                 'root_radius_T1':root_radius_T1,
                 'root_radius_T2':root_radius_T2,
                 'root_radius_R1':root_radius_R1,
                 'root_radius_R2':root_radius_R2,
                 'transverse_pressure_angle_rack_R1':transverse_pressure_angle_rack_R1,
                 'transverse_pressure_angle_rack_R2':transverse_pressure_angle_rack_R2}
        
        self.InitX0()
        self.save=self.Xu[:]
        self.solutions=[]
        
        self.DefXU()
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.GearAssembly=GearAssembly(**self.xt)
        
    def Update(self,Xu):
        
        self.DefXU(Xu)
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.DefVar(**self.xt)
        self.GearAssembly.Update(**self.xt)
        self.save=Xu[:]
        
    def GearAssemblyParam(self,Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,
           circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,
           gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.Z1=Z1
        self.Z2=Z2
        self.maximum_torque=maximum_torque
        self.center_distance=[]
        self.transverse_pressure_angle=[]
        self.helix_angle=[]
        self.coefficient_profile_shift1=[]
        self.coefficient_profile_shift2=[]
        self.gear_width=[]
        
    def DefVar(self,Z1=None,Z2=None,center_distance=None,transverse_pressure_angle=None,
           helix_angle=None,coefficient_profile_shift1=None,
           coefficient_profile_shift2=None,gear_width=None,maximum_torque=None,
           transverse_pressure_angle_rack_T1=None,
           transverse_pressure_angle_rack_T2=None,
           circular_tooth_thickness_rack1=None,circular_tooth_thickness_rack2=None,
           gear_addendum_rack1=None,
           gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
           root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
           root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        self.Z1=Z1
        self.Z2=Z2
        self.maximum_torque=maximum_torque
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.helix_angle=helix_angle
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2
        self.gear_width=gear_width
    
    def DefXU(self,Xu=None):
        #Définition des inconnues utilisées dans l'optimiseur
        if type(Xu)==type(npy.array(())):
            self.center_distance=Xu[0,0]
            self.transverse_pressure_angle=Xu[1,0]
            self.helix_angle=Xu[2,0]
            self.coefficient_profile_shift1=Xu[3,0]
            self.coefficient_profile_shift2=Xu[4,0]
            self.gear_width=Xu[5,0]
            
        self.xu={'center_distance':self.center_distance,
                 'transverse_pressure_angle':self.transverse_pressure_angle,
                 'helix_angle':self.helix_angle,
                 'coefficient_profile_shift1':self.coefficient_profile_shift1,
                 'coefficient_profile_shift2':self.coefficient_profile_shift2,
                 'gear_width':self.gear_width}
        self.Xu=npy.array([[(self.center_distance),
                 (self.transverse_pressure_angle),
                 (self.helix_angle),
                 (self.coefficient_profile_shift1),
                 (self.coefficient_profile_shift2),
                 (self.gear_width)]])
        self.Xu=npy.transpose(self.Xu)
        
    def InitX0(self):
        
        self.DefXU()
        for i in self.bounds:
            self.i=((self.bounds[i][1]-self.bounds[i][0])*npy.random.random(1)+self.bounds[i][0])[0]
            self.xu[i]=self.i
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.DefVar(**self.xt)
        self.DefXU()
        
    def CriteriaEq(self):
        
#        crit=[self.GearAssembly.Rack1.transverse_pressure_angle_T-self.GearAssembly.Rack1.transverse_pressure_angle_R,
#              self.GearAssembly.Rack1.root_radius_T-self.GearAssembly.Rack1.root_radius_R,
#              self.GearAssembly.Rack2.transverse_pressure_angle_T-self.GearAssembly.Rack2.transverse_pressure_angle_R,
#              self.GearAssembly.Rack2.root_radius_T-self.GearAssembly.Rack2.root_radius_R]
#        crit=[self.transverse_radial_pitch_rack1-self.transverse_radial_pitch_rack2]
        crit=[0]
        return crit
        
    def CriteriaIneq(self):
        
        crit=[0]
        return crit
        
    def Objective(self,x):

        FEQ=self.feq(x)
        FINEQ=self.fineq(x)
        obj=0
#        for i in FINEQ:
#            obj=obj+(i**2)
        obj+=0.01*((self.GearAssembly.radial_contact_ratio-1.2)**2)
        obj+=0.01*((self.GearAssembly.linear_backlash-0.1)**2)
        for i in FINEQ:
            if i < 0:
                obj+=-10*i
            else:
                obj+=0.0001*i
        
        return obj
             
    def fineq(self,x):
        
#        x=npy.transpose([x])
#        self.GearAssembly.Update(x)
        if False in (npy.transpose([x])==self.save):
            x=npy.transpose([x])
            self.Update(x)
        self.save=self.Xu[:]
        
        ineq=[]
        #ineq.extend(self.GearAssembly.Gear1.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Gear2.CriteriaIneq())
#        ineq.extend(self.CriteriaIneq())
        ineq.extend(self.GearAssembly.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Rack1.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Rack2.CriteriaIneq())
        
        return ineq
    
        
    def feq(self,x):
        
        eq=[0]
#        eq.extend(self.CriteriaEq())
#        eq.extend(self.GearAssembly.CriteriaEq())
        return eq
        
    def Optimize(self):
        
        boucle=20
        i=0
        arret=0
        while i<boucle and arret==0:
#            print('Boucle d\'itération locale {}'.format(i))
#            dim=npy.shape(self.save)[0]
#            sol=npy.random.random(dim)
#            x0=(self.GearAssembly.bounds[:,1]-self.GearAssembly.bounds[:,0])*sol+self.GearAssembly.bounds[:,0]
            self.InitX0()
            x0=npy.array(self.Xu[:])
#            self.GearAssembly.Update(x0)
            self.Update(x0)

            cons = ({'type': 'eq','fun' : self.feq},{'type': 'ineq','fun' : self.fineq})
            cons = {'type': 'ineq','fun' : self.fineq}
            opt = {'maxiter':10000,'xtol': 1e-8, 'disp': False}
            #[x, _, _, imode,_] = fmin_slsqp(objective, x0,bounds=bounds, ieqcons=[const2], full_output=1,args=(p0,dR),iter=1000,iprint=0,acc=1e-1)
            #[x, _, imode,_] = root(const1, x0, method='lm')
#                cx = minimize(self.Objective, x0, bounds=self.GearAssembly.Bounds,constraints=cons,options=opt)
            cx = minimize(self.Objective, x0, bounds=self.Bounds, options=opt)
            FEQ=self.feq(cx.x)
            FINEQ=self.fineq(cx.x)
#            print('Status de convergence {}, Valeur de la fonctionnelle {}'.format(cx.status,cx.fun))
#                print(cx)
            if cx.fun<10 and min(FINEQ)>-1e-4:
                xsol=cx.x
                self.solutions.append(xsol)
                arret=1
#                print('Convergence atteinte avec le statut {}, Valeur de la fonctionnelle {}'.format(cx.status,cx.fun))
            i=i+1
            
            
class GearAssemblyOptimizer:
    def __init__(self,ratio,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,coefficient_profile_shift2,gear_width,maximum_torque):
        
        self.ratio=ratio
        self.Z1=Z1
        self.Z2=Z2
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.helix_angle=helix_angle
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2
        self.gear_width=gear_width
        self.maximum_torque=maximum_torque
        
        self.AnalyzeCombination()
        self.solutions=[]
        
    def AnalyzeCombination(self):
        
        mat1=list(npy.arange(self.Z1['min'],self.Z1['max']+1))
        mat2=list(npy.arange(self.Z2['min'],self.Z2['max']+1))
        plex=itertools.product(mat1,mat2)
        self.plex_calcul=[]
        for i in plex:
            drap=0
            ratio=i[0]/i[1]
            if self.ratio['min']==None:
                if not self.ratio['max']==None:
                    if ratio<=self.ratio['max']:
                        drap=1
                else:
                    drap=1
            else:
                if self.ratio['max']==None:
                    if ratio>=self.ratio['min']:
                        drap=1
                else:
                    if ratio>=self.ratio['min']:
                        if ratio<=self.ratio['max']:
                            drap=1
            if drap==1:
                Temp=[(i[0],i[0])]
                Temp1={}
                Temp1['Z1']=i[0]
                Temp.append((i[1],i[1]))
                Temp1['Z2']=i[1]
                data=[self.center_distance,self.transverse_pressure_angle,self.helix_angle,self.coefficient_profile_shift1,self.coefficient_profile_shift2,self.gear_width,self.maximum_torque]
                for j in data:
                    Temp.append((j['min'],j['max']))
                Temp1['center_distance']=(self.center_distance['min'],self.center_distance['max'])
                Temp1['transverse_pressure_angle']=(self.transverse_pressure_angle['min'],self.transverse_pressure_angle['max'])
                Temp1['helix_angle']=(self.helix_angle['min'],self.helix_angle['max'])
                Temp1['coefficient_profile_shift1']=(self.coefficient_profile_shift1['min'],self.coefficient_profile_shift1['max'])
                Temp1['coefficient_profile_shift2']=(self.coefficient_profile_shift2['min'],self.coefficient_profile_shift2['max'])
                Temp1['gear_width']=(self.gear_width['min'],self.gear_width['max'])
                Temp1['maximum_torque']=self.maximum_torque['min']
                self.plex_calcul.append(Temp1)
                
        
    def Optimize(self):
        lpx=len(self.plex_calcul)
        for ii,i in enumerate(self.plex_calcul):
            print('{}%'.format(ii/lpx*100))
            A1=ContinuousGearAssemblyOptimizer(**i)
            A1.Optimize()
            try:
                xsol=npy.transpose([A1.solutions[-1]])
                A1.DefXU(xsol)
                xt=dict(list(A1.xk.items())+list(A1.xu.items())+list(A1.xo.items()))
                self.solutions.append(GearAssembly(**xt))
            except:
                pass    
            
            
class GearAssemblyOptimizerWizard:
    
    def __init__(self,data):
        
        def DefZ(type,nom=None,min=None,max=None):
            if type=='Z1':
                self.Z1={'nom':nom,'min':min,'max':max}
            elif type=='Z2':
                self.Z2={'nom':nom,'min':min,'max':max}
                
        def DefRatio(type,nom=None,min=None,max=None,err=None,prog=None,prem=None,pgcd=None):
            self.ratio={'nom':nom,'min':min,'max':max,'err':err,'prog':prog,'prem':prem,'pgcd':pgcd}
            
        def DefGeneral(type,nom=None,min=None,max=None,err=None,prog=None,prem=None,pgcd=None):
            if type=='center_distance':
                self.center_distance={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='transverse_pressure_angle':
                self.transverse_pressure_angle={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='helix_angle':
                self.helix_angle={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='coefficient_profile_shift1':
                self.coefficient_profile_shift1={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='coefficient_profile_shift2':
                self.coefficient_profile_shift2={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='gear_width':
                self.gear_width={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
            elif type=='maximum_torque':
                self.maximum_torque={'nom':nom,'min':min,'max':max,'err':err,'prog':prog}
                
        #Initialisation des dictionnaires (passage à None l'ensemble des arg)
        DefZ('Z1')
        DefZ('Z2')
        DefRatio('ratio')
        DefGeneral('center_distance')
        DefGeneral('transverse_pressure_angle')
        DefGeneral('helix_angle')
        DefGeneral('coefficient_profile_shift1')
        DefGeneral('coefficient_profile_shift2')
        DefGeneral('gear_width')
        DefGeneral('maximum_torque')
        #Prise en compte des datas renseignees par la variable d'entree
        for i in data:
            if i['type'] in ['Z1','Z2']:
                DefZ(**i)
            elif i['type'] in ['ratio']:
                DefRatio(**i)
            else:
                DefGeneral(**i)
        
        #Analyse conformite des datas et lancement des calculs
        self.error=self.AnalyzeDataSet()
        if self.error==True:
            self.DefaultDataSet()
            
            
    def Optimize(self):
        M1=GearAssemblyOptimizer({'min':self.ratio['min'],'max':self.ratio['max']},
           {'min':self.Z1['min'],'max':self.Z1['max']},
           {'min':self.Z2['min'],'max':self.Z2['max']},
           {'min':self.center_distance['min'],'max':self.center_distance['max']},
           {'min':self.transverse_pressure_angle['min'],'max':self.transverse_pressure_angle['max']},
           {'min':self.helix_angle['min'],'max':self.helix_angle['max']},
           {'min':self.coefficient_profile_shift1['min'],'max':self.coefficient_profile_shift1['max']},
           {'min':self.coefficient_profile_shift2['min'],'max':self.coefficient_profile_shift2['max']},
           {'min':self.gear_width['min'],'max':self.gear_width['max']},
           {'min':self.maximum_torque['min'],'max':self.maximum_torque['max']})
        M1.Optimize()
        self.solutions=M1.solutions
        
            
    def AnalyzeDataSet(self):
        
        if self.ratio['nom']==None:
            if self.ratio['min']==None:
                if self.Z1['nom']==None:
                    if self.Z1['min']==None:
                        return False
                    elif self.Z1['max']==None:
                        return False
                elif self.Z2['nom']==None:
                    if self.Z2['min']==None:
                        return False
                    elif self.Z2['max']==None:
                        return False
            elif self.ratio['max']==None:
                if self.Z1['nom']==None:
                    if self.Z1['min']==None:
                        return False
                    elif self.Z1['max']==None:
                        return False
                elif self.Z2['nom']==None:
                    if self.Z2['min']==None:
                        return False
                    elif self.Z2['max']==None:
                        return False
        else:
            if not self.ratio['min']==None:
                return False
            elif not self.ratio['max']==None:
                return False
                    
        data=[self.Z1,self.Z2,self.center_distance,self.transverse_pressure_angle,self.helix_angle,self.coefficient_profile_shift1,self.coefficient_profile_shift2,self.gear_width,self.maximum_torque]
        for i in data:
            if not i['nom']==None:
                if not i['min']==None:
                    return False
                elif not i['max']==None:
                    return False
            
        return True
        
    def DefaultDataSet(self):
        
        self.bounds={}
        if not self.ratio['nom']==None:
            if self.ratio['err']==None:
                self.ratio['err']=0.05
            self.ratio['min']=self.ratio['nom']*(1-self.ratio['err'])
            self.ratio['max']=self.ratio['nom']*(1+self.ratio['err'])
            
        if self.Z1['nom']==None:
            if self.Z1['min']==None:
                self.Z1['min']=10
            if self.Z1['max']==None:
                self.Z1['max']=100
        else:
            self.Z1['min']=self.Z1['nom']
            self.Z1['max']=self.Z1['nom']
            
        if self.Z2['nom']==None:
            if self.Z2['min']==None:
                self.Z2['min']=10
            if self.Z2['max']==None:
                self.Z2['max']=100
        else:
            self.Z2['min']=self.Z2['nom']
            self.Z2['max']=self.Z2['nom']

        def analyze(i,err,mini,maxi):
            if not i['nom']==None:
                if i['err']==None:
                    i['err']=err
                i['min']=i['nom']*(1-i['err'])
                i['max']=i['nom']*(1+i['err'])
            elif i['min']==None:
                i['min']=mini
                if i['max']==None:
                    i['max']=maxi
            elif i['max']==None:
                i['max']=maxi
                
        analyze(self.helix_angle,0.05,0.25,0.4)
        analyze(self.transverse_pressure_angle,0.05,0.25,0.4)
        analyze(self.center_distance,0.05,40,100)
        analyze(self.coefficient_profile_shift1,0.05,-1,1)
        analyze(self.coefficient_profile_shift2,0.05,-1,1)
        analyze(self.gear_width,0.05,15,25)
        analyze(self.maximum_torque,0.05,100,150)
        

    
class GearAssemblyOptimizationResults(persistent.Persistent):
    
    def __init__(self,list_solutions,bounds):
        
        self.solutions=list_solutions
        self.input_data=bounds
#        self.solutions[family]={}
#        self.solutions[family]['obj']=[]
#        self.solutions[family]['bnds']=[]
#        self.Add(list_solutions,bounds,family)
#        self.type='mc_gears_assembly'
#            
#    def Add(self,list_solutions,bounds,family):
#        
#        for i in list_solutions:
#            self.solutions[family]['obj'].append(i)
#            self.solutions[family]['bnds'].append(bounds)
    
    def CSVExport(self,name,opt='w',family='Famille_A'):
        if self.solutions!=[]:
            (temp1,temp2)=self.solutions[family]['obj'][0].CSVExport()
            temp=temp1[0]
            for i in temp1[1::]:
                temp+=','+i
            if opt=='a':
                fichier=open(name,'r')
                temp1=fichier.read()
                temp1=temp1.split('\n')[0]
                temp1=temp1.split(',')
                fichier.close()
            fichier=open(name,opt)
            if not opt=='a':
                fichier.write(temp+'\n')
            for GA in self.solutions[family]['obj']:
                (temp3,temp4)=GA.CSVExport()
                temp=''
                for i in temp1:
                    add=temp3.index(i)
                    temp+=str(temp4[add])+','
                fichier.write(temp[0:-1]+'\n')
            fichier.close()
            
    def Dict(self):
        d=[]
        for ga in self.solutions:
            d.append(ga.Dict())
        return d
        
#class GearAssemblyDBClient(ResultsDBClient):
#    def __init__(self,address):
#        ResultsDBClient.__init__(address,'mc_gear_assembly')