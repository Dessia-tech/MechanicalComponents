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

#import pyDOE

class Rack:
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
    
class Gear:
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
        print(DistX,DistY,angle)
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
        
class AssemblyGear:
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle=0,coefficient_profile_shift1=0,
                 coefficient_profile_shift2=0,gear_width=20,maximum_torque=100,
                 transverse_pressure_angle_rack_T1=None,
                 transverse_pressure_angle_rack_T2=None,circular_tooth_thickness_rack1=None,
                 circular_tooth_thickness_rack2=None,
                 gear_addendum_rack1=None,gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
                 root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
                 root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.AssemblyGearParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
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
        
        self.AssemblyGearParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
    def AssemblyGearParam(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
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
        
    def AssemblyGearParam2(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.transverse_base_pitch=self.transverse_radial_pitch*npy.cos(self.transverse_pressure_angle)
        
        #self.transverse_radial_pitch_rack1=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T1)
        #self.transverse_radial_pitch_rack2=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T2)
        self.circular_tooth_thickness1=self.Rack1.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T1)/npy.cos(self.transverse_pressure_angle)
        self.circular_tooth_thickness2=self.Rack1.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T2)/npy.cos(self.transverse_pressure_angle)
        self.space_width1=self.transverse_radial_pitch-self.circular_tooth_thickness1
        self.space_width2=self.transverse_radial_pitch-self.circular_tooth_thickness2

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
        
        self.AssemblyGearParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
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
        
        self.AssemblyGearParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
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
              self.radial_contact_ratio-1.2,
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
        bnds = (0,1)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle1=xsol[0]-npy.arctan(xsol[0])
        fun = (lambda tan_alpha : (norm(self.Gear2.Involute(tan_alpha))-(self.DF2/2))**2)
        bnds = (0,1)
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
        
#class MasterAssemblyGear:
#    def __init__(self,dico_bounds):
#            
#        self.ratio=dico_bounds['ratio']
#        self.Z1=dico_bounds['Z1']
#        self.Z2=dico_bounds['Z2']
#        self.center_distance=dico_bounds['center_distance']
#        self.transverse_pressure_angle=dico_bounds['transverse_pressure_angle']
#        self.helix_angle=dico_bounds['helix_angle']
#        self.coefficient_profile_shift1=dico_bounds['coefficient_profile_shift1']
#        self.coefficient_profile_shift2=dico_bounds['coefficient_profile_shift2']
#        self.gear_width=dico_bounds['gear_width']
#        self.maximum_torque=dico_bounds['maximum_torque']
#        
#    def AnalyzeCombination(self):
#        
#        case=[]
#        for i in [self.ratio,self.Z1,self.Z2]:
#            if i==None:
#                case.append(0)
#            if type(i)==int:
#                case.append(1)
#            if type(i)==tuple:
#                case.append(2)
#                
#        if case in [[0,1,2],[0,2,1],[0,2,2]]]:
#            itertool.product
#        if case in [[0,1,1]]:
#            
#        if case in [[2,2,2],[2,2,1],[2,1,2]]
#        
#    def Optimize(self):
#        
#        list_bounds=[(30,30),(56,56),(60,70),(0.01,0.5),(0.3,0.3),(-1,1),(-1,1),(20,20),(200,200)]
#        A1=AssemblyGearOptimize1(list_bounds)
#        O1=Optimization(A1)
#        O1.Optimize()
#        
#    def Storage(self):
        
        
#class ComplexeAssemblyGear:
#    def __init__(self):
#        list_gear=[(13,(0,0.3)),(56,(-0.4,0.3)),((30,35),(0,0.3)),(82,(-0.5,0.3))]
#        list_assembly=[(0,1),(1,2),(2,3)]
#        list_axis=[{0:[(0,0),(0,0)]},{1:[(0,100),(0,100)]},{2:[(0,100),(0,100)]},{3:[(150,150),(150,150)]}]
#        list_entraxe=[(40,100),(40,100),(40,100)]
#        list_pression=[()]
        
        
class AssemblyGearOptimize1:
    def __init__(self,list_bounds,x=None):
        
        self.bounds=npy.array(list_bounds)
        if x==None:
            x=self.InitX0()
        self.AssemblyGearParam(x)
        self.save=x[:]
        self.AssemblyGear=AssemblyGear(self.Z1,self.Z2,self.center_distance,self.transverse_pressure_angle,
                                       self.helix_angle,self.coefficient_profile_shift1,
                                       self.coefficient_profile_shift2,self.gear_width,self.maximum_torque,
                                       self.transverse_pressure_angle_rack_T1,
                                       self.transverse_pressure_angle_rack_T2,
                                       self.circular_tooth_thickness_rack1,self.circular_tooth_thickness_rack2,
                                       self.gear_addendum_rack1,
                                       self.gear_addendum_rack2,self.gear_dedendum_rack1,self.gear_dedendum_rack2,
                                       self.root_radius_T1,self.root_radius_T2,self.root_radius_R1,
                                       self.root_radius_R2,self.transverse_pressure_angle_rack_R1,self.transverse_pressure_angle_rack_R2)
    
    def InitX0(self):
        
        dim=npy.shape(self.bounds)[0]
        sol=npy.random.random(dim)
        x0=(self.bounds[:,1]-self.bounds[:,0])*sol+self.bounds[:,0]
        return x0
        
    def AssemblyGearParam(self,x):
        
        self.Z1=x[0]
        self.Z2=x[1]
        self.center_distance=x[2]
        self.transverse_pressure_angle=x[3]
        self.helix_angle=x[4]
        self.coefficient_profile_shift1=x[5]
        self.coefficient_profile_shift2=x[6]
        self.gear_width=x[7]
        self.maximum_torque=x[8]
        self.transverse_pressure_angle_rack_T1=None
        self.transverse_pressure_angle_rack_T2=None
        self.circular_tooth_thickness_rack1=None
        self.circular_tooth_thickness_rack2=None
        self.gear_addendum_rack1=None
        self.gear_addendum_rack2=None
        self.gear_dedendum_rack1=None
        self.gear_dedendum_rack2=None
        self.root_radius_T1=None
        self.root_radius_T2=None
        self.root_radius_R1=None
        self.root_radius_R2=None
        self.transverse_pressure_angle_rack_R1=None
        self.transverse_pressure_angle_rack_R2=None
        
    def Update(self,x):
        
        self.AssemblyGearParam(x)
        self.AssemblyGear.Update(self.Z1,self.Z2,self.center_distance,self.transverse_pressure_angle,
                                       self.helix_angle,self.coefficient_profile_shift1,
                                       self.coefficient_profile_shift2,self.gear_width,self.maximum_torque,
                                       self.transverse_pressure_angle_rack_T1,
                                       self.transverse_pressure_angle_rack_T2,
                                       self.circular_tooth_thickness_rack1,self.circular_tooth_thickness_rack2,
                                       self.gear_addendum_rack1,
                                       self.gear_addendum_rack2,self.gear_dedendum_rack1,self.gear_dedendum_rack2,
                                       self.root_radius_T1,self.root_radius_T2,self.root_radius_R1,
                                       self.root_radius_R2,self.transverse_pressure_angle_rack_R1,self.transverse_pressure_angle_rack_R2)
        self.save=x[:]

    def CriteriaEq(self):
        
#        crit=[self.AssemblyGear.Rack1.transverse_pressure_angle_T-self.AssemblyGear.Rack1.transverse_pressure_angle_R,
#              self.AssemblyGear.Rack1.root_radius_T-self.AssemblyGear.Rack1.root_radius_R,
#              self.AssemblyGear.Rack2.transverse_pressure_angle_T-self.AssemblyGear.Rack2.transverse_pressure_angle_R,
#              self.AssemblyGear.Rack2.root_radius_T-self.AssemblyGear.Rack2.root_radius_R]
        crit=[self.transverse_radial_pitch_rack1-self.transverse_radial_pitch_rack2]
        return crit
        
    def CriteriaIneq(self):
        
        #crit=[self.AssemblyGear.radial_contact_ratio-0.6]
        crit=[0]
        return crit
        
class Optimization:
    def __init__(self,Assembly):
        
        self.AssemblyGear=Assembly
        self.solution=[]
        self.solution2=[]
        
    def Objective(self,x):

        self.AssemblyGear.Update(x)
        FEQ=self.feq(x)
        FINEQ=self.fineq(x)
        obj=0
        for i in FEQ:
            obj=obj+(i**2)

        return obj
             
    def fineq(self,x):
        
        self.AssemblyGear.Update(x)
        ineq=[]
        #ineq.extend(self.AssemblyGear.AssemblyGear.Gear1.CriteriaIneq())
        #ineq.extend(self.AssemblyGear.AssemblyGear.Gear2.CriteriaIneq())
        ineq.extend(self.AssemblyGear.CriteriaIneq())
        ineq.extend(self.AssemblyGear.AssemblyGear.CriteriaIneq())
        #ineq.extend(self.AssemblyGear.AssemblyGear.Rack1.CriteriaIneq())
        #ineq.extend(self.AssemblyGear.AssemblyGear.Rack2.CriteriaIneq())
        return ineq
    
        
    def feq(self,x):
        
        self.AssemblyGear.Update(x)
        eq=[]
#        eq.extend(self.AssemblyGear.CriteriaEq())
        eq.extend(self.AssemblyGear.AssemblyGear.CriteriaEq())
        
        return eq
        
    def Optimize(self):
        
        boucle=100
        i=0
        arret=0
        while i<boucle and arret==0:
            print(i)
            dim=npy.shape(self.AssemblyGear.save)[0]
            sol=npy.random.random(dim)
            x0=(self.AssemblyGear.bounds[:,1]-self.AssemblyGear.bounds[:,0])*sol+self.AssemblyGear.bounds[:,0]
#            self.AssemblyGear.Update(x0)
            self.AssemblyGear.Update(x0)
            cons = ({'type': 'eq','fun' : self.feq},{'type': 'ineq','fun' : self.fineq})
            opt = {'maxiter':1000}
            try:
                #[x, _, _, imode,_] = fmin_slsqp(objective, x0,bounds=bounds, ieqcons=[const2], full_output=1,args=(p0,dR),iter=1000,iprint=0,acc=1e-1)
                #[x, _, imode,_] = root(const1, x0, method='lm')
                cx = minimize(self.Objective, x0, bounds=self.AssemblyGear.bounds,constraints=cons,options=opt)
                FEQ=self.feq(cx.x)
                FINEQ=self.fineq(cx.x)
                #print(FEQ,FINEQ)
                if max(FEQ)<1e-6 and min(FEQ)>-1e-6 and min(FINEQ)>-1e-6:
                    xsol=cx.x
                    self.solution.append(xsol)
                    arret=1
                    print(self.AssemblyGear.AssemblyGear.Gear1.root_diameter_active,22)
                    print(FINEQ)
            except:
                print('erreur de nom')
            i=i+1
    
class SVGTrace:
    def __init__(self,scale=1,swift=0,angle=0,decal_x=0,decal_y=0):
        
        self.scale=scale
        self.swift=swift
        self.angle=angle
        self.decal_x=decal_x
        self.decal_y=decal_y
        self.data=[]
        self.box=[]
        self.zpd=[]
        self.begin=[]
        self.end=[]
        self.list_name={}
        self.indice=1
        self.InitBorne()
        
    def InitBorne(self):
        
        self.boxX_min=[npy.inf]
        self.boxX_max=[-npy.inf]
        self.boxY_min=[npy.inf]
        self.boxY_max=[-npy.inf]
        
    def UpdateBox(self,x1,y1,r1=None):
        
        X1=x1
        Y1=y1
        if not r1==None:
            print(x1,y1,r1)
            X1=[x1[0]-r1,x1[0],x1[0]+r1,x1[0]]
            Y1=[y1[0],y1[0]+r1,y1[0],y1[0]-r1]
            print(x1,y1,r1)
        self.boxX_min=min(self.boxX_min,[min(X1)])
        self.boxX_max=max(self.boxX_max,[max(X1)])
        self.boxY_min=min(self.boxY_min,[min(Y1)])
        self.boxY_max=max(self.boxY_max,[max(Y1)])
        
    def CreateBegin(self):
        
        scale_html=100-100/(1500/(1000/(self.boxX_max[0]-self.boxX_min[0])*(self.boxY_max[0]-self.boxY_min[0])))
        self.begin.append('''<html>
<head>
<meta charset="UTF-8"> 
<script src="snap.svg-min.js"></script>
    <script src="snap.svg.zpd.js"></script>

</script>
</head>

<body>
<svg id="Gear" width="100%" height="100%">

</svg>

<script>''')
        
    def CreateEnd(self,animate):
        
        if not animate==None:
            master=list(animate.keys())[0]
            i=0
            Temp=[]
            for (j,k) in animate.items():
                self.end.append('var myMatrix'+str(int(i))+' = new Snap.Matrix();')
                self.end.append('myMatrix'+str(int(i))+'.add(1,0,0,1,'+str(k['R'][1])+', '+str(k['R'][2])+');')
                self.end.append(j+'.transform(myMatrix'+str(int(i))+');')
                Temp.append(j+'.attr({ transform: myMatrix'+str(int(i))+'});')
                i+=1
                self.end.append('var myMatrix'+str(int(i))+' = new Snap.Matrix();')
                self.end.append('myMatrix'+str(int(i))+'.add('+str(npy.cos(k['R'][0]))+','+str(npy.sin(k['R'][0]))+','+str(-npy.sin(k['R'][0]))+','+str(npy.cos(k['R'][0]))+','+str(k['R'][1])+', '+str(k['R'][2])+');')
                i+=1
                
            Temp2=[]
            i=0
            for (j,k) in animate.items():
                if not j==master:
                    Temp2.append(str(j)+'''.animate(
			{ transform: myMatrix'''+str(int(i+1))+'''}, 
			1000)''')
                if j==master:
                    indice=i
                i+=2
                
            self.end.append('''function anim(){
		'''+str(master)+'''.animate(
			{ transform: myMatrix'''+str(int(indice+1))+'''}, 
			1000,
			function(){ 
				'''+Temp[0]+'''
                  '''+Temp[1]+'''
				anim();
			});
             '''+Temp2[0]+'''
	}''')
    
            self.end.append('s.click( anim);') 
        self.end.append('</script></body></html>')
        
    def CreateBox(self):
        
        self.view_x=(self.boxX_max[0]-self.boxX_min[0])
        self.view_y=(self.boxY_max[0]-self.boxY_min[0])
        self.box.append('var s = Snap("#Gear");')
#        self.box.append('var s = Snap('+str(self.view_x)+','+str(self.view_y)+');')
#        self.box.append('s.attr({ viewBox: "'+str(self.boxX_min[0]+self.decal_x)+' '+str(self.boxY_min[0]+self.decal_y)+' '+str(self.boxX_max[0]-self.boxX_min[0]+selfe.decal_x)+' '+str(self.boxY_max[0]-self.boxY_min[0]+self.decal_y)+'" });')

    def Transform(self,data,init=None):
        
        alpha=self.angle
        MatRot=npy.array([[npy.cos(alpha),-npy.sin(alpha),0],[npy.sin(alpha),npy.cos(alpha),0],[0,0,1]])
        if self.swift==1:
            MatPermut=npy.array([[0,1,0],[1,0,0],[0,0,1]])
            Mat=npy.matmul(MatRot,MatPermut)
        else:
            Mat=MatRot
        dataS=npy.matmul(Mat,data)
        
        if not init==None:
            self.InitBorne()
        self.UpdateBox(list(dataS[0,:]),list(dataS[1,:]))
        
        scaleX=self.scale/(self.boxX_max[0]-self.boxX_min[0])
        scaleY=self.scale/(self.boxY_max[0]-self.boxY_min[0])
        scaleZPD=min(scaleX,scaleY)
        MatScale=npy.array([[scaleZPD,0,-scaleZPD*self.boxX_min[0]],[0,scaleZPD,-scaleZPD*self.boxY_min[0]],[0,0,1]])
        
#        Mat=npy.matmul(Mat1,npy.array([[1,0,-self.boxX_min[0]],[0,1,-5],[0,0,1]]))
        Mat=npy.matmul(MatScale,Mat)
        return Mat
    
    def CreateZPD(self):
        
        data=npy.array([[self.boxX_min[0],self.boxX_max[0],self.boxX_max[0],self.boxX_min[0]],[self.boxY_max[0],self.boxY_max[0],self.boxY_min[0],self.boxY_min[0]],[1,1,1,1]])
        Mat=self.Transform(data,'ok')
        
#        self.zpd.append('s.zpd({ load: {a:'+str(scaleZPD)+',b:0,c:0,d:'+str(scaleZPD)+',e:'+str(-scaleZPD*self.boxX_min[0])+',f:'+str(-scaleZPD*self.boxY_min[0])+'}});')
        
        for n,m in self.list_name.items():
            Temp='var '+n+' = s.group('
            for i in m:
                Temp+=i+','
            self.zpd.append(Temp[0:-1]+');')
            
        self.zpd.append('s.zpd({ load: {a:'+str(Mat[0,0])+',b:'+str(Mat[1,0])+',c:'+str(Mat[0,1])+',d:'+str(Mat[1,1])+',e:'+str(Mat[0,2])+',f:'+str(Mat[1,2])+'}});')
        self.zpd.append('s.zpd({ drag: false}, function (err, paper) {console.log(paper)});')
                
    def ConvertPrimitive2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'primitives2D' in str(j.__class__):
                for n,m in enumerate(j.primitives):
                    if 'Line2D' in str(m.__class__):
                        temp='var primitive2D'+str(self.indice)+' = s.polyline(['
                        self.list_name[group_name].append('primitive2D'+str(self.indice))
                        self.indice+=1
                        for k in m.points[0:-1]:
                            if impact_box==0:
                                self.UpdateBox([k.vector[0]],[k.vector[1]])
                            temp+=str(k.vector[0])+','+str(k.vector[1])+','
                        temp+=str(m.points[-1].vector[0])+','+str(m.points[-1].vector[1])
                        temp+=']).attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"});'
                        self.data.append(temp)
                    if 'Arc2D' in str(m.__class__):
                        temp='var path'+str(self.indice)+' = s.path('
                        self.list_name[group_name].append('path'+str(self.indice))
                        self.indice+=1
                        temp+='\'M '+str(m.start.vector[0])+','+str(m.start.vector[1])
                        temp+=' A '+str(m.radius)+','+str(m.radius)+' '+str((m.angle1+m.angle2)/2/npy.pi*180+90)+' 0 '+' 1 '+str(m.end.vector[0])+' '+str(m.end.vector[1])
                        temp+='\').attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"});'
                        self.data.append(temp)
                        
            
    def ConvertCircle2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Circle2D' in str(j.__class__):
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]],j.radius)
                temp='var circle2D'+str(self.indice)+' = s.circle('
                self.list_name[group_name].append('circle2D'+str(self.indice))
                self.indice+=1
                temp+=str(j.center.vector[0])+','+str(j.center.vector[1])+','
                temp+=str(j.radius)
                temp+=').attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"});'
                self.data.append(temp)
                
    def ConvertLine2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Line2D' in str(j.__class__):
                temp='var primitive2D'+str(self.indice)+' = s.polyline(['
                self.list_name[group_name].append('primitive2D'+str(self.indice))
                self.indice+=1
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]])
                for k in j.points[0:-1]:
                    if impact_box==0:
                        self.UpdateBox([k.vector[0]],[k.vector[1]])
                    temp+=str(k.vector[0])+','+str(k.vector[1])+','
                temp+=str(j.points[-1].vector[0])+','+str(j.points[-1].vector[1])
                temp+=']).attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"});'
                self.data.append(temp)
    
    def Convert(self,L,group_name,stroke,strokeWidth,impact_box=0,strokeDasharray="none"):
        
        self.list_name[group_name]=[]
        self.ConvertPrimitive2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertCircle2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertLine2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
    
    def Export(self,name='export.html',animate=None):
        fichier=open(name,'w')
        self.CreateBegin()
        for i in self.begin:
            fichier.write(i+'\n')
        self.CreateZPD()
        self.CreateBox()
        for i in self.box:
            fichier.write(i+'\n')
        for i in self.data:
            fichier.write(i+'\n')
        for i in self.zpd:
            fichier.write(i+'\n')
        self.CreateEnd(animate)
        for i in self.end:
            fichier.write(i+'\n')
        fichier.close()

##Initialisation
#list_bounds=[(30,30),(56,56),(60,70),(0.01,0.5),(0.3,0.3),(-1,1),(-1,1),(20,20),(200,200)]
#A1=AssemblyGearOptimize1(list_bounds)
#O1=Optimization(A1)
#O1.Optimize()
##O1.Optimize2(list(O1.solution[-1]))

xsol=list(npy.array([30.,50.,115.41878235,0.38523225,-0.94639435,0.26922009,0.,20.]))
#xsol=list(npy.array([3.00000000e+01,5.10000000e+01,9.74954025e+01,4.25119842e-01,3.00000000e-01,9.95153575e-01,7.35600926e-02,2.00000000e+01,2.00000000e+02]))
#AG1=AssemblyGear(*list(O1.solution[-1]))
AG1=AssemblyGear(*xsol)

### Pied de dent
L=AG1.Gear1.TrochoideTrace(500,1,'T')
L.extend(AG1.Gear1.GearContours(10,[1]))
#G1=vm.Contour2D(L)
#G1.MPLPlot()
print(AG1.Gear1.root_diameter_active)

LL=AG1.Rack1.RackComplete([1])
SVG1=SVGTrace(700,0)
SVG1.Convert(LL,'G1','black',0.02,0)
SVG1.Export('Essai2.html')

### Trace cremaillere Gear1
L1=AG1.Gear1.GearContours(20,[int(AG1.Z1)-1,0,1])
L2=[]
for i in npy.linspace(AG1.Gear1.alpha_root_diameter_active,AG1.Gear1.alpha_outside_diameter,5):
    L=AG1.Gear1.PosRack(i,[-1,0,1])
    L2.extend(L[1])
L3=AG1.Gear1.TrochoideSecondary([0,AG1.Z1-1],'T',0.8,1000)
Temp=AG1.Gear1.TrochoideSecondary([0,AG1.Z1-1],'R',0.8,1000)
L3.extend(Temp)
#G1=vm.Contour2D(L1)
#G1.MPLPlot()
SVG1=SVGTrace(700,0,-npy.pi/2)
SVG1.Convert(L1,'Rack','black',0.02,0)
SVG1.Convert(L2,'Rack','black',0.01,0,'0.01px, 0.08px')
SVG1.Convert(L3,'Rack','blue',0.03,0)
SVG1.Export('Cremaillere_Z1.html')

### Détail dent 1
def DetailDent(assembly,dent,diam_base,diam_fonct,nom_html):
    Ldev=[dent.InvoluteTrace(20,0,'T')]
    Ltroc=[dent.TrochoideTrace(40,0,'T')[0]]
    Lpied=[dent.TrochoideTrace(40,0,'T')[1]]
    Lout=[dent.OutsideTrace(0)]
    Lcomplet=dent.GearContours(10,[-1,0,1])
    sol=dent.Involute(npy.linspace(0,npy.tan(dent.alpha_outside_diameter),20))
    Lconst=[dent.Trace(sol[0],sol[1])]
    Lconst.extend(dent.TrochoideSecondary([0],'T',0.8,1000))
    L2=[vm.Circle2D(vm.Point2D((0,0)),diam_base/2)]
    L2.append(vm.Circle2D(vm.Point2D((0,0)),diam_fonct/2))
    L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter_active/2))
    L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter/2))
    L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.outside_diameter/2))
    L3=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2,0)))]
    L4=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2),dent.root_diameter/2*npy.sin(-dent.root_angle/2))))]
    L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2-AG1.Gear1.phi0),dent.root_diameter/2*npy.sin(-dent.root_angle/2-dent.phi0)))))
    L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2*npy.cos(-dent.root_angle/2-dent.phi_trochoide),dent.outside_diameter/2*npy.sin(-dent.root_angle/2-dent.phi_trochoide)))))
    SVG1=SVGTrace(700,1,-npy.pi/2)
    SVG1.Convert(Ldev,'Ldev','black',0.02,0)
    SVG1.Convert(Ltroc,'Ltroc','blue',0.02,0)
    SVG1.Convert(Lpied,'Lpied','red',0.02,0)
    SVG1.Convert(Lout,'Lout','green',0.02,0)
    SVG1.Convert(Lcomplet,'Lcomplet','black',0.01,1,'0.01px, 0.08px')
    SVG1.Convert(Lconst,'Gc','blue',0.01,1,'0.01px, 0.08px')
    SVG1.Convert(L2,'L2','black',0.01,1,'0.01px, 0.08px')
    SVG1.Convert(L3,'L3','black',0.03,1)
    SVG1.Convert(L4,'L4','black',0.01,1)
    SVG1.Export(nom_html)

DetailDent(AG1,AG1.Gear1,AG1.DB1,AG1.DF1,'Creation_Dent1.html')
DetailDent(AG1,AG1.Gear2,AG1.DB2,AG1.DF2,'Creation_Dent2.html')

### Trace cremaillere Gear2
L1=AG1.Gear2.GearContours(20,[int(AG1.Z2)-1,0,1])
L2=[]
for i in npy.linspace(AG1.Gear2.alpha_root_diameter_active,AG1.Gear2.alpha_outside_diameter,5):
    L=AG1.Gear2.PosRack(i,[-1,0,1])
    L2.extend(L[1])
L3=AG1.Gear2.TrochoideSecondary([0,AG1.Z2-1],'T',0.8,1000)
Temp=AG1.Gear2.TrochoideSecondary([0,AG1.Z2-1],'R',0.8,1000)
L3.extend(Temp)
#G1=vm.Contour2D(L1)
#G1.MPLPlot()
SVG1=SVGTrace(700,0,-npy.pi/2)
SVG1.Convert(L1,'Rack','black',0.02,0)
SVG1.Convert(L2,'Rack','black',0.01,0,'0.01px, 0.08px')
SVG1.Convert(L3,'Rack','blue',0.03,0)
SVG1.Export('Cremaillere_Z2.html')

######################
### Trace engrenenemnt
#xsol=list(O1.solution[-1])
TG1=AG1.Gear1.GearContours(10)
TG2=AG1.Gear2.GearContours(10)
list_rot=AG1.InitialPosition()
L1=AG1.GearAssemblyTrace([TG1,TG2],[(0,0),(0,0)],list_rot)
L2=[]
L2.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DF1/2))
L2.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.DF2/2))
L2.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DB1/2))
L2.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.DB2/2))
L2.append(vm.Circle2D(vm.Point2D((0,0)),AG1.pitch_diameter_factory1/2))
L2.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.pitch_diameter_factory2/2))
L3=[]
L3.append(vm.Circle2D(vm.Point2D((0,0)),AG1.Gear1.root_diameter_active/2))
L3.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.Gear2.root_diameter_active/2))
#G1=vm.Contour2D(LR)
#G1.MPLPlot()
SVG1=SVGTrace(700)
SVG1.Convert(L1[0],'G1','black',0.04,0)
SVG1.Convert(L1[1],'G2','red',0.04,0)
SVG1.Convert(L2,'Construction','blue',0.06,0,'0.1px, 0.3px')
SVG1.Convert(L3,'Construction','red',0.03,0,'0.1px, 0.4px')
data=npy.array([[0,0],[0,AG1.center_distance],[1,1]])
SVG2=SVG1
dataP=npy.matmul(SVG2.Transform(data),data)
SVG1.Export('Gear.html',{'G1':{'R':[2*npy.pi/AG1.Z1,0,0]},'G2':{'R':[-2*npy.pi/AG1.Z2,AG1.center_distance,0]}})

