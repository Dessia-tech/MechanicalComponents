import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *
from scipy.interpolate import splprep, splev

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
        
    def RackContours(self,number):
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*npy.tan(self.transverse_pressure_angle_R),self.gear_addendum))
        p4=p1.Translation((self.circular_tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*npy.tan(self.transverse_pressure_angle_T),self.gear_addendum))
        p5=p4.Translation((self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T),-self.gear_dedendum))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_R),-self.gear_dedendum))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.root_radius_T,5:self.root_radius_R},False)
        Rack_Elem=[L]
        for i in range(number):
            Rack_Elem.append(L.Translation(((i+1)*(p7.vector-p1.vector))))
        return Rack_Elem
        
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
        
    def GearContours(self,discret=10,number=None):
        #Analytical involute profil
        self.RootDiameterActive()
        if number==None:
            number=int(self.tooth_number)
        L=[self.InvoluteTrace(discret,0,'T')]
        L.extend(self.TrochoideTrace(2*discret,0,'T'))
        L.append(self.InvoluteTrace(discret,0,'R'))
        L.extend(self.TrochoideTrace(2*discret,0,'R'))
        L.append(self.OutsideTrace(0))
        #print(L)
        for i in range(1,number):
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
    
    def TrochoideTrace(self,discret,number,ind='T'):
        if ind=='T':
            drap=1
        else:
            drap=-1
        
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        theta0=npy.arctan(abs(a)/(self.pitch_diameter_factory-b))
        
        ref=[]
        self.RootDiameterActive()
        for t in npy.linspace(drap*theta0,self.phi_trochoide,discret):
            ref.append(vm.Point2D((self.Trochoide(drap*t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        
        if ind=='T':
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            
        theta4=-self.root_angle/2
        p1=vm.Point2D((self.root_diameter/2*npy.cos(theta4),self.root_diameter/2*npy.sin(theta4)))
        p2=vm.Point2D((self.Trochoide(theta0,ind)))
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
        
    def PosRack(self,alpha):
        
        angle=npy.tan(alpha)-self.rack.transverse_pressure_angle_T
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle)
        ptT=vm.Point2D(self.Involute(npy.tan(alpha)))
        DistY=npy.dot(ptT.vector-repere[0].vector,repere[2].vector-repere[0].vector)
        DistX=npy.dot(ptT.vector-repere[0].vector,repere[1].vector-repere[0].vector)
        decal=DistY+self.rack.circular_tooth_thickness-npy.tan(self.rack.transverse_pressure_angle_T)*DistX
        temp=self.rack.RackContours(0)[0]
        temp=temp.Rotation(vm.Point2D((0,0)),-npy.pi/2)
        temp=temp.Translation((self.pitch_diameter_factory/2,decal))
        rack_profil=temp.Rotation(vm.Point2D((0,0)),angle)
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
    
    def DistRackInvolute(self):
        
        repere,rack_profil=self.PosRack(self.alpha_root_diameter_active)
        
        l1=rack_profil.primitives[5]
        p1=l1.points[0]
        p2=l1.points[1]
        p3=p2.Translation((p2.vector-p1.vector))
        l1=vm.Line2D(p1,p3)
        
        L=self.InvoluteTrace(10,1,'R')
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])

        tck1,pt1=self.BSpline(PT)
        p8=vm.Point2D.MiddlePoint(p1,p2)
        p9=vm.Point2D.MiddlePoint(p2,p3)
        tck2,pt2=self.BSpline([p1,p8,p2,p9,p3])        
        fun = (lambda t : norm(npy.array([splev(t[0],tck1)[0],splev(t[0],tck1)[1]])-npy.array([splev(t[1],tck2)[0],splev(t[1],tck2)[1]]))**2)
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,0], bounds=bnds)
        xsol=sol.x
        dist=fun(xsol)
        
#        L=[primitives2D.RoundedLines2D(PT,{},False)]
#        L.append(l1)
#        L.append(vm.Point2D(splev(xsol[0],tck1)))
#        L.append(vm.Point2D(splev(xsol[1],tck2)))
#        G2=vm.Contour2D(L)
#        G2.MPLPlot()
        
        return dist
        
    def TrochoidalRackEstimation(self,discret=10):
        
        #Initialisation rack position
        alpha=self.alpha_root_diameter_active
        repere,self.rack_profilf=self.PosRack(alpha)
        alpha=self.alpha_outside_diameter
        repere,self.rack_profil0=self.PosRack(alpha)
        
        alpha=2.5*npy.pi/self.tooth_number
        repere,rack_profil=self.EvolRack(alpha,repere[:],self.rack_profil0)
        
        #Primary trochoide
        nb=discret
        liste_ligneT=[rack_profil.primitives[2]]
        liste_pointT=[liste_ligneT[0].points[0]]
        liste_ligneR=[rack_profil.primitives[5]]
        liste_pointR=[liste_ligneR[0].points[0]]
        liste_circleT=[rack_profil.primitives[-2]]
        liste_circleR=[rack_profil.primitives[-1]]
        Ltrajet=[rack_profil,rack_profil.Translation(-(self.circular_tooth_thickness+self.tooth_space)*(repere[2].vector-repere[0].vector))]
        for i in range(nb-1):
            Alpha=(-5*alpha)/nb
            repere,rack_profil=self.EvolRack(Alpha,repere[:],rack_profil)
            ligne_temp=rack_profil.primitives[2]
            liste_pointT.append(vm.Point2D.LinesIntersection(liste_ligneT[-1],ligne_temp))
            liste_ligneT.append(ligne_temp)
            liste_circleT.append(rack_profil.primitives[-2])
            ligne_temp=rack_profil.primitives[5]
            liste_pointR.append(vm.Point2D.LinesIntersection(liste_ligneR[-1],ligne_temp))
            liste_ligneR.append(ligne_temp)
            liste_circleR.append(rack_profil.primitives[-1])
            Ltrajet.append(rack_profil)
            Ltrajet.append(rack_profil.Translation(-(self.circular_tooth_thickness+self.tooth_space)*(repere[2].vector-repere[0].vector)))
        TrochoidalPrimaryT=[]
        TrochoidalPrimaryR=[]
        for i in liste_circleT:
            TrochoidalPrimaryT.append(i.center)
        for i in liste_circleR:
            TrochoidalPrimaryR.append(i.center)
        self.TrochoidalPrimaryControlT,pt1=self.BSpline(TrochoidalPrimaryT)
        self.TrochoidalPrimaryControlR,pt1=self.BSpline(TrochoidalPrimaryR)
        self.Ltrajet=Ltrajet
        
        #Secondary trochoide of the drive coast 
        cm=liste_circleT[-1].center
        cp=liste_circleT[-2].center
        u=cp.vector-cm.vector
        v=npy.array([u[1],-u[0]])/norm(u)
        liste_ptT=[vm.Point2D(cm.vector+v*self.rack.root_radius_T)]
        liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.root_radius_T))
        for i in liste_circleT[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptT.append(vm.Point2D(cm.vector+v*self.rack.root_radius_T))
            liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.root_radius_T))
            cm=cp
        self.TrochoidalSecondaryControlT,pt1=self.BSpline(liste_ptT[2::])
        
        #Secondary trochoide of the rear coast 
        cm=liste_circleR[-1].center
        cp=liste_circleR[-2].center
        u=cp.vector-cm.vector
        v=npy.array([u[1],-u[0]])/norm(u)
        liste_ptR=[vm.Point2D(cm.vector+v*self.rack.root_radius_T)]
        liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.root_radius_T))
        for i in liste_circleR[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptR.append(vm.Point2D(cm.vector+v*self.rack.root_radius_T))
            liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.root_radius_T))
            cm=cp
        self.TrochoidalSecondaryControlR,pt1=self.BSpline(liste_ptR[2::])
        
        #Construction arc cercle pied et tete de dent
        angle=self.base_circular_tooth_thickness*2/self.base_diameter-(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)
        liste_pointDE=[vm.Point2D((self.outside_diameter/2*npy.cos(angle),self.outside_diameter/2*npy.sin(angle)))]
        nb=discret
        for i in range(nb):
            liste_pointDE.append(liste_pointDE[-1].Rotation(vm.Point2D((0,0)),-angle/nb))
        self.DEControl,pt1=self.BSpline(liste_pointDE)
        
        #Estimation new internal radius
        angle=npy.pi/self.tooth_number
        Control1=self.TrochoidalSecondaryControlT
        fun = (lambda t : norm(npy.array([splev(t,Control1)[0],splev(t,Control1)[1]])))
        sol=minimize(fun,0)
        self.root_diameter_factory=2*fun(sol.x)
        liste_pointDI=[vm.Point2D((self.root_diameter_factory/2*npy.cos(angle),self.root_diameter_factory/2*npy.sin(angle)))]
        delta_angle=angle+(2*npy.pi/self.tooth_number)
        for i in range(nb):
            liste_pointDI.append(liste_pointDI[-1].Rotation(vm.Point2D((0,0)),-delta_angle/nb))
        self.DIControl,pt1=self.BSpline(liste_pointDI)
        
        return liste_ptT
        
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
        
    def GearContoursRack(self,discret=10):
        
        discret=10
        post=self.TrochoidalRackEstimation(discret)
        
        L=self.InvoluteTrace(10,0,'T')
        InvoluteL=L
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])
        PT.append(i.points[1])
        temp=PT
        self.InvoluteControlT,InvolutePT=self.BSpline(PT)
        L=self.InvoluteTrace(10,1,'R')
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])
        PT.append(i.points[1])
        self.InvoluteControlR,InvolutePR=self.BSpline(PT)
        
        Control1=self.TrochoidalSecondaryControlT
        Control2=self.DIControl
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]])))
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,0], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolI1=sol.x
        Control1=self.TrochoidalSecondaryControlR
        Control2=self.DIControl
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]])))
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,1], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolI2=sol.x
        
        Control1=self.TrochoidalSecondaryControlT
        Control2=self.InvoluteControlT
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((-0.1,xsolI1[0]), (0,1))
        sol=minimize(fun,[0,0.5],bounds=bnds, method='SLSQP', tol=1e-20)
        xsolT=sol.x
        Control1=self.TrochoidalSecondaryControlR
        Control2=self.InvoluteControlR
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((xsolI2[0],1.1), (0,1))
        sol=minimize(fun,[1,0.5], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolR=sol.x
        
        Control1=self.DEControl
        Control2=self.InvoluteControlT
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((0,1), (0,1))
        sol=minimize(fun,[0,1],bounds=bnds, method='SLSQP', tol=1e-20)
        xsolE1=sol.x
        
        PT=[]
        dis=20
        pas=(xsolE1[0])/dis
        for j in range(dis):
            i=j*pas
            PT.append(vm.Point2D((splev(i,self.DEControl))))       
        pas=(xsolE1[1]-xsolT[1])/dis
        for j in range(dis):
            i=(xsolT[1]+dis*pas)-(j*pas)
            PT.append(vm.Point2D((splev(i,self.InvoluteControlT))))
        pas=(xsolT[0]-xsolI1[0])/dis
        for j in range(dis):
            i=(xsolI1[0]+dis*pas)-(j*pas)
            PT.append(vm.Point2D((splev(i,self.TrochoidalSecondaryControlT))))
        pas=(xsolI2[1]-xsolI1[1])/dis
        for j in range(dis):
            i=xsolI1[1]+j*pas
            PT.append(vm.Point2D((splev(i,self.DIControl))))
        pas=(xsolI2[0]-xsolR[0])/dis
        for j in range(dis):
            i=(xsolR[0]+dis*pas)-(j*pas)
            PT.append(vm.Point2D((splev(i,self.TrochoidalSecondaryControlR))))
        pas=(1-xsolR[1])/dis
        for j in range(dis+1):
            i=xsolR[1]+j*pas
            PT.append(vm.Point2D((splev(i,self.InvoluteControlR))))
        L=[primitives2D.RoundedLines2D(PT,{},False)]
        
        Lintersect=[vm.Point2D(splev(xsolT[0],self.TrochoidalSecondaryControlT))]
        Lintersect.append(vm.Point2D(splev(xsolT[1],self.InvoluteControlT)))
        Lintersect.append(vm.Point2D(splev(xsolR[0],self.TrochoidalSecondaryControlR)))
        Lintersect.append(vm.Point2D(splev(xsolR[1],self.InvoluteControlR)))
        Lintersect.append(vm.Point2D(splev(xsolI1[0],self.TrochoidalSecondaryControlT)))
        Lintersect.append(vm.Point2D(splev(xsolI1[1],self.DIControl)))
        Lintersect.append(vm.Point2D(splev(xsolI2[0],self.TrochoidalSecondaryControlR)))
        Lintersect.append(vm.Point2D(splev(xsolI2[1],self.DIControl)))
        Lintersect.append(vm.Point2D(splev(xsolE1[1],self.InvoluteControlT)))
        Lintersect.append(vm.Point2D(splev(xsolE1[0],self.DEControl)))
        Lintersect.extend(L)
        G2=vm.Contour2D(Lintersect)
        G2.MPLPlot()

        return L
    
    def WheelContours(self,pos=0,number=None):
        
        if number==None:
            number=range(int(self.tooth_number))
        L=self.GearContoursRack()
        Lint=[]
        for i in number:
            angle=i*2*npy.pi/self.tooth_number+pos
            temp=[]
            for j in L:
                Lint.append(j.Rotation(vm.Point2D((0,0)),angle))
        return Lint
    
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
                 coefficient_profile_shift2=0,gear_width=20,transverse_pressure_angle_rack_T1=None,
                 transverse_pressure_angle_rack_T2=None,circular_tooth_thickness_rack1=None,circular_tooth_thickness_rack2=None,
                 gear_addendum_rack1=None,gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
                 root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
                 root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.AssemblyGearParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
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
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
    def AssemblyGearParam(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
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
        
    def AssemblyGearParam2(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
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
        
        
        
    def Update(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.AssemblyGearParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
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
                       coefficient_profile_shift2,gear_width,transverse_pressure_angle_rack_T1,
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
            TG.extend(temp)
        return TG
        
        
class AssemblyGearOptimize1:
    def __init__(self,x=None):
        
        if x==None:
            x=self.Bounds()
        self.AssemblyGearParam(x)
        self.save=x[:]
        self.AssemblyGear=AssemblyGear(self.Z1,self.Z2,self.center_distance,self.transverse_pressure_angle,
                                       self.helix_angle,self.coefficient_profile_shift1,
                                       self.coefficient_profile_shift2,self.gear_width,self.transverse_pressure_angle_rack_T1,
                                       self.transverse_pressure_angle_rack_T2,
                                       self.circular_tooth_thickness_rack1,self.circular_tooth_thickness_rack2,
                                       self.gear_addendum_rack1,
                                       self.gear_addendum_rack2,self.gear_dedendum_rack1,self.gear_dedendum_rack2,
                                       self.root_radius_T1,self.root_radius_T2,self.root_radius_R1,
                                       self.root_radius_R2,self.transverse_pressure_angle_rack_R1,self.transverse_pressure_angle_rack_R2)
        
    def Bounds(self):
        
        bounds=[(30,30)] #bound p1 x
        bounds.append((50,56))
        bounds.append((60,150))
        bounds.append((0.01,0.5))
        bounds.append((-1,1))
        bounds.append((-1,1))
        bounds.append((0,0))
        bounds.append((20,20))
        self.bounds=npy.array(bounds)
        dim=npy.shape(bounds)[0]
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
                                       self.coefficient_profile_shift2,self.gear_width,self.transverse_pressure_angle_rack_T1,
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
            
            
    def MPLPlot(self,x,args):

        if args==0:
            #Trochoidal primary with two teeth
            self.AssemblyGear.Update(x)
            TG1=self.AssemblyGear.AssemblyGear.Gear1.GearContours(10,2)
            TG2=self.AssemblyGear.AssemblyGear.Gear2.GearContours(10,2)
            TrochoideT1,TrochoideR1=self.AssemblyGear.AssemblyGear.Gear1.TrochoidalRackEstimation()
            TrochoideT2,TrochoideR2=self.AssemblyGear.AssemblyGear.Gear2.TrochoidalRackEstimation()
            L=TrochoideT1
            L.extend(TrochoideR1)
            L.extend(TG1)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profil0)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profili)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilj)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilf)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.pti)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.ptj)
            G2=vm.Contour2D(L)
            G2.MPLPlot()
            
        if args==1:
        
            TG1=self.AssemblyGear.AssemblyGear.Gear1.GearContours(10,2)
            TG2=self.AssemblyGear.AssemblyGear.Gear2.GearContours(10,2)
            TG=self.AssemblyGear.AssemblyGear.WheelAssembly()
            c1=vm.Contour2D(TG)
            c1.MPLPlot()
            
            G2=[self.AssemblyGear.AssemblyGear.Gear1.rack_profil0]
            G2.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilf)
            G2.extend(TG1)
            G2.extend(self.AssemblyGear.AssemblyGear.Gear1.Ltrajet)
            G2=vm.Contour2D(G2)
            G2.MPLPlot()
            
            G2=[self.AssemblyGear.AssemblyGear.Gear2.rack_profil0]
            G2.append(self.AssemblyGear.AssemblyGear.Gear2.rack_profilf)
            G2.extend(TG2)
            G2.extend(self.AssemblyGear.AssemblyGear.Gear2.Ltrajet)
            G2=vm.Contour2D(G2)
            G2.MPLPlot()
            
            LR=self.AssemblyGear.AssemblyGear.Rack1.RackContours(10)
            R1=vm.Contour2D(LR)
            R1.MPLPlot()
            
            LR=self.AssemblyGear.AssemblyGear.Rack2.RackContours(10)
            R1=vm.Contour2D(LR)
            R1.MPLPlot()
            
            DistRI=self.AssemblyGear.AssemblyGear.Gear1.DistRackInvolute()
    
        
#Initialisation
A1=AssemblyGearOptimize1()
O1=Optimization(A1)
O1.Optimize()
#O1.Optimize2(list(O1.solution[-1]))

AG1=AssemblyGear(*list(O1.solution[-1]))
print(AG1.Gear1.root_diameter_active)
L=AG1.Gear1.TrochoideTrace(500,1,'T')
L.extend(AG1.Gear1.GearContours(10,1))
G1=vm.Contour2D(L)
G1.MPLPlot()
print(AG1.Gear1.root_diameter_active)

L=AG1.Rack1.RackContours(5)
G1=vm.Contour2D(L)
G1.MPLPlot()

xsol=list(O1.solution[-1])
TG1=AG1.Gear1.GearContours(10)
TG2=AG1.Gear2.GearContours(10)
list_rot=AG1.InitialPosition()
LR=AG1.GearAssemblyTrace([TG1,TG2],[(0,0),(AG1.center_distance,0)],list_rot)
LR.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DF1/2))
LR.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.DF2/2))
LR.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DB1/2))
LR.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.DB2/2))
LR.append(vm.Circle2D(vm.Point2D((0,0)),AG1.pitch_diameter_factory1/2))
LR.append(vm.Circle2D(vm.Point2D((AG1.center_distance,0)),AG1.pitch_diameter_factory2/2))
G1=vm.Contour2D(LR)
G1.MPLPlot()


#TG1=A1.AssemblyGear.Gear1.GearContours(10)
#TG2=A1.AssemblyGear.Gear2.GearContours(10)
#LR=A1.AssemblyGear.GearAssemblyTrace([TG1,TG2],[(0,0),(A1.AssemblyGear.center_distance,0)],[0.5,1])
#LR.append(vm.Circle2D(vm.Point2D((0,0)),A1.AssemblyGear.DF1/2))
#LR.append(vm.Circle2D(vm.Point2D((A1.AssemblyGear.center_distance,0)),A1.AssemblyGear.DF2/2))
#G1=vm.Contour2D(LR)
#G1.MPLPlot()

#O1.Optimize()
#print('calcul finalisé')
#O1.MPLPlot(O1.solution[-1],1)
            

    
#sol=npy.array([ 20.        ,  49.        ,  70.55799373,   6.20167685,
#         0.34906585,   0.26179939,   0.33565841,   3.17991653,
#         2.5373728 ,   3.00000002,   3.        ,   3.00000008,
#         3.00000001,   0.30000001,   1.08524723,   0.30000001,
#         1.08524723,   0.26179939,   0.33565841])
#A1=AssemblyGearOptimize1(sol)
#O1=Optimization(A1)
#O1.MPLPlot(sol,1)

#center_distance=100
#Z1=12
#Z2=46
#transverse_radial_pitch=6 #circular pitch on the pitch circle
#transverse_pressure_angle=18/180*npy.pi #real presure angle
#
#transverse_base_pitch=transverse_radial_pitch*npy.cos(transverse_pressure_angle)
#DF1=2*center_distance*Z1/Z2/(1+Z1/Z2)
#DF2=2*center_distance-DF1
#DB1=DF1*npy.cos(transverse_pressure_angle)
#DB2=DF2*npy.cos(transverse_pressure_angle)
#
##base_diameter=RackPitchDiameter*npy.cos(Racktransverse_pressure_angle)       
#        
##Rack definition
#transverse_pressure_angle_rack_T1=21/180*npy.pi
#transverse_pressure_angle_rack_T2=20/180*npy.pi
#transverse_radial_pitch_rack1=transverse_base_pitch/npy.cos(transverse_pressure_angle_rack_T1)
#transverse_radial_pitch_rack2=transverse_base_pitch/npy.cos(transverse_pressure_angle_rack_T2)
#Rack1=Rack(transverse_pressure_angle_rack_T1,transverse_radial_pitch_rack1,transverse_radial_pitch_rack1*0.6,1*transverse_radial_pitch_rack1/npy.pi,0.5*transverse_radial_pitch_rack1/npy.pi)
#Rack2=Rack(transverse_radial_pitch_rack2)
##LR=Rack2.RackContours(10)
##R1=vm.Contour2D(LR)
##R1.MPLPlot()
#
#Gear1=Gear(Z1,Rack2,0)
##LR=Gear1.GearContours(20)
##LR.append(vm.Circle2D(vm.Point2D((0,0)),Gear1.root_diameter/2))
##LR.append(vm.Circle2D(vm.Point2D((0,0)),Gear1.outside_diameter/2))
##LR.append(vm.Circle2D(vm.Point2D((0,0)),Gear1.pitch_diameter_factory/2))
##LR.append(vm.Circle2D(vm.Point2D((0,0)),Gear1.root_diameter_active/2))
##G1=vm.Contour2D(LR)
##G1.MPLPlot()
##
##Gear1=Gear(Z1,Rack1)
##Gear2=Gear(Z2,Rack2)
##Gear1.GearContoursRack(40)
#center_distance=60
#AG1=AssemblyGear(13,56,center_distance,0.34,0,0)
#TG1=AG1.Gear1.GearContours(10)
#TG2=AG1.Gear2.GearContours(10)
#LR=AG1.GearAssemblyTrace([TG1,TG2],[(0,0),(center_distance,0)],[0.5,1])
#LR.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DF1/2))
#LR.append(vm.Circle2D(vm.Point2D((center_distance,0)),AG1.DF2/2))
#LR.append(vm.Circle2D(vm.Point2D((0,0)),AG1.DB1/2))
#LR.append(vm.Circle2D(vm.Point2D((center_distance,0)),AG1.DB2/2))
#G1=vm.Contour2D(LR)
#G1.MPLPlot()